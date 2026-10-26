"""Full-transcriptome, patient-replicated CD8 analysis. No outcome-driven tuning."""
from pathlib import Path
import os,sys,json,gzip,hashlib,time,shutil
ROOT=Path(r'D:\Downloads\Virtual_Biotech_Study_20261025\biology_rebuild_20261026')
OLD=Path(r'D:\CD8_Immunotherapy_20261023');PREV=ROOT.parent;EXT=OLD/'analysis_extensions'
sys.path[:0]=[str(ROOT/'packages'),str(EXT/'packages'),str(PREV/'scripts')]
for k in ['TEMP','TMP','MPLCONFIGDIR','NUMBA_CACHE_DIR']:os.environ[k]=str(ROOT/'work')
os.environ['OMP_NUM_THREADS']='4';os.environ['OPENBLAS_NUM_THREADS']='4';os.environ['NUMBA_NUM_THREADS']='4'
import numpy as np,pandas as pd
from scipy import sparse,stats
import anndata as ad
from evidence_core import scoped_clone,paired_inference,bh_adjust
def bootstrap_mean(v):
 v=np.asarray(v,float);v=v[np.isfinite(v)];rng=np.random.default_rng(20261026)
 bs=rng.choice(v,(5000,len(v)),replace=True).mean(1)
 return v.mean(),*np.quantile(bs,[.025,.975])
for p in ['data','results','work','references','figures','manuscript']:(ROOT/p).mkdir(exist_ok=True)
TARGETS=json.loads((PREV/'work/panel.json').read_text())['targets']
PROGRAMS={
 'Memory':['TCF7','LEF1','CCR7','SELL','MAL','LTB'],
 'Cytotoxic':['NKG7','GNLY','PRF1','GZMB','GZMH','FGFBP2','CTSW'],
 'Dysfunction':['CXCL13','TOX','TOX2','LAYN','PRDM1','BATF'],
 'Cycling':['MKI67','TOP2A','STMN1','TYMS','UBE2C','PCNA'],
 'Interferon':['ISG15','IFIT1','IFIT2','IFIT3','MX1','OAS1','IFI6','STAT1'],
 'Residency':['ZNF683','ITGA1','CD69','RUNX3','RGS1']}
assert not set(sum(PROGRAMS.values(),[]))&set(TARGETS)
MARKERS=list(dict.fromkeys(sum(PROGRAMS.values(),[])+TARGETS+['CD3D','CD3E','CD8A','CD8B','CD4','FOXP3','MS4A1','LYZ','LST1','HSPA1A','HSPA1B','DNAJB1']))
(ROOT/'work/programme_definitions.json').write_text(json.dumps({'programmes':PROGRAMS,'marker_exclusion':TARGETS,'score':'Mean of within-cohort gene z scores of log1p(CP10k), clipped to [-3,3]; at least half markers measured','classification':'largest score >=0.25 and top-second margin >=0.10, otherwise Mixed','neighbours':15,'PCs':30,'HVGs':2000,'Leiden_resolution':0.6,'seed':20261026,'exploratory':True},indent=2),encoding='utf-8')
SOURCES=[]
def record(p,role):
 p=Path(p);SOURCES.append({'path':str(p),'bytes':p.stat().st_size,'role':role})
 (ROOT/'work/source_manifest.json').write_text(json.dumps(SOURCES,indent=2),encoding='utf-8')
def eligible(m):
 n=m.groupby(['patient','time']).size();ok=n[n>=20].reset_index()[['patient','time']];m=m.merge(ok,on=['patient','time'],validate='many_to_one');pairs=m.groupby('patient').time.nunique();m=m[m.patient.isin(pairs[pairs==2].index)].copy();return m
def save_counts(x,m,genes,c):
 genes=['VSIR' if g=='C10orf54' else g for g in genes]
 if pd.Index(genes).duplicated().any():
  unique=list(dict.fromkeys(genes));lookup={g:i for i,g in enumerate(unique)}
  merge=sparse.csr_matrix((np.ones(len(genes),dtype=np.int32),(np.arange(len(genes)),[lookup[g] for g in genes])),shape=(len(genes),len(unique)))
  x=x@merge;genes=unique
 assert len(genes)==x.shape[1] and len(m)==x.shape[0] and not pd.Index(genes).duplicated().any() and m.cell_id.is_unique
 x=x.astype(np.int32).tocsr();assert x.min()>=0
 a=ad.AnnData(x,obs=m.set_index('cell_id'),var=pd.DataFrame(index=pd.Index(genes,name='gene')))
 a.obs['cohort']=c;a.write_h5ad(ROOT/'data'/(c+'_counts.h5ad'),compression='gzip');print('COUNTS',c,a.shape,'patients',a.obs.patient.nunique(),flush=True)
def extract_bcc():
 ip=OLD/'data/processed/BCC_matrix.json';info=json.loads(ip.read_text());mp=OLD/'data/processed/BCC_CD8_cell_scores.csv.gz';m=pd.read_csv(mp)
 # Refined author T-cell annotation is used, not the conflicting broad full-atlas label.
 m=m[['cell_id','patient','time','author_cluster','clone_beta','detected_genes','total_counts','mitochondrial_fraction']].rename(columns={'author_cluster':'author_state','clone_beta':'clone'})
 m=eligible(m);ci={c:i for i,c in enumerate(info['cells'])};cols=np.array([ci[c] for c in m.cell_id]);p=OLD/'data/processed/BCC_expression.float32.memmap';capacity=p.stat().st_size//4//len(info['cells']);arr=np.memmap(p,dtype='float32',mode='r',shape=(capacity,len(info['cells'])))
 blocks=[]
 for start in range(0,len(info['genes']),512):
  z=np.asarray(arr[start:min(start+512,len(info['genes'])),:][:,cols]);assert np.isfinite(z).all() and (z==np.floor(z)).all();blocks.append(sparse.csr_matrix(z.T))
 x=sparse.hstack(blocks,format='csr');save_counts(x,m,info['genes'],'BCC');record(ip,'Refined T-cell CD8 membership and full UMI gene names');record(mp,'Source refined CD8 metadata');record(p,'Full raw UMI expression memmap')
def extract_nsclc():
 base=EXT/'independent_GSE179994';ip=base/'matrix/matrix.json';info=json.loads(ip.read_text());p=base/'primary_paired_CD8_continuous_scores.csv.gz'
 m=pd.read_csv(p,usecols=['cell_id','patient','sample','cluster','clone_beta','time','available_matrix_UMI','detected_genes','mitochondrial_fraction']).rename(columns={'cluster':'author_state','clone_beta':'clone','available_matrix_UMI':'total_counts'});m=eligible(m)
 ci={c:i for i,c in enumerate(info['cells'])};cols=np.array([ci[c] for c in m.cell_id]);a=info['sparse_arrays'];ix=np.load(a['i']['path'],mmap_mode='r');v=np.load(a['x']['path'],mmap_mode='r');ptr=np.load(a['p']['path'],mmap_mode='r');mat=sparse.csc_matrix((v,ix,ptr),shape=tuple(info['shape']));x=mat[:,cols].T.tocsr();assert np.allclose(np.asarray(x.sum(1)).ravel(),m.total_counts);save_counts(x,m,info['genes'],'NSCLC');record(p,'Paired author CD8 metadata');record(ip,'Original reconstructed dgCMatrix');[record(d['path'],'Raw sparse count array') for d in a.values()]
def read_h5(p):
 import h5py
 with h5py.File(p) as f:
  h=f['matrix'];dec=lambda a:np.asarray([x.decode() if isinstance(x,bytes) else str(x) for x in a]);genes=dec(h['features/name'][:]);ft=dec(h['features/feature_type'][:]);bars=dec(h['barcodes'][:]);x=sparse.csc_matrix((h['data'][:],h['indices'][:],h['indptr'][:]),shape=tuple(h['shape'][:]))
 return x,genes,bars,ft
def extract_oral():
 base=EXT/'GSE200996_paired_validation';mp=base/'tables/paired_tumor_CD8_scores_TCR.csv.gz';m=pd.read_csv(mp);m=eligible(m);manifest=json.loads((base/'download_manifest.json').read_text());parts=[];meta=[];common=None
 for rec in [r for r in manifest if r['kind']=='GEX']:
  p=Path(rec['path']);u=m[m.patient.eq(rec['patient'])&m.time.eq(rec['time'])].copy();x,g,b,ft=read_h5(p);keep=ft=='Gene Expression';genes=g[keep].tolist();lookup={v:i for i,v in enumerate(b)};xx=x[keep,:][:,[lookup[v] for v in u.barcode_raw]].T.tocsr()
  assert common is None or common==genes;common=genes;assert np.allclose(np.asarray(xx.sum(1)).ravel(),u.matrix_UMI);parts.append(xx);meta.append(u);record(p,'Original oral tumor RNA H5');print('ORAL_LIBRARY',rec['patient'],rec['time'],len(u),flush=True)
 m=pd.concat(meta,ignore_index=True).rename(columns={'clone_beta':'clone','matrix_UMI':'total_counts'});save_counts(sparse.vstack(parts),m,common,'ORAL');record(mp,'Tumor CD8 author annotation and TCR')
def extract_breast():
 m=pd.read_csv(PREV/'data/protein_primary_cells.csv.gz');base=EXT/'GSE267552_protein_validation/raw';parts=[];meta=[];common=None
 for lane in ['A','B','C']:
  p=base/f'GSE267552_lane{lane}-CITE-GEX.h5';u=m[m.lane.eq(lane)].copy();x,g,b,ft=read_h5(p);keep=ft=='Gene Expression';genes=g[keep].tolist();assert common is None or common==genes;common=genes;ci={v:i for i,v in enumerate(b)};xx=x[keep,:][:,[ci[v] for v in u.barcode]].T.tocsr();assert np.allclose(np.asarray(xx.sum(1)).ravel(),u.rna_total_umi);parts.append(xx);meta.append(u);record(p,'Independent breast GEX/ADT original H5')
 m=pd.concat(meta,ignore_index=True);m['patient']=m.donor;m['time']=m.tissue;m['clone']='';save_counts(sparse.vstack(parts),m,common,'BREAST')
def process(c):
 import scanpy as sc
 print('PROCESS',c,flush=True);a=ad.read_h5ad(ROOT/'data'/(c+'_counts.h5ad'));m=a.obs.copy().reset_index();m['cohort']=c
 # Model raw whole-gene pseudobulk before any cell-level scaling.
 if c!='BREAST':
  group=m.patient.astype(str)+'__'+m.time.astype(str);lev=sorted(group.unique());design=pd.DataFrame({'sample':lev});design['patient']=design['sample'].str.split('__').str[0];design['condition']=design['sample'].str.split('__').str[1]
  membership=sparse.csr_matrix((np.ones(len(m),dtype=np.int32),([lev.index(v) for v in group],np.arange(len(m)))),shape=(len(lev),len(m)))
  pb=np.asarray((membership@a.X).todense());pd.DataFrame(pb,index=lev,columns=a.var_names).to_csv(ROOT/'data'/(c+'_pseudobulk_counts.csv.gz'));design.to_csv(ROOT/'data'/(c+'_pseudobulk_design.csv'),index=False)
  assert int(pb.sum())==int(a.X.sum())
 # Store unscaled per-cell data only for designated marker and target genes.
 idx=[g for g in MARKERS if g in a.var_names];raw=a[:,idx].X.toarray();tot=np.asarray(a.X.sum(1)).ravel();logged=np.log1p(raw*10000/np.maximum(tot[:,None],1));m['library_UMI']=tot
 for j,g in enumerate(idx):m['raw_'+g]=raw[:,j];m['log_'+g]=logged[:,j]
 scores=[];coverage=[]
 for prog,gs in PROGRAMS.items():
  use=[g for g in gs if g in idx];assert len(use)>=len(gs)/2;z=logged[:,[idx.index(g) for g in use]];z=(z-z.mean(0))/np.maximum(z.std(0),1e-5);score=np.clip(z,-3,3).mean(1);m['score_'+prog]=score;scores.append(score);coverage.append({'cohort':c,'programme':prog,'n_markers':len(use),'markers':';'.join(use)})
 s=np.asarray(scores).T;rank=np.sort(s,axis=1);labels=np.array(list(PROGRAMS))[s.argmax(1)];labels[(rank[:,-1]<.25)|((rank[:,-1]-rank[:,-2])<.1)]='Mixed';m['programme']=labels
 pd.DataFrame(coverage).to_csv(ROOT/'results'/(c+'_programme_coverage.csv'),index=False)
 # Cohort-specific graph, no mixing of cancer biology with cohort correction.
 a.X=a.X.astype(np.float32);sc.pp.normalize_total(a,target_sum=10000);sc.pp.log1p(a)
 sc.pp.highly_variable_genes(a,n_top_genes=2000,flavor='seurat',batch_key='patient')
 hv=a.var.highly_variable.to_numpy()&~a.var_names.str.match(r'^(MT-|RPL|RPS|TRAV|TRBV|IG[HKL])')
 b=a[:,hv].copy();sc.pp.scale(b,zero_center=False,max_value=10);sc.pp.pca(b,n_comps=30,zero_center=True,svd_solver='arpack',random_state=20261026)
 sc.pp.neighbors(b,n_neighbors=15,n_pcs=30,random_state=20261026);sc.tl.umap(b,min_dist=.35,random_state=20261026);sc.tl.leiden(b,resolution=.6,key_added='leiden',flavor='igraph',n_iterations=2,directed=False,random_state=20261026)
 m['UMAP1']=b.obsm['X_umap'][:,0];m['UMAP2']=b.obsm['X_umap'][:,1];m['leiden']=b.obs.leiden.to_numpy();np.save(ROOT/'data'/(c+'_PCA.npy'),b.obsm['X_pca']);pd.Series(b.var_names).to_csv(ROOT/'results'/(c+'_HVGs.csv'),index=False)
 m.to_csv(ROOT/'data'/(c+'_cells.csv.gz'),index=False);m[['cell_id','patient','time','cohort','programme','leiden','UMAP1','UMAP2']].to_csv(ROOT/'results'/(c+'_atlas_coordinates.csv.gz'),index=False)
 print('ATLAS_DONE',c,len(m),m.programme.value_counts().to_dict(),'Leiden',m.leiden.nunique(),flush=True)
def de(c):
 from pydeseq2.dds import DeseqDataSet
 from pydeseq2.ds import DeseqStats
 from pydeseq2.default_inference import DefaultInference
 import gseapy as gp
 x=pd.read_csv(ROOT/'data'/(c+'_pseudobulk_counts.csv.gz'),index_col=0);meta=pd.read_csv(ROOT/'data'/(c+'_pseudobulk_design.csv')).set_index('sample').loc[x.index]
 n=max(3,int(np.ceil(len(x)/4)));keep=(x>=10).sum(0)>=n;x=x.loc[:,keep].astype(int)
 assert x.index.equals(meta.index) and meta.groupby('patient').condition.nunique().eq(2).all()
 inference=DefaultInference(n_cpus=4);dds=DeseqDataSet(counts=x,metadata=meta,design='~ patient + condition',refit_cooks=True,inference=inference,quiet=True);dds.deseq2();ds=DeseqStats(dds,contrast=['condition','post','pre'],inference=inference,quiet=True);ds.summary();r=ds.results_df;r.to_csv(ROOT/'results'/(c+'_paired_DESeq2.csv'))
 norm=pd.DataFrame(dds.layers['normed_counts'],index=x.index,columns=x.columns);norm.to_csv(ROOT/'data'/(c+'_normalized_pseudobulk.csv.gz'))
 rank=r['stat'].replace([np.inf,-np.inf],np.nan).dropna().sort_values(ascending=False);rank.to_csv(ROOT/'results'/(c+'_Wald_ranking.csv'))
 gmt=ROOT/'references/MSigDB_Hallmark_2020.gmt'
 if not gmt.exists():shutil.copy2(EXT/'references/MSigDB_Hallmark_2020.gmt',gmt)
 pre=gp.prerank(rnk=rank,gene_sets=str(gmt),threads=4,permutation_num=1000,min_size=15,max_size=500,seed=20261026,outdir=None,verbose=False)
 import re
 def canonical(s):
  s=re.sub(r'[^A-Z0-9]+','_',s.upper()).strip('_')
  for a,b in [('TNF_ALPHA','TNFA'),('NF_KB','NFKB'),('IL_2','IL2'),('IL_6','IL6'),('G2_M','G2M')]:s=s.replace(a,b)
  return s if s.startswith('HALLMARK_') else 'HALLMARK_'+s
 pre.res2d.to_csv(ROOT/'results'/(c+'_Hallmark_GSEA_original_names.csv'),index=False)
 enrichment=pre.res2d.copy();enrichment['source_term']=enrichment.Term;enrichment['Term']=enrichment.Term.map(canonical);enrichment.to_csv(ROOT/'results'/(c+'_Hallmark_GSEA.csv'),index=False)
 curves=[]
 for term,res in pre.results.items():
  term=canonical(term)
  if term in ['HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_INTERFERON_ALPHA_RESPONSE','HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT','HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_OXIDATIVE_PHOSPHORYLATION']:
   curve=np.asarray(res['RES']);hits=set(res['hits']);curves.extend({'cohort':c,'term':term,'rank':j,'running_ES':v,'hit':j in hits} for j,v in enumerate(curve))
 pd.DataFrame(curves).to_csv(ROOT/'results'/(c+'_GSEA_curves.csv.gz'),index=False)
 (ROOT/'results'/(c+'_model_audit.json')).write_text(json.dumps({'n_samples':len(x),'patients':meta.patient.nunique(),'genes_tested':len(r),'filter_count':10,'filter_min_samples':n,'design':'~ patient + condition','contrast':'post versus pre','gene_FDR_005':int((r.padj<.05).sum()),'up':int(((r.padj<.05)&(r.log2FoldChange>0)).sum()),'down':int(((r.padj<.05)&(r.log2FoldChange<0)).sum()),'GSEA_permutations':1000},indent=2));print('DE_DONE',c,'genes',len(r),'FDR',int((r.padj<.05).sum()),flush=True)
def infer_delta(vals):
 v=np.asarray(vals,float);v=v[np.isfinite(v)];n=len(v)
 if n==0:return {'n':0,'mean':np.nan,'low':np.nan,'high':np.nan,'p':np.nan}
 signs=1-2*((np.arange(2**n)[:,None]>>np.arange(n))&1);p=np.mean(np.abs(signs@v/n)>=abs(v.mean())-1e-12);mean,lo,hi=bootstrap_mean(v)
 return {'n':n,'mean':mean,'low':lo,'high':hi,'p':p}
def summarise():
 prop=[];effects=[];dep=[];depeff=[];rhos=[];ct=[];flow=[];clonemet=[];support=[]
 for c in ['BCC','NSCLC','ORAL']:
  m=pd.read_csv(ROOT/'data'/(c+'_cells.csv.gz'));programmes=list(PROGRAMS)+['Mixed']
  for (p,t),u in m.groupby(['patient','time']):
   for s in programmes:prop.append({'cohort':c,'patient':p,'time':t,'programme':s,'fraction':u.programme.eq(s).mean(),'n_cells':len(u),'mean_score':u['score_'+s].mean() if s!='Mixed' else np.nan})
  for s in programmes:
   pp=pd.DataFrame(prop);pp=pp[pp.cohort.eq(c)&pp.programme.eq(s)];w=pp.pivot(index='patient',columns='time',values='fraction');effects.append({'cohort':c,'programme':s,'endpoint':'fraction',**infer_delta(w.post-w.pre)})
   if s!='Mixed':
    w=pp.pivot(index='patient',columns='time',values='mean_score');effects.append({'cohort':c,'programme':s,'endpoint':'score',**infer_delta(w.post-w.pre)})
  # Symmetric Kitagawa decomposition on the explicitly supported subpopulation.
  for p,u in m.groupby('patient'):
   counts=pd.crosstab(u.programme,u.time).reindex(index=programmes,columns=['pre','post'],fill_value=0).fillna(0);shared=counts.index[counts.min(axis=1)>=5];a=u[u.programme.isin(shared)];coverage=a.groupby('time').size()/u.groupby('time').size()
   support.append({'cohort':c,'patient':p,'pre_coverage':coverage.get('pre',0),'post_coverage':coverage.get('post',0),'n_supported_programmes':len(shared)})
   if len(shared)==0 or coverage.min()<.5:continue
   freq=pd.crosstab(a.programme,a.time).reindex(index=shared,columns=['pre','post']).astype(float);freq=freq/freq.sum(0)
   for g in TARGETS:
    mu=a.assign(value=a['raw_'+g]>0).groupby(['programme','time']).value.mean().unstack().loc[shared,['pre','post']]
    composition=((freq.post-freq.pre)*(mu.post+mu.pre)/2).sum();within=((mu.post-mu.pre)*(freq.post+freq.pre)/2).sum();total=(freq.post*mu.post-freq.pre*mu.pre).sum();full=u.assign(value=u['raw_'+g]>0).groupby('time').value.mean();assert abs(composition+within-total)<1e-10
    dep.append({'cohort':c,'patient':p,'target':g,'composition':composition,'within_programme':within,'supported_total':total,'full_total':full['post']-full['pre'],'pre_coverage':coverage.pre,'post_coverage':coverage.post})
  for p,u in m.groupby('patient'):
   for s in PROGRAMS:
    for g in TARGETS:
     rho=stats.spearmanr(u['score_'+s],u['log_'+g]).statistic if u['log_'+g].std()>0 else np.nan;rhos.append({'cohort':c,'patient':p,'programme':s,'target':g,'rho':rho})
  if c=='ORAL':m=m[~m.patient.eq('P23')]
  assigned=m[m.clone.notna()&m.clone.astype(str).ne('')].copy()
  # All clone identities remain scoped to both cohort and patient.
  for (p,t),u in assigned.groupby(['patient','time']):
   n=u.groupby('clone').size().to_numpy();q=n/n.sum();gini=np.abs(n[:,None]-n[None,:]).sum()/(2*len(n)*n.sum()) if len(n)<5000 else np.nan
   clonemet.append({'cohort':c,'patient':p,'time':t,'assigned_cells':len(u),'n_clones':len(n),'clonality':1+np.sum(q*np.log(q))/np.log(len(q)) if len(q)>1 else np.nan,'gini':gini,'top10_fraction':np.sort(q)[-10:].sum()})
  for p,u in assigned.groupby('patient'):
   baseline=set(u.loc[u.time.eq('pre'),'clone']);post=u[u.time.eq('post')].copy();post['history']=np.where(post.clone.isin(baseline),'Previously detected','Newly detected')
   for h,v in post.groupby('history'):
    for s in programmes:ct.append({'cohort':c,'patient':p,'history':h,'programme':s,'fraction':v.programme.eq(s).mean(),'n_cells':len(v)})
   dom=[]
   for (cl,t),v in u.groupby(['clone','time']):
    if len(v)<3:continue
    z=v.programme.value_counts();label=z.index[0] if len(z)==1 or z.iloc[0]>z.iloc[1] else 'Mixed';dom.append({'clone':cl,'time':t,'programme':label,'n_cells':len(v)})
   if not dom:continue
   w=pd.DataFrame(dom).pivot(index='clone',columns='time',values='programme')
   if not {'pre','post'}.issubset(w):continue
   w=w.dropna();
   for (s0,s1),n in w.groupby(['pre','post']).size().items():flow.append({'cohort':c,'patient':p,'pre_programme':s0,'post_programme':s1,'n_clones':int(n),'shared_clones':len(w),'patient_fraction':n/len(w)})
 for name,rows in [('programme_patient',prop),('programme_effects',effects),('decomposition_patient',dep),('decomposition_support',support),('receptor_programme_patient_rho',rhos),('clone_history_programme',ct),('shared_clone_links',flow),('clonality_patient',clonemet)]:
  d=pd.DataFrame(rows)
  if name=='programme_effects':d['q']=d.groupby(['cohort','endpoint']).p.transform(lambda x:bh_adjust(x.to_numpy(),family_size=len(x)))
  d.to_csv(ROOT/'results'/(name+'.csv'),index=False)
 dp=pd.DataFrame(dep)
 for (c,g),u in dp.groupby(['cohort','target']):
  for term in ['composition','within_programme','supported_total','full_total']:depeff.append({'cohort':c,'target':g,'term':term,**infer_delta(u[term])})
 pd.DataFrame(depeff).to_csv(ROOT/'results/decomposition_effects.csv',index=False)
 # Correlations across cells are descriptive; donor summaries are independent units.
 rr=pd.DataFrame(rhos);summary=rr.groupby(['cohort','programme','target']).rho.agg(['median','mean','count']).reset_index();summary.to_csv(ROOT/'results/receptor_programme_rho_summary.csv',index=False)
 b=pd.read_csv(ROOT/'data/BREAST_cells.csv.gz');p=pd.read_csv(PREV/'data/protein_primary_cells.csv.gz');missing=[g+'_adt_corrected' for g in TARGETS if g+'_adt_corrected' in p and g+'_adt_corrected' not in b]
 if missing:b=b.merge(p[['cell_id']+missing],on='cell_id',validate='one_to_one')
 rows=[]
 for donor,u in b.groupby('patient'):
  for s in PROGRAMS:
   for g in TARGETS:
    col=g+'_adt_corrected'
    if col in u:rows.append({'donor':donor,'programme':s,'target':g,'rho':stats.spearmanr(u['score_'+s],u[col]).statistic})
 pd.DataFrame(rows).to_csv(ROOT/'results/protein_programme_patient_rho.csv',index=False)
 # External functionally labelled clones, clone-equal within each donor.
 ext=pd.read_csv(EXT/'GSE176021_antigen_validation/tables/tumor_CD8_expression_TCR_functional_labels.csv.gz');out=[]
 for (p,cl,label),u in ext[ext.antigen.isin(['MANA','viral'])&ext.clone_beta.notna()].groupby(['patient','clone_beta','antigen']):
  if len(u)<3:continue
  row={'patient':p,'clone':cl,'antigen':label,'n_cells':len(u)}
  for g in ['IL7R','PDCD1','LAG3','ENTPD1','CXCL13','TCF7','GZMB','TOX']:row[g]=u['detected_'+g].mean()
  out.append(row)
 pd.DataFrame(out).to_csv(ROOT/'results/external_antigen_clone_detection.csv',index=False)
 summary={'cells':{},'programme_effects_q005':{},'decomposition_max_identity_error':float(abs(dp.composition+dp.within_programme-dp.supported_total).max())}
 for c in ['BCC','NSCLC','ORAL','BREAST']:
  d=pd.read_csv(ROOT/'results'/(c+'_atlas_coordinates.csv.gz'));summary['cells'][c]={'cells':len(d),'patients':d.patient.nunique(),'clusters':d.leiden.nunique(),'programmes':d.programme.value_counts().to_dict()}
 for c,u in pd.DataFrame(effects).groupby('cohort'):summary['programme_effects_q005'][c]=pd.read_csv(ROOT/'results/programme_effects.csv').query('cohort==@c and q<0.05')[['programme','endpoint','mean','q']].to_dict('records')
 (ROOT/'results/biological_analysis_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print('SUMMARY',json.dumps(summary),flush=True)

if __name__=='__main__':
 mode=sys.argv[1]
 if mode=='extract':
  for c,f in [('BCC',extract_bcc),('NSCLC',extract_nsclc),('ORAL',extract_oral),('BREAST',extract_breast)]:
   if not (ROOT/'data'/(c+'_counts.h5ad')).exists():f()
 elif mode=='process':
  for c in sys.argv[2:] or ['BCC','NSCLC','ORAL','BREAST']:process(c)
 elif mode=='de':
  for c in sys.argv[2:] or ['BCC','NSCLC','ORAL']:de(c)
 elif mode=='summarise':summarise()
