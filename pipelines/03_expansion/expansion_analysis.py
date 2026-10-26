"""Patient-replicated reference-guided expansion; no new outcome labels."""
from pathlib import Path
import os,sys,json,hashlib,itertools
ROOT=Path(r'D:\Downloads\Virtual_Biotech_Study_20261025\biology_rebuild_20261026')
OUT=ROOT/'reference_expansion_20261026';EXT=Path(r'D:\CD8_Immunotherapy_20261023\analysis_extensions')
sys.path.insert(0,str(ROOT/'packages'))
os.environ['OMP_NUM_THREADS']='4';os.environ['OPENBLAS_NUM_THREADS']='4'
import numpy as np,pandas as pd
from scipy import stats
from sklearn.metrics import roc_curve,roc_auc_score
CO=['BCC','NSCLC','ORAL'];ALL=CO+['BREAST'];P=['Memory','Cytotoxic','Dysfunction','Cycling','Interferon','Residency','Mixed']
DEFS=json.loads((ROOT/'work/programme_definitions.json').read_text())['programmes']
TARGETS=json.loads((ROOT.parent/'work/panel.json').read_text())['targets']
FOCAL=['IL7R','PDCD1','LAG3','TIGIT','ENTPD1','CTLA4'];GENES=list(dict.fromkeys(TARGETS+sum(DEFS.values(),[])+['CD8A','CD8B','CD3D','CD4','HSPA1A']))
def save(d,name):
 d.to_csv(OUT/'results'/(name+'.csv'),index=False);return d
def bh(p):
 p=np.asarray(p,float);q=np.full(len(p),np.nan);idx=np.flatnonzero(np.isfinite(p));order=idx[np.argsort(p[idx])];q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(order)+1))[::-1])[::-1]);return q
def paired(v):
 v=np.asarray(v,float);v=v[np.isfinite(v)];n=len(v)
 if not n:return {'n':0,'mean':np.nan,'low':np.nan,'high':np.nan,'p':np.nan}
 signs=1-2*((np.arange(2**n)[:,None]>>np.arange(n))&1)
 p=float(np.mean(np.abs(signs@v/n)>=abs(v.mean())-1e-12));rng=np.random.default_rng(20261026)
 b=rng.choice(v,(5000,n),replace=True).mean(1);lo,hi=np.quantile(b,[.025,.975])
 return {'n':n,'mean':v.mean(),'low':lo,'high':hi,'p':p}
def assigned(d,c):
 d=d[d.clone.notna()&d.clone.astype(str).ne('')].copy()
 if c=='ORAL':d=d[~d.patient.eq('P23')]
 d['scoped_clone']=c+'|'+d.patient.astype(str)+'|'+d.clone.astype(str)
 return d
def source_hash(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return {'path':str(p),'bytes':p.stat().st_size,'sha256':h.hexdigest()}
def run():
 for n in ['results','data','figure_source_data','scripts','references','work']:(OUT/n).mkdir(parents=True,exist_ok=True)
 markers=[];clones=[];repertoire=[];sharing=[];replacement=[];coexpr=[];crosswalk=[];qc=[];threshold=[];hetero=[];clinical=[];sources=[]
 for c in ALL:
  path=ROOT/'data'/(c+'_cells.csv.gz');d=pd.read_csv(path);assert d.cell_id.is_unique;sources.append(source_hash(path));genes=[g for g in GENES if 'log_'+g in d]
  if c=='NSCLC':
   cp=EXT/'independent_GSE179994/primary_paired_sample_clinical_metadata.csv';cm=pd.read_csv(cp);assert not cm.duplicated(['patient','time']).any()
   d=d.merge(cm[['patient','time','Response','Treatment','Biospy Site']],on=['patient','time'],validate='many_to_one',indicator=True);assert d._merge.eq('both').all();d=d.drop(columns='_merge');sources.append(source_hash(cp))
  for (p,t),u in d.groupby(['patient','time'],sort=True):
   for g in genes:markers.append({'cohort':c,'patient':p,'time':t,'gene':g,'n_cells':len(u),'mean_log':u['log_'+g].mean(),'detected_fraction':u['raw_'+g].gt(0).mean()})
   for prog in P:hetero.append({'cohort':c,'patient':p,'time':t,'programme':prog,'fraction':u.programme.eq(prog).mean(),'n_cells':len(u)})
   for field in ['library_UMI','detected_genes','mitochondrial_fraction','rna_n_genes','rna_mito_fraction']:
    if field in u:qc.append({'cohort':c,'patient':p,'time':t,'metric':field,'median':u[field].median(),'n_cells':len(u)})
   if c in CO:
    for g in ['CXCL13','CD8A','IL7R']:
     for k in [1,2,3]:threshold.append({'cohort':c,'patient':p,'time':t,'gene':g,'min_UMI':k,'fraction':u['raw_'+g].ge(k).mean(),'n_cells':len(u)})
   elif c=='BREAST':pass
   if c in ['NSCLC','ORAL']:
    fields=['Response'] if c=='NSCLC' else ['Path_response','regimen']
    for f in fields:assert u[f].nunique()==1
    clinical.append({'cohort':c,'patient':p,'time':t,'n_cells':len(u),**{f:u[f].iloc[0] for f in fields}})
  for p,u in d.groupby('patient'):
   corr=u[['log_'+g for g in TARGETS]].corr(method='spearman')
   for a,g in enumerate(TARGETS):
    for b,h in enumerate(TARGETS):coexpr.append({'cohort':c,'patient':p,'gene1':g,'gene2':h,'rho':corr.iloc[a,b]})
  cw=pd.crosstab(d.author_state.astype(str),d.programme).reindex(columns=P,fill_value=0)
  for state,row in cw.iterrows():
   for prog in P:crosswalk.append({'cohort':c,'author_state':state,'programme':prog,'n_cells':int(row[prog]),'fraction':row[prog]/row.sum()})
  if c=='BREAST':continue
  a=assigned(d,c)
  for (p,t),u in a.groupby(['patient','time']):
   total=len(u);ct=u.groupby('clone').size();q=ct/ct.sum();n=len(q)
   repertoire.append({'cohort':c,'patient':p,'time':t,'assigned_cells':total,'n_clones':n,'hill_D1':float(np.exp(-(q*np.log(q)).sum())),'inverse_simpson':1/float((q*q).sum()),'expanded_cell_fraction':ct[ct>=2].sum()/total,'large_clone_cell_fraction':ct[ct>=20].sum()/total,'top10_fraction':ct.nlargest(10).sum()/total,'tcr_coverage':total/len(d[d.patient.eq(p)&d.time.eq(t)])})
   for cl,v in u.groupby('clone',sort=False):
    row={'cohort':c,'patient':p,'time':t,'clone':str(cl),'scoped_clone':c+'|'+str(p)+'|'+str(cl),'n_cells':len(v),'assigned_cells':total,'frequency':len(v)/total}
    for g in ['IL7R','TCF7','CXCL13','PDCD1','LAG3','ENTPD1','TIGIT','CTLA4','GZMB','PRF1','NKG7','TOX','MKI67','STMN1','ZNF683','HSPA1A']:
     row[g]=v['log_'+g].mean()
    for prog in P:row['fraction_'+prog]=v.programme.eq(prog).mean()
    clones.append(row)
   # Network edges require >=2 cells from a clone in a programme, within patient and visit.
   tab=pd.crosstab(u.clone,u.programme).reindex(columns=P,fill_value=0)
   sets={s:set(tab.index[tab[s]>=2]) for s in P}
   for s0,s1 in itertools.product(P,P):
    union=sets[s0]|sets[s1];inter=sets[s0]&sets[s1]
    sharing.append({'cohort':c,'patient':p,'time':t,'programme1':s0,'programme2':s1,'shared':len(inter),'union':len(union),'jaccard':len(inter)/len(union) if union else np.nan,'n1':len(sets[s0]),'n2':len(sets[s1])})
  for p,u in a.groupby('patient'):
   pre=u[u.time.eq('pre')];post=u[u.time.eq('post')];baseline=set(pre.clone);postset=set(post.clone)
   replacement.append({'cohort':c,'patient':p,'pre_cells':len(pre),'post_cells':len(post),'pre_clones':len(baseline),'post_clones':len(postset),'new_post_clones':len(postset-baseline),'new_post_clone_fraction':len(postset-baseline)/len(postset) if postset else np.nan,'new_post_cell_fraction':(~post.clone.isin(baseline)).mean(),'shared_clones':len(postset&baseline)})
  print('ANALYSED',c,len(d),len(a),flush=True)
 md=save(pd.DataFrame(markers),'gene_patient_visit');cd=save(pd.DataFrame(clones),'clone_patient_visit');rep=save(pd.DataFrame(repertoire),'repertoire_patient_visit')
 save(pd.DataFrame(sharing),'programme_clone_sharing');save(pd.DataFrame(replacement),'clone_replacement');save(pd.DataFrame(coexpr),'target_coexpression_patient');save(pd.DataFrame(crosswalk),'author_state_crosswalk');save(pd.DataFrame(qc),'qc_patient_visit');save(pd.DataFrame(threshold),'detection_threshold_patient');save(pd.DataFrame(hetero),'programme_heterogeneity');save(pd.DataFrame(clinical),'clinical_patient_visit')
 effects=[]
 for (c,g),u in md[md.cohort.isin(CO)&md.gene.isin(TARGETS)].groupby(['cohort','gene']):
  for endpoint in ['detected_fraction','mean_log']:
   w=u.pivot(index='patient',columns='time',values=endpoint).dropna();effects.append({'cohort':c,'gene':g,'endpoint':endpoint,**paired(w.post-w.pre)})
 ef=pd.DataFrame(effects);ef['q']=ef.groupby(['cohort','endpoint']).p.transform(bh);save(ef,'target_paired_effects')
 re=[]
 for (c,),u in rep.groupby(['cohort']):
  for endpoint in ['hill_D1','inverse_simpson','expanded_cell_fraction','large_clone_cell_fraction','top10_fraction','tcr_coverage']:
   w=u.pivot(index='patient',columns='time',values=endpoint).dropna();re.append({'cohort':c,'endpoint':endpoint,**paired(w.post-w.pre)})
 re=pd.DataFrame(re);re['q']=re.groupby('cohort').p.transform(bh);save(re,'repertoire_paired_effects')
 # Independent source functional labels: never inferred from expression.
 ep=EXT/'GSE176021_antigen_validation/tables/tumor_CD8_expression_TCR_functional_labels.csv.gz';e=pd.read_csv(ep);sources.append(source_hash(ep))
 eligible=e[e.antigen.isin(['MANA','viral'])&e.clone_beta.notna()];rows=[]
 eg=['CXCL13','IL7R','TCF7','PDCD1','LAG3','ENTPD1','GZMB','TOX']
 for (p,cl,antigen),u in eligible.groupby(['patient','clone_beta','antigen']):
  if len(u)<3:continue
  rows.append({'patient':p,'clone':str(cl),'antigen':antigen,'n_cells':len(u),**{g:u['log_'+g].mean() for g in eg},'memory_score':u.memory_no_IL7R_rank.mean(),'exhaustion_score':u.exhaustion_rank.mean()})
 ex=pd.DataFrame(rows);qualified=[p for p,u in ex.groupby('patient') if all(u.antigen.value_counts().get(k,0)>=2 for k in ['MANA','viral'])];ex=ex[ex.patient.isin(qualified)].sort_values(['patient','antigen','clone']);save(ex,'external_functional_clone_logexpr')
 roc=[];rocpoints=[]
 for p,u in ex.groupby('patient'):
  y=u.antigen.eq('MANA').astype(int).to_numpy()
  for g in eg:
   s=u[g].to_numpy();auc=roc_auc_score(y,s);pos=s[y==1];neg=s[y==0];direct=((pos[:,None]>neg).sum()+.5*(pos[:,None]==neg).sum())/(len(pos)*len(neg));assert abs(auc-direct)<1e-12
   fpr,tpr,cut=roc_curve(y,s);roc.append({'patient':p,'gene':g,'AUC':auc,'direct_pairwise_AUC':direct,'MANA_clones':int(y.sum()),'viral_clones':int((1-y).sum())})
   for x,z,k in zip(fpr,tpr,cut):rocpoints.append({'patient':p,'gene':g,'fpr':x,'tpr':z,'threshold':k})
 save(pd.DataFrame(roc),'external_within_donor_auc');save(pd.DataFrame(rocpoints),'external_within_donor_roc')
 # Retain the source row identity for all numerical dotplots and original state labels.
 exid=e[e.patient.isin(qualified)&e.clone_beta.astype(str).isin(ex.clone)].copy();exid[['cell_id','patient','clone_beta','antigen']+[f'log_{g}' for g in eg]].to_csv(OUT/'data/external_qualified_cells.csv.gz',index=False)
 coverage=[
  ('Expression embedding / gene densities','Supported','Four frozen full-transcriptome cohorts; within-cohort graph coordinates'),
  ('Cell-level and clone-level expression heatmaps','Supported','Patient-balanced display subsets; full source tables retained'),
  ('Paired programmes and receptor fractions','Supported','25 independent patient pairs; no cell-level significance tests'),
  ('Clone sizes, replacement and sharing networks','Supported','Patient-scoped TCR; oral P23 excluded from paired TCR analyses'),
  ('Functional-label marker ROC','Descriptive only','Two qualified donors, 15 clones; within-donor ranking, no classifier fitting'),
  ('Clinical response comparisons in primary lung pairs','Unavailable contrast','All eight patients have original Response=Yes'),
  ('Oral pathology/regimen context','Descriptive only','Original Low/Medium/High and Mono/Combo retained; High has one patient'),
  ('Protein and RNA concordance','Supported','Five paired tissue donors, 17 assayed proteins; seven proteins remain unassayed'),
  ('ATAC accessibility tracks','Additional source required','No linked fragment/peak object in the four frozen RNA/ADT cohorts'),
  ('Tumour versus blood clone tracking','Additional source required','Frozen study populations are tumour/paired peritumour tissue; no paired blood compartment'),
  ('Residual tumour size regression','Additional source required','No verified continuous residual-tumour measurement joined to these cohorts'),
  ('CD4 subtype response panels','Different population','Current frozen biological study and programmes concern author-annotated CD8 cells')]
 save(pd.DataFrame(coverage,columns=['reference_analysis','status','evidence_basis']),'reference_feasibility')
 (OUT/'Source_Manifest.json').write_text(json.dumps(sources,indent=2),encoding='utf-8')
 summary={'cohorts':ALL,'cells':{c:len(pd.read_csv(ROOT/'results'/(c+'_atlas_coordinates.csv.gz'))) for c in ALL},'new_statistical_families':'24 target endpoints within each cohort and measurement type; six repertoire endpoints within each cohort','target_effects_FDR005':ef[ef.q<.05].to_dict('records'),'functional_qualified_donors':qualified,'functional_clones':len(ex),'functional_cells':int(ex.n_cells.sum()),'univariate_roc_direction':'Higher marker expression ranks MANA positive; direction never reversed to maximize AUC','NSCLC_primary_Response':'All Yes; no response contrast','density':'binned expression-weighted spatial density; within-cohort normalization, not a hypothesis test','models_refitted':False}
 (OUT/'Analysis_Summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print('DONE',json.dumps(summary),flush=True)
if __name__=='__main__':run()
