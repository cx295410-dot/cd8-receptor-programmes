"""Independent checks of count provenance, biological summaries and figure exports."""
from pathlib import Path
import os,sys,json,hashlib,re
ROOT=Path(r'D:\Downloads\Virtual_Biotech_Study_20261025\biology_rebuild_20261026');PREV=ROOT.parent
sys.path[:0]=[str(ROOT/'packages'),str(Path(r'D:\CD8_Immunotherapy_20261023\analysis_extensions\packages'))]
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from PIL import Image
from pypdf import PdfReader
from xml.etree import ElementTree as ET
OUT={};checks=[]
def check(name,condition,details=None):
 assert bool(condition),(name,details)
 checks.append({'check':name,'passed':True,'details':details})
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for z in iter(lambda:f.read(2**22),b''):h.update(z)
 return h.hexdigest()
defs=json.loads((ROOT/'work/programme_definitions.json').read_text());target=json.loads((PREV/'work/panel.json').read_text())['targets'];markers=sum(defs['programmes'].values(),[])
check('38 receptor-independent programme markers',len(markers)==38 and not set(markers)&set(target))
def canonical(s):
 s=re.sub(r'[^A-Z0-9]+','_',s.upper()).strip('_')
 for a,b in [('TNF_ALPHA','TNFA'),('NF_KB','NFKB'),('IL_2','IL2'),('IL_6','IL6'),('G2_M','G2M')]:s=s.replace(a,b)
 return s if s.startswith('HALLMARK_') else 'HALLMARK_'+s
gene_sets={}
for line in (ROOT/'references/MSigDB_Hallmark_2020.gmt').read_text().splitlines():
 q=line.split('\t');gene_sets[canonical(q[0])]=set(q[2:])
eligibility=[]
expected={'BCC':(12364,23309,11,13),'NSCLC':(36172,19790,8,12),'ORAL':(8767,36591,6,12),'BREAST':(5061,36591,5,9)}
matrix_hashes=[]
for c,(n,g,npairs,ncluster) in expected.items():
 a=ad.read_h5ad(ROOT/'data'/(c+'_counts.h5ad'));m=pd.read_csv(ROOT/'data'/(c+'_cells.csv.gz'));a=a[m.cell_id.astype(str),:]
 check(c+' matrix shape and unique IDs',a.shape==(n,g) and m.cell_id.is_unique and a.var_names.is_unique)
 check(c+' nonnegative integer counts',a.X.data.min()>=0 and np.equal(a.X.data,np.floor(a.X.data)).all())
 check(c+' donor count and cluster count',m.patient.nunique()==npairs and m.leiden.nunique()==ncluster)
 check(c+' finite embedding and scores',np.isfinite(m[['UMAP1','UMAP2']+['score_'+p for p in defs['programmes']]].to_numpy()).all())
 check(c+' exact recorded library depths',np.array_equal(np.asarray(a.X.sum(1)).ravel(),m.library_UMI.to_numpy()))
 for t in target:check(c+' raw target '+t,np.array_equal(a[:,t].X.toarray().ravel(),m['raw_'+t].to_numpy()))
 if c!='BREAST':
  old=pd.read_csv(PREV/'data'/(c+'_target_cells.csv.gz')).set_index('cell_id').loc[m.cell_id]
  # The fixed-panel extraction stores gene symbols directly as raw-count columns.
  for t in target:
   col='raw_'+t if 'raw_'+t in old else t+'_count' if t+'_count' in old else t
   check(c+' previous independent panel count '+t,np.array_equal(old[col].to_numpy(),m['raw_'+t].to_numpy()))
  visit=m.groupby('patient').time.nunique();check(c+' complete pairs',visit.eq(2).all() and len(visit)==npairs)
  pb=pd.read_csv(ROOT/'data'/(c+'_pseudobulk_counts.csv.gz'),index_col=0);group=m.patient.astype(str)+'__'+m.time.astype(str)
  for label in pb.index:
   calc=np.asarray(a[group.eq(label).to_numpy(),:].X.sum(0)).ravel();check(c+' pseudobulk '+label,np.array_equal(calc,pb.loc[label].to_numpy()))
  de=pd.read_csv(ROOT/'results'/(c+'_paired_DESeq2.csv'),index_col=0);minimum=max(3,int(np.ceil(len(pb)/4)));retained=pb.columns[(pb>=10).sum(0)>=minimum]
  check(c+' expression filter and tested family',set(retained)==set(de.index))
  ns={'BCC':1,'NSCLC':424,'ORAL':241};check(c+' significant genes',int((de.padj<.05).sum())==ns[c])
  gsea=pd.read_csv(ROOT/'results'/(c+'_Hallmark_GSEA.csv'));rank=pd.read_csv(ROOT/'results'/(c+'_Wald_ranking.csv'),index_col=0);sizes={k:len(v&set(rank.index)) for k,v in gene_sets.items()};eligible={k for k,v in sizes.items() if 15<=v<=500}
  check(c+' all eligible Hallmarks from 50-set collection',set(gsea.Term)==eligible and gsea.Term.is_unique,{'tested':len(gsea),'candidate':len(gene_sets)})
  eligibility.extend({'cohort':c,'term':k,'matched_genes':v,'eligible':k in eligible} for k,v in sizes.items())
  before=pd.read_csv(ROOT/'results'/(c+'_Hallmark_GSEA_original_names.csv'))
  for field in ['NES','ES','FDR q-val']:check(c+' unchanged original enrichment '+field,np.allclose(gsea[field],before[field],equal_nan=True))
  OUT[c]={'cells':n,'patients':npairs,'genes':g,'clusters':ncluster,'genes_tested':len(de),'FDR005_genes':ns[c],'Hallmark_sets_tested':len(gsea)}
 else:
  check('Breast tissue counts',m.tissue.value_counts().sort_values().tolist()==[1621,3440])
 matrix_hashes.append({'file':str(ROOT/'data'/(c+'_counts.h5ad')),'sha256':sha(ROOT/'data'/(c+'_counts.h5ad'))})
 print('COUNTS VERIFIED',c,flush=True)
decomp=pd.read_csv(ROOT/'results/decomposition_patient.csv');err=float(abs(decomp.composition+decomp.within_programme-decomp.supported_total).max())
check('All 600 primary decomposition identities',len(decomp)==25*24 and err<1e-12,err)
check('Primary support threshold',decomp[['pre_coverage','post_coverage']].min().min()>=.5)
sen=pd.read_csv(ROOT/'results/decomposition_sensitivity_patient.csv');check('Six state and support specifications',len(sen[['state_definition','minimum_state_cells']].drop_duplicates())==6)
check('All sensitivity identities',abs(sen.composition+sen.within_programme-sen.supported_total).max()<1e-12)
prog=pd.read_csv(ROOT/'results/programme_effects.csv');check('Programme fraction/score test families',len(prog)==39 and not (prog.q<.05).any())
flows=pd.read_csv(ROOT/'results/shared_clone_links.csv')
for c,ncl,npat in [('BCC',97,9),('NSCLC',214,8),('ORAL',111,5)]:
 u=flows[flows.cohort.eq(c)];check(c+' shared clone counts',u.n_clones.sum()==ncl and u.patient.nunique()==npat)
 check(c+' clone-equal patient weights',np.allclose(u.groupby('patient').patient_fraction.sum(),1))
protein=pd.read_csv(ROOT/'results/protein_programme_patient_rho.csv');check('Independent protein programme associations',len(protein)==5*6*17 and protein.donor.nunique()==5 and protein.target.nunique()==17)
partial=pd.read_csv(ROOT/'results/selected_partial_associations.csv');rho=float(partial.loc[partial.cohort.eq('BREAST')&partial.programme.eq('Memory')&partial.target.eq('IL7R'),'rho'].iloc[0]);check('Adjusted independent memory-CD127 association',np.isclose(rho,.245057,atol=1e-6),rho)
ant=pd.read_csv(ROOT/'results/external_antigen_clone_detection.csv');qual=[p for p,u in ant.groupby('patient') if all(u.antigen.value_counts().get(k,0)>=2 for k in ['MANA','viral'])];u=ant[ant.patient.isin(qual)];check('Source antigen context 2 donors 15 clones 404 cells',len(qual)==2 and len(u)==15 and int(u.n_cells.sum())==404)
curves=json.loads((ROOT/'results/GSEA_curve_verification.json').read_text());check('18 independently reconstructed GSEA curves',len(curves)==18 and max(q['independent_running_ES_error'] for q in curves)<.002)
OUT['maximum_decomposition_error']=err;OUT['maximum_GSEA_ES_error']=max(q['independent_running_ES_error'] for q in curves)
figs=['Fig'+str(i) for i in range(1,8)]+['FigS'+str(i) for i in range(1,9)]
for f in figs:
 p=ROOT/'figures'/(f+'.tif');im=Image.open(p);check(f+' 600 dpi RGB TIFF',im.mode=='RGB' and im.width>=4320 and all(abs(d-600)<1 for d in im.info['dpi']))
 tree=ET.parse(ROOT/'figures'/(f+'.svg'));ns={'s':'http://www.w3.org/2000/svg'};check(f+' native editable SVG text',len(tree.findall('.//s:text',ns))>0 and len(tree.findall('.//s:image',ns))==0)
 reader=PdfReader(ROOT/'figures'/(f+'.pdf'));check(f+' one-page vector PDF',len(reader.pages)==1 and float(reader.pages[0].mediabox.width)>510)
pd.DataFrame(eligibility).to_csv(ROOT/'results/Hallmark_eligibility.csv',index=False)
OUT['figure_files']=len(figs)*4;OUT['checks']=checks;OUT['status']='passed';OUT['full_count_checksums']=matrix_hashes
(ROOT/'Verification_Report.json').write_text(json.dumps(OUT,indent=2),encoding='utf-8');print('PASSED',len(checks),'checks',flush=True)
