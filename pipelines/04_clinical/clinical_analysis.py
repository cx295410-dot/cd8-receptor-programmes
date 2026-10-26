"""Prespecified, patient-replicated pathology extension; no fitted predictors."""
from pathlib import Path
import sys,os,json,tarfile,gzip,hashlib
sys.stdout.reconfigure(encoding='utf-8');os.environ['OMP_NUM_THREADS']='2';os.environ['OPENBLAS_NUM_THREADS']='2'
ROOT=Path('D:/Downloads/Virtual_Biotech_Study_20261025/biology_rebuild_20261026');O=ROOT/'clinical_validation_20261026'
sys.path.insert(0,str(ROOT/'packages'))
import numpy as np,pandas as pd
from scipy.io import mmread
from scipy import sparse,stats
from threadpoolctl import threadpool_limits
SPEC=json.loads((O/'analysis_protocol_frozen.json').read_text());DEFS=SPEC['programmes'];TARGETS=SPEC['targets'];PROGS=list(DEFS)
GENES=sorted(set(TARGETS+sum(DEFS.values(),[])+['CD3D','CD3E','CD8A','CD8B']))
assert not set(TARGETS)&set(sum(DEFS.values(),[]))
def save(d,name,folder='results'):d.to_csv(O/folder/(name+'.csv'),index=False);return d
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def member(tar,name):
 m=[m for m in tar.getmembers() if m.isfile() and Path(m.name).name==name];assert len(m)==1
 return tar.extractfile(m[0])
def extract():
 if '--reuse-scored' in sys.argv:
  d=pd.read_csv(O/'data/CAUSHI_CD8_cells_frozen_scores.csv.gz');return patient_summary(d)
 anno=pd.read_csv(O/'references/GSE176021_CD8_annotations.csv.gz',index_col=0).rename(columns={'barcode':'cell_id','imid':'patient','CellType':'author_state'})
 assert anno.cell_id.is_unique
 manifest=json.loads((O/'download_manifest.json').read_text());cells=[];audits=[]
 cache=O/'work/library_features';cache.mkdir(exist_ok=True)
 for record in manifest:
  p=Path(record['path']);assert p.stat().st_size==record['bytes']
  dest=cache/(record['gsm']+'.csv.gz');qa=cache/(record['gsm']+'.json')
  if dest.exists() and qa.exists():
   cells.append(pd.read_csv(dest));audits.append(json.loads(qa.read_text()));print('CACHED',record['gsm'],flush=True);continue
  prefix=f"{record['patient']}:tumor-{record['lane']}_";a=anno[anno.cell_id.str.startswith(prefix)].copy()
  if len(a)==0:
   audit=dict(library=record['gsm'],patient=record['patient'],lane=record['lane'],author_CD8_cells=0,excluded_reason='No author CD8 cell with established patient ID; source annotation has NA:tumor-25 cells with missing imid, retained separately and not assigned by filename',sha256=sha(p),source=str(p))
   assert audit['sha256']==record['sha256'];audits.append(audit);print('ZERO_AUTHOR_CD8',record['gsm'],flush=True);continue
  with tarfile.open(p,mode='r:gz') as tar:
   features=pd.read_csv(gzip.GzipFile(fileobj=member(tar,'features.tsv.gz')),sep='\t',header=None)
   barcodes=pd.read_csv(gzip.GzipFile(fileobj=member(tar,'barcodes.tsv.gz')),sep='\t',header=None)[0].astype(str).tolist()
   with threadpool_limits(limits=2):mat=mmread(gzip.GzipFile(fileobj=member(tar,'matrix.mtx.gz')),spmatrix=True).tocsc()
  assert mat.shape==(len(features),len(barcodes));assert (mat.data>=0).all() and np.equal(mat.data,np.floor(mat.data)).all()
  lookup={prefix+b:i for i,b in enumerate(barcodes)};assert a.cell_id.isin(lookup).all()
  mat=mat[:,[lookup[c] for c in a.cell_id]];non_gex=0
  if features.shape[1]>=3:
   keep=features[2].eq('Gene Expression').to_numpy();non_gex=int((~keep).sum());mat=mat[keep,:];features=features[keep]
  symbols=features[1].astype(str).to_numpy();codes,unique=pd.factorize(symbols,sort=False);duplicates=len(symbols)-len(unique)
  if duplicates:mat=sparse.csr_matrix((np.ones(len(codes)),(codes,np.arange(len(codes)))),shape=(len(unique),len(codes)))@mat;symbols=np.asarray(unique)
  gi={g:i for i,g in enumerate(symbols)};missing=[g for g in GENES if g not in gi];assert not missing,missing
  total=np.asarray(mat.sum(axis=0)).ravel();assert (total>0).all();raw=mat[[gi[g] for g in GENES],:].T.toarray();log=np.log1p(raw*10000/total[:,None])
  values={}
  for j,g in enumerate(GENES):values['raw_'+g]=raw[:,j].astype(np.int32);values['log_'+g]=log[:,j]
  a=pd.concat([a.reset_index(drop=True),pd.DataFrame(values)],axis=1)
  a['library_UMI']=total;a['library']=record['gsm'];a['lane']=record['lane']
  # The original annotation stores the author's original full-count UMI total.
  relative=np.abs(a.nCount_RNA.to_numpy()-total)/np.maximum(total,1)
  audit=dict(library=record['gsm'],patient=record['patient'],lane=record['lane'],author_CD8_cells=len(a),raw_matrix_cells=len(barcodes),features=len(features),non_GEX_features_excluded=non_gex,duplicate_symbols_merged=duplicates,genes_requested=len(GENES),genes_measured=len(GENES),original_UMI_exact_matches=int(np.sum(relative<1e-10)),original_UMI_max_relative_difference=float(relative.max()),sha256=sha(p),source=str(p))
  assert audit['sha256']==record['sha256']
  a.to_csv(dest,index=False);qa.write_text(json.dumps(audit,indent=2));cells.append(a);audits.append(audit)
  print('EXTRACTED',record['gsm'],len(a),'CD8 cells',flush=True);del mat,raw,log
 d=pd.concat(cells,ignore_index=True);assert d.cell_id.is_unique
 tumor_anno=anno[anno.cell_id.str.contains(':tumor-',regex=False)].copy();unlinked=tumor_anno[~tumor_anno.patient.isin(d.patient.unique())].copy()
 save(unlinked,'unlinked_author_CD8_missing_patient','data')
 expected=tumor_anno[tumor_anno.patient.isin(d.patient.unique())];assert set(d.cell_id)==set(expected.cell_id)
 assert d.groupby('patient').size().ge(SPEC['minimum_cells']).all()
 mu=d[['log_'+g for g in GENES]].mean();sd=d[['log_'+g for g in GENES]].std(ddof=0)
 norm=pd.DataFrame({'gene':GENES,'mean_logCP10k':mu.values,'sd_logCP10k':sd.values});save(norm,'marker_normalization')
 z=(d[['log_'+g for g in GENES]].to_numpy()-mu.to_numpy())/np.maximum(sd.to_numpy(),1e-10);z=np.clip(z,-3,3);gi={g:i for i,g in enumerate(GENES)}
 ranks=d[['log_'+g for g in GENES]].rank(pct=True,method='average')
 for p,gs in DEFS.items():
  d['score_'+p]=z[:,[gi[g] for g in gs]].mean(1)
  d['rank_'+p]=ranks[['log_'+g for g in gs]].mean(1)
 scores=d[['score_'+p for p in PROGS]].to_numpy();order=np.argsort(scores,axis=1);best=scores[np.arange(len(d)),order[:,-1]];gap=best-scores[np.arange(len(d)),order[:,-2]]
 d['programme']=np.where((best>=.25)&(gap>=.1),np.asarray(PROGS)[order[:,-1]],'Mixed')
 d.to_csv(O/'data/CAUSHI_CD8_cells_frozen_scores.csv.gz',index=False);save(pd.DataFrame(audits),'library_matching_audit')
 return patient_summary(d)
def patient_summary(d):
 fields=['score_'+p for p in PROGS]+['rank_'+p for p in PROGS]+['log_'+g for g in GENES]
 pm=d.groupby('patient')[fields].mean();pm['n_cells']=d.groupby('patient').size();pm['n_libraries']=d.groupby('patient').library.nunique()
 lm=d.groupby(['patient','library'])[fields].mean().groupby('patient').mean().add_prefix('equal_library_');pm=pm.join(lm)
 cl=pd.read_csv(O/'data/clinical_original_Table1.csv');matched=pm.reset_index().merge(cl,on='patient',how='left',validate='one_to_one');assert matched.residual_tumor_pct.notna().all()
 matched['MPR']=matched.residual_tumor_pct.le(10);assert (matched.MPR==matched.author_MPR.str.strip().eq('MPR')).all()
 save(matched,'patient_clinical_features');save(d.groupby(['patient','author_state','programme']).size().rename('n_cells').reset_index(),'author_programme_crosswalk')
 return d,matched
def association(x,y,seed=20261026):
 x=np.asarray(x,float);y=np.asarray(y,float);n=len(x);rx=stats.rankdata(x);ry=stats.rankdata(y);rx-=rx.mean();ry-=ry.mean();den=np.linalg.norm(rx)*np.linalg.norm(ry)
 if den==0:return dict(n=n,rho=np.nan,p=np.nan,low=np.nan,high=np.nan,bootstrap_valid=0)
 observed=float(rx@ry/den);rng=np.random.default_rng(seed);perms=rng.permuted(np.broadcast_to(ry,(99999,n)),axis=1);rr=perms@rx/den
 p=(1+np.sum(np.abs(rr)>=abs(observed)-1e-12))/100000
 ix=rng.integers(0,n,(10000,n));bx=stats.rankdata(x[ix],axis=1);by=stats.rankdata(y[ix],axis=1);bx-=bx.mean(1,keepdims=True);by-=by.mean(1,keepdims=True)
 bn=np.linalg.norm(bx,axis=1)*np.linalg.norm(by,axis=1);valid=bn>0;br=np.sum(bx[valid]*by[valid],axis=1)/bn[valid];lo,hi=np.quantile(br,[.025,.975])
 return dict(n=n,rho=observed,p=float(p),low=float(lo),high=float(hi),bootstrap_valid=int(valid.sum()))
def bh(p):
 p=np.asarray(p,float);order=np.argsort(p);q=np.empty_like(p);q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q
def analyse(d,m):
 rows=[]
 for subset in ['primary','all_patients']:
  a=m[m.primary_analysis] if subset=='primary' else m
  for kind,features in [('programme',PROGS),('target',TARGETS)]:
   for f in features:
    col=('score_' if kind=='programme' else 'log_')+f;r=association(a[col],a.residual_tumor_pct)
    rows.append(dict(subset=subset,kind=kind,feature=f,**r))
 res=pd.DataFrame(rows);res['q_secondary']=np.nan
 for subset in ['primary','all_patients']:
  mask=res.subset.eq(subset)&~(res.kind.eq('programme')&res.feature.eq('Memory'));res.loc[mask,'q_secondary']=bh(res.loc[mask,'p'])
 save(res,'clinical_associations_prespecified')
 a=m[m.primary_analysis];sensitivity=[]
 for label,col in [('equal_library','equal_library_score_Memory'),('rank_scoring','rank_Memory')]:sensitivity.append(dict(analysis=label,**association(a[col],a.residual_tumor_pct)))
 save(pd.DataFrame(sensitivity),'primary_sensitivities')
 loo=[]
 for p in a.patient:
  u=a[a.patient.ne(p)];loo.append(dict(omitted_patient=p,n=len(u),rho=stats.spearmanr(u.score_Memory,u.residual_tumor_pct).statistic))
 save(pd.DataFrame(loo),'primary_leave_one_out')
 replic=[]
 for patient,u in d.groupby('patient'):
  replic.append(dict(patient=patient,n_cells=len(u),IL7R_Memory_rho=float(stats.spearmanr(u.log_IL7R,u.score_Memory).statistic),residual_tumor_pct=float(m.set_index('patient').loc[patient,'residual_tumor_pct'])))
 rep=save(pd.DataFrame(replic),'IL7R_Memory_patient_correlations');rng=np.random.default_rng(20261026);v=rep.IL7R_Memory_rho.to_numpy();b=rng.choice(v,(10000,len(v)),replace=True).mean(1)
 summary=dict(cohort='GSE176021',libraries=int(d.library.nunique()),downloaded_libraries=52,unlinked_author_CD8_cells=2602,patients=len(m),primary_patients=int(m.primary_analysis.sum()),CD8_cells=len(d),primary_CD8_cells=int(m.loc[m.primary_analysis,'n_cells'].sum()),MPR_patients=int(m.MPR.sum()),primary_MPR_patients=int(m.loc[m.primary_analysis,'MPR'].sum()),residual_range=[float(m.residual_tumor_pct.min()),float(m.residual_tumor_pct.max())],primary_memory=res[(res.subset=='primary')&(res.feature=='Memory')].iloc[0].to_dict(),all_patient_memory=res[(res.subset=='all_patients')&(res.feature=='Memory')].iloc[0].to_dict(),secondary_BH_significant_primary=res[(res.subset=='primary')&res.q_secondary.lt(.05)][['feature','kind','rho','p','q_secondary']].to_dict('records'),IL7R_Memory_mean=float(v.mean()),IL7R_Memory_CI=np.quantile(b,[.025,.975]).tolist(),IL7R_Memory_positive_patients=int((v>0).sum()),independent_replicate='patient',endpoint='post-treatment pathological residual viable tumor percent',protocol_sha256=sha(O/'analysis_protocol_frozen.json'))
 for field in ['primary_memory','all_patient_memory']:summary[field]['q_secondary']=None
 (O/'results/clinical_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False));print('SUMMARY',json.dumps(summary),flush=True)
 sources=[O/'analysis_protocol_frozen.json',O/'references/41586_2021_3752_MOESM8_ESM.xlsx',O/'references/GSE176021_CD8_annotations.rds.gz',O/'references/GSE176021_CD8_annotations.csv.gz']
 (O/'source_provenance.json').write_text(json.dumps([dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sources],indent=2));shutil_copy()
def shutil_copy():
 import shutil
 shutil.copy2(__file__,O/'scripts'/Path(__file__).name)
if __name__=='__main__':
 d,m=extract();analyse(d,m)
