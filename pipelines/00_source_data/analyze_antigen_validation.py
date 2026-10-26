"""Read actual public UMI/VDJ files; join functional labels; infer at patient level."""
from pathlib import Path
import json,sys,os,io,gzip,tarfile,itertools,time
ROOT=Path('D:/CD8_Immunotherapy_20261023');OUT=ROOT/'analysis_extensions/GSE176021_antigen_validation'
os.environ['MPLCONFIGDIR']=str(ROOT/'cache');sys.path.insert(0,str(ROOT/'scripts'))
sys.stdout.reconfigure(encoding='utf-8',errors='replace')
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits
from analysis_core import SIGNATURES,bh_adjust
T=OUT/'tables';C=OUT/'scored_libraries';T.mkdir(exist_ok=True);C.mkdir(exist_ok=True)
MARKERS=['CD3D','CD3E','CD8A','CD8B','IL7R','TCF7','LEF1','CCR7','SELL','SLAMF6','ENTPD1','CXCL13','PDCD1','LAG3','HAVCR2','TIGIT','TOX','NKG7','PRF1','GZMB','IFNG','ITGAE','ZNF683','EOMES','TOX2','KLRG1']
ENDPOINTS=['log_IL7R','log_TCF7','memory_no_IL7R_rank','exhaustion_rank','cytotoxic_rank','IL7_induced_rank']
def get_vdj(record):
 b=Path(record['path']).read_bytes();layers=0
 while b.startswith(b'\x1f\x8b'):b=gzip.decompress(b);layers+=1
 assert 1<=layers<=3
 with tarfile.open(fileobj=io.BytesIO(b),mode='r:') as tar:
  names=[m for m in tar.getmembers() if m.isfile() and Path(m.name).name=='filtered_contig_annotations.csv']
  assert len(names)==1
  d=pd.read_csv(tar.extractfile(names[0]))
 return d
def member(tar,basename):
 matches=[m for m in tar.getmembers() if m.isfile() and Path(m.name).name==basename]
 assert len(matches)==1,(basename,[m.name for m in matches])
 return tar.extractfile(matches[0])
def truth_table():
 d=pd.read_csv(OUT/'functional_antigen_clonotypes_original_table8.csv')
 rows=[]
 for (patient,aa),g in d[~d.compound_dual_TRB_entry].groupby(['Patient','CDR3_AA_clean']):
  group=set('MANA' if s=='MANA' else 'viral' for s in g['Type of antigen recognized (viral or MANA)'])
  rows.append(dict(patient=patient,TRB_AA=aa,antigen=next(iter(group)) if len(group)==1 else 'ambiguous',source_antigen=';'.join(sorted(g['Type of antigen recognized (viral or MANA)'].unique())),source_epitope=';'.join(sorted(g['Epitope recognized in MANA/ViraFEST assay'].astype(str).unique())),source_double_star=bool(g['TCR clonotype (AA)'].str.contains(r'\*\*').any()),VDJdb_assignment=';'.join(sorted(g['Antigen specificity matched to VDJdb'].dropna().unique()))))
 result=pd.DataFrame(rows);assert not result[['patient','TRB_AA']].duplicated().any()
 result.to_csv(T/'functional_TRB_AA_label_dictionary.csv',index=False)
 return result
def main():
 manifest=json.loads((OUT/'download_manifest.json').read_text());assert all('error' not in r for r in manifest)
 assert len(manifest)==58
 anno=pd.read_csv(ROOT/'analysis_extensions/oncology_letters_extension/new_data/GSE176021_CD8_annotations.csv.gz',index_col=0)
 assert anno.barcode.is_unique;anno=anno.rename(columns={'barcode':'cell_id','imid':'patient','CellType':'author_state'})
 labels=truth_table()
 il7=pd.read_excel(OUT/'supplements/41586_2021_3752_MOESM4_ESM.xlsx',sheet_name='2.3 IL7 induced genes',header=None).iloc[2:,0].dropna().astype(str).str.strip().tolist()
 il7=sorted(set(g for g in il7 if g!='IL7R'))
 programs={'memory_no_IL7R_rank':[g for g in SIGNATURES['memory_progenitor'] if g!='IL7R'],'exhaustion_rank':SIGNATURES['exhaustion'],'cytotoxic_rank':SIGNATURES['cytotoxic'],'IL7_induced_rank':il7}
 (OUT/'analysis_specification.json').write_text(json.dumps(dict(endpoints=ENDPOINTS,programs=programs,normalization='raw UMI / sum of Gene Expression matrix UMI * 10000; log1p',scoring='custom positive-only average-tie rank capped at 1500, not exact UCell',primary='MANA minus viral within patient; clone equal; >=3 cells per unique TRB nt clone; >=1 eligible clone per category',tissue='original tumor only',functional_match='patient + original Table8 TRB CDR3 amino acid; compound dual-TRB source entry excluded; beta match lacks direct per-cell paired-alpha functional confirmation',statistical_family='6 exploratory endpoint tests, patient exact sign flip and 5000 bootstrap',E_access_performed=False),indent=2))
 allgenes=sorted(set(MARKERS+[g for v in programs.values() for g in v]));audits=[];cells=[];coverage=[]
 vdjs={(r['patient'],r['lane']):r for r in manifest if r['kind']=='VDJ'}
 gex=[r for r in manifest if r['kind']=='GEX'];assert len(gex)==29 and len(vdjs)==29
 for record in gex:
  patient=record['patient'];lane=record['lane'];prefix=f'{patient}:tumor-{lane}_';name=record['gsm'];cache=C/(name+'_CD8_scores.csv.gz');qa=C/(name+'_audit.json')
  if cache.exists() and qa.exists():
   cells.append(pd.read_csv(cache));audits.append(json.loads(qa.read_text()));print('CACHED',name,flush=True);continue
  a=anno[anno.cell_id.str.startswith(prefix)].copy();assert (a.patient==patient).all()
  with tarfile.open(record['path'],mode='r:gz') as tar:
   features=pd.read_csv(gzip.GzipFile(fileobj=member(tar,'features.tsv.gz')),sep='\t',header=None)
   barcodes=pd.read_csv(gzip.GzipFile(fileobj=member(tar,'barcodes.tsv.gz')),sep='\t',header=None)[0].astype(str).tolist()
   assert features[0].is_unique and len(barcodes)==len(set(barcodes))
   with threadpool_limits(limits=2):mat=mmread(gzip.GzipFile(fileobj=member(tar,'matrix.mtx.gz')),spmatrix=True).tocsc()
  assert mat.shape==(len(features),len(barcodes))
  assert np.isfinite(mat.data).all() and (mat.data>=0).all() and np.equal(mat.data,np.floor(mat.data)).all()
  raw_shape=mat.shape
  lookup={prefix+b:i for i,b in enumerate(barcodes)}
  present=a.cell_id.isin(lookup);missing_ids=a.loc[~present,'cell_id'].tolist();assert not missing_ids,(name,missing_ids[:5])
  mat=mat[:,[lookup[c] for c in a.cell_id]].copy()
  symbols=features[1].astype(str).tolist();feature_excluded=0
  if features.shape[1]>=3:
   keep=features[2].eq('Gene Expression').to_numpy();feature_excluded=int((~keep).sum());mat=mat[keep,:];symbols=np.array(symbols)[keep].tolist()
  codes,unique=pd.factorize(np.asarray(symbols),sort=False);duplicate_symbols=len(symbols)-len(unique)
  if duplicate_symbols:
   merge=sparse.csr_matrix((np.ones(len(codes)),(codes,np.arange(len(codes)))),shape=(len(unique),len(codes)))
   mat=(merge@mat).tocsc();symbols=list(unique)
  mat.sum_duplicates();mat.eliminate_zeros();mat.sort_indices()
  gi={g:i for i,g in enumerate(symbols)};available=[g for g in allgenes if g in gi];target=np.array([gi[g] for g in available]);li={g:i for i,g in enumerate(available)}
  for program,genes in programs.items():
   existing=[g for g in genes if g in gi];assert len(existing)>=3
   coverage.append(dict(library=name,program=program,present=len(existing),requested=len(genes),missing=';'.join(g for g in genes if g not in gi)))
  values=np.zeros((len(a),len(available)));ranks=np.zeros_like(values);total=np.asarray(mat.sum(axis=0)).ravel();assert (total>0).all()
  for k in range(len(a)):
   start,stop=mat.indptr[k:k+2];rows=mat.indices[start:stop];vals=mat.data[start:stop]
   ix=np.searchsorted(rows,target);valid=ix<len(rows);valid[valid]=rows[ix[valid]]==target[valid]
   ranks[k,valid]=((1501-np.minimum(rankdata(-vals,method='average'),1501))/1500)[ix[valid]];values[k,valid]=vals[ix[valid]]
  norm=np.log1p(values*10000/total[:,None])
  for g in MARKERS:
   assert g in li,(g,name)
   a['log_'+g]=norm[:,li[g]];a['detected_'+g]=values[:,li[g]]>0
  for program,genes in programs.items():a[program]=ranks[:,[li[g] for g in genes if g in li]].mean(axis=1)
  a['available_matrix_UMI']=total;a['detected_matrix_genes']=np.diff(mat.indptr);a['library']=name;a['lane']=lane
  v=get_vdj(vdjs[(patient,lane)])
  use=v.chain.eq('TRB')
  for c in ['is_cell','high_confidence','full_length','productive']:use &= v[c].astype(str).str.lower().eq('true')
  beta=v[use].dropna(subset=['cdr3','cdr3_nt']).drop_duplicates(['barcode','cdr3','cdr3_nt'])
  number=beta.groupby('barcode').size();single=number.index[number==1]
  beta=beta[beta.barcode.isin(single)][['barcode','cdr3','cdr3_nt']].rename(columns={'cdr3':'TRB_AA','cdr3_nt':'TRB_NT'})
  beta['cell_id']=prefix+beta.barcode;assert beta.cell_id.is_unique
  a=a.merge(beta.drop(columns='barcode'),on='cell_id',how='left',validate='one_to_one')
  a=a.merge(labels,on=['patient','TRB_AA'],how='left',validate='many_to_one');a['antigen']=a.antigen.fillna('unassigned');a['clone_beta']=np.where(a.TRB_NT.notna(),patient+'|TRB:'+a.TRB_NT.fillna(''),'')
  a['source_double_star']=a.source_double_star.fillna(False).astype(bool)
  audit=dict(library=name,patient=patient,lane=lane,source_matrix_shape=list(raw_shape),author_CD8_cells=len(a),matched_author_cells=int(present.sum()),unique_TRB_cells=int(a.TRB_NT.notna().sum()),multiple_TRB_barcodes=int((number>1).sum()),MANA_cells=int(a.antigen.eq('MANA').sum()),viral_cells=int(a.antigen.eq('viral').sum()),duplicate_symbols_merged=duplicate_symbols,non_GEX_features_excluded=feature_excluded)
  a.to_csv(cache,index=False);qa.write_text(json.dumps(audit,indent=2));cells.append(a);audits.append(audit)
  print('SCORED',json.dumps(audit),flush=True)
  del mat,values,ranks,norm
 d=pd.concat(cells,ignore_index=True);assert d.cell_id.is_unique
 d.to_csv(T/'tumor_CD8_expression_TCR_functional_labels.csv.gz',index=False)
 pd.DataFrame(audits).to_csv(T/'library_matching_audit.csv',index=False)
 # Rebuild coverage consistently even if library-scoring caches were reused.
 (OUT/'program_gene_sets.json').write_text(json.dumps(programs,indent=2))
 if coverage:pd.DataFrame(coverage).to_csv(T/'signature_coverage.csv',index=False)
 d.groupby(['patient','antigen']).agg(n_cells=('cell_id','size'),n_beta_nt=('clone_beta',lambda v:v[v!=''].nunique())).to_csv(T/'functional_label_patient_coverage.csv')
 summary={'libraries':len(audits),'patients':int(d.patient.nunique()),'author_CD8_cells':len(d),'unique_TRB_cells':int(d.TRB_NT.notna().sum()),'MANA_cells':int(d.antigen.eq('MANA').sum()),'viral_cells':int(d.antigen.eq('viral').sum()),'E_access_performed':False}
 (OUT/'expression_linkage_audit.json').write_text(json.dumps(summary,indent=2));print('LINKAGE_COMPLETE',summary,flush=True)
if __name__=='__main__':main()
