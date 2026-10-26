from pathlib import Path
import os,sys,json
ROOT=Path(r'D:\Downloads\Virtual_Biotech_Study_20261025\biology_rebuild_20261026');OUT=ROOT/'reference_expansion_20261026'
sys.path[:0]=[str(ROOT/'packages'),r'D:\CD8_Immunotherapy_20261023\analysis_extensions\packages']
for k in ['TEMP','TMP','MPLCONFIGDIR','NUMBA_CACHE_DIR']:os.environ[k]=str(OUT/'work')
os.environ['OMP_NUM_THREADS']='4';os.environ['OPENBLAS_NUM_THREADS']='4';os.environ['NUMBA_NUM_THREADS']='4'
import numpy as np,pandas as pd,scanpy as sc,anndata as ad
rows=[];audits=[]
P=['Memory','Cytotoxic','Dysfunction','Cycling','Interferon','Residency','Mixed']
for c in ['BCC','NSCLC','ORAL','BREAST']:
 d=pd.read_csv(ROOT/'results'/(c+'_atlas_coordinates.csv.gz'));x=np.load(ROOT/'data'/(c+'_PCA.npy'));assert len(d)==len(x)
 a=ad.AnnData(np.zeros((len(d),1),dtype=np.float32),obs=d.set_index('cell_id'));a.obs['programme']=pd.Categorical(a.obs.programme,categories=P);a.obsm['X_pca']=x
 sc.pp.neighbors(a,n_neighbors=15,n_pcs=30,use_rep='X_pca',random_state=20261026)
 sc.tl.paga(a,groups='programme');conn=a.uns['paga']['connectivities'].toarray();assert np.allclose(conn,conn.T) and (conn>=0).all() and (conn<=1+1e-8).all()
 for i,p in enumerate(P):
  for j,q in enumerate(P):rows.append({'cohort':c,'programme1':p,'programme2':q,'connectivity':conn[i,j]})
 audits.append({'cohort':c,'cells':len(a),'graph_edges':int(a.obsp['connectivities'].nnz),'n_neighbors':15,'n_pcs':30,'random_state':20261026,'groups':P,'interpretation':'undirected transcriptional neighbourhood connectivity; no lineage or temporal direction'})
 print('PAGA',c,len(a),flush=True)
pd.DataFrame(rows).to_csv(OUT/'results/paga_programme_connectivity.csv',index=False);(OUT/'Graph_Analysis_Audit.json').write_text(json.dumps(audits,indent=2),encoding='utf-8')
