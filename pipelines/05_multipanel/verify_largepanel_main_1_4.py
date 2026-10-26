from largepanel_core import *
from pypdf import PdfReader
import xml.etree.ElementTree as ET
report={'figures':{},'source_checks':{}}
manifest=json.loads((OUT/'manifest_main_1_4.json').read_text())
assert [r['panels'] for r in manifest]==[7,7,8,8]
for r in manifest:
 stem=r['figure'];files={}
 assert len(r['panel_details'])==r['panels']
 for k,v in r['panel_details'].items():
  assert k in r['panel_ids']
  for p in v['sources']:assert Path(p).exists(),p
 for ext in ['png','tif']:
  im=Image.open(OUT/'figures'/(stem+'.'+ext));dpi=im.info.get('dpi');assert im.mode in ['RGB','RGBA'];assert abs(dpi[0]-(600 if ext=='tif' else 300))<1;assert im.width==(4322 if ext=='tif' else 2161);files[ext]={'mode':im.mode,'pixels':im.size,'dpi':tuple(float(v) for v in dpi)}
 svg=ET.parse(OUT/'figures'/(stem+'.svg'));nt=len(svg.findall('.//{http://www.w3.org/2000/svg}text'));assert nt>30;files['svg_editable_text_nodes']=nt
 page=PdfReader(OUT/'figures'/(stem+'.pdf')).pages[0];width=float(page.mediabox.width)*25.4/72;assert abs(width-183)<.01;files['pdf_width_mm']=width;report['figures'][stem]=files
w=pd.read_csv(OUT/'figure_source_data/Fig1_patient_composition.csv');assert np.allclose(w[P].sum(axis=1),1);assert len(w)==50;report['source_checks']['composition_50_paired_patient_visits_sum_to_one']=True
for c in CO:
 d=cells(c);genes=sum(DEFS.values(),[])
 for fname,gs in [('Fig2_'+c+'_marker_z',genes),('Fig4_'+c+'_receptor_profile_z',TARGETS)]:
  obs=pd.read_csv(OUT/'figure_source_data'/(fname+'.csv'),index_col=0);means=[]
  for prog in P:
   values=[]
   for patient,u in d[d.programme.eq(prog)].groupby('patient'):
    if len(u)>=5:values.append(u[['log_'+g for g in gs]].to_numpy().mean(0))
   means.append(np.mean(values,axis=0))
  raw=np.asarray(means).T;z=(raw-raw.mean(1,keepdims=True))/raw.std(1,keepdims=True);expected=pd.DataFrame(z,index=gs,columns=AB);assert np.allclose(obs,expected.loc[obs.index],atol=1e-10,equal_nan=True)
 report['source_checks'][c+'_patient_equal_marker_and_target_profiles_recomputed']=True
 edges=pd.read_csv(OUT/'figure_source_data'/('Fig4_'+c+'_coexpression_edges.csv'));assert np.all(edges.median_within_patient_rho.abs()>=.15)
 if len(edges):
  e=edges.iloc[0];rr=[stats.spearmanr(u['log_'+e.gene1],u['log_'+e.gene2]).statistic for _,u in d.groupby('patient')];assert np.isclose(np.nanmedian(rr),e.median_within_patient_rho,atol=1e-10);report['source_checks'][c+'_network_edge_recomputed_from_cells']={'gene1':e.gene1,'gene2':e.gene2,'median_rho':float(np.nanmedian(rr)),'n_edges':len(edges)}
for c in CO:
 z=np.load(OUT/'figure_source_data'/('Fig1_'+c+'_CXCL13_density.npz'))['relative_density'];assert z.shape==(110,110) and np.isclose(z.max(),1)
report['source_checks']['all_density_grids_110x110_and_unit_maximum']=True
for term in ['IFNG','E2F']:
 w=pd.read_csv(OUT/'figure_source_data'/('Fig3_'+term+'_leadingedge_z.csv'),index_col=0);assert w.shape==(14,50);assert np.isfinite(w).all().all()
report['source_checks']['leadingedge_matrices_two_14_genes_50_samples']=True
(OUT/'work/verification_main_1_4.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
