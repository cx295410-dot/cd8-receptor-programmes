from largepanel_core import *
from pypdf import PdfReader
import hashlib,xml.etree.ElementTree as ET
from scipy.stats import false_discovery_control
man=json.loads((OUT/'manifest_supp_1_15.json').read_text());assert len(man)==15;counts=[9,6,6,7,8,8,8,8,6,6,6,6,8,8,8];assert [r['panels'] for r in man]==counts and sum(counts)==108
report={'figures':{},'numeric_checks':{}}
for r in man:
 s=r['figure'];assert len(r['panel_ids'])==len(r['panel_details'])==r['panels'];assert set(r['panel_ids'])==set(string.ascii_lowercase[:r['panels']]);assert len(r['caption'])>200
 for v in r['panel_details'].values():
  assert v['identity'] and v['quantitative_meaning']
  for p in v['sources']:assert Path(p).exists(),p
 for ext,dpi in [('png',300),('tif',600)]:
  im=Image.open(OUT/'figures'/(s+'.'+ext));assert abs(float(im.info['dpi'][0])-dpi)<1;assert im.width==int(183/25.4*dpi)
  if ext=='tif':assert im.mode=='RGB'
 width=float(PdfReader(OUT/'figures'/(s+'.pdf')).pages[0].mediabox.width)*25.4/72;assert abs(width-183)<.01
 nt=len(ET.parse(OUT/'figures'/(s+'.svg')).findall('.//{http://www.w3.org/2000/svg}text'));assert nt>10;report['figures'][s]={'panels':r['panels'],'pdf_width_mm':width,'svg_text_nodes':nt,'source_paths_exist':True,'RGB_TIFF_600dpi_PNG_300dpi':True}
for j,c in enumerate(ALL):
 srcstem='FigS'+str(9+j);deststem='FigS'+str(5+j)
 for ext in ['png','svg','pdf','tif']:assert hashlib.sha256((REV/'figures'/(srcstem+'.'+ext)).read_bytes()).digest()==hashlib.sha256((OUT/'figures'/(deststem+'.'+ext)).read_bytes()).digest()
report['numeric_checks']['four_density_atlases_identical_all_four_formats']=True
for c in CO:
 a=pd.read_csv(OUT/'figure_source_data'/('FigS2_'+c+'_marker_z.csv'),index_col=0);b=pd.read_csv(OUT/'figure_source_data'/('Fig2_'+c+'_marker_z.csv'),index_col=0);assert a.shape==(38,7) and np.allclose(a,b)
 a=pd.read_csv(OUT/'figure_source_data'/('FigS9_'+c+'_target_programme_z.csv'),index_col=0);b=pd.read_csv(OUT/'figure_source_data'/('Fig4_'+c+'_receptor_profile_z.csv'),index_col=0);assert a.shape==(24,7) and np.allclose(a,b)
 ids=pd.read_csv(OUT/'figure_source_data'/('FigS9_'+c+'_display_cells.csv'));assert ids.cell_id.is_unique;assert ids.groupby(['patient','time','programme']).size().max()<=30
 mat=pd.read_csv(OUT/'figure_source_data'/('FigS9_'+c+'_balanced_cell_z.csv'),index_col=0);d=cells(c).set_index('cell_id').loc[ids.cell_id];z=d.log_CXCL13.to_numpy();z=(z-z.mean())/z.std();assert np.allclose(mat.loc['CXCL13'].to_numpy(),z)
 tab=pd.read_csv(OUT/'figure_source_data'/('FigS3_'+c+'_IFNA_running.csv'));es=tab.running_ES.iloc[np.argmax(np.abs(tab.running_ES))];original=old(c+'_Hallmark_GSEA').set_index('Term').loc['HALLMARK_INTERFERON_ALPHA_RESPONSE','ES'];assert np.isclose(es,original,atol=1e-10),(c,es,original)
 report['numeric_checks'][c+'_marker_target_cell_display_and_IFNA_ES_verified']={'display_cells':len(ids),'IFNA_running_ES':float(es),'frozen_ES':float(original)}
 # Recompute the original six-score BH family.
 ef=old('programme_effects');u=ef[ef.cohort.eq(c)&ef.endpoint.eq('score')];assert len(u)==6 and np.allclose(false_discovery_control(u.p.to_numpy()),u.q.to_numpy())
report['numeric_checks']['six_score_endpoint_BH_families_unchanged']=True
te=read('target_paired_effects')
for (c,endpoint),u in te.groupby(['cohort','endpoint']):assert len(u)==24 and np.allclose(false_discovery_control(u.p.to_numpy()),u.q.to_numpy())
report['numeric_checks']['all24_target_BH_families_unchanged']=True
seen=[]
for s in ['FigS13','FigS14','FigS15']:
 u=pd.read_csv(OUT/'figure_source_data'/(s+'_patient_target_detection.csv'));assert set(u.cohort)==set(CO);assert u.gene.nunique()==8;assert len(u)==400;seen+=list(u.gene.unique());original=read('gene_patient_visit');expected=original[original.cohort.isin(CO)&original.gene.isin(u.gene.unique())];pd.testing.assert_frame_equal(u.reset_index(drop=True),expected.reset_index(drop=True),check_dtype=False)
assert set(seen)==set(TARGETS) and len(seen)==24;report['numeric_checks']['24_targets_once_each_and_all1200_patient_target_visits_exact']=True
w=old('decomposition_patient').pivot(index=['cohort','patient','target'],columns='unused',values='composition') if False else old('decomposition_patient');assert np.allclose(w.composition+w.within_programme,w.supported_total)
report['numeric_checks']['frozen_decomposition_components_conserve_total']=True
(OUT/'work/verification_supp_1_15.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
