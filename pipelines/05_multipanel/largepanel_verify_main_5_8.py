from largepanel_core import *
from pypdf import PdfReader
report={'figures':[],'scientific_checks':{}}
manifest=json.loads((OUT/'manifest_main_5_8.json').read_text())
assert [r['figure'] for r in manifest]==['Fig5','Fig6','Fig7','Fig8']
for r in manifest:
 stem=r['figure'];assert r['panels']==8 and r['panel_ids']==list('abcdefgh')
 for p in r['source']:assert Path(p).exists(),p
 png=Image.open(OUT/'figures'/(stem+'.png'));tif=Image.open(OUT/'figures'/(stem+'.tif'))
 assert tif.mode=='RGB' and min(tif.info['dpi'])==600
 assert 2150<png.width<2170 and 4310<tif.width<4340
 svg=(OUT/'figures'/(stem+'.svg')).read_text(encoding='utf-8');assert svg.count('<text')>30
 text=' '.join(p.extract_text() for p in PdfReader(OUT/'figures'/(stem+'.pdf')).pages);assert len(text)>400
 report['figures'].append({'figure':stem,'data_panels':8,'sources_exist':True,'png_px':png.size,'png_dpi':list(map(float,png.info['dpi'])),'tiff_px':tif.size,'tiff_dpi':list(map(float,tif.info['dpi'])),'tiff_mode':tif.mode,'editable_svg_text_nodes':svg.count('<text'),'pdf_extractable_text_chars':len(text)})
ex=read('external_functional_clone_logexpr');auc=read('external_within_donor_auc')
for r in auc.itertuples():
 u=ex[ex.patient.eq(r.patient)];v=u[u.antigen.eq('MANA')][r.gene].to_numpy();n=u[u.antigen.eq('viral')][r.gene].to_numpy();a=((v[:,None]>n).sum()+.5*(v[:,None]==n).sum())/(len(v)*len(n));assert np.isclose(a,r.AUC)
assert len(ex)==15 and ex.n_cells.sum()==404
for c in CO:
 m=pd.read_csv(OUT/'figure_source_data'/f'Fig5_{c}_flow.csv').iloc[:,1:];assert np.isclose(m.to_numpy().sum(),1)
spec=pd.read_csv(OUT/'figure_source_data/Fig5_repertoire_size_patient.csv');assert np.allclose(spec[['one','two_four','five_nineteen','twenty_plus']].sum(axis=1),1)
q=pd.read_csv(OUT/'figure_source_data/Fig7_qualified_functional_cells.csv');assert len(q)==404 and q.cell_id.is_unique
m=pd.read_csv(CLIN/'results/patient_clinical_features.csv');p=m[m.primary_analysis];rho=stats.spearmanr(p.score_Memory,p.residual_tumor_pct).statistic;assert abs(rho-.6703832596185881)<1e-12
ass=pd.read_csv(CLIN/'results/clinical_associations_prespecified.csv');primary=ass[ass.subset.eq('primary')&ass.feature.eq('Memory')].iloc[0];assert primary.p==.0152 and not ass[ass.subset.eq('primary')].q_secondary.lt(.05).any()
report['scientific_checks']={'functional_clones':15,'functional_cells':404,'functional_AUC_direct_concordance_verified':True,'flow_and_repertoire_weights_conserved':True,'clinical_primary_n':len(p),'primary_rho':float(rho),'primary_P':float(primary.p),'secondary_BH_significant':0,'new_genomewide_exploratory':True}
report['visual_QA']='Every main figure inspected; corrected flow labels, colorbar margin, density key, volcano labels, clinical heading/statistic overlaps.'
(OUT/'work/verification_main_5_8.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
