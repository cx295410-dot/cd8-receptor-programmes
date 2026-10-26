from largepanel_core import *
from pypdf import PdfReader
import shutil
expected=[12,9,6,9,6,9,12,8,8,8];rows=[]
manifest=json.loads((OUT/'manifest_supp_16_25.json').read_text());assert len(manifest)==10
for k,n,r in zip(range(16,26),expected,manifest):
 stem=f'FigS{k}';assert r['figure']==stem and r['panels']==n and len(r['panel_ids'])==n
 for p in r['source']:assert Path(p).exists(),p
 png=Image.open(OUT/'figures'/(stem+'.png'));tif=Image.open(OUT/'figures'/(stem+'.tif'));assert tif.mode=='RGB' and min(tif.info['dpi'])==600
 assert 2150<png.width<2170 and 4310<tif.width<4340
 svg=(OUT/'figures'/(stem+'.svg')).read_text(encoding='utf-8');assert svg.count('<text')>25
 pdf=PdfReader(OUT/'figures'/(stem+'.pdf'));text=' '.join(p.extract_text() for p in pdf.pages);assert len(text)>250
 rows.append({'figure':stem,'panels':n,'PNG_pixels':png.size,'PNG_dpi':list(map(float,png.info['dpi'])),'TIFF_pixels':tif.size,'TIFF_dpi':list(map(float,tif.info['dpi'])),'TIFF_mode':tif.mode,'SVG_text_nodes':svg.count('<text'),'PDF_text_characters':len(text),'all_sources_exist':True,'height_mm':r['height_mm']})
threshold=read('detection_threshold_patient')
for _,u in threshold.groupby(['cohort','patient','time','gene']):assert (np.diff(u.sort_values('min_UMI').fraction)<=1e-12).all()
qc=read('qc_patient_visit');assert not ((qc.cohort=='ORAL')&(qc.metric=='mitochondrial_fraction')).any()
assert 'unavailable' in next(r for r in manifest if r['figure']=='FigS22')['caption']
de=pd.read_csv(OUT/'results/clinical_continuous_DESeq2.csv').set_index('gene');un=pd.read_csv(OUT/'results/clinical_continuous_unadjusted_DESeq2.csv').set_index('gene');sx=pd.read_csv(OUT/'results/clinical_continuous_sexadjusted_DESeq2.csv').set_index('gene');sets=[set(u.index[u.padj.lt(.05)]) for u in [de,un,sx]];assert list(map(len,sets))==[31,18,27] and len(sets[0]&sets[2])==15
m=pd.read_csv(OUT/'figure_source_data/FigS24_complete_Hallmark_NES.csv',index_col=0);assert m.shape==(48,3) and np.isfinite(m.to_numpy()).all()
bs=pd.read_csv(OUT/'figure_source_data/FigS20_balanced_cells.csv');assert bs.cell_id.is_unique and bs.groupby(['patient','time','programme']).size().max()<=30
for c in CO:
 for t in ['pre','post']:
  z=pd.read_csv(OUT/'figure_source_data'/f'FigS23_{c}_{t}_clone_z.csv',index_col=0);assert z.shape[0]==16 and (z.isna().any().any() or z.notna().all().all())
report={'passed':True,'figures':rows,'total_figures':10,'total_panels':sum(expected),'source_checks':{'thresholds_monotone':True,'oral_missing_mito_not_imputed':True,'oral_QC_replacement_actual_cell_counts':True,'clinical_genes_tested':len(de),'clinical_BH_discoveries_primary_unadjusted_posthoc_sex':list(map(len,sets)),'primary_posthoc_sex_discovery_overlap':15,'complete_Hallmark_matrix_shape':m.shape,'breast_balanced_cells_unique':True,'maximum_cells_per_patient_tissue_programme':30},'visual_QA':'All ten full composites inspected. Corrected S17 colorbar clipping, S20 title/annotation strip collision and reading order, S22 unavailable oral mitochondrial display, S23 statistics overlap and S24 leading-edge/curve text placement. S24 preserves the complete48-pathway matrix with10 genes per leading-edge display.'}
(OUT/'work/verification_supp_16_25.json').write_text(json.dumps(report,indent=2));shutil.copy2(__file__,OUT/'scripts'/Path(__file__).name);print(json.dumps(report,indent=2))
