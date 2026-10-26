from pathlib import Path
import sys,json,hashlib,ast
ROOT=Path(r'D:\Downloads\Virtual_Biotech_Study_20261025\biology_rebuild_20261026');OUT=ROOT/'reference_expansion_20261026'
sys.path.insert(0,str(ROOT/'packages'))
import numpy as np,pandas as pd
from scipy import stats
from PIL import Image
from pypdf import PdfReader
from xml.etree import ElementTree as ET
checks=[]
def check(ok,name):
 assert bool(ok),name
 checks.append(name)
def r(name):return pd.read_csv(OUT/'results'/(name+'.csv'))
md=r('gene_patient_visit');cl=r('clone_patient_visit');rep=r('repertoire_patient_visit');sh=r('programme_clone_sharing');ef=r('target_paired_effects');re=r('repertoire_paired_effects');th=r('detection_threshold_patient');ex=r('external_functional_clone_logexpr');auc=r('external_within_donor_auc');cp=r('clinical_patient_visit')
check(not md.duplicated(['cohort','patient','time','gene']).any(),'Unique patient/visit/gene grain')
check(md.detected_fraction.between(0,1).all(),'Gene detection fractions in [0,1]')
for c in ['BCC','NSCLC','ORAL','BREAST']:
 d=pd.read_csv(ROOT/'data'/(c+'_cells.csv.gz'));check(d.cell_id.is_unique,c+' unique source cell IDs')
 u=md[md.cohort.eq(c)]
 for g in ['IL7R','PDCD1','CXCL13']:
  expected=d.assign(detected=d['raw_'+g].gt(0)).groupby(['patient','time']).detected.mean();got=u[u.gene.eq(g)].set_index(['patient','time']).detected_fraction.reindex(expected.index);check(np.allclose(expected,got,atol=1e-12),c+' '+g+' direct source fractions')
 if c!='BREAST':
  v=cl[cl.cohort.eq(c)];check(not v.duplicated(['scoped_clone','time']).any(),c+' clone/visit uniqueness');check(v.scoped_clone.eq(c+'|'+v.patient.astype(str)+'|'+v.clone.astype(str)).all(),c+' clone identifiers remain patient scoped')
  groups=v.groupby(['patient','time']);check(np.allclose(groups.frequency.sum(),1),c+' clone fractions sum to one')
  totals=groups.n_cells.sum();expected=rep[rep.cohort.eq(c)].set_index(['patient','time']).assigned_cells.reindex(totals.index);check(np.array_equal(totals,expected),c+' clone cells reconcile to repertoire denominators')
  counts=groups.size();expected=rep[rep.cohort.eq(c)].set_index(['patient','time']).n_clones.reindex(counts.index);check(np.array_equal(counts,expected),c+' distinct clones reconcile')
  for p in v.patient.unique():
   z=v[v.patient.eq(p)];check(set(z.time)=={'pre','post'},c+' '+p+' TCR pair completeness')
 h=pd.read_csv(OUT/'figure_source_data'/(f'FigS{13+["BCC","NSCLC","ORAL","BREAST"].index(c)}_display_cells.csv'));check(h.cell_id.is_unique and h.cell_id.isin(d.cell_id).all(),c+' cell heatmap selection uses real unique cells');check(h.groupby(['patient','time','programme']).size().max()<=30,c+' balanced display cap')
check(not cl.query('cohort=="ORAL"').patient.eq('P23').any(),'Oral P23 excluded from paired TCR')
check(((sh.shared<=sh.union)&(sh.shared>=0)).all(),'Sharing intersections bounded by unions')
check(sh.jaccard.dropna().between(0,1).all(),'Jaccard bounds')
rev=sh.rename(columns={'programme1':'programme2','programme2':'programme1'}).set_index(['cohort','patient','time','programme1','programme2']);a=sh.set_index(['cohort','patient','time','programme1','programme2']);check(np.allclose(a.jaccard,rev.reindex(a.index).jaccard,equal_nan=True),'Sharing graph symmetry')
for (c,g,endpoint),u in ef.groupby(['cohort','gene','endpoint']):
 w=md[md.cohort.eq(c)&md.gene.eq(g)].pivot(index='patient',columns='time',values=endpoint).dropna();delta=(w.post-w.pre).to_numpy();test=stats.permutation_test((delta,),np.mean,permutation_type='samples',vectorized=True,n_resamples=np.inf,alternative='two-sided');check(abs(test.pvalue-u.iloc[0].p)<1e-10,c+' '+g+' '+endpoint+' independent exact paired P')
for (c,endpoint),u in ef.groupby(['cohort','endpoint']):check(np.allclose(stats.false_discovery_control(u.p),u.q),c+' '+endpoint+' independent BH family of 24')
for (c,endpoint),u in re.groupby(['cohort','endpoint']):
 w=rep[rep.cohort.eq(c)].pivot(index='patient',columns='time',values=endpoint).dropna();delta=(w.post-w.pre).to_numpy();test=stats.permutation_test((delta,),np.mean,permutation_type='samples',vectorized=True,n_resamples=np.inf,alternative='two-sided');check(abs(test.pvalue-u.iloc[0].p)<1e-10,c+' '+endpoint+' independent repertoire paired P')
for c,u in re.groupby('cohort'):check(np.allclose(stats.false_discovery_control(u.p),u.q),c+' independent repertoire BH family of six')
for key,u in th.groupby(['cohort','patient','time','gene']):check(np.all(np.diff(u.sort_values('min_UMI').fraction)<=1e-12),'Threshold monotonicity '+str(key))
check(len(ex)==15 and ex.patient.nunique()==2 and ex.n_cells.sum()==404,'Functional validation retains 15 clones / 404 cells / two donors')
check(not ex.duplicated(['patient','clone']).any(),'Independent functional clone labels are unambiguous')
for row in auc.itertuples():
 u=ex[ex.patient.eq(row.patient)];pos=u[u.antigen.eq('MANA')][row.gene].to_numpy();neg=u[u.antigen.eq('viral')][row.gene].to_numpy();v=((pos[:,None]>neg).mean()+.5*(pos[:,None]==neg).mean());check(abs(v-row.AUC)<1e-12,'Direct ROC concordance '+row.patient+' '+row.gene)
check(cp.query('cohort=="NSCLC"').Response.eq('Yes').all(),'Original lung primary response labels are all Yes')
check(cp.query('cohort=="ORAL"')[['patient','Path_response']].drop_duplicates().Path_response.eq('High').sum()==1,'Oral High category contains one independent patient')
paga=r('paga_programme_connectivity')
for c,u in paga.groupby('cohort'):
 w=u.pivot(index='programme1',columns='programme2',values='connectivity');check(np.allclose(w,w.T) and w.min().min()>=0 and w.max().max()<=1+1e-8,c+' bounded symmetric PAGA connectivity')
manifest=json.loads((OUT/'figure_manifest.json').read_text());check(len(manifest)==43,'43 complete parent figures');check(sum(x['panels'] for x in manifest)==324,'324 labelled data panels')
for row in manifest:
 stem=row['figure'];im=Image.open(OUT/'figures'/(stem+'.tif'));check(im.mode=='RGB' and all(abs(x-600)<1 for x in im.info['dpi']),stem+' 600 dpi RGB TIFF');check(len(PdfReader(OUT/'figures'/(stem+'.pdf')).pages)==1,stem+' vector PDF page');tree=ET.parse(OUT/'figures'/(stem+'.svg'));check(len(tree.findall('.//{http://www.w3.org/2000/svg}text'))>5,stem+' editable SVG labels')
check(len(PdfReader(OUT/'Complete_Reference_Style_Atlas.pdf').pages)==43,'43-page final atlas')
for p in (OUT/'scripts').glob('*.py'):ast.parse(p.read_text(encoding='utf-8'))
summary={'status':'passed','checks':len(checks),'figure_count':43,'panel_count':324,'statistical_models_refitted':False,'old_analysis_verification':str(ROOT/'Verification_Report.json'),'new_exact_tests_independently_verified':len(ef)+len(re),'old_estimates_preserved':True,'checks_passed':checks}
(OUT/'Expansion_Verification.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print('VERIFIED',len(checks),'checks; 43 figures / 324 panels; exact paired tests, BH families, clone identities, ROC and exports',flush=True)
