from pathlib import Path
import sys,json,hashlib,shutil
ROOT=Path('D:/Downloads/Virtual_Biotech_Study_20261025/biology_rebuild_20261026');O=ROOT/'clinical_validation_20261026';sys.path.insert(0,str(ROOT/'packages'));sys.stdout.reconfigure(encoding='utf-8')
import numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from PIL import Image
from pypdf import PdfReader
CHECK=[]
def check(name,condition):
 CHECK.append(dict(check=name,passed=bool(condition)))
 if not condition:raise AssertionError(name)
def run():
 spec=json.loads((O/'analysis_protocol_frozen.json').read_text());m=pd.read_csv(O/'results/patient_clinical_features.csv');d=pd.read_csv(O/'data/CAUSHI_CD8_cells_frozen_scores.csv.gz');a=pd.read_csv(O/'results/clinical_associations_prespecified.csv');orig=pd.read_excel(O/'references/41586_2021_3752_MOESM8_ESM.xlsx',sheet_name='Supplementary Table 1',header=1);orig=orig.set_index('Patient study ID')
 ann=pd.read_csv(O/'references/GSE176021_CD8_annotations.csv.gz');expected=set(ann.loc[ann.barcode.str.contains(':tumor-',regex=False)&ann.imid.isin(m.patient),'barcode']);check('Exact original author tumor CD8 membership with known patient',set(d.cell_id)==expected);unlinked=pd.read_csv(O/'data/unlinked_author_CD8_missing_patient.csv');check('2602 missing-patient author cells separately retained',len(unlinked)==2602 and unlinked.patient.isna().all());check('Unique cell identifiers',d.cell_id.is_unique);check('15 pathology patients and 13 primary patients',len(m)==15 and m.primary_analysis.sum()==13);check('Exact primary overlap exclusions',set(m.loc[~m.primary_analysis,'patient'])=={'MD01-004','MD01-005'})
 check('Primary and sensitivity endpoint source values',np.allclose(m.residual_tumor_pct,orig.loc[m.patient,'% residual tumor']));check('Original MPR classification unchanged after whitespace normalization',(m.MPR==orig.loc[m.patient,'MPR status'].str.strip().eq('MPR').to_numpy()).all());check('Pathology range is percent',m.residual_tumor_pct.between(0,100).all());check('Cell counts sum correctly',m.n_cells.sum()==len(d));check('All matched patients satisfy frozen minimum',m.n_cells.ge(spec['minimum_cells']).all())
 audit=pd.read_csv(O/'results/library_matching_audit.csv');check('52 GEX libraries audited',len(audit)==52);check('51 contribute author CD8 cells',audit.author_CD8_cells.gt(0).sum()==51);check('Exactly one zero-author-CD8 library retained in audit',set(audit.loc[audit.author_CD8_cells.eq(0),'library'])=={'GSM5352922'});check('Full gene coverage in every retained library',audit.loc[audit.author_CD8_cells.gt(0),'genes_measured'].eq(66).all())
 marker=set(sum(spec['programmes'].values(),[]));check('38 programme markers exclude all 24 target genes',len(marker)==38 and len(spec['targets'])==24 and not marker&set(spec['targets']))
 means=pd.read_csv(O/'results/marker_normalization.csv').set_index('gene')
 for p,gs in spec['programmes'].items():
  raw=d[['log_'+g for g in gs]].to_numpy();mu=means.loc[gs,'mean_logCP10k'].to_numpy();sd=means.loc[gs,'sd_logCP10k'].to_numpy();v=np.clip((raw-mu)/np.maximum(sd,1e-10),-3,3).mean(1);check(p+' frozen score independent reconstruction',np.allclose(v,d['score_'+p],atol=1e-10))
 for g in ['IL7R','TCF7','CD8A','CXCL13']:
  check(g+' raw-UMI normalization',np.allclose(np.log1p(d['raw_'+g]*10000/d.library_UMI),d['log_'+g],atol=1e-10))
 for _,r in a.iterrows():
  u=m[m.primary_analysis] if r.subset=='primary' else m;col=('score_' if r.kind=='programme' else 'log_')+r.feature;rho=stats.spearmanr(u[col],u.residual_tumor_pct).statistic;check(f'{r.subset} {r.feature} Spearman independently recomputed',np.isclose(r.rho,rho,atol=1e-12));check(f'{r.subset} {r.feature} interval bounds',-1.000001<=r.low<=r.high<=1.000001)
 for subset in ['primary','all_patients']:
  u=a[(a.subset==subset)&~((a.kind=='programme')&(a.feature=='Memory'))];check(subset+' 29-hypothesis secondary family',len(u)==29);check(subset+' BH independently recomputed',np.allclose(multipletests(u.p,method='fdr_bh')[1],u.q_secondary,atol=1e-12))
 u=m[m.primary_analysis];obs=a[(a.subset=='primary')&(a.feature=='Memory')].iloc[0];rx=stats.rankdata(u.score_Memory);ry=stats.rankdata(u.residual_tumor_pct);rng=np.random.default_rng(90261026);null=np.array([np.corrcoef(rx,rng.permutation(ry))[0,1] for _ in range(10000)]);p=(1+(np.abs(null)>=abs(obs.rho)-1e-12).sum())/10001;check('Primary permutation P independently approximated',abs(obs.p-p)<.025)
 for r in json.loads((O/'figure_manifest.json').read_text()):
  for ext in ['png','svg','pdf','tif']:check(r['figure']+' '+ext+' exists',(O/'figures'/(r['figure']+'.'+ext)).stat().st_size>1000)
  im=Image.open(O/'figures'/(r['figure']+'.tif'));check(r['figure']+' TIFF RGB 600 dpi',im.mode=='RGB' and all(abs(v-600)<1 for v in im.info['dpi']));check(r['figure']+' single-page vector PDF',len(PdfReader(O/'figures'/(r['figure']+'.pdf')).pages)==1)
  svg=(O/'figures'/(r['figure']+'.svg')).read_text(encoding='utf-8');check(r['figure']+' live text SVG','<text' in svg)
 check('Clinical atlas contains four parent figures',len(PdfReader(O/'Clinical_Validation_Figure_Atlas.pdf').pages)==4)
 report=dict(passed=all(x['passed'] for x in CHECK),n_checks=len(CHECK),primary_permutation_independent_approximation=p,checks=CHECK)
 (O/'Verification_Clinical.json').write_text(json.dumps(report,indent=2));shutil.copy2(__file__,O/'scripts'/Path(__file__).name);print(json.dumps({k:v for k,v in report.items() if k!='checks'}),flush=True)
if __name__=='__main__':run()
