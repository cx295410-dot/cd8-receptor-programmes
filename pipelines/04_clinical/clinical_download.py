"""Freeze an analysis plan and complete all public primary-resection tumor GEX libraries."""
from pathlib import Path
import sys,json,re,hashlib,urllib.request,concurrent.futures,gzip,time,shutil
sys.stdout.reconfigure(encoding='utf-8')
import pandas as pd
ROOT=Path('D:/Downloads/Virtual_Biotech_Study_20261025/biology_rebuild_20261026')
O=ROOT/'clinical_validation_20261026'; B=Path('D:/CD8_Immunotherapy_20261023/analysis_extensions')
OLD=B/'GSE176021_antigen_validation'
for n in ['raw','references','data','results','scripts','figures','work','manuscript','delivery']:(O/n).mkdir(exist_ok=True,parents=True)
definitions=json.loads((ROOT/'work/programme_definitions.json').read_text())
targets=json.loads((ROOT.parent/'work/panel.json').read_text())['targets']
spec=dict(version='1.0',frozen_UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),cohort='GSE176021, Caushi et al. Nature 2021',
 eligibility='All author-annotated CD8 cells with raw GEX from primary resected tumor; exact patient IDs in original Supplementary Table 1; no selection on functional assay availability or observed outcomes. Exclude normal lung, lymph node, postoperative metastatic recurrence and stimulated cultures.',
 endpoint='Original Supplementary Table 1 % residual tumor (pathological residual viable tumor percentage, 0–100); cross-sectional post-treatment association, not tumor diameter and not baseline response prediction.',
 primary='Two-sided patient-level Spearman association of mean frozen Memory programme score with residual tumor percentage; primary set excludes MD01-004 and MD01-005 already shown in functional-clone analysis. No outcome-based exclusions or direction changes.',
 primary_exclusions=['MD01-004','MD01-005'],minimum_cells=30,programmes=definitions['programmes'],targets=targets,
 normalization='Gene Expression raw UMI / all Gene Expression UMI * 10000; log1p. For each marker, cohort-wide cell z score clipped to [-3,3]; frozen programme is mean of its markers. All target genes excluded from programme definitions.',
 secondary='Remaining five programmes and all 24 fixed target genes: one BH family of 29 hypotheses in primary patient set. Mean log1pCP10K for target genes. MPR <=10% descriptive only.',
 inference='Patient, not cell, is the independent replicate. Two-sided Spearman permutation p with 99,999 permutations (plus observed) and 10,000 patient bootstrap resamples for percentile 95% CI. Fixed seed 20261026. No fitted classifier or expression cutpoint.',
 sensitivities=['All pathology-matched primary-resection patients, including the two previous functional donors','Leave-one-patient-out Spearman correlations for primary Memory endpoint','Equal-library rather than equal-cell averaging within patient','Patient rank-based Memory scoring without cohort-wide gene scaling'],
 programme_association_replication='Within each patient, Spearman correlation of IL7R log expression and Memory score. Display all patient correlations and patient bootstrap mean CI; no cell-level significance tests.',
 prior_exposure='This source was previously used for two-donor functional-clone context and some downloaded libraries. This is a new patient-level pathology reanalysis, not a blinded untouched validation cohort.',
 negative_results='Retain and report all prespecified results, regardless of significance. Do not turn exploratory gene associations into response-prediction claims.')
sp=O/'analysis_protocol_frozen.json'
if not sp.exists():sp.write_text(json.dumps(spec,indent=2),encoding='utf-8')
else:assert json.loads(sp.read_text())['version']=='1.0'
x=OLD/'supplements/41586_2021_3752_MOESM8_ESM.xlsx';s=pd.read_excel(x,sheet_name='Supplementary Table 1',header=1)
s=s[s['Patient study ID'].astype(str).str.fullmatch(r'(?:MD01|MD043|NY016)-\d{3}')].copy()
s=s.rename(columns={'Patient study ID':'patient','% residual tumor':'residual_tumor_pct','MPR status':'author_MPR'})
assert s.patient.is_unique and s.residual_tumor_pct.between(0,100).all()
s['primary_analysis']=~s.patient.isin(spec['primary_exclusions']);s.to_csv(O/'data/clinical_original_Table1.csv',index=False)
z=pd.read_excel(x,sheet_name='Supplementary Table 2',header=3);z['Patient study ID']=z['Patient study ID'].ffill()
z=z[z['Assays performed'].astype(str).str.contains('Single cell',case=False)].copy()
assert z.loc[z['Patient study ID'].isin(s.patient),'Resected tumor (fresh/viable)'].notna().sum()==15
z.to_csv(O/'data/clinical_specimen_eligibility_Table2.csv',index=False)
existing=json.loads((OLD/'download_manifest.json').read_text());existing={r['name']:r for r in existing if r['kind']=='GEX'}
items=[]
for line in (B/'oncology_letters_extension/new_data/GSE176021_filelist.txt').read_text().splitlines():
 p=line.split('\t')
 if len(p)<5 or p[0]!='File':continue
 m=re.match(r'(GSM\d+)_((?:MD01|MD043|NY016)-\d{3})_tumor_(\d+)\.tar\.gz$',p[1])
 if not m:continue
 gsm,patient,lane=m.groups();assert patient in set(s.patient)
 items.append(dict(name=p[1],bytes=int(p[3]),gsm=gsm,patient=patient,lane=int(lane),kind='GEX',url=f'https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/{p[1]}'))
plan=dict(n_libraries=len(items),n_patients=len(set(i['patient'] for i in items)),new_files=sum(i['name'] not in existing for i in items),new_bytes=sum(i['bytes'] for i in items if i['name'] not in existing),items=items)
(O/'download_plan.json').write_text(json.dumps(plan,indent=2));print('PLAN',{k:v for k,v in plan.items() if k!='items'},flush=True)
assert shutil.disk_usage(O).free>plan['new_bytes']*2+2**30
def fetch(item):
 if item['name'] in existing:
  r=existing[item['name']];assert Path(r['path']).stat().st_size==item['bytes'];return dict(r,reused=True)
 dest=O/'raw'/item['name'];record=dest.with_suffix('.verified.json')
 if dest.exists() and record.exists():
  r=json.loads(record.read_text());assert dest.stat().st_size==item['bytes'];return r
 for retry in range(3):
  try:
   part=dest.with_suffix(dest.suffix+'.part')
   with urllib.request.urlopen(urllib.request.Request(item['url'],headers={'User-Agent':'Mozilla/5.0'}),timeout=90) as f,part.open('wb') as w:
    while b:=f.read(2**20):w.write(b)
   assert part.stat().st_size==item['bytes']
   with gzip.open(part,'rb') as f:
    while f.read(4*2**20):pass
   sha=hashlib.sha256(part.read_bytes()).hexdigest();part.replace(dest)
   r=dict(item,path=str(dest),sha256=sha,gzip_CRC_verified=True,reused=False);record.write_text(json.dumps(r,indent=2));print('DOWNLOADED',dest.name,flush=True);return r
  except Exception as e:
   print('RETRY',item['name'],retry,str(e),flush=True)
 raise RuntimeError(item['name'])
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as p:records=list(p.map(fetch,items))
(O/'download_manifest.json').write_text(json.dumps(records,indent=2))
for source in [x,ROOT/'work/programme_definitions.json',ROOT.parent/'work/panel.json',B/'oncology_letters_extension/new_data/GSE176021_CD8_annotations.csv.gz',B/'oncology_letters_extension/new_data/GSE176021_CD8_annotations.rds.gz',B/'oncology_letters_extension/new_data/GSE176021_filelist.txt']:
 shutil.copy2(source,O/'references'/source.name)
shutil.copy2(__file__,O/'scripts'/Path(__file__).name)
print('COMPLETE',len(records),'libraries',flush=True)
