from pathlib import Path
import sys,json,re,hashlib,os
ROOT=Path(os.environ.get('CD8_STUDY_ROOT',r'D:\Downloads\Virtual_Biotech_Study_20261025\biology_rebuild_20261026'))
OUT=ROOT/'clinical_outcomes_20261027'
sys.path.insert(0,str(OUT/'packages'))
import pandas as pd,numpy as np
import pdfplumber
from lifelines import KaplanMeierFitter
from statsmodels.duration.hazard_regression import PHReg

results=pd.read_csv(OUT/'tables/all_programme_outcome_statistics.csv')
ds={n:pd.read_csv(OUT/'tables'/f'{n}_patient_scores_and_outcomes.csv',index_col=0) for n in ['GSE78220','GSE91061','GSE135222']}
checks=[]
def check(label,value):
    checks.append({'check':label,'passed':bool(value)})
    assert value,label
check('One clinical record per patient',all(d.index.is_unique for d in ds.values()))
check('Hugo baseline patient count 26',len(ds['GSE78220'])==26)
check('Hugo on-treatment Pt16 excluded','Pt16' not in ds['GSE78220'].index)
check('Hugo Pt27 biopsies averaged',ds['GSE78220'].loc['Pt27','biopsies_averaged']==2)
check('Hugo OS missing Pt8 excluded only from OS',pd.isna(ds['GSE78220'].loc['Pt8','time']) and ds['GSE78220'].loc['Pt8','response']==1)
check('Riaz 51 baseline and 49 response-evaluable',len(ds['GSE91061'])==51 and ds['GSE91061'].response.notna().sum()==49)
check('All 24 endpoint tests disclosed',len(results)==24)
check('All eight primary-programme tests disclosed',results.q_primary_8_tests.notna().sum()==8)
cov=pd.read_csv(OUT/'tables/gene_coverage.csv')
check('All original markers measured in each programme/cohort',(cov.measured==cov.defined).all())
pg=pdfplumber.open(OUT/'raw/41467_2019_12159_MOESM1_ESM.pdf').pages[11]
table=[]
for line in pg.extract_text().splitlines():
    vals=line.split()
    if len(vals)==7 and vals[0].isdigit():
        table.append({'patient':'NSCLC'+vals[0],'pfs_source_months':pd.to_numeric(vals[4],errors='coerce'),
                      'event_source':pd.to_numeric(vals[5],errors='coerce'),'clinical_benefit':vals[6]})
tab=pd.DataFrame(table).set_index('patient');check('Source supplement Table 3 has 60 patients',len(tab)==60)
match=ds['GSE135222'][['time','event']].join(tab,how='left')
check('All 27 lung patients match original source table',match.pfs_source_months.notna().all())
check('PFS event coding agrees for all 27 patients',np.allclose(match.event,match.event_source))
check('GEO PFS time equals source months times 30 for all patients',np.allclose(match.time,match.pfs_source_months*30,atol=1e-5))
match.to_csv(OUT/'tables/GSE135222_source_clinical_crosscheck.csv')
for _,r in results[results.endpoint!='Response'].iterrows():
    d=ds[r.cohort].dropna(subset=['time','event']);x=(d[r.programme]-d[r.programme].mean())/d[r.programme].std(ddof=0)
    fit=PHReg(d.time,x.to_numpy()[:,None],status=d.event,ties='efron').fit()
    check(f'Independent Cox implementation {r.cohort} {r.programme}',np.isclose(np.exp(fit.params[0]),r.HR,rtol=2e-4))
follow=[]
for name,d in ds.items():
    if 'time' not in d:continue
    s=d.dropna(subset=['time','event']);rev=KaplanMeierFitter().fit(s.time,event_observed=1-s.event)
    follow.append({'cohort':name,'patients':len(s),'events':int(s.event.sum()),'maximum_recorded_days':s.time.max(),
      'reverse_KM_median_followup_days':rev.median_survival_time_, 'patients_at_risk_1year':int((s.time>=365.25).sum()),
      'patients_at_risk_3years':int((s.time>=1095.75).sum())})
pd.DataFrame(follow).to_csv(OUT/'tables/survival_followup_summary.csv',index=False)
check('No BH-adjusted positive result asserted',not (results.q_all_24_tests<.05).any())
report={'checks':checks,'all_passed':all(c['passed'] for c in checks),'source_verification':{
 'GSE135222_baseline_and_time_origin':'Jung et al. 2019 Methods: tumour obtained before therapy; PFS from treatment start to progression/death.',
 'GSE135222_units':'GEO days match original Table 3 PFS months multiplied by 30 across all 27 patients.',
 'GSE135222_event_coding':'Original Table 3 header PD_Event.1_Censoring.0; all 27 match GEO.',
 'GSE78220_OS':'GEO author-recorded OS in days and Alive/Dead status; precise OS time origin not specified in downloaded GEO fields.'}}
(OUT/'Verification_Report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({'checks':len(checks),'all_passed':report['all_passed'],'followup':follow},indent=2,default=float))
