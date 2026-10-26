from pathlib import Path
import re,json,shutil
import pandas as pd,numpy as np
R=Path(r'D:\Downloads\Virtual_Biotech_Study_20261025\biology_rebuild_20261026')
def canonical(s):
 s=re.sub(r'[^A-Z0-9]+','_',s.upper()).strip('_')
 for a,b in [('TNF_ALPHA','TNFA'),('NF_KB','NFKB'),('IL_2','IL2'),('IL_6','IL6'),('G2_M','G2M')]:s=s.replace(a,b)
 return s if s.startswith('HALLMARK_') else 'HALLMARK_'+s
gene_sets={}
for line in (R/'references/MSigDB_Hallmark_2020.gmt').read_text().splitlines():
 a=line.split('\t');gene_sets[canonical(a[0])]=set(a[2:])
aud=[]
for c in ['BCC','NSCLC','ORAL']:
 p=R/'results'/(c+'_Hallmark_GSEA.csv');u=pd.read_csv(p)
 if not (R/'results'/(c+'_Hallmark_GSEA_original_names.csv')).exists():shutil.copy2(p,R/'results'/(c+'_Hallmark_GSEA_original_names.csv'))
 if 'source_term' not in u:u['source_term']=u.Term
 u['Term']=u.Term.map(canonical);u.to_csv(p,index=False)
 rank=pd.read_csv(R/'results'/(c+'_Wald_ranking.csv'),index_col=0).iloc[:,0];curve=[]
 for term in ['HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_INTERFERON_ALPHA_RESPONSE','HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT','HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_OXIDATIVE_PHOSPHORYLATION']:
  hits=rank.index.isin(gene_sets[term]);weights=np.abs(rank.to_numpy())*hits;step=np.where(hits,weights/weights.sum(),-1/(len(rank)-hits.sum()));running=np.cumsum(step);es=running[np.abs(running).argmax()];expected=u.set_index('Term').loc[term,'ES'];error=abs(es-expected);assert error<.002,(c,term,error)
  aud.append({'cohort':c,'term':term,'independent_running_ES_error':error});curve.extend({'cohort':c,'term':term,'rank':i,'running_ES':v,'hit':bool(h)} for i,(v,h) in enumerate(zip(running,hits)))
 pd.DataFrame(curve).to_csv(R/'results'/(c+'_GSEA_curves.csv.gz'),index=False)
 print(c,'GSEA q<.05',int((u['FDR q-val']<.05).sum()),u.loc[u.Term.isin(['HALLMARK_INTERFERON_ALPHA_RESPONSE','HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_E2F_TARGETS']),['Term','NES','FDR q-val']].to_dict('records'))
(R/'results/GSEA_curve_verification.json').write_text(json.dumps(aud,indent=2))
