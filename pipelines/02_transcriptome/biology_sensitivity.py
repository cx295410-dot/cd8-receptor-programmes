from pathlib import Path
import sys,json
import numpy as np,pandas as pd
from scipy import stats
ROOT=Path(r'D:\Downloads\Virtual_Biotech_Study_20261025\biology_rebuild_20261026');sys.path.insert(0,str(ROOT.parent/'scripts'))
T=json.loads((ROOT.parent/'work/panel.json').read_text())['targets'];P=list(json.loads((ROOT/'work/programme_definitions.json').read_text())['programmes']);CO=['BCC','NSCLC','ORAL'];decs=[];coverage=[];rr=[]
def partial(x,nuisance):
 z=stats.rankdata(x,axis=0);design=np.column_stack([np.ones(len(x)),nuisance]);r=z-design@np.linalg.lstsq(design,z,rcond=None)[0];return np.corrcoef(r,rowvar=False)[:6,6:]
for c in CO:
 m=pd.read_csv(ROOT/'data'/(c+'_cells.csv.gz'))
 for p,u in m.groupby('patient'):
  x=u[['score_'+g for g in P]+['log_'+g for g in T]].to_numpy();nuis=np.column_stack([stats.rankdata(np.log1p(u.library_UMI)),u.time.eq('post')]);rho=partial(x,nuis)
  rr.extend({'cohort':c,'patient':p,'programme':s,'target':g,'rho':rho[i,j]} for i,s in enumerate(P) for j,g in enumerate(T))
 for col,minn in [('programme',3),('programme',5),('programme',10),('leiden',3),('leiden',5),('leiden',10)]:
  for p,u in m.groupby('patient'):
   n=pd.crosstab(u[col],u.time).reindex(columns=['pre','post'],fill_value=0);common=n.index[n.min(axis=1)>=minn];a=u[u[col].isin(common)];cover=a.groupby('time').size()/u.groupby('time').size();coverage.append({'cohort':c,'patient':p,'state_definition':col,'minimum_state_cells':minn,'pre_coverage':cover.get('pre',0),'post_coverage':cover.get('post',0),'included':len(common)>0 and cover.min()>=.5})
   if len(common)==0 or cover.min()<.5:continue
   f=pd.crosstab(a[col],a.time).reindex(columns=['pre','post']).astype(float);f=f/f.sum(axis=0);det=pd.DataFrame(a[['patient','time',col]].copy());det[T]=np.column_stack([(a['raw_'+g]>0).to_numpy() for g in T]);mu=det.groupby([col,'time'])[T].mean();pre=mu.xs('pre',level='time').reindex(common);post=mu.xs('post',level='time').reindex(common);f=f.reindex(common)
   comp=((post+pre).mul((f.post-f.pre)/2,axis=0)).sum();within=((post-pre).mul((f.post+f.pre)/2,axis=0)).sum();total=post.mul(f.post,axis=0).sum()-pre.mul(f.pre,axis=0).sum();assert np.allclose(comp+within,total)
   decs.extend({'cohort':c,'patient':p,'state_definition':col,'minimum_state_cells':minn,'target':g,'composition':comp[g],'within_programme':within[g],'supported_total':total[g]} for g in T)
pd.DataFrame(decs).to_csv(ROOT/'results/decomposition_sensitivity_patient.csv',index=False);pd.DataFrame(coverage).to_csv(ROOT/'results/decomposition_sensitivity_coverage.csv',index=False);pd.DataFrame(rr).to_csv(ROOT/'results/RNA_partial_programme_rho.csv',index=False)
b=pd.read_csv(ROOT/'data/BREAST_cells.csv.gz');prot=[g for g in T if g+'_adt_corrected' in b];r=[]
for p,u in b.groupby('patient'):
 nuisance=pd.get_dummies(u[['tissue','lane']].astype(str),drop_first=True).to_numpy(dtype=float);nuisance=np.column_stack([stats.rankdata(np.log1p(u.library_UMI)),stats.rankdata(np.log1p(u.adt_total)),nuisance]);x=u[['score_'+s for s in P]+[g+'_adt_corrected' for g in prot]].to_numpy();rho=partial(x,nuisance);r.extend({'cohort':'BREAST','patient':p,'programme':s,'target':g,'rho':rho[i,j]} for i,s in enumerate(P) for j,g in enumerate(prot))
pd.DataFrame(r).to_csv(ROOT/'results/protein_partial_programme_rho.csv',index=False)
d=pd.DataFrame(decs).groupby(['state_definition','minimum_state_cells','cohort','target'])[['composition','within_programme']].mean();summary=(d.within_programme.abs()>d.composition.abs()).groupby(['state_definition','minimum_state_cells','cohort']).sum().reset_index(name='within_larger_targets');summary.to_csv(ROOT/'results/decomposition_sensitivity_summary.csv',index=False);print(summary.to_string(index=False))
examples=[('Memory','IL7R'),('Dysfunction','PDCD1'),('Dysfunction','LAG3')]
pr=pd.concat([pd.DataFrame(rr),pd.DataFrame(r)]);out=pr.groupby(['cohort','programme','target']).rho.median().reset_index();out=out[[((s,g) in examples) for s,g in zip(out.programme,out.target)]];out.to_csv(ROOT/'results/selected_partial_associations.csv',index=False);print(out.to_string(index=False))
