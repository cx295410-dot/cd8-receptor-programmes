from pathlib import Path
import pandas as pd,json
R=Path(r'D:\Downloads\Virtual_Biotech_Study_20261025\biology_rebuild_20261026\results')
d=pd.read_csv(R/'receptor_programme_rho_summary.csv');selected=[('Memory','IL7R'),('Dysfunction','PDCD1'),('Dysfunction','LAG3')];q=d[[((p,g) in selected) for p,g in zip(d.programme,d.target)]]
p=pd.read_csv(R/'protein_programme_patient_rho.csv');ps=p.groupby(['programme','target']).rho.agg(['median','min','max']).reset_index();ps=ps[[((p,g) in selected) for p,g in zip(ps.programme,ps.target)]]
s=pd.read_csv(R/'decomposition_support.csv');f=pd.read_csv(R/'shared_clone_links.csv');e=pd.read_csv(R/'decomposition_effects.csv').pivot(index=['cohort','target'],columns='term',values='mean')
out={'RNA_examples':q.to_dict('records'),'protein_examples':ps.to_dict('records'),'support':s.groupby('cohort')[['pre_coverage','post_coverage']].agg(['min','median']).round(5).to_string(),'flow_patients':f.groupby('cohort').patient.nunique().to_dict(),'shared_clones':f.groupby('cohort').n_clones.sum().to_dict(),'within_bigger':(e.within_programme.abs()>e.composition.abs()).groupby('cohort').sum().to_dict(),'BCC_top':pd.read_csv(R/'BCC_paired_DESeq2.csv').sort_values('padj').head(3).to_dict('records')}
(R/'key_findings.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out,indent=2))
