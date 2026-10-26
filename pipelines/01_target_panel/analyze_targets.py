"""Fixed-panel transcript availability, patient replication, TCR history, normal-cell context."""
from pathlib import Path
import os,json,itertools,zipfile
ROOT=Path(r'D:\Downloads\Virtual_Biotech_Study_20261025');OLD=Path(r'D:\CD8_Immunotherapy_20261023')
for key in ['TEMP','TMP','MPLCONFIGDIR']:os.environ[key]=str(ROOT/'work')
import numpy as np,pandas as pd
from scipy import stats
from evidence_core import donor_summary,paired_inference,bh_adjust,post_clone_history,pareto_front
TARGETS=json.loads((ROOT/'work/panel.json').read_text())['targets'];SEED=20261025
RESULT=ROOT/'results';DATA=ROOT/'data';rng=np.random.default_rng(SEED)

def save(df,name):df.to_csv(RESULT/(name+'.csv'),index=False)
def bootstrap_mean(values,nboot=5000):
    v=np.asarray(values,float);v=v[np.isfinite(v)]
    if not len(v):return np.nan,np.nan,np.nan
    ci=np.quantile(np.random.default_rng(SEED).choice(v,(nboot,len(v))).mean(1),[.025,.975])
    return float(v.mean()),float(ci[0]),float(ci[1])

def load(cohort):
    m=pd.read_csv(DATA/f'{cohort}_target_cells.csv.gz')
    if cohort in ['BCC','SCC']:
        frozen=pd.read_csv(DATA/f'{cohort}_refined_CD8_identifiers_and_library.csv.gz')
        assert set(frozen.cell_id).issubset(set(m.cell_id))
        m=m[m.cell_id.isin(frozen.cell_id)].copy()
        totals=frozen.set_index('cell_id').library_total;m['library_total']=m.cell_id.map(totals)
        assert m.library_total.notna().all()
    elif cohort=='NSCLC':
        m=m[m.primary_paired].copy();m['library_total']=m.available_matrix_UMI
    elif cohort=='ORAL':m['library_total']=m.matrix_UMI
    else:m['library_total']=np.nan
    assert m.cell_id.is_unique
    m['cohort']=cohort
    if 'clone' in m:m['clone']=m['clone'].replace('',np.nan)
    derived={}
    for g in TARGETS:
        if cohort!='MEL':
            derived['log_'+g]=np.log1p(m['raw_'+g]*10000/np.maximum(m.library_total,1))
            for depth in [500,1000,2000]:
                keep=np.minimum(1,depth/m.library_total.to_numpy(float))
                derived[f'depth{depth}_{g}']=1-np.power(1-keep,m['raw_'+g].to_numpy(float))
    if derived:
        m=m.drop(columns=[c for c in derived if c in m]).copy()
        m=pd.concat([m,pd.DataFrame(derived,index=m.index)],axis=1)
    return m

def primary_analysis(cohort,m):
    levels=[20,10,50];patient_parts=[];stats_parts=[];availability=[];loo=[];summary=[];depthparts=[];normalparts=[]
    for threshold in levels:
        for variant in ['all_donors']+(['exclude_shared_su010'] if cohort=='SCC' else []):
            mm=m[~m.patient.eq('su010')].copy() if variant!='all_donors' else m
            pm=donor_summary(mm,TARGETS,threshold);pm['cohort']=cohort;pm['minimum_cells']=threshold;pm['variant']=variant
            patient_parts.append(pm)
            eligible=pm[['patient','time']];cm=mm.merge(eligible,on=['patient','time'],validate='many_to_one')
            pairs=set(pm.groupby('patient').time.nunique().loc[lambda s:s==2].index)
            paired_cells=cm[cm.patient.isin(pairs)]
            rows=[]
            for g in TARGETS:
                r=paired_inference(pm[['patient','time',g]],g);r.update(cohort=cohort,target=g,minimum_cells=threshold,variant=variant)
                pre=paired_cells.loc[paired_cells.time.eq('pre'),g].dropna();post=paired_cells.loc[paired_cells.time.eq('post'),g].dropna()
                p0=pre.mean();p1=post.mean();n0=len(pre);n1=len(post)
                if min(n0,n1)>0:
                    pooled=(pre.sum()+post.sum())/(n0+n1);se=np.sqrt(pooled*(1-pooled)*(1/n0+1/n1))
                    r['p_cell_independent']=float(2*stats.norm.sf(abs((p1-p0)/se))) if se>0 else 1.
                else:r['p_cell_independent']=np.nan
                r.update(cell_weighted_delta=p1-p0,n_pre_cells=n0,n_post_cells=n1)
                rows.append(r)
                for tm,sub in pm.groupby('time'):
                    mean,lo,hi=bootstrap_mean(sub[g]);availability.append(dict(cohort=cohort,target=g,time=tm,minimum_cells=threshold,
                        variant=variant,n_patients=sub[g].notna().sum(),mean_detection=mean,ci_low=lo,ci_high=hi,
                        donors_above10pct=int((sub[g]>=.1).sum()),donors_above20pct=int((sub[g]>=.2).sum())))
            s=pd.DataFrame(rows);s['q_patient']=bh_adjust(s.p_signflip,family_size=24);s['q_cell_independent']=bh_adjust(s.p_cell_independent,family_size=24);stats_parts.append(s)
            if threshold==20 and variant=='all_donors':
                summary.append(dict(cohort=cohort,n_cd8_cells=len(m),n_primary_cells=len(cm),n_patients=len(pairs),
                    n_pre_cells=int(paired_cells.time.eq('pre').sum()),n_post_cells=int(paired_cells.time.eq('post').sum()),
                    patient_significant=int((s.q_patient<.05).sum()),cell_independent_significant=int((s.q_cell_independent<.05).sum()),
                    cell_significant_without_patient_support=int(((s.q_cell_independent<.05)&(s.q_patient>=.05)).sum()),
                    opposite_delta_direction=int(((s.mean_delta*s.cell_weighted_delta)<0).sum())))
                for donor in sorted(pairs):
                    one=pm[pm.patient.ne(donor)]
                    post=one[one.time.eq('post')][TARGETS].mean();ranking=post.rank(ascending=False,method='min')
                    for g in TARGETS:
                        w=one.pivot(index='patient',columns='time',values=g);d=(w.post-w.pre).dropna()
                        loo.append(dict(cohort=cohort,omitted_patient=donor,target=g,mean_delta=d.mean(),post_rank=ranking[g]))
                if cohort!='MEL':
                    for depth in [500,1000,2000]:
                        cols=[f'depth{depth}_{g}' for g in TARGETS]
                        dm=donor_summary(mm,cols,threshold)
                        rr=[]
                        for g in TARGETS:
                            x=f'depth{depth}_{g}';r=paired_inference(dm[['patient','time',x]],x)
                            r.update(cohort=cohort,target=g,depth=depth);rr.append(r)
                        sd=pd.DataFrame(rr);sd['q_bh']=bh_adjust(sd.p_signflip,family_size=24);depthparts.append(sd)
                    cols=['log_'+g for g in TARGETS];nm=donor_summary(mm,cols,threshold)
                    rr=[]
                    for g in TARGETS:
                        x='log_'+g;r=paired_inference(nm[['patient','time',x]],x);r.update(cohort=cohort,target=g);rr.append(r)
                    nmstat=pd.DataFrame(rr);nmstat['q_bh']=bh_adjust(nmstat.p_signflip,family_size=24);normalparts.append(nmstat)
    return dict(patients=pd.concat(patient_parts),paired=pd.concat(stats_parts),availability=pd.DataFrame(availability),
                loo=pd.DataFrame(loo),summary=pd.DataFrame(summary),depth=pd.concat(depthparts) if depthparts else pd.DataFrame(),
                normalized=pd.concat(normalparts) if normalparts else pd.DataFrame())

def compartment_analysis():
    m=pd.read_csv(DATA/'BCC_target_cells.csv.gz')
    def group(s):
        if s.startswith('CD8_'):return 'CD8 T cells'
        if s.startswith('CD4_'):return 'CD4 T cells'
        if 'Treg' in s:return 'Regulatory T cells'
        if s.startswith('Tumor'):return 'Tumor'
        return s
    m['compartment']=m.author_state.map(group)
    rows=[]
    for (patient,time,comp),sub in m.groupby(['patient','time','compartment']):
        if len(sub)<20:continue
        for g in TARGETS:rows.append(dict(patient=patient,time=time,compartment=comp,target=g,n_cells=len(sub),fraction=sub[g].mean()))
    detail=pd.DataFrame(rows);save(detail,'BCC_compartment_patient_detection')
    donor=detail.groupby(['patient','compartment','target'],as_index=False).fraction.mean()
    save(donor,'BCC_compartment_donor_detection')
    agg=donor.groupby(['compartment','target'],as_index=False).agg(mean_detection=('fraction','mean'),n_donors=('fraction','size'))
    save(agg,'BCC_compartment_detection')

def clone_analysis(cohort,m):
    if 'clone' not in m:return {}
    if cohort=='ORAL':m=m[~m.patient.eq('P23')].copy()
    counts=m.groupby(['patient','time']).size().unstack(fill_value=0)
    good=counts.index[(counts.get('pre',0)>=20)&(counts.get('post',0)>=20)];m=m[m.patient.isin(good)].copy()
    m['history']=post_clone_history(m)
    # If baseline contains no assigned TCR, post clones cannot be newly classified.
    assigned=m[m['clone'].notna()];observed=assigned.groupby(['patient','time']).size().unstack(fill_value=0)
    valid=observed.index[(observed.get('pre',0)>0)&(observed.get('post',0)>0)];m=m[m.patient.isin(valid)].copy()
    post=m[m.time.eq('post')];cov=[];detail=[];clonedetail=[];sampling=[]
    for (patient,time),sub in m.groupby(['patient','time']):
        cov.append(dict(cohort=cohort,patient=patient,time=time,n_cd8=len(sub),n_TCR=int(sub['clone'].notna().sum()),
            fraction_TCR=sub['clone'].notna().mean(),n_clones=sub['clone'].nunique(),
            fraction_new_among_assigned=(sub.history.eq('newly_detected').sum()/sub['clone'].notna().sum()) if time=='post' else np.nan))
    for (patient,history),sub in post.groupby(['patient','history']):
        if history not in ['newly_detected','previously_detected']:continue
        if len(sub)>=10:
            for g in TARGETS:detail.append(dict(cohort=cohort,patient=patient,history=history,target=g,n_cells=len(sub),fraction=sub[g].mean()))
        cl=sub.groupby('clone')[TARGETS].mean()
        if len(cl)>=3 and len(sub)>=10:
            for g in TARGETS:clonedetail.append(dict(cohort=cohort,patient=patient,history=history,target=g,n_clones=len(cl),fraction=cl[g].mean()))
    def contrasts(df):
        rows=[]
        for g in TARGETS:
            sub=df[df.target.eq(g)].rename(columns={'history':'time','fraction':g})
            sub['time']=sub.time.map({'previously_detected':'pre','newly_detected':'post'})
            r=paired_inference(sub[['patient','time',g]],g);r.update(cohort=cohort,target=g);rows.append(r)
        s=pd.DataFrame(rows);s['q_bh']=bh_adjust(s.p_signflip,family_size=24);return s
    detail=pd.DataFrame(detail);clonedetail=pd.DataFrame(clonedetail)
    # Sampling exercise estimates observability under reduced baseline TCR sampling.
    for patient,sub in m.groupby('patient'):
        pre=sub[(sub.time=='pre')&sub['clone'].notna()];po=sub[(sub.time=='post')&sub['clone'].notna()]
        if min(len(pre),len(po))<20:continue
        choices=pre['clone'].to_numpy();postcl=po['clone'].to_numpy()
        for fraction in [.25,.5,1.]:
            reps=1 if fraction==1 else 200
            vals=[]
            for rep in range(reps):
                subset=choices if fraction==1 else rng.choice(choices,max(1,int(len(choices)*fraction)),replace=False)
                vals.append((~np.isin(postcl,np.unique(subset))).mean())
            sampling.append(dict(cohort=cohort,patient=patient,baseline_fraction=fraction,mean_new_fraction=np.mean(vals),
                sampling_low=np.quantile(vals,.025),sampling_high=np.quantile(vals,.975),n_baseline_TCR=len(pre),n_post_TCR=len(po)))
    return {'clone_coverage':pd.DataFrame(cov),'clone_history_patient':detail,'clone_history_contrasts':contrasts(detail),
            'clone_equal_patient':clonedetail,'clone_equal_contrasts':contrasts(clonedetail),'clone_sampling':pd.DataFrame(sampling)}

def atlas_analysis():
    def subset_zip(name):
        with zipfile.ZipFile(DATA/name) as z:
            with z.open(z.namelist()[0]) as f:
                parts=[]
                for chunk in pd.read_csv(f,sep='\t',chunksize=150000):
                    parts.append(chunk[chunk['Gene name'].isin(TARGETS)])
        return pd.concat(parts,ignore_index=True)
    cells=subset_zip('rna_single_cell_type.tsv.zip');tissue=subset_zip('rna_tissue_consensus.tsv.zip')
    mapping=pd.read_csv(DATA/'rna_single_cell_type_cell_types.tsv.zip',sep='\t')
    print('HPA_SCHEMA',cells.columns.tolist(),mapping.columns.tolist(),flush=True)
    # Map through official atlas classes, not a handpicked list of vital organs.
    mapping.columns=[str(c).strip() for c in mapping.columns]
    cellcol=next(c for c in mapping.columns if c.lower()=='cell type')
    classcol=next(c for c in mapping.columns if c.lower()=='cell type class')
    cells=cells.merge(mapping[[cellcol,classcol]],left_on='Cell type',right_on=cellcol,how='left',validate='many_to_one')
    assert cells[classcol].notna().all()
    assert not cells.duplicated(['Gene name','Cell type']).any()
    cells['outside_blood_immune_class']=cells[classcol].str.lower().ne('blood and immune cells')
    save(cells,'HPA_target_normal_cell_expression');save(tissue,'HPA_target_normal_tissue_expression')
    rows=[]
    for gene,sub in cells.groupby('Gene name'):
        non=sub[sub.outside_blood_immune_class];immune=sub[~sub.outside_blood_immune_class]
        maximum=non.loc[non.nCPM.idxmax()]
        rows.append(dict(target=gene,n_normal_types=len(sub),n_outside_blood_immune_types=len(non),
            normal_outside_blood_immune_max_nCPM=non.nCPM.max(),normal_outside_blood_immune_max_cell=maximum['Cell type'],
            normal_outside_blood_immune_breadth_ge1=int((non.nCPM>=1).sum()),
            normal_outside_blood_immune_breadth_ge10=int((non.nCPM>=10).sum()),
            normal_outside_blood_immune_breadth_fraction_ge10=(non.nCPM>=10).mean(),
            normal_immune_max_nCPM=immune.nCPM.max(),n_immune_types=len(immune)))
    result=pd.DataFrame(rows);save(result,'normal_expression_context');return result

def evidence_board(availability,normal):
    a=availability[(availability.minimum_cells==20)&(availability.variant=='all_donors')&availability.time.eq('post')]
    main=a[a.cohort.isin(['BCC','NSCLC','ORAL'])].pivot(index='target',columns='cohort',values='mean_detection').reindex(TARGETS)
    main.columns=['post_'+c for c in main.columns];out=main.reset_index().merge(normal,on='target',how='left',validate='one_to_one')
    assert len(out)==24
    for threshold in [.05,.1,.2]:out[f'cohorts_mean_detection_ge{int(threshold*100)}']=(main.to_numpy()>=threshold).sum(1)
    values=out[['post_BCC','post_NSCLC','post_ORAL']].to_numpy()
    out['availability_pareto']=pareto_front(values)
    out['availability_normal_context_pareto']=pareto_front(np.column_stack([values,-out.normal_outside_blood_immune_breadth_fraction_ge10]))
    save(out,'target_evidence_board');return out

def main():
    outputs={k:[] for k in ['patients','paired','availability','loo','summary','depth','normalized']};clones={};cohorts={}
    for c in ['BCC','SCC','NSCLC','ORAL','MEL']:
        m=load(c);cohorts[c]=m
        result=primary_analysis(c,m)
        for key,val in result.items():
            if len(val):outputs[key].append(val)
        if c!='MEL':
            for key,val in clone_analysis(c,m).items():clones.setdefault(key,[]).append(val)
        print('ANALYZED',c,len(m),flush=True)
    outputs={k:pd.concat(v,ignore_index=True) for k,v in outputs.items()}
    for key,val in outputs.items():save(val,key)
    for key,val in clones.items():save(pd.concat(val,ignore_index=True),key)
    compartment_analysis();normal=atlas_analysis();board=evidence_board(outputs['availability'],normal)
    depthrows=[]
    for c,m in cohorts.items():
        if c=='MEL':continue
        for (patient,time),sub in m.groupby(['patient','time']):
            depthrows.append(dict(cohort=c,patient=patient,time=time,n_cd8=len(sub),median_UMI=sub.library_total.median(),
                mean_UMI=sub.library_total.mean(),fraction_below500=(sub.library_total<500).mean(),
                fraction_below1000=(sub.library_total<1000).mean(),fraction_below2000=(sub.library_total<2000).mean()))
    save(pd.DataFrame(depthrows),'patient_sample_depth')
    # Alternative author annotation and common-sort sensitivity are kept separate.
    broad=pd.read_csv(DATA/'BCC_target_cells.csv.gz');sensitivity=[]
    for label,sub in [('all_metadata_CD8',broad[broad.author_cd8]),
                      ('refined_CD8_common_CD45_CD3_sort',cohorts['BCC'][cohorts['BCC']['sort'].eq('CD45+ CD3+')])]:
        pm=donor_summary(sub,TARGETS,20)
        for g in TARGETS:
            r=paired_inference(pm[['patient','time',g]],g);r.update(annotation=label,target=g)
            r['post_mean_detection']=pm.loc[pm.time.eq('post'),g].mean();sensitivity.append(r)
    sens=pd.DataFrame(sensitivity);sens['q_bh']=sens.groupby('annotation').p_signflip.transform(lambda s:bh_adjust(s,family_size=24));save(sens,'BCC_annotation_sort_sensitivity')
    metric={'run_date':'2026-10-25','primary_cohorts':['BCC','NSCLC','ORAL'],'fixed_targets':TARGETS,
        'summary':outputs['summary'].to_dict('records'),
        'multi_cohort_available_10pct':board.loc[board.cohorts_mean_detection_ge10.eq(3),'target'].tolist(),
        'availability_pareto':board.loc[board.availability_pareto,'target'].tolist(),
        'context_pareto':board.loc[board.availability_normal_context_pareto,'target'].tolist(),
        'interpretation':'Retrospective transcript evidence; no treatment efficacy/causal effect or AI performance benchmark.'}
    (RESULT/'analysis_metrics.json').write_text(json.dumps(metric,indent=2),encoding='utf-8')
    print(json.dumps(metric,indent=2),flush=True)

if __name__=='__main__':main()
