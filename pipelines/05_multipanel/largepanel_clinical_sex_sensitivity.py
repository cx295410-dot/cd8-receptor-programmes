"""POST HOC recorded-Gender adjustment after observed sex-linked gene findings.

This sensitivity never overwrites the primary or original unadjusted estimates.
"""
from pathlib import Path
import sys,os,json,datetime,shutil,importlib.metadata
import largepanel_clinical_genome as base
from largepanel_clinical_genome import ROOT,CLIN,OUT,sha,write_json,bh_all,canonical,SEED
import numpy as np,pandas as pd
from scipy.stats import norm,false_discovery_control

def main():
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from pydeseq2.default_inference import DefaultInference
    import gseapy as gp
    primary_files=[OUT/'results/clinical_continuous_DESeq2.csv',OUT/'results/clinical_continuous_unadjusted_DESeq2.csv',OUT/'results/clinical_Hallmark_GSEA.csv',OUT/'results/clinical_Hallmark_GSEA_unadjusted.csv',OUT/'data/clinical_normalized_pseudobulk.csv.gz',OUT/'data/clinical_pseudobulk_counts.csv.gz',OUT/'data/clinical_pseudobulk_design.csv',OUT/'results/clinical_genome_summary.json',OUT/'results/clinical_genome_source_hashes.json',OUT/'results/clinical_genome_independent_verification.json']
    preserved={str(p):sha(p) for p in primary_files}
    counts=pd.read_csv(OUT/'data/clinical_pseudobulk_counts.csv.gz',index_col=0)
    design=pd.read_csv(OUT/'data/clinical_pseudobulk_design.csv',index_col=0)
    clinical=pd.read_csv(CLIN/'results/patient_clinical_features.csv').set_index('patient')
    meta=design.loc[design.primary_analysis].copy();meta['Gender']=clinical.loc[meta.index,'Gender'].str.strip()
    assert meta.Gender.notna().all() and set(meta.Gender)=={'M','F'}
    primary=pd.read_csv(OUT/'results/clinical_continuous_DESeq2.csv',index_col=0)
    x=counts.loc[meta.index,primary.index].astype(np.int64)
    assert len(meta)==13 and set(design.index[~design.primary_analysis])=={'MD01-004','MD01-005'}
    assert np.array_equal(counts.loc[meta.index].ge(10).sum(axis=0).ge(4),counts.columns.isin(primary.index))
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
    addendum=dict(analysis_type='POST HOC confounding sensitivity',frozen_before_new_fit=True,frozen_at_Asia_Shanghai=now,prompt_date='2026-10-26',trigger='The already completed primary prefix-adjusted genome-wide analysis included sex-linked RPS4Y1 and DDX3Y among31 BH-FDR<0.05 genes. The user-directed coordinator then requested an additional Gender-adjusted sensitivity.',rationale='Assess measured demographic confounding after this observed finding, without altering the primary analysis or optimizing significance.',formula='~ source_prefix + Gender + residual_fraction',Gender_source='Source-recorded Gender column from original clinical Table1, carried unchanged in clinical_validation_20261026/results/patient_clinical_features.csv; no inference from expression',patients='Same13 primary patients, excluding previously functional-context MD01-004 and MD01-005; no Gender-based exclusions',counts='Exactly the primary13 x14,129 raw UMI pseudobulk matrix and fixed >=10 UMI in >=4 primary-patient filter',estimability='Require full design rank and residual degrees of freedom >=5; if not met, audit and stop this sensitivity only',endpoint='Residual viable tumour fraction0 to1; positive coefficient indicates higher residual tumour',BH='Full14,129 tested-gene family, independent filtering disabled, nonfinite p count in family and retain NA',cooks_filter=True,refit_cooks=True,GSEA='Same existing Hallmark GMT, complete finite Wald ranking,1000 gene-set permutations, min15/max500 overlap, seed20261026; GSEA FDR estimates are not BH-adjusted p values',no_patient_exclusion_by_Gender=True,no_selected_gene_tests=True,no_model_tuning=True,preserved_primary_file_sha256=preserved,clinical_source_sha256=sha(CLIN/'results/patient_clinical_features.csv'))
    addendum_path=OUT/'references/clinical_sexadjusted_posthoc_addendum.json'
    if addendum_path.exists():
        previous=json.loads(addendum_path.read_text());assert previous['formula']==addendum['formula'];addendum=previous
    else:write_json(addendum_path,addendum)
    print('FROZEN_POSTHOC_ADDENDUM',addendum['frozen_at_Asia_Shanghai'],flush=True)
    meta.to_csv(OUT/'data/clinical_sexadjusted_pseudobulk_design.csv')
    distribution=meta.groupby(['source_prefix','Gender'],observed=True).residual_tumor_pct.agg(['count','min','median','max']).reset_index()
    distribution.to_csv(OUT/'results/clinical_sexadjusted_source_Gender_distribution.csv',index=False)
    inference=DefaultInference(n_cpus=2)
    dds=DeseqDataSet(counts=x,metadata=meta,design='~ source_prefix + Gender + residual_fraction',refit_cooks=True,inference=inference,quiet=True)
    matrix=dds.obsm['design_matrix'];rank=int(np.linalg.matrix_rank(matrix.to_numpy()));df=len(meta)-rank
    audit=dict(analysis_type='POST HOC',design='~ source_prefix + Gender + residual_fraction',n_patients=len(meta),n_genes_tested=x.shape[1],design_columns=list(matrix.columns),design_rank=rank,design_columns_count=matrix.shape[1],residual_degrees_of_freedom=df,full_rank=rank==matrix.shape[1],minimum_residual_df=5,Gender_counts={str(k):int(v) for k,v in meta.Gender.value_counts().items()},addendum_sha256=sha(addendum_path))
    matrix.to_csv(OUT/'data/clinical_sexadjusted_model_matrix.csv')
    if rank!=matrix.shape[1] or df<5:
        audit['status']='not_estimable_sensitivity_stopped';write_json(OUT/'results/clinical_sexadjusted_model_audit.json',audit);write_json(OUT/'results/clinical_sexadjusted_summary_addendum.json',audit);print('NOT_ESTIMABLE',audit,flush=True);return
    print('FIT_POSTHOC',rank,df,x.shape,meta.Gender.value_counts().to_dict(),flush=True)
    dds.deseq2();contrast=np.zeros(matrix.shape[1]);contrast[list(matrix.columns).index('residual_fraction')]=1.
    ds=DeseqStats(dds,contrast=contrast,inference=inference,quiet=True,independent_filter=False,cooks_filter=True);ds.summary()
    r=ds.results_df.copy();r.index.name='gene';r['padj_package']=r.padj;r['padj']=bh_all(r.pvalue)
    r['log2FoldChange_per_10pct']=r.log2FoldChange/10.;r['lfcSE_per_10pct']=r.lfcSE/10.
    r.to_csv(OUT/'results/clinical_continuous_sexadjusted_DESeq2.csv')
    normalized=pd.DataFrame(dds.layers['normed_counts'],index=x.index,columns=x.columns)
    normalized.to_csv(OUT/'data/clinical_sexadjusted_normalized_pseudobulk.csv.gz')
    factors=pd.Series(dds.obs['size_factors'],index=x.index,name='size_factor');factors.to_csv(OUT/'data/clinical_sexadjusted_size_factors.csv')
    audit.update(status='fitted',BH_family_size=len(r),finite_p_values=int(r.pvalue.notna().sum()),nonfinite_or_Cook_filtered_p_values=int(r.pvalue.isna().sum()),BH_FDR_005=int(r.padj.lt(.05).sum()),positive_BH_FDR_005=int((r.padj.lt(.05)&r.log2FoldChange.gt(0)).sum()),negative_BH_FDR_005=int((r.padj.lt(.05)&r.log2FoldChange.lt(0)).sum()),coefficient_scale='per residual viable tumour fraction0 to1 (100 percentage points)',independent_filter=False,cooks_filter=True,refit_cooks=True)
    write_json(OUT/'results/clinical_sexadjusted_model_audit.json',audit)
    ranking=r.stat.replace([np.inf,-np.inf],np.nan).dropna().sort_values(ascending=False,kind='stable');ranking.name='Wald_stat';ranking.to_csv(OUT/'results/clinical_sexadjusted_Wald_ranking.csv')
    gmt=ROOT/'references/MSigDB_Hallmark_2020.gmt'
    pre=gp.prerank(rnk=ranking,gene_sets=str(gmt),threads=2,permutation_num=1000,min_size=15,max_size=500,seed=SEED,outdir=None,verbose=False)
    enrichment=pre.res2d.copy();enrichment['source_term']=enrichment.Term;enrichment['Term']=enrichment.Term.map(canonical);enrichment.to_csv(OUT/'results/clinical_sexadjusted_Hallmark_GSEA.csv',index=False)
    curves=[];leading=[];membership=[]
    for term,res in pre.results.items():
        term=canonical(term);hits=set(res['hits'])
        for gene in res.get('lead_genes','').split(';'):
            if gene:leading.append(dict(model='posthoc_sexadjusted',term=term,gene=gene,Wald_stat=float(ranking[gene]),log2FoldChange=float(r.loc[gene,'log2FoldChange']),padj=r.loc[gene,'padj']))
        curves.extend(dict(model='posthoc_sexadjusted',term=term,rank=j+1,gene=ranking.index[j],Wald_stat=float(ranking.iloc[j]),running_ES=float(v),hit=j in hits) for j,v in enumerate(res['RES']))
        membership.extend(dict(model='posthoc_sexadjusted',term=term,rank=j+1,gene=ranking.index[j]) for j in hits)
    pd.DataFrame(curves).to_csv(OUT/'results/clinical_sexadjusted_GSEA_curves.csv.gz',index=False)
    pd.DataFrame(leading).to_csv(OUT/'results/clinical_sexadjusted_GSEA_leading_edge.csv',index=False)
    pd.DataFrame(membership).to_csv(OUT/'results/clinical_sexadjusted_GSEA_membership.csv',index=False)
    gene_comparison=primary[['log2FoldChange','stat','pvalue','padj']].add_prefix('primary_').join(r[['log2FoldChange','stat','pvalue','padj']].add_prefix('posthoc_sexadjusted_'))
    gene_comparison['same_effect_direction']=gene_comparison.primary_log2FoldChange*gene_comparison.posthoc_sexadjusted_log2FoldChange>0
    gene_comparison.to_csv(OUT/'results/clinical_sexadjusted_full_gene_stability.csv')
    eprimary=pd.read_csv(OUT/'results/clinical_Hallmark_GSEA.csv').set_index('Term')
    esex=enrichment.set_index('Term')
    stability=eprimary[['NES','NOM p-val','FDR q-val']].add_prefix('primary_').join(esex[['NES','NOM p-val','FDR q-val']].add_prefix('posthoc_sexadjusted_'))
    stability['same_NES_direction']=stability.primary_NES*stability.posthoc_sexadjusted_NES>0
    stability.to_csv(OUT/'results/clinical_sexadjusted_Hallmark_stability.csv')
    source_audit=pd.read_csv(OUT/'results/clinical_library_pseudobulk_audit.csv')
    sources=[CLIN/'results/patient_clinical_features.csv',CLIN/'data/clinical_original_Table1.csv',OUT/'data/clinical_pseudobulk_counts.csv.gz',OUT/'data/clinical_pseudobulk_design.csv',OUT/'results/clinical_continuous_DESeq2.csv',gmt,addendum_path,Path(__file__)]
    write_json(OUT/'results/clinical_sexadjusted_source_hashes.json',[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sources])
    original_table=pd.read_csv(CLIN/'data/clinical_original_Table1.csv').set_index('patient')
    assert clinical.loc[meta.index,'Gender'].str.strip().equals(original_table.loc[meta.index,'Gender'].str.strip())
    assert np.allclose(normalized,x/factors.to_numpy()[:,None])
    assert np.allclose(normalized,pd.read_csv(OUT/'data/clinical_normalized_pseudobulk.csv.gz',index_col=0))
    p=r.pvalue.to_numpy();q=false_discovery_control(np.where(np.isfinite(p),p,1.),method='bh');q[~np.isfinite(p)]=np.nan
    assert np.allclose(r.padj,q,equal_nan=True)
    assert np.allclose(r.log2FoldChange/r.lfcSE,r.stat,rtol=1e-6,atol=1e-9)
    finite=np.isfinite(p);assert np.allclose(p[finite],2*norm.sf(np.abs(r.stat[finite])),rtol=1e-6,atol=1e-12)
    gmt_sets={canonical(v[0]):set(v[2:]) for v in [line.split('\t') for line in gmt.read_text().splitlines()]}
    curve_table=pd.read_csv(OUT/'results/clinical_sexadjusted_GSEA_curves.csv.gz');max_error=0.
    for term,u in curve_table.groupby('term'):
        u=u.sort_values('rank');assert u.gene.tolist()==ranking.index.tolist()
        hits=ranking.index.isin(gmt_sets[term]);weights=np.abs(ranking.to_numpy())*hits
        es=np.cumsum(weights/weights.sum()-(~hits)/(~hits).sum());max_error=max(max_error,float(np.max(np.abs(es-u.running_ES))))
        assert np.isclose(es[np.argmax(np.abs(es))],esex.loc[term,'ES'])
    assert max_error<1e-10
    for p,digest in preserved.items():assert sha(p)==digest,p
    shutil.copy2(__file__,OUT/'scripts'/Path(__file__).name)
    significant_primary=primary.padj.lt(.05);significant_sex=r.padj.lt(.05)
    primary_sig_hallmark=eprimary['FDR q-val'].lt(.05);sex_sig_hallmark=esex['FDR q-val'].lt(.05)
    summary=dict(analysis_type='POST HOC',trigger=addendum['trigger'],addendum_frozen_at=addendum['frozen_at_Asia_Shanghai'],model=audit,primary_results_unchanged=True,primary_gene_FDR005=int(significant_primary.sum()),posthoc_gene_FDR005=int(significant_sex.sum()),primary_significant_genes_also_posthoc_significant=int((significant_primary&significant_sex).sum()),full_gene_direction_agreement_fraction=float(gene_comparison.same_effect_direction.mean()),full_gene_Wald_Spearman=float(gene_comparison[['primary_stat','posthoc_sexadjusted_stat']].corr(method='spearman').iloc[0,1]),GSEA=dict(gmt_sha256=sha(gmt),seed=SEED,permutations=1000,min_size=15,max_size=500,n_pathways=len(enrichment),ranking_genes=len(ranking),primary_FDR005=int(primary_sig_hallmark.sum()),posthoc_FDR005=int(sex_sig_hallmark.sum()),primary_FDR025=int(eprimary['FDR q-val'].lt(.25).sum()),posthoc_FDR025=int(esex['FDR q-val'].lt(.25).sum()),same_NES_direction_pathways=int(stability.same_NES_direction.sum()),NES_Spearman=float(stability[['primary_NES','posthoc_sexadjusted_NES']].corr(method='spearman').iloc[0,1]),primary_significant_pathways_retained=int((primary_sig_hallmark&sex_sig_hallmark).sum()),primary_FDR005_pathway_comparisons=stability.loc[primary_sig_hallmark].reset_index().to_dict('records')),source_quality=dict(Gender_matches_original_clinical_Table1=True,source_recorded_Gender_counts=audit['Gender_counts'],raw_counts_exactly_primary=True,all_frozen_raw_cell_UMI_exact=True,author_nCount_RNA_discrepancy_preserved_in_original_library_audit=True,author_nCount_exact_matches=int(source_audit.author_nCount_exact_matches.fillna(0).sum()),raw_GEX_used_without_annotation_count_rescaling=True),verification=dict(full_family_BH_independently_recomputed=True,Wald_statistics_pvalues_verified=True,normalization_matches_primary=True,all_GSEA_curves_independently_verified=True,max_GSEA_curve_error=max_error,preserved_primary_files=len(preserved),primary_file_hashes_unchanged=True),versions={p:importlib.metadata.version(p) for p in ['pydeseq2','gseapy','numpy','pandas','scipy']})
    write_json(OUT/'results/clinical_sexadjusted_summary_addendum.json',summary)
    print('POSTHOC_COMPLETE',json.dumps(summary),flush=True)

if __name__=='__main__':main()
