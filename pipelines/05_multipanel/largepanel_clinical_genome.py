"""Fixed-protocol exploratory whole-transcriptome clinical pseudobulk analysis.

Independent replicate is the patient. Positive coefficients indicate higher residual
viable tumour fraction. No outcome-selected model, feature set or classifier.
"""
from pathlib import Path
import os, sys, json, tarfile, gzip, hashlib, shutil, re, importlib.metadata
ROOT=Path('D:/Downloads/Virtual_Biotech_Study_20261025/biology_rebuild_20261026')
CLIN=ROOT/'clinical_validation_20261026'; OUT=ROOT/'multipanel_rebuild_20261026'
sys.path[:0]=[str(ROOT/'packages'), 'D:/CD8_Immunotherapy_20261023/analysis_extensions/packages']
os.environ['OMP_NUM_THREADS']='2';os.environ['OPENBLAS_NUM_THREADS']='2'
for key in ['TEMP','TMP','NUMBA_CACHE_DIR']:os.environ[key]=str(OUT/'work')
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np, pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.stats import norm as normal
from threadpoolctl import threadpool_limits
SEED=20261026

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()

def write_json(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8')

def member(tar,name):
    matches=[m for m in tar.getmembers() if m.isfile() and Path(m.name).name==name]
    assert len(matches)==1,(name,len(matches))
    return tar.extractfile(matches[0])

def bh_all(p):
    p=np.asarray(p,float); finite=np.isfinite(p); values=np.where(finite,p,1.)
    order=np.argsort(values); adjusted=np.empty(len(p)); adjusted[order]=np.minimum(1,np.minimum.accumulate((values[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    adjusted[~finite]=np.nan
    return adjusted

def aggregate():
    manifest=json.loads((CLIN/'download_manifest.json').read_text())
    anno=pd.read_csv(CLIN/'references/GSE176021_CD8_annotations.csv.gz',index_col=0).rename(columns={'barcode':'cell_id','imid':'patient','CellType':'author_state'})
    cells=pd.read_csv(CLIN/'data/CAUSHI_CD8_cells_frozen_scores.csv.gz')
    clinical=pd.read_csv(CLIN/'results/patient_clinical_features.csv').set_index('patient')
    assert cells.cell_id.is_unique and len(cells)==151601 and clinical.index.is_unique
    assert set(cells.patient)==set(clinical.index) and clinical.primary_analysis.sum()==13
    unassigned=anno[anno.patient.isna()&anno.cell_id.str.contains(':tumor-',regex=False)]
    assert len(unassigned)==2602
    cache=OUT/'work/clinical_library_pseudobulk';cache.mkdir(parents=True,exist_ok=True)
    audits=[];parts=[];feature_maps=[];seen=[]
    for record in manifest:
        p=Path(record['path']);assert p.stat().st_size==record['bytes']
        prefix=f"{record['patient']}:tumor-{record['lane']}_"
        a=anno[anno.cell_id.str.startswith(prefix)].copy()
        frozen=cells[cells.library.eq(record['gsm'])]
        assert set(a.cell_id)==set(frozen.cell_id),(record['gsm'],len(a),len(frozen))
        base=dict(library=record['gsm'],patient=record['patient'],lane=record['lane'],source=str(p),source_sha256=record['sha256'],author_CD8_cells=len(a))
        dest=cache/(record['gsm']+'.csv.gz');qa=cache/(record['gsm']+'.json')
        if dest.exists() and qa.exists():
            audit=json.loads(qa.read_text());assert audit['source_sha256']==record['sha256']
            if len(a):parts.append(pd.read_csv(dest,index_col=0).iloc[:,0].rename(record['patient']));seen.extend(a.cell_id)
            audits.append(audit);print('CACHED',record['gsm'],len(a),flush=True);continue
        actual_sha=sha(p);assert actual_sha==record['sha256']
        if len(a)==0:
            base.update(excluded_reason='No author CD8 annotation with established patient ID; unassigned NA:tumor-25 cells excluded')
            write_json(qa,base);audits.append(base);print('SKIP_ZERO_AUTHOR',record['gsm'],flush=True);continue
        with tarfile.open(p,mode='r:gz') as tar:
            features=pd.read_csv(gzip.GzipFile(fileobj=member(tar,'features.tsv.gz')),sep='\t',header=None)
            barcodes=pd.read_csv(gzip.GzipFile(fileobj=member(tar,'barcodes.tsv.gz')),sep='\t',header=None)[0].astype(str).tolist()
            with threadpool_limits(limits=2):mat=mmread(gzip.GzipFile(fileobj=member(tar,'matrix.mtx.gz')),spmatrix=True).tocsc()
        assert mat.shape==(len(features),len(barcodes)) and (mat.data>=0).all() and np.equal(mat.data,np.floor(mat.data)).all()
        lookup={prefix+b:i for i,b in enumerate(barcodes)};assert a.cell_id.isin(lookup).all()
        mat=mat[:,[lookup[c] for c in a.cell_id]]
        keep=features[2].eq('Gene Expression').to_numpy() if features.shape[1]>=3 else np.ones(len(features),bool)
        non_gex=int((~keep).sum());mat=mat[keep,:];features=features.loc[keep]
        feature_map=pd.DataFrame({'feature_id':features[0].astype(str),'gene':features[1].astype(str),'library':record['gsm']});feature_maps.append(feature_map)
        symbols=features[1].astype(str).to_numpy();codes,unique=pd.factorize(symbols,sort=False)
        if len(unique)<len(symbols):mat=sparse.csr_matrix((np.ones(len(codes),dtype=np.int64),(codes,np.arange(len(codes)))),shape=(len(unique),len(codes)))@mat
        symbols=np.asarray(unique)
        totals=np.asarray(mat.sum(0)).ravel().astype(np.int64)
        assert np.array_equal(totals,frozen.set_index('cell_id').loc[a.cell_id,'library_UMI'].astype(np.int64))
        author_difference=np.abs(a.nCount_RNA.to_numpy()-totals)/np.maximum(totals,1)
        raw=np.asarray(mat.sum(1)).ravel().astype(np.int64);assert raw.sum()==totals.sum()
        s=pd.Series(raw,index=pd.Index(symbols,name='gene'),name=record['patient'])
        markers=[c[4:] for c in frozen if c.startswith('raw_')]
        marker_errors={g:int(s[g]-frozen['raw_'+g].sum()) for g in markers}
        assert all(v==0 for v in marker_errors.values()),marker_errors
        assert int(totals.sum())==int(frozen.library_UMI.sum())
        s.to_frame('counts').to_csv(dest)
        base.update(raw_matrix_cells=len(barcodes),genes=len(s),non_GEX_features_excluded=non_gex,duplicate_symbols_merged=len(features)-len(unique),CD8_UMI=int(totals.sum()),all_frozen_raw_cell_UMI_exact=True,author_nCount_exact_matches=int((author_difference<1e-10).sum()),author_nCount_max_relative_difference=float(author_difference.max()),marker_count_checks=len(markers),marker_max_absolute_error=max(abs(v) for v in marker_errors.values()),source_sha256_verified=actual_sha)
        write_json(qa,base);parts.append(s);seen.extend(a.cell_id);audits.append(base)
        print('AGGREGATED',record['gsm'],len(a),'cells',int(totals.sum()),'UMI',flush=True)
        del mat
    assert len(seen)==len(set(seen))==len(cells) and set(seen)==set(cells.cell_id)
    counts=pd.concat(parts,axis=1).fillna(0).astype(np.int64).T.groupby(level=0).sum().sort_index();counts.index.name='patient'
    assert int(counts.to_numpy().sum())==int(cells.library_UMI.sum())
    counts.to_csv(OUT/'data/clinical_pseudobulk_counts.csv.gz')
    design=clinical[['primary_analysis','residual_tumor_pct','MPR','n_cells','n_libraries']].copy().loc[counts.index]
    design['residual_fraction']=design.residual_tumor_pct/100.;design['source_prefix']=design.index.str.split('-').str[0]
    assert set(design.index[~design.primary_analysis])=={'MD01-004','MD01-005'}
    design.to_csv(OUT/'data/clinical_pseudobulk_design.csv')
    pd.DataFrame(audits).to_csv(OUT/'results/clinical_library_pseudobulk_audit.csv',index=False)
    if feature_maps:pd.concat(feature_maps).drop_duplicates().to_csv(OUT/'references/clinical_raw_gene_feature_map.csv.gz',index=False)
    return counts,design,audits

def fit_model(counts,meta,label,formula):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from pydeseq2.default_inference import DefaultInference
    inference=DefaultInference(n_cpus=2)
    dds=DeseqDataSet(counts=counts,metadata=meta,design=formula,refit_cooks=True,inference=inference,quiet=True)
    matrix=dds.obsm['design_matrix'];rank=int(np.linalg.matrix_rank(matrix.to_numpy()));df=len(meta)-rank
    assert rank==matrix.shape[1] and df>=5
    matrix.to_csv(OUT/f'data/clinical_{label}_model_matrix.csv')
    contrast=np.zeros(matrix.shape[1]);contrast[list(matrix.columns).index('residual_fraction')]=1.
    print('FIT',label,'patients',len(meta),'genes',counts.shape[1],'design rank',rank,'residual df',df,flush=True)
    dds.deseq2();ds=DeseqStats(dds,contrast=contrast,inference=inference,quiet=True,independent_filter=False,cooks_filter=True);ds.summary()
    r=ds.results_df.copy();r.index.name='gene';r['padj_package']=r.padj
    # Fixed all-tested-gene family: nonfinite tests retain NA and count in m.
    r['padj']=bh_all(r.pvalue);r['log2FoldChange_per_10pct']=r.log2FoldChange/10.;r['lfcSE_per_10pct']=r.lfcSE/10.
    r.to_csv(OUT/('results/clinical_continuous_DESeq2.csv' if label=='adjusted' else 'results/clinical_continuous_unadjusted_DESeq2.csv'))
    normalized=pd.DataFrame(dds.layers['normed_counts'],index=counts.index,columns=counts.columns)
    normalized.to_csv(OUT/('data/clinical_normalized_pseudobulk.csv.gz' if label=='adjusted' else 'data/clinical_unadjusted_normalized_pseudobulk.csv.gz'))
    size_factors=pd.Series(dds.obs['size_factors'],index=counts.index,name='size_factor');size_factors.to_csv(OUT/f'data/clinical_{label}_size_factors.csv')
    assert np.allclose(normalized.to_numpy(),counts.to_numpy()/size_factors.to_numpy()[:,None])
    p=r.pvalue.to_numpy();stat=r.stat.to_numpy();finite=np.isfinite(p)&np.isfinite(stat)
    assert np.allclose(p[finite],2*normal.sf(np.abs(stat[finite])),rtol=1e-6,atol=1e-12)
    audit=dict(label=label,design=formula,n_patients=len(meta),design_columns=list(matrix.columns),design_rank=rank,residual_degrees_of_freedom=df,n_genes_tested=len(r),n_finite_Wald_tests=int(np.isfinite(p).sum()),n_cooks_filtered_or_nonfinite=int((~np.isfinite(p)).sum()),gene_FDR_005=int(r.padj.lt(.05).sum()),positive_FDR_005=int((r.padj.lt(.05)&r.log2FoldChange.gt(0)).sum()),negative_FDR_005=int((r.padj.lt(.05)&r.log2FoldChange.lt(0)).sum()),BH_family_size=len(r),independent_filter=False,cooks_filter=True,refit_cooks=True,coefficient_scale='per 1.0 residual viable tumor fraction (100 percentage points); per_10pct columns divide by 10',size_factors={str(k):float(v) for k,v in size_factors.items()})
    write_json(OUT/f'results/clinical_{label}_model_audit.json',audit)
    print('FIT_DONE',label,audit['gene_FDR_005'],flush=True)
    return r,audit

def canonical(term):
    term=re.sub(r'[^A-Z0-9]+','_',term.upper()).strip('_')
    for a,b in [('TNF_ALPHA','TNFA'),('NF_KB','NFKB'),('IL_2','IL2'),('IL_6','IL6'),('G2_M','G2M')]:term=term.replace(a,b)
    return term if term.startswith('HALLMARK_') else 'HALLMARK_'+term

def gsea(result,label='adjusted'):
    import gseapy as gp
    ranking=result.stat.replace([np.inf,-np.inf],np.nan).dropna().sort_values(ascending=False,kind='stable');ranking.name='Wald_stat'
    ranking.to_csv(OUT/f'results/clinical_{label}_Wald_ranking.csv')
    gmt=ROOT/'references/MSigDB_Hallmark_2020.gmt';shutil.copy2(gmt,OUT/'references/MSigDB_Hallmark_2020.gmt')
    pre=gp.prerank(rnk=ranking,gene_sets=str(gmt),threads=2,permutation_num=1000,min_size=15,max_size=500,seed=SEED,outdir=None,verbose=False)
    enrichment=pre.res2d.copy();enrichment['source_term']=enrichment.Term;enrichment['Term']=enrichment.Term.map(canonical)
    suffix='' if label=='adjusted' else '_unadjusted'
    enrichment.to_csv(OUT/f'results/clinical_Hallmark_GSEA{suffix}.csv',index=False)
    curves=[];leading=[];membership=[]
    for term,res in pre.results.items():
        term=canonical(term);hits=set(res['hits']);lead=res.get('lead_genes','').split(';')
        for gene in lead:
            if gene:leading.append(dict(model=label,term=term,gene=gene,Wald_stat=float(ranking[gene]),log2FoldChange=float(result.loc[gene,'log2FoldChange']),padj=None if pd.isna(result.loc[gene,'padj']) else float(result.loc[gene,'padj'])))
        for j,v in enumerate(res['RES']):curves.append(dict(model=label,term=term,rank=j+1,gene=ranking.index[j],Wald_stat=float(ranking.iloc[j]),running_ES=float(v),hit=j in hits))
        membership.extend(dict(model=label,term=term,rank=j+1,gene=ranking.index[j]) for j in hits)
    pd.DataFrame(curves).to_csv(OUT/f'results/clinical_GSEA_curves{suffix}.csv.gz',index=False)
    pd.DataFrame(leading).to_csv(OUT/f'results/clinical_GSEA_leading_edge{suffix}.csv',index=False)
    pd.DataFrame(membership).to_csv(OUT/f'results/clinical_GSEA_membership{suffix}.csv',index=False)
    return dict(gmt=str(gmt),gmt_sha256=sha(gmt),seed=SEED,permutations=1000,min_size=15,max_size=500,ranking_genes=len(ranking),tested_pathways=len(enrichment),FDR_005=int(pd.to_numeric(enrichment['FDR q-val']).lt(.05).sum()),FDR_025=int(pd.to_numeric(enrichment['FDR q-val']).lt(.25).sum()),rank_statistic='complete finite Wald statistics from fixed model; positive = higher residual fraction')

def main():
    for folder in ['data','results','references','work','scripts']:(OUT/folder).mkdir(parents=True,exist_ok=True)
    protocol=dict(status='fixed before fitting',primary_patients=13,excluded_patients=['MD01-004','MD01-005'],unassigned_author_tumor_CD8=2602,endpoint='residual viable tumour fraction 0 to 1',planned_design='~ source_prefix + residual_fraction',fallback='unadjusted only if adjusted design not full rank or residual df <5',count_filter='raw UMI >=10 in at least max(3,ceil(n_primary/4)) primary patients',BH='all genes passing fixed count filter; nonfinite p count in denominator and retain NA',GSEA_seed=SEED,GSEA_permutations=1000,model_tuning=False)
    write_json(OUT/'references/clinical_genome_protocol.json',protocol)
    counts,design,audits=aggregate();meta=design[design.primary_analysis].copy();x=counts.loc[meta.index]
    min_samples=max(3,int(np.ceil(len(meta)/4)));keep=x.ge(10).sum(0).ge(min_samples)
    pd.DataFrame({'gene':x.columns,'primary_patients_count_ge10':x.ge(10).sum(0).values,'passes_fixed_filter':keep.values,'total_UMI_all15':counts.sum(0).values}).to_csv(OUT/'results/clinical_gene_filter_audit.csv',index=False)
    x=x.loc[:,keep].astype(np.int64)
    matrix=pd.concat([pd.Series(1,index=meta.index,name='Intercept'),pd.get_dummies(meta.source_prefix,drop_first=True,dtype=float),meta.residual_fraction],axis=1)
    rank=int(np.linalg.matrix_rank(matrix.to_numpy()));valid=rank==matrix.shape[1] and len(meta)-rank>=5
    assert valid,'Fixed source-prefix model is not estimable; explicit audit required before unadjusted fallback'
    adjusted,aa=fit_model(x,meta,'adjusted','~ source_prefix + residual_fraction')
    unadjusted,ua=fit_model(x,meta,'unadjusted','~ residual_fraction')
    ga=gsea(adjusted);gu=gsea(unadjusted,'unadjusted')
    distribution=meta.groupby('source_prefix').residual_tumor_pct.agg(['count','min','median','max']).reset_index();distribution.to_csv(OUT/'results/clinical_source_prefix_distribution.csv',index=False)
    sources=[CLIN/'download_manifest.json',CLIN/'references/GSE176021_CD8_annotations.csv.gz',CLIN/'data/CAUSHI_CD8_cells_frozen_scores.csv.gz',CLIN/'results/patient_clinical_features.csv',ROOT/'references/MSigDB_Hallmark_2020.gmt',OUT/'references/clinical_genome_protocol.json',Path(__file__)]
    provenance=[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sources]
    provenance.extend(dict(path=a['source'],sha256=a['source_sha256'],role='raw GEX library',verified=a.get('source_sha256_verified')==a['source_sha256'] or a['author_CD8_cells']==0) for a in audits)
    write_json(OUT/'results/clinical_genome_source_hashes.json',provenance)
    summary=dict(exploratory=True,patient_replicates=True,primary_patients=len(meta),all_patients=len(design),all_CD8_cells=int(design.n_cells.sum()),primary_CD8_cells=int(meta.n_cells.sum()),source_libraries=len(audits),retained_libraries=sum(a['author_CD8_cells']>0 for a in audits),unassigned_CD8_excluded=2602,total_GEX_UMI=int(counts.to_numpy().sum()),full_genes=counts.shape[1],genes_tested=x.shape[1],count_filter_min_UMI=10,count_filter_min_primary_patients=min_samples,adjusted=aa,unadjusted=ua,GSEA=ga,GSEA_unadjusted=gu,source_prefix_distribution=distribution.to_dict('records'),verification={'exact_author_membership':True,'all_frozen_raw_cell_UMI_exact':True,'author_nCount_RNA_discrepancy_retained_in_library_audit':True,'all_frozen_marker_raw_sums_exact':True,'count_conservation':True,'normalization_recomputed':True,'Wald_pvalues_recomputed':True,'BH_recomputed_all_tested_genes':True},versions={p:importlib.metadata.version(p) for p in ['pydeseq2','gseapy','numpy','pandas','scipy']})
    write_json(OUT/'results/clinical_genome_summary.json',summary)
    shutil.copy2(__file__,OUT/'scripts'/Path(__file__).name)
    print('COMPLETE',json.dumps(summary),flush=True)

def verify_outputs():
    """Read delivered files back and independently reconstruct key calculations."""
    from scipy.stats import false_discovery_control
    counts=pd.read_csv(OUT/'data/clinical_pseudobulk_counts.csv.gz',index_col=0)
    design=pd.read_csv(OUT/'data/clinical_pseudobulk_design.csv',index_col=0)
    summary=json.loads((OUT/'results/clinical_genome_summary.json').read_text())
    assert counts.shape==(15,33514) and counts.index.equals(design.index)
    cache=OUT/'work/clinical_library_pseudobulk'
    audits=pd.read_csv(OUT/'results/clinical_library_pseudobulk_audit.csv')
    for patient in counts.index:
        series=[pd.read_csv(cache/(r.library+'.csv.gz'),index_col=0).iloc[:,0] for r in audits.itertuples() if r.patient==patient and r.author_CD8_cells>0]
        independently_summed=pd.concat(series,axis=1).sum(axis=1).reindex(counts.columns).fillna(0)
        assert np.array_equal(counts.loc[patient].to_numpy(),independently_summed.to_numpy())
    assert int(counts.to_numpy().sum())==summary['total_GEX_UMI']==490215866
    primary=counts.loc[design.primary_analysis];passes=primary.ge(10).sum(axis=0).ge(4)
    tested=primary.columns[passes]
    verifications=[];coverage=[]
    sets={}
    for line in (OUT/'references/MSigDB_Hallmark_2020.gmt').read_text().splitlines():
        values=line.split('\t');sets[canonical(values[0])]=set(values[2:])
    for label,suffix in [('adjusted',''),('unadjusted','_unadjusted')]:
        path=OUT/('results/clinical_continuous_DESeq2.csv' if label=='adjusted' else 'results/clinical_continuous_unadjusted_DESeq2.csv')
        result=pd.read_csv(path,index_col=0);assert result.index.equals(pd.Index(tested,name='gene'))
        p=result.pvalue.to_numpy();q=false_discovery_control(np.where(np.isfinite(p),p,1.),method='bh');q[~np.isfinite(p)]=np.nan
        assert np.allclose(result.padj,q,equal_nan=True)
        assert np.allclose(result.log2FoldChange/result.lfcSE,result.stat,rtol=1e-6,atol=1e-9)
        p_finite=np.isfinite(p);assert np.allclose(p[p_finite],2*normal.sf(np.abs(result.stat[p_finite])),rtol=1e-6,atol=1e-12)
        rank=pd.read_csv(OUT/f'results/clinical_{label}_Wald_ranking.csv',index_col=0).iloc[:,0]
        assert len(rank)==len(result) and rank.index.is_unique and rank.is_monotonic_decreasing
        enrichment=pd.read_csv(OUT/f'results/clinical_Hallmark_GSEA{suffix}.csv').set_index('Term')
        curves=pd.read_csv(OUT/f'results/clinical_GSEA_curves{suffix}.csv.gz')
        lead=pd.read_csv(OUT/f'results/clinical_GSEA_leading_edge{suffix}.csv')
        max_error=0.
        for term,genes in sets.items():
            overlap=genes&set(rank.index);eligible=15<=len(overlap)<=500
            coverage.append(dict(model=label,term=term,source_gene_set_size=len(genes),ranking_overlap=len(overlap),eligible=eligible,reason='tested' if eligible else 'overlap below 15' if len(overlap)<15 else 'overlap above 500'))
            assert (term in enrichment.index)==eligible
            if not eligible:continue
            u=curves[curves.term.eq(term)].sort_values('rank')
            assert u.gene.tolist()==rank.index.tolist()
            hits=rank.index.isin(genes);weights=np.abs(rank.to_numpy())*hits
            independent=np.cumsum(weights/weights.sum()-(~hits)/(~hits).sum())
            error=float(np.max(np.abs(independent-u.running_ES.to_numpy())));max_error=max(max_error,error);assert error<1e-10
            es=independent[np.argmax(np.abs(independent))];assert np.isclose(es,enrichment.loc[term,'ES'])
            assert set(lead.loc[lead.term.eq(term),'gene'])==set(enrichment.loc[term,'Lead_genes'].split(';'))
        normalization=pd.read_csv(OUT/('data/clinical_normalized_pseudobulk.csv.gz' if label=='adjusted' else 'data/clinical_unadjusted_normalized_pseudobulk.csv.gz'),index_col=0)
        factors=pd.read_csv(OUT/f'data/clinical_{label}_size_factors.csv',index_col=0).iloc[:,0]
        assert normalization.index.equals(primary.index) and normalization.columns.equals(tested)
        assert np.allclose(normalization,primary.loc[:,tested]/factors.to_numpy()[:,None])
        verifications.append(dict(model=label,genes=len(result),BH_FDR_005=int(result.padj.lt(.05).sum()),all_tested_BH_independently_verified=True,Wald_statistics_and_pvalues_verified=True,normalization_verified=True,pathway_count=len(enrichment),all_running_enrichment_curves_independently_verified=True,max_running_ES_error=max_error,all_leading_edges_verified=True))
    pd.DataFrame(coverage).to_csv(OUT/'results/clinical_Hallmark_pathway_coverage.csv',index=False)
    sources=json.loads((OUT/'results/clinical_genome_source_hashes.json').read_text())
    for entry in sources:
        p=Path(entry['path'])
        if p.resolve()==Path(__file__).resolve():entry['bytes']=p.stat().st_size;entry['sha256']=sha(p)
        assert sha(p)==entry['sha256'],str(p)
    write_json(OUT/'results/clinical_genome_source_hashes.json',sources)
    shutil.copy2(__file__,OUT/'scripts'/Path(__file__).name)
    report=dict(passed=True,source_hashes_verified=len(sources),patient_raw_counts_independently_reconstructed=True,patient_counts=15,full_genes=33514,all_source_UMI=490215866,models=verifications)
    write_json(OUT/'results/clinical_genome_independent_verification.json',report)
    print('VERIFIED',json.dumps(report),flush=True)

if __name__=='__main__':
    if '--verify' not in sys.argv:main()
    verify_outputs()
