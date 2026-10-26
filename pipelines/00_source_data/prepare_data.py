"""Memory-bounded import and predefined rank scores; no full dense RAM matrix."""
import os
from pathlib import Path
ROOT=Path(r'D:\CD8_Immunotherapy_20261023')
for k in ['TEMP','TMP','MPLCONFIGDIR']: os.environ[k]=str(ROOT/'cache')
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[k]='2'
import gzip,json,time,sys
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from analysis_core import SIGNATURES,clone_key,validated_join
RAW=ROOT/'data'/'raw'; OUT=ROOT/'data'/'processed'; TABLE=ROOT/'results'/'tables'
MARKERS=['CD3D','CD3E','TRAC','CD8A','CD8B','CD4','MS4A1','CD79A','LYZ','LST1','KLRD1','GNLY','PTPRC']

def metadata(cohort):
    if cohort in ['BCC','SCC']:
        tag=cohort.lower(); base=RAW/'GSE123813'
        name=f'GSE123813_{tag}_'+('tcell_metadata' if tag=='bcc' else 'metadata')+'.txt.gz'
        m=pd.read_csv(base/name,sep='\t').rename(columns={'cell.id':'cell_id','treatment':'time','cluster':'author_cluster'})
        m['sample']=m.cell_id.str.split('_').str[0]
        t=pd.read_csv(base/f'GSE123813_{tag}_tcr.txt.gz',sep='\t',index_col=0).rename_axis('cell_id').reset_index()
        m=validated_join(m,t); m['is_cd8']=m.author_cluster.str.startswith('CD8_')
        m['annotation_source']='original author T-cell cluster'
        for mode in ['beta','paired']:
            m['clone_'+mode]=[clone_key(cohort+':'+p,s,mode=mode) for p,s in zip(m.patient,m.cdr3s_nt)]
        path=base/f'GSE123813_{tag}_scRNA_counts.txt.gz'; skip=1; exprmode='UMI counts'
    else:
        source=RAW/'GSE120575'/'GSE120575_patient_ID_single_cells.txt.gz'
        m=pd.read_csv(source,sep='\t',skiprows=19,nrows=16291,encoding='latin1').iloc[:,:7]
        m.columns=['sample_name','cell_id','source','organism','sample','response_recorded','therapy']
        m['patient']=m['sample'].str.extract(r'(P\d+)')[0]; m['time']=m['sample'].str.split('_').str[0].str.lower()
        m['author_cluster']='not_available'; m['annotation_source']='candidate CD8: >=2 T-lineage genes and CD8B detected'
        path=RAW/'GSE120575'/'GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz'; skip=2
        exprmode='author processed TPM-derived values; rank invariant to monotonic transformations'
    m['cohort']=cohort
    if m.cell_id.duplicated().any(): raise ValueError('Nonunique metadata cells')
    return m,path,skip,exprmode

def prepare(cohort):
    print('START',cohort,flush=True); m,path,skip,exprmode=metadata(cohort)
    with gzip.open(path,'rt',encoding='utf-8') as f: rawhead=f.readline().rstrip('\n\r').split('\t')
    header=[s.strip('"') for s in (rawhead[1:] if rawhead[0]=='' else rawhead)]
    if len(set(header))!=len(header): raise ValueError('Repeated matrix cells')
    wanted=m.cell_id[m.is_cd8].tolist() if 'is_cd8' in m else m.cell_id.tolist()
    hmap={s:i for i,s in enumerate(header)}; missing=[s for s in wanted if s not in hmap]
    if missing: raise ValueError(f'{cohort}: {len(missing)} requested cells absent')
    indices=np.array([hmap[s] for s in wanted]); n=len(indices)
    mempath=OUT/f'{cohort}_expression.float32.memmap'; cache=OUT/f'{cohort}_matrix.json'
    if cache.exists() and mempath.exists():
        info=json.loads(cache.read_text()); genes=info['genes']; totals=np.array(info['totals']); detected=np.array(info['detected']); mito=np.array(info['mito'])
        capacity=info.get('capacity',40000)
        mat=np.memmap(mempath,mode='r',dtype='float32',shape=(capacity,n)); g=len(genes)
        print('REUSE_MATRIX',cohort,g,n,flush=True)
    else:
        with gzip.open(path,'rt',encoding='utf-8') as counter: capacity=sum(1 for _ in counter)-skip
        print('GENE_ROWS',cohort,capacity,flush=True)
        mat=np.memmap(mempath,mode='w+',dtype='float32',shape=(capacity,n))
        totals=np.zeros(n); detected=np.zeros(n,dtype=np.int32); mito=np.zeros(n); genes=[]; seen=set(); tick=time.monotonic()
        with gzip.open(path,'rt',encoding='utf-8') as f:
            for _ in range(skip): next(f)
            for line in f:
                gene,sep,rest=line.rstrip('\r\n').partition('\t'); gene=gene.strip('"')
                if not sep: continue
                allvalues=np.fromstring(rest,sep='\t',dtype=np.float32)
                if len(allvalues)!=len(header): raise ValueError(f'{cohort} {gene}: wrong row width {len(allvalues)}')
                if np.any(~np.isfinite(allvalues)) or np.any(allvalues<0): raise ValueError('Invalid expression')
                if gene in seen: raise ValueError('Duplicate gene '+gene)
                if len(genes)>=capacity: raise ValueError('Unexpected matrix row count')
                seen.add(gene); vals=allvalues[indices]; mat[len(genes),:]=vals; genes.append(gene)
                totals+=vals; detected+=vals>0
                if gene.upper().startswith('MT-'): mito+=vals
                if time.monotonic()-tick>25: print('IMPORT',cohort,len(genes),'genes',flush=True); tick=time.monotonic()
        mat.flush(); g=len(genes)
        cache.write_text(json.dumps({'genes':genes,'cells':wanted,'totals':totals.tolist(),'detected':detected.tolist(),'mito':mito.tolist(),'expression_mode':exprmode,'capacity':capacity}),encoding='utf-8')
    m=m.set_index('cell_id').loc[wanted].reset_index(); m['detected_genes']=detected
    if cohort!='MEL':
        m['total_counts']=totals; m['mitochondrial_fraction']=np.divide(mito,totals,out=np.zeros(n),where=totals>0)
    gi={s:i for i,s in enumerate(genes)}
    marker=pd.DataFrame({s:np.asarray(mat[gi[s],:]) if s in gi else np.zeros(n) for s in MARKERS})
    if cohort=='MEL':
        lineage=(marker[['CD3D','CD3E','TRAC']]>0).sum(axis=1)>=2
        m['is_cd8']=(lineage & (marker.CD8B>0)).to_numpy()
        m['is_cd8_broad']=(lineage & ((marker.CD8A>0)|(marker.CD8B>0))).to_numpy()
    use=m.is_cd8.to_numpy(); idsall=np.flatnonzero(use); cover=[]
    for name,sig in SIGNATURES.items():
        present=[s for s in sig if s in gi]
        cover.append({'cohort':cohort,'program':name,'requested':len(sig),'present':len(present),'genes':';'.join(present),'missing':';'.join(s for s in sig if s not in gi)})
        if len(present)<3: raise ValueError('Insufficient signature coverage')
        m[name]=np.nan; m[name+'_logmean']=np.nan
    six={name:np.array([gi[s] for s in sig if s in gi]) for name,sig in SIGNATURES.items()}
    scores={name:np.empty(len(idsall)) for name in SIGNATURES}; means={name:np.empty(len(idsall)) for name in SIGNATURES}
    for start in range(0,len(idsall),96):
        ids=idsall[start:start+96]; values=np.array(mat[:g,ids],dtype=np.float32)
        ranks=rankdata(-values,axis=0,method='average'); ranks[values<=0]=1501
        contribution=(1501-np.minimum(ranks,1501))/1500
        norm=np.log1p(values*10000/np.maximum(totals[ids],1)) if cohort!='MEL' else values
        for name,ix in six.items():
            scores[name][start:start+len(ids)]=contribution[ix].mean(axis=0)
            means[name][start:start+len(ids)]=norm[ix].mean(axis=0)
        if start%1920==0: print('SCORE',cohort,start,'/',len(idsall),flush=True)
    for name in SIGNATURES:
        m.loc[use,name]=scores[name]; m.loc[use,name+'_logmean']=means[name]
    for name in MARKERS: m['expr_'+name]=marker[name].to_numpy()
    m.to_csv(OUT/f'{cohort}_all_cell_metadata.csv.gz',index=False)
    cd8=m.loc[use].copy(); cd8.to_csv(OUT/f'{cohort}_CD8_cell_scores.csv.gz',index=False)
    pd.DataFrame(cover).to_csv(TABLE/f'{cohort}_signature_coverage.csv',index=False)
    count=cd8.groupby(['patient','time']).size().unstack(fill_value=0); count.to_csv(TABLE/f'{cohort}_CD8_patient_time_counts.csv')
    cd8.groupby(['patient','sample','time']).agg(n_cd8=('cell_id','size'),median_genes=('detected_genes','median')).to_csv(TABLE/f'{cohort}_sample_audit.csv')
    audit={'cohort':cohort,'source_matrix_cells':len(header),'metadata_selected_cells':len(m),'CD8_cells':len(cd8),'patients':int(cd8.patient.nunique()),'genes':g,'mode':exprmode,'n_pairs_20':int(((count.get('pre',0)>=20)&(count.get('post',0)>=20)).sum()),'matrix_metadata_missing':len(missing)}
    if 'clone_beta' in cd8:
        audit['CD8_beta_TCR_coverage']=float(cd8.clone_beta.notna().mean()); audit['CD8_paired_TCR_coverage']=float(cd8.clone_paired.notna().mean())
    (TABLE/f'{cohort}_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print('AUDIT',json.dumps(audit),flush=True); del mat

if __name__=='__main__':
    for cohort in (sys.argv[1:] or ['BCC','SCC','MEL']): prepare(cohort)
