"""Extract the fixed surface-target panel from original public expression files.

The previous projects are read-only. New data and hashes are saved on D:.
"""
import os
from pathlib import Path
ROOT=Path(r'D:\Downloads\Virtual_Biotech_Study_20261025')
OLD=Path(r'D:\CD8_Immunotherapy_20261023')
for key in ['TEMP','TMP','MPLCONFIGDIR']: os.environ[key]=str(ROOT/'work')
import gzip,json,hashlib,time,sys
import numpy as np
import pandas as pd
import h5py
from scipy import sparse
from evidence_core import scoped_clone
TARGETS=json.loads((ROOT/'work'/'panel.json').read_text())['targets']
ALIASES={'VSIR':['C10orf54']}
SOURCES=[]

def source(path,role,hash_file=True):
    p=Path(path); rec={'path':str(p),'bytes':p.stat().st_size,'role':role}
    if hash_file:
        h=hashlib.sha256()
        with p.open('rb') as f:
            for chunk in iter(lambda:f.read(2**22),b''):h.update(chunk)
        rec['sha256']=h.hexdigest()
    SOURCES.append(rec)

def emit(frame,name,audit):
    if frame.cell_id.duplicated().any():raise ValueError('Duplicated cells '+name)
    frame.to_csv(ROOT/'data'/f'{name}_target_cells.csv.gz',index=False)
    (ROOT/'data'/f'{name}_extraction.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print('EXTRACTED',name,len(frame),'cells',frame.patient.nunique(),'patients',flush=True)

def original_yost(cohort):
    tag=cohort.lower(); base=OLD/'data/raw/GSE123813'
    meta=base/(f'GSE123813_{tag}_all_metadata.txt.gz' if cohort=='BCC' else f'GSE123813_{tag}_metadata.txt.gz')
    path=base/f'GSE123813_{tag}_scRNA_counts.txt.gz'
    m=pd.read_csv(meta,sep='\t').rename(columns={'cell.id':'cell_id','treatment':'time','cluster':'author_state'})
    tpath=base/f'GSE123813_{tag}_tcr.txt.gz'
    t=pd.read_csv(tpath,sep='\t',index_col=0).rename_axis('cell_id').reset_index()
    assert m.cell_id.is_unique and t.cell_id.is_unique
    m=m.merge(t[['cell_id','cdr3s_nt']],on='cell_id',how='left',validate='one_to_one')
    m['clone']=[scoped_clone(p,s) for p,s in zip(m.patient,m.cdr3s_nt)]
    m['clone_paired']=[scoped_clone(p,s,True) for p,s in zip(m.patient,m.cdr3s_nt)]
    m['author_cd8']=m.author_state.str.startswith('CD8_')
    info_path=OLD/f'data/processed/{cohort}_matrix.json';info=json.loads(info_path.read_text())
    frozen=pd.DataFrame({'cell_id':info['cells'],'library_total':info['totals']})
    frozen.to_csv(ROOT/'data'/f'{cohort}_refined_CD8_identifiers_and_library.csv.gz',index=False)
    source(info_path,'frozen refined author-CD8 membership and library totals')
    values={}; tick=time.monotonic(); readrows=0
    with gzip.open(path,'rt',encoding='utf-8') as f:
        header=[v.strip('"') for v in next(f).rstrip('\r\n').split('\t')]
        if header[0]=='':header=header[1:]
        for line in f:
            gene,sep,rest=line.rstrip('\r\n').partition('\t');gene=gene.strip('"')
            source_gene=gene
            if gene not in TARGETS:
                gene=next((official for official,aliases in ALIASES.items() if gene in aliases),gene)
            if gene not in TARGETS:continue
            v=np.fromstring(rest,sep='\t',dtype=np.float32)
            if len(v)!=len(header):raise ValueError('Wrong raw row width '+gene)
            if not np.isfinite(v).all() or (v<0).any() or np.any(v!=np.floor(v)):raise ValueError('Invalid UMI')
            if gene in values:raise ValueError('Duplicate target symbol '+gene)
            values[gene]=v;readrows+=1
    raw=pd.DataFrame({g:values.get(g,np.full(len(header),np.nan)) for g in TARGETS},index=header)
    assert raw.index.is_unique and set(m.cell_id).issubset(raw.index)
    vv=raw.loc[m.cell_id].reset_index(drop=True)
    for g in TARGETS:m['raw_'+g]=vv[g].to_numpy();m[g]=np.where(vv[g].isna(),np.nan,(vv[g]>0).astype(float))
    m['cohort']=cohort;m['sample']=m.patient.astype(str)+'.'+m.time
    for p in [meta,path,tpath]:source(p,'original author GSE123813 source')
    emit(m,cohort,dict(source_matrix_cells=len(header),selected_cells=len(m),CD8_cells=int(m.author_cd8.sum()),
         patients=m.patient.nunique(),genes_found=sorted(values),genes_missing=sorted(set(TARGETS)-set(values)),
         measurement='raw UMI >0; no batch harmonization',symbol_aliases=ALIASES,seconds=time.monotonic()-tick))

def nsclc():
    base=OLD/'analysis_extensions/independent_GSE179994'
    ipath=base/'matrix/matrix.json';info=json.loads(ipath.read_text())
    mpath=base/'all_CD8_continuous_scores.csv.gz'
    m=pd.read_csv(mpath,usecols=['cell_id','patient','sample','cluster','clone_beta','time','available_matrix_UMI'])
    m=m.rename(columns={'cluster':'author_state','clone_beta':'clone'});m['author_cd8']=True;m['cohort']='NSCLC'
    primary=pd.read_csv(base/'primary_paired_CD8_continuous_scores.csv.gz',usecols=['cell_id'])
    m['primary_paired']=m.cell_id.isin(primary.cell_id)
    genes=info['genes'];ci={c:i for i,c in enumerate(info['cells'])};gi={g:i for i,g in enumerate(genes)}
    ix=np.load(info['sparse_arrays']['i']['path'],mmap_mode='r')
    xp=np.load(info['sparse_arrays']['x']['path'],mmap_mode='r')
    ptr=np.load(info['sparse_arrays']['p']['path'],mmap_mode='r')
    lookup=np.array([next((gi[s] for s in [g]+ALIASES.get(g,[]) if s in gi),-1) for g in TARGETS]);a=np.zeros((len(m),len(TARGETS)),np.float32)
    for k,cell in enumerate(m.cell_id):
        c=ci[cell];start,end=ptr[c:c+2];rows=ix[start:end];vals=xp[start:end]
        pos=np.searchsorted(rows,lookup);ok=(pos<len(rows))&(lookup>=0)
        ok[ok]=rows[pos[ok]]==lookup[ok]
        a[k,ok]=vals[pos[ok]]
    for j,g in enumerate(TARGETS):
        if lookup[j]<0:a[:,j]=np.nan
        m['raw_'+g]=a[:,j];m[g]=np.where(np.isnan(a[:,j]),np.nan,a[:,j]>0)
        m['log_'+g]=np.log1p(a[:,j]*10000/np.maximum(m.available_matrix_UMI,1))
    for g in ['PDCD1','LAG3','IL7R']:
        old=pd.read_csv(mpath,usecols=['cell_id','expr_'+g])
        assert np.array_equal(m['raw_'+g].to_numpy(),old['expr_'+g].to_numpy()),g
    for p in [ipath,mpath,base/'primary_paired_CD8_continuous_scores.csv.gz']:source(p,'author-annotated CD8 with reconstructed original sparse matrix')
    for d in info['sparse_arrays'].values():source(d['path'],'raw matrix extraction array')
    emit(m,'NSCLC',dict(source_matrix_shape=info['shape'],selected_cells=len(m),primary_cells=int(m.primary_paired.sum()),
         primary_patients=int(m.loc[m.primary_paired,'patient'].nunique()),genes_missing=[g for g in TARGETS if g not in gi],
         measurement='raw UMI >0; only first post.1 in paired analyses',raw_source=info['source']))

def oral():
    base=OLD/'analysis_extensions/GSE200996_paired_validation'
    manifest=json.loads((base/'download_manifest.json').read_text())
    mm=base/'tables/paired_tumor_CD8_scores_TCR.csv.gz';m=pd.read_csv(mm)
    m=m.rename(columns={'author_state':'old_author_state','clone_beta':'clone'})
    m['author_state']=m.old_author_state.astype(str);m['author_cd8']=True;m['cohort']='ORAL'
    parts=[];audits=[]
    for rec in [r for r in manifest if r['kind']=='GEX']:
        p=Path(rec['path']);a=m[(m.patient==rec['patient'])&(m.time==rec['time'])].copy()
        with h5py.File(p) as f:
            h=f['matrix'];decode=lambda v:[x.decode() if isinstance(x,bytes) else str(x) for x in v]
            genes=decode(h['features/name'][:]);barcodes=decode(h['barcodes'][:]);ft=decode(h['features/feature_type'][:])
            x=sparse.csc_matrix((h['data'][:],h['indices'][:],h['indptr'][:]),shape=tuple(h['shape'][:]))
        ci={b:i for i,b in enumerate(barcodes)};cols=np.array([ci[b] for b in a.barcode_raw]);gene_mask=np.array(ft)=='Gene Expression'
        totals=np.asarray(x[gene_mask,:][:,cols].sum(axis=0)).ravel()
        assert np.allclose(totals,a.matrix_UMI)
        for g in TARGETS:
            rows=[j for j,name in enumerate(genes) if name in [g]+ALIASES.get(g,[]) and gene_mask[j]]
            values=np.asarray(x[rows,:][:,cols].sum(axis=0)).ravel() if rows else np.full(len(a),np.nan)
            a['raw_'+g]=values;a[g]=np.where(np.isnan(values),np.nan,values>0)
            a['log_'+g]=np.log1p(values*10000/np.maximum(totals,1))
        for g in ['IL7R','SLAMF6']:
            assert np.allclose(a['log_'+g],m.loc[a.index,'log_'+g]),g
        parts.append(a);source(p,'original 10x oral tumor Gene Expression H5')
        audits.append({'patient':rec['patient'],'time':rec['time'],'n_cells':len(a),'file':str(p)})
        print('ORAL_LIBRARY',rec['patient'],rec['time'],len(a),flush=True)
    m=pd.concat(parts,ignore_index=True);source(mm,'original annotated oral CD8 and quality-filtered TCR mapping')
    emit(m,'ORAL',dict(cells=len(m),patients=m.patient.nunique(),libraries=audits,
        TCR_exception='P23 has no baseline TCR; exclude from clone-history comparisons',measurement='raw UMI >0; duplicates summed'))

def melanoma():
    base=OLD/'data/processed';ipath=base/'MEL_matrix.json';info=json.loads(ipath.read_text())
    mpath=base/'MEL_CD8_cell_scores.csv.gz';m=pd.read_csv(mpath)
    m=m[['cell_id','patient','sample','time']].copy();m['author_cd8']=False;m['author_state']='lineage-filtered candidate CD8';m['cohort']='MEL'
    ci={c:i for i,c in enumerate(info['cells'])};gi={g:i for i,g in enumerate(info['genes'])}
    path=base/'MEL_expression.float32.memmap';x=np.memmap(path,dtype='float32',mode='r',shape=(info['capacity'],len(info['cells'])))
    cols=np.array([ci[c] for c in m.cell_id])
    for g in TARGETS:
        symbol=next((s for s in [g]+ALIASES.get(g,[]) if s in gi),None)
        v=np.asarray(x[gi[symbol],cols]) if symbol else np.full(len(m),np.nan)
        m['raw_'+g]=v;m[g]=np.where(np.isnan(v),np.nan,v>0)
    for p in [ipath,mpath]:source(p,'lineage-filtered independent melanoma support')
    emit(m,'MEL',dict(cells=len(m),patients=m.patient.nunique(),measurement='author TPM-derived >0; not UMI; no author CD8 cell annotation'))

if __name__=='__main__':
    for task in (sys.argv[1:] or ['BCC','SCC','NSCLC','ORAL','MEL']):
        if task in ['BCC','SCC']:original_yost(task)
        else:{'NSCLC':nsclc,'ORAL':oral,'MEL':melanoma}[task]()
        (ROOT/'data'/'target_extraction_sources.json').write_text(json.dumps(SOURCES,indent=2),encoding='utf-8')
