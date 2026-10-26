"""Patient-scoped TCR definitions and patient-level inference."""
import itertools
import numpy as np
import pandas as pd
from scipy.stats import rankdata

SIGNATURES={
 'exhaustion':['PDCD1','LAG3','HAVCR2','TIGIT','TOX','ENTPD1','CXCL13','LAYN'],
 'cytotoxic':['NKG7','PRF1','GZMB','GZMH','CTSW','IFNG','CCL5'],
 'memory_progenitor':['TCF7','LEF1','IL7R','CCR7','SELL','SLAMF6'],
 'proliferation':['MKI67','TOP2A','STMN1','TYMS','TUBB'],
 'stress':['HSPA1A','HSPA1B','HSP90AA1','HSPH1','DNAJB1']}

def clone_key(patient,nt,mode='beta'):
    if not isinstance(nt,str) or not nt.strip(): return None
    chains=sorted(set(s.strip() for s in nt.split(';') if s.strip()))
    beta=[s for s in chains if s.startswith('TRB:')]; alpha=[s for s in chains if s.startswith('TRA:')]
    if len(beta)!=1: return None
    if mode=='beta': sequence=beta[0]
    elif mode=='paired':
        if len(alpha)!=1: return None
        sequence=';'.join(sorted(alpha+beta))
    else: raise ValueError(mode)
    return str(patient)+'|'+sequence

def bh_adjust(pvalues):
    p=np.asarray(pvalues,dtype=float); out=np.full(len(p),np.nan)
    valid=np.flatnonzero(np.isfinite(p)); order=valid[np.argsort(p[valid])]
    if len(order):
        q=p[order]*len(order)/np.arange(1,len(order)+1)
        out[order]=np.minimum(1,np.minimum.accumulate(q[::-1])[::-1])
    return out

def paired_summary(df,value,seed=20261023):
    wide=df.groupby(['patient','time'])[value].mean().unstack('time')
    delta=(wide['post']-wide['pre']).dropna().to_numpy(float) if {'pre','post'}.issubset(wide.columns) else np.array([])
    n=len(delta)
    if n==0: return {'n_patients':0,'mean_delta':np.nan,'median_delta':np.nan,'ci_low':np.nan,'ci_high':np.nan,'p_signflip':np.nan}
    rng=np.random.default_rng(seed); observed=abs(delta.mean())
    if n<=16:
        signs=np.array(list(itertools.product([-1,1],repeat=n)),dtype=float)
        pv=float((np.abs((signs*delta).mean(axis=1))>=observed-1e-12).mean())
    else:
        signs=rng.choice([-1,1],(20000,n)); perm=np.abs((signs*delta).mean(axis=1))
        pv=float((1+(perm>=observed-1e-12).sum())/(len(perm)+1))
    boot=rng.choice(delta,(5000,n),replace=True).mean(axis=1)
    ci=np.quantile(boot,[.025,.975]) if n>1 else [np.nan,np.nan]
    return {'n_patients':n,'mean_delta':float(delta.mean()),'median_delta':float(np.median(delta)),
            'ci_low':float(ci[0]),'ci_high':float(ci[1]),'p_signflip':pv}

def validated_join(left,right):
    if left['cell_id'].duplicated().any() or right['cell_id'].duplicated().any(): raise ValueError('Duplicate cell IDs')
    return left.merge(right,on='cell_id',how='left',validate='one_to_one')

def rank_enrichment(values,signature_indices,max_rank=1500):
    values=np.asarray(values); r=rankdata(-values,method='average'); r[values<=0]=max_rank+1
    contribution=(max_rank+1-np.minimum(r,max_rank+1))/max_rank
    return float(contribution[np.asarray(signature_indices,dtype=int)].mean())

def decompose_observed_clones(table):
    n0=table.pre_n.to_numpy(float); n1=table.post_n.to_numpy(float)
    if min(n0.sum(),n1.sum())<=0: raise ValueError('Both observed samples are required')
    p0=n0/n0.sum(); p1=n1/n1.sum(); mu0=table.pre_score.to_numpy(float); mu1=table.post_score.to_numpy(float)
    a=n0>0; b=n1>0; shared=a&b; pre_only=a&~b; post_only=b&~a
    if not np.all(np.isfinite(mu0[a])) or not np.all(np.isfinite(mu1[b])): raise ValueError('Missing observed clone score')
    mean0=np.sum(p0[a]*mu0[a]); mean1=np.sum(p1[b]*mu1[b]); anchor=(mean0+mean1)/2
    state=np.sum((p0[shared]+p1[shared])/2*(mu1[shared]-mu0[shared]))
    reweight=np.sum((p1[shared]-p0[shared])*((mu0[shared]+mu1[shared])/2-anchor))
    detection=np.sum(p1[post_only]*(mu1[post_only]-anchor))-np.sum(p0[pre_only]*(mu0[pre_only]-anchor))
    return {'observed_delta':float(mean1-mean0),'within_shared_state':float(state),'shared_reweighting':float(reweight),
            'detection_composition':float(detection),'reconciliation_error':float(mean1-mean0-state-reweight-detection)}
