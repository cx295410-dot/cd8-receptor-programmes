"""Patient-level evidence operations. No model-derived therapeutic efficacy score."""
import itertools
import numpy as np
import pandas as pd

def donor_summary(frame,targets,min_cells=20):
    grouped=frame.groupby(['patient','time'],sort=True,observed=True)
    means=grouped[targets].mean()
    means['n_cells']=grouped.size()
    return means.loc[means.n_cells>=min_cells].reset_index()

def bh_adjust(values,family_size=None):
    p=np.asarray(values,float); q=np.full(len(p),np.nan)
    valid=np.flatnonzero(np.isfinite(p)); order=valid[np.argsort(p[valid])]
    if len(order):
        family=len(order) if family_size is None else family_size
        if family<len(order):raise ValueError('Family smaller than observed p-values')
        raw=p[order]*family/np.arange(1,len(order)+1)
        q[order]=np.minimum(1,np.minimum.accumulate(raw[::-1])[::-1])
    return q

def paired_inference(frame,target,seed=20261025,n_boot=5000):
    if frame.duplicated(['patient','time']).any():
        raise ValueError('Input must contain one record per patient and time')
    wide=frame.pivot(index='patient',columns='time',values=target)
    if {'pre','post'}.issubset(wide.columns):
        delta=(wide.post-wide.pre).dropna().to_numpy(float)
    else: delta=np.array([],float)
    n=len(delta)
    out={'n_pairs':n,'mean_delta':np.nan,'ci_low':np.nan,'ci_high':np.nan,'p_signflip':np.nan,
         'positive_pairs':0,'negative_pairs':0,'zero_pairs':0}
    if not n: return out
    rng=np.random.default_rng(seed); obs=abs(delta.mean())
    if n<=18:
        signs=np.array(list(itertools.product([-1.,1.],repeat=n)))
        p=float((np.abs(signs@delta/n)>=obs-1e-12).mean())
    else:
        signs=rng.choice([-1.,1.],(20000,n))
        p=float((1+(np.abs(signs@delta/n)>=obs-1e-12).sum())/20001)
    boot=rng.choice(delta,(n_boot,n),replace=True).mean(1)
    ci=np.quantile(boot,[.025,.975]) if n>1 else [np.nan,np.nan]
    out.update(mean_delta=float(delta.mean()),ci_low=float(ci[0]),ci_high=float(ci[1]),p_signflip=p,
               positive_pairs=int((delta>0).sum()),negative_pairs=int((delta<0).sum()),zero_pairs=int((delta==0).sum()))
    return out

def scoped_clone(patient,sequence,paired=False):
    if not isinstance(sequence,str) or not sequence.strip(): return None
    chains=sorted(set(s.strip() for s in sequence.split(';') if s.strip()))
    beta=[s for s in chains if s.startswith('TRB:')]; alpha=[s for s in chains if s.startswith('TRA:')]
    if len(beta)!=1 or (paired and len(alpha)!=1): return None
    return str(patient)+'|'+(';'.join(sorted(alpha+beta)) if paired else beta[0])

def post_clone_history(frame,clone_column='clone'):
    baseline=frame.loc[frame.time.eq('pre')]
    pre_clones=set(
        baseline.loc[baseline[clone_column].notna(),['patient',clone_column]].itertuples(index=False,name=None))
    paired=set(baseline.patient)
    def label(row):
        if row.time=='pre': return 'pre'
        if row.patient not in paired: return 'no_paired_baseline'
        clone=getattr(row,clone_column)
        if pd.isna(clone) or not str(clone).strip(): return 'unassigned'
        return 'previously_detected' if (row.patient,clone) in pre_clones else 'newly_detected'
    return pd.Series([label(r) for r in frame.itertuples()],index=frame.index)

def pareto_front(values):
    """All columns maximize. Rows with missing evidence do not enter this front."""
    x=np.asarray(values,float); good=np.isfinite(x).all(1); front=np.zeros(len(x),bool)
    for i in np.flatnonzero(good):
        front[i]=not np.any(np.all(x[good]>=x[i],axis=1)&np.any(x[good]>x[i],axis=1))
    return front
