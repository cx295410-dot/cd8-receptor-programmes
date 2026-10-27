"""Patient-level post hoc outcome extension. All data/output paths are on D:."""
from pathlib import Path
import sys, csv, gzip, json, re, hashlib, shutil, platform, os
ROOT = Path(os.environ.get('CD8_STUDY_ROOT',r'D:\Downloads\Virtual_Biotech_Study_20261025\biology_rebuild_20261026'))
OUT = ROOT / 'clinical_outcomes_20261027'
sys.path.insert(0, str(OUT / 'packages'))
import numpy as np
import pandas as pd
import scipy
from scipy.stats import mannwhitneyu, spearmanr
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
from sklearn.metrics import roc_auc_score, roc_curve
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.patches import Patch
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import proportional_hazard_test, logrank_test

RAW = OUT / 'raw'
for sub in ['tables','figures','source_data','scripts']:
    (OUT/sub).mkdir(parents=True, exist_ok=True)
PROGRAMMES=json.loads((ROOT/'work/programme_definitions.json').read_text())['programmes']
PRIMARY=['Memory','Dysfunction']
COLOURS={'Memory':'#5587B8','Dysfunction':'#D75564','Cytotoxic':'#40A68F',
         'Cycling':'#A87BC0','Interferon':'#DDA354','Residency':'#7C9C63'}
RESPONSE={0:'#63AAA3',1:'#CF5262'}
CMAP=LinearSegmentedColormap.from_list('expression',['#5478AE','#F9F8F5','#D84C58'])
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,
 'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'axes.spines.top':False,
 'axes.spines.right':False,'axes.linewidth':.6,'pdf.fonttype':42,'ps.fonttype':42,
 'svg.fonttype':'none','savefig.facecolor':'white'})
RNG=np.random.default_rng(20261027)

def metadata(name):
    rows=list(csv.reader(gzip.open(RAW/f'{name}_series_matrix.txt.gz','rt'),delimiter='\t'))
    ids=next(r[1:] for r in rows if r and r[0]=='!Sample_geo_accession')
    titles=next(r[1:] for r in rows if r and r[0]=='!Sample_title')
    ds=[{'gsm':i,'title':t} for i,t in zip(ids,titles)]
    for row in rows:
        if row and row[0]=='!Sample_characteristics_ch1':
            for d,v in zip(ds,row[1:]):
                if ': ' in v:
                    k,v=v.split(': ',1); d[k]=v
    df=pd.DataFrame(ds)
    df.to_csv(OUT/'source_data'/f'{name}_metadata_all_samples.csv',index=False)
    return df

def score(expr,meta,name):
    assert meta.index.is_unique and set(meta.index)==set(expr.columns)
    expr=expr.loc[:,meta.index].astype(float)
    assert np.isfinite(expr.to_numpy()).all() and (expr.to_numpy()>=0).all()
    logged=np.log1p(expr)
    pts=list(dict.fromkeys(meta.patient))
    out=[]; cols=[]
    for pt in pts:
        mm=meta[meta.patient==pt]
        for field in ['response_label','time','event']:
            if field in mm:
                assert mm[field].fillna('NA').nunique()==1,(name,pt,field)
        rr=mm.iloc[0].to_dict(); rr['biopsies_averaged']=len(mm)
        rr['sample_ids']=';'.join(mm.index.astype(str))
        rr['patient']=pt; out.append(rr)
        cols.append(logged.loc[:,mm.index].mean(axis=1).rename(pt))
    patient=pd.DataFrame(out).set_index('patient')
    lg=pd.concat(cols,axis=1)
    genez=lg.sub(lg.mean(axis=1),axis=0).div(lg.std(axis=1,ddof=0).replace(0,np.nan),axis=0).fillna(0).clip(-3,3)
    coverage=[]
    for prog,genes in PROGRAMMES.items():
        measured=[g for g in genes if g in genez.index]
        assert len(measured)>=np.ceil(len(genes)/2),(name,prog,measured)
        patient[prog]=genez.loc[measured].mean(axis=0)
        coverage.append({'cohort':name,'programme':prog,'measured':len(measured),
                         'defined':len(genes),'genes':';'.join(measured)})
    patient.to_csv(OUT/'tables'/f'{name}_patient_scores_and_outcomes.csv')
    genez.to_csv(OUT/'source_data'/f'{name}_marker_zscores.csv')
    return patient,genez,coverage

def prepare():
    datasets={}; genes={}; coverage=[]
    m=metadata('GSE78220')
    m['patient']=m['patient id']; m['response_label']=m['anti-pd-1 response']
    m['time']=pd.to_numeric(m['overall survival (days)'],errors='coerce')
    m['event']=m['vital status'].map({'Dead':1,'Alive':0})
    m['response']=m.response_label.map({'Complete Response':1,'Partial Response':1,'Stable Disease':0,'Progressive Disease':0})
    m=m[m['biopsy time']=='pre-treatment'].copy()
    m['expression_id']=m.title+'.baseline'; m=m.set_index('expression_id')
    e=pd.read_excel(RAW/'GSE78220_PatientFPKM.xlsx',sheet_name='FPKM',index_col=0)
    e=e.groupby(level=0,sort=False).mean()
    datasets['GSE78220'],genes['GSE78220'],c=score(e.loc[:,m.index],m,'GSE78220'); coverage+=c

    m=metadata('GSE91061'); m['patient']=m.title.str.extract(r'^(Pt\d+)',expand=False)
    m['response_label']=m.response
    m['response']=m.response_label.map({'PRCR':1,'SD':0,'PD':0})
    m=m[m['visit (pre or on treatment)']=='Pre'].set_index('title')
    e=pd.read_csv(RAW/'GSE91061_fpkm.csv.gz',index_col=0)
    info=pd.read_csv(RAW/'Homo_sapiens.gene_info.gz',sep='\t',dtype=str,low_memory=False)
    lookup=info.drop_duplicates('GeneID').set_index('GeneID').Symbol
    missing=set(e.index.astype(str))-set(lookup.index)
    e.index=e.index.astype(str).map(lookup)
    e=e[e.index.notna()].groupby(level=0,sort=False).mean()
    datasets['GSE91061'],genes['GSE91061'],c=score(e.loc[:,m.index],m,'GSE91061'); coverage+=c

    m=metadata('GSE135222'); m['patient']=m.title.str.replace(' ','',regex=False)
    m['time']=pd.to_numeric(m['pfs.time'],errors='coerce')
    m['event']=pd.to_numeric(m['progression-free survival (pfs)'],errors='coerce')
    m=m.set_index('patient',drop=False)
    e=pd.read_csv(RAW/'GSE135222_expression.tsv.gz',sep='\t',index_col=0)
    gm=pd.read_csv(ROOT/'multipanel_rebuild_20261026/references/clinical_raw_gene_feature_map.csv.gz')
    gm=gm[['feature_id','gene']].drop_duplicates()
    ambiguous=gm.groupby('feature_id').gene.nunique(); assert not (ambiguous>1).any()
    gm=gm.drop_duplicates('feature_id').set_index('feature_id').gene
    gm.to_csv(OUT/'source_data'/'ensembl_symbol_map.csv')
    e.index=e.index.str.replace(r'\.\d+$','',regex=True).map(gm)
    e=e[e.index.notna()].groupby(level=0,sort=False).mean()
    datasets['GSE135222'],genes['GSE135222'],c=score(e.loc[:,m.index],m,'GSE135222'); coverage+=c
    pd.DataFrame(coverage).to_csv(OUT/'tables'/'gene_coverage.csv',index=False)
    return datasets,genes

def auc_ci(y,x):
    y=np.asarray(y,dtype=int); x=np.asarray(x); ix=[np.flatnonzero(y==v) for v in [0,1]]
    vals=[]
    for _ in range(3000):
        b=np.concatenate([RNG.choice(v,len(v),replace=True) for v in ix])
        vals.append(roc_auc_score(y[b],x[b]))
    return np.percentile(vals,[2.5,97.5])

def fit_cox(d,programme):
    z=(d[programme]-d[programme].mean())/d[programme].std(ddof=0)
    a=pd.DataFrame({'time':d.time,'event':d.event,'score':z})
    model=CoxPHFitter().fit(a,'time','event')
    r=model.summary.loc['score']
    ph=proportional_hazard_test(model,a,time_transform='rank').summary.loc['score','p']
    return model,{'HR':r['exp(coef)'],'CI_low':r['exp(coef) lower 95%'],
      'CI_high':r['exp(coef) upper 95%'],'p':r['p'],'PH_rank_p':ph,
      'C_index_apparent':model.concordance_index_}

def analyse(datasets):
    rows=[]; counts=[]
    for name,d in datasets.items():
        counts.append({'cohort':name,'baseline_patients':len(d),'endpoint':('PFS' if name=='GSE135222' else 'ORR/OS' if name=='GSE78220' else 'ORR'),
         'response_evaluable':int(d.response.notna().sum()) if 'response' in d else 0,
         'responders':int(d.response.sum()) if 'response' in d else 0,
         'survival_evaluable':int(d[['time','event']].notna().all(axis=1).sum()) if 'time' in d else 0,
         'events':int(d.event.sum()) if 'event' in d else 0})
        if 'response' in d:
            r=d[d.response.notna()]
            for prog in PROGRAMMES:
                _,p=mannwhitneyu(r.loc[r.response==1,prog],r.loc[r.response==0,prog],alternative='two-sided',method='asymptotic')
                low,high=auc_ci(r.response,r[prog])
                rows.append({'cohort':name,'endpoint':'Response','programme':prog,'n':len(r),
                 'responders':int(r.response.sum()),'AUC':roc_auc_score(r.response,r[prog]),'AUC_low':low,'AUC_high':high,'p':p})
        if 'time' in d:
            s=d.dropna(subset=['time','event']); assert (s.time>0).all() and set(s.event.unique())<={0,1}
            for prog in PROGRAMMES:
                model,r=fit_cox(s,prog)
                rows.append({'cohort':name,'endpoint':'OS' if name=='GSE78220' else 'PFS','programme':prog,
                             'n':len(s),'events':int(s.event.sum()),**r})
    res=pd.DataFrame(rows)
    res['q_all_24_tests']=multipletests(res.p,method='fdr_bh')[1]
    prim=res.programme.isin(PRIMARY)
    res.loc[prim,'q_primary_8_tests']=multipletests(res.loc[prim,'p'],method='fdr_bh')[1]
    res.to_csv(OUT/'tables'/'all_programme_outcome_statistics.csv',index=False)
    pd.DataFrame(counts).to_csv(OUT/'tables'/'patient_flow.csv',index=False)
    return res,pd.DataFrame(counts)

def panel(ax,letter,title):
    ax.set_title(title,loc='left',pad=9,fontweight='normal')
    ax.text(-.15,1.13,letter,transform=ax.transAxes,fontsize=13,fontweight='bold',va='top')

def save(fig,name):
    for ext in ['png','pdf','svg','tiff']:
        kw={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}
        fig.savefig(OUT/'figures'/f'{name}.{ext}',dpi=300,bbox_inches='tight',**kw)
    plt.close(fig)

def marker_heatmap(fig,spec,d,z,letter,title):
    # Expression clustering is outcome-blind; annotation follows the same column order.
    order=leaves_list(linkage(z.T,method='average',metric='correlation'))
    sub=spec.subgridspec(3,2,height_ratios=[.65,.16,4],width_ratios=[20,.65],hspace=.03,wspace=.08)
    tree=fig.add_subplot(sub[0,0]); dendrogram(linkage(z.T,method='average',metric='correlation'),ax=tree,no_labels=True,color_threshold=0,above_threshold_color='#555555');tree.axis('off')
    tree.text(-.11,1.43,letter,transform=tree.transAxes,fontsize=13,fontweight='bold')
    tree.set_title(title,loc='left',fontsize=9,pad=12)
    anno=fig.add_subplot(sub[1,0]); vals=d.response.reindex(z.columns[order]).to_numpy() if 'response' in d else np.full(len(d),np.nan)
    anno.imshow(np.array([np.where(np.isnan(vals),2,vals)]),cmap=ListedColormap([RESPONSE[0],RESPONSE[1],'#D4D4D4']),vmin=0,vmax=2,aspect='auto');anno.set_axis_off()
    ax=fig.add_subplot(sub[2,0]); image=ax.imshow(z.iloc[:,order],cmap=CMAP,vmin=-2,vmax=2,aspect='auto',interpolation='nearest')
    ax.set_yticks(range(len(z)));ax.set_yticklabels(z.index,fontsize=6,fontstyle='italic')
    tick_indices=np.arange(0,len(d),2 if len(d)>35 else 1)
    ax.set_xticks(tick_indices);ax.set_xticklabels(z.columns[order][tick_indices],rotation=90,fontsize=5)
    for sp in ax.spines.values():sp.set_visible(False)
    cb=fig.add_subplot(sub[2,1]);fig.colorbar(image,cax=cb,ticks=[-2,0,2]);cb.set_title('Z',fontsize=7,pad=4)

def km(fig,spec,d,prog,letter,title):
    s=d.dropna(subset=['time','event']); threshold=d[prog].median()
    inner=spec.subgridspec(2,1,height_ratios=[3.8,1],hspace=.12)
    ax=fig.add_subplot(inner[0]);rt=fig.add_subplot(inner[1]);rt.axis('off')
    group=s[prog]>=threshold
    cols=['#578FB4','#D45868']; times=np.array([0,6,12,18,24,30]) if s.time.max()>750 else np.array([0,3,6,9,12,15,18])
    riskrows=[]
    for value,label,col in [(False,'Lower',cols[0]),(True,'Higher',cols[1])]:
        sub=s[group==value]; t=sub.time/30.4375
        f=KaplanMeierFitter().fit(t,sub.event,label=f'{label} (n={len(sub)})')
        f.plot_survival_function(ax=ax,color=col,ci_show=True,ci_alpha=.10,show_censors=True,lw=1.25,censor_styles={'marker':'+','ms':4})
        riskrows.append([int((t>=v).sum()) for v in times])
    lr=logrank_test(s.loc[group,'time'],s.loc[~group,'time'],s.loc[group,'event'],s.loc[~group,'event'])
    panel(ax,letter,title);ax.set(xlabel='',ylabel='Survival probability',ylim=(0,1.03),xticks=times,xlim=(0,s.time.max()/30.4375*1.02))
    ax.text(.04,.08,f'Descriptive log-rank P={lr.p_value:.3g}',transform=ax.transAxes,fontsize=6.8)
    ax.legend(frameon=False,fontsize=6.5,loc='upper right');ax.set_title(title,loc='left',pad=10)
    rt.set_xlim(ax.get_xlim());rt.set_ylim(-.7,2.3)
    rt.text(-.12,1.03,'At risk',transform=rt.transAxes,fontsize=6.5)
    for j,(label,col) in enumerate(zip(['Lower','Higher'],cols)):
        rt.text(-.08,1-j,label,ha='right',va='center',fontsize=6.5,color=col,transform=rt.get_yaxis_transform())
        for t,v in zip(times,riskrows[j]):rt.text(t,1-j,str(v),ha='center',va='center',fontsize=6.5)
    rt.text(.5,-.28,'Recorded follow-up (months)',ha='center',transform=rt.transAxes,fontsize=7)

def make_figures(ds,genes,res):
    markerlist=sum(PROGRAMMES.values(),[])
    fig=plt.figure(figsize=(12.1,14.5)); gs=fig.add_gridspec(3,3,height_ratios=[1.9,1.1,1.2],hspace=.43,wspace=.65)
    for j,name in enumerate(['GSE78220','GSE91061']):
        z=genes[name].loc[[g for g in markerlist if g in genes[name].index]]
        marker_heatmap(fig,gs[0,j],ds[name],z,'ab'[j],f'{name}: baseline marker expression')
    ax=fig.add_subplot(gs[0,2]); rr=res[res.endpoint=='Response']; names=['Memory','Dysfunction','Cytotoxic','Cycling','Interferon','Residency']
    a=rr.pivot(index='programme',columns='cohort',values='AUC').reindex(index=names,columns=['GSE78220','GSE91061'])
    im=ax.imshow(a,vmin=0,vmax=1,cmap=CMAP,aspect='auto');ax.set(xticks=[0,1],xticklabels=['Hugo','Riaz'],yticks=range(6),yticklabels=names)
    for i in range(6):
        for j in range(2):ax.text(j,i,f'{a.iloc[i,j]:.2f}',ha='center',va='center',fontsize=8)
    fig.colorbar(im,ax=ax,fraction=.055,pad=.06,label='AUC (higher score = response)'); panel(ax,'c','All programmes: response discrimination')
    ax=fig.add_subplot(gs[1,0]);panel(ax,'d','Fixed-score ROC: Memory programme')
    for name,col in [('GSE78220','#D75564'),('GSE91061','#5587B8')]:
        d=ds[name].dropna(subset=['response']);fpr,tpr,_=roc_curve(d.response,d.Memory); r=res[(res.cohort==name)&(res.endpoint=='Response')&(res.programme=='Memory')].iloc[0]
        ax.plot(fpr,tpr,lw=1.5,color=col,label=f'{name}: {r.AUC:.2f} [{r.AUC_low:.2f}, {r.AUC_high:.2f}]')
    ax.plot([0,1],[0,1],'--',color='#999999',lw=.7);ax.set(xlabel='False positive rate',ylabel='True positive rate',xlim=(0,1),ylim=(0,1.02));ax.legend(frameon=False,fontsize=6.5,loc='lower right')
    ax=fig.add_subplot(gs[1,1]);panel(ax,'e','Programme expression by response')
    for i,(name,prog) in enumerate([(n,p) for n in ['GSE78220','GSE91061'] for p in PRIMARY]):
        d=ds[name].dropna(subset=['response'])
        for val,offset in [(0,-.18),(1,.18)]:
            x=d.loc[d.response==val,prog].to_numpy(); pos=i+offset
            parts=ax.violinplot(x,[pos],widths=.3,showextrema=False)
            for body in parts['bodies']:body.set_facecolor(RESPONSE[val]);body.set_alpha(.23)
            ax.scatter(pos+RNG.uniform(-.08,.08,len(x)),x,s=7,color=RESPONSE[val],alpha=.8,edgecolors='none')
            ax.plot([pos-.10,pos+.10],[np.median(x)]*2,color=RESPONSE[val],lw=1.4)
    ax.set(xticks=range(4),xticklabels=['Hugo\nMemory','Hugo\nDysfunction','Riaz\nMemory','Riaz\nDysfunction'],ylabel='Bulk-adapted programme score');ax.axhline(0,c='#CCCCCC',lw=.6)
    ax.legend(handles=[Patch(color=RESPONSE[0],label='SD/PD'),Patch(color=RESPONSE[1],label='CR/PR')],frameon=False,fontsize=6.5,loc='upper left')
    ax=fig.add_subplot(gs[1,2]);panel(ax,'f','Continuous survival associations')
    rr=res[(res.endpoint!='Response')&res.programme.isin(PRIMARY)].reset_index(drop=True)
    for i,r in rr.iterrows():
        ax.errorbar(r.HR,i,xerr=[[r.HR-r.CI_low],[r.CI_high-r.HR]],fmt='o',ms=4,color=COLOURS[r.programme],capsize=2,lw=1)
        ax.text(1.02,i,f'{r.HR:.2f} [{r.CI_low:.2f}, {r.CI_high:.2f}]\nP={r.p:.3g}; q={r.q_primary_8_tests:.3g}',transform=ax.get_yaxis_transform(),va='center',fontsize=6)
    ax.axvline(1,color='#999999',ls='--',lw=.8);ax.set_xscale('log');ax.set(yticks=range(len(rr)),yticklabels=[f'{r.cohort}\n{r.endpoint}: {r.programme}' for _,r in rr.iterrows()],xlabel='HR per programme-score SD');ax.invert_yaxis()
    km(fig,gs[2,0],ds['GSE78220'],'Memory','g','GSE78220: OS by Memory score')
    km(fig,gs[2,1],ds['GSE78220'],'Dysfunction','h','GSE78220: OS by Dysfunction score')
    km(fig,gs[2,2],ds['GSE135222'],'Memory','i','GSE135222: recorded PFS by Memory')
    save(fig,'Figure_Clinical_Outcomes_9panels')

    fig=plt.figure(figsize=(11.8,14));gs=fig.add_gridspec(3,3,height_ratios=[1.8,1.1,1.1],hspace=.53,wspace=.65)
    z=genes['GSE135222'].loc[[g for g in markerlist if g in genes['GSE135222'].index]]
    marker_heatmap(fig,gs[0,0],ds['GSE135222'],z,'a','GSE135222: baseline TPM markers')
    ax=fig.add_subplot(gs[0,1]);panel(ax,'b','Continuous Cox results: all programmes')
    rr=res[res.endpoint!='Response'];aa=rr.pivot(index='programme',columns='endpoint',values='HR').reindex(index=list(PROGRAMMES),columns=['OS','PFS'])
    im=ax.imshow(np.log2(aa),cmap=CMAP,vmin=-1,vmax=1,aspect='auto');ax.set(yticks=range(6),yticklabels=list(PROGRAMMES),xticks=[0,1],xticklabels=['Hugo OS','Kim PFS'])
    for i in range(6):
        for j in range(2):ax.text(j,i,f'{aa.iloc[i,j]:.2f}',ha='center',va='center',fontsize=8)
    fig.colorbar(im,ax=ax,fraction=.06,label='log2(HR)')
    ax=fig.add_subplot(gs[0,2]);panel(ax,'c','PH assumption diagnostics')
    a=rr.pivot(index='programme',columns='endpoint',values='PH_rank_p').reindex(index=list(PROGRAMMES),columns=['OS','PFS'])
    ax.imshow(-np.log10(a),cmap='Purples',vmin=0,vmax=2,aspect='auto');ax.set(yticks=range(6),yticklabels=list(PROGRAMMES),xticks=[0,1],xticklabels=['Hugo OS','Kim PFS'])
    for i in range(6):
        for j in range(2):ax.text(j,i,f'{a.iloc[i,j]:.3f}',ha='center',va='center',fontsize=8)
    for i,(name,prog) in enumerate([('GSE78220','Memory'),('GSE78220','Dysfunction'),('GSE135222','Memory')]):
        ax=fig.add_subplot(gs[1,i]);panel(ax,'def'[i],f'{name}: {prog} and follow-up')
        d=ds[name].dropna(subset=['time','event']).sort_values(prog)
        for k,(_,r) in enumerate(d.iterrows()):
            col=COLOURS[prog];ax.plot([0,r.time/30.4375],[k,k],color=col,alpha=.25,lw=2)
            ax.scatter(r.time/30.4375,k,marker='x' if r.event else 'o',s=13,color=col,facecolors=None if r.event else 'white',linewidth=.8)
        ax.set(xlabel='Recorded follow-up (months)',ylabel=f'Patients ordered by {prog}',yticks=[])
        ax.text(.02,.97,'x event; open circle censored',transform=ax.transAxes,fontsize=6.5,va='top')
    for i,name in enumerate(ds):
        ax=fig.add_subplot(gs[2,i]);panel(ax,'ghi'[i],f'{name}: programme correlation')
        corr=ds[name][list(PROGRAMMES)].corr(method='spearman');im=ax.imshow(corr,cmap=CMAP,vmin=-1,vmax=1)
        ax.set(xticks=range(6),xticklabels=list(PROGRAMMES),yticks=range(6),yticklabels=list(PROGRAMMES));ax.tick_params(axis='x',rotation=70)
        for a in range(6):
            for b in range(6):ax.text(b,a,f'{corr.iloc[a,b]:.1f}',ha='center',va='center',fontsize=6)
    save(fig,'Figure_Supplement_Outcome_Audit_9panels')

def write_report(ds,res,flow):
    def md(df):
        columns=list(df.columns)
        lines=['| '+' | '.join(map(str,columns))+' |','| '+' | '.join(['---']*len(columns))+' |']
        lines+=['| '+' | '.join('' if pd.isna(x) else str(x) for x in row)+' |' for row in df.itertuples(index=False,name=None)]
        return '\n'.join(lines)
    lines=['# External clinical outcome extension','',
      'Post hoc exploratory analyses use fixed programmes from the single-cell study. Higher scores always define the positive direction in response AUCs. The same patients contribute to response and OS in GSE78220; they are not independent endpoint cohorts.','',
      md(flow),'', '## All observed results', '',md(res.round(4)),'',
      '## Endpoint provenance and interpretation', '',
      '- GSE78220: author-recorded complete/partial response versus stable/progressive disease; one on-treatment sample excluded, Pt27 biopsies averaged on log expression. OS uses author-recorded overall survival (days) and vital status. Missing OS time is excluded. Exact time origin is not specified in the GEO characteristic field; refer to the original publication before equating this with treatment-start OS.',
      '- GSE91061: baseline samples only, PRCR versus SD/PD; UNK excluded from response analysis, retained for unsupervised expression standardisation. No survival data were inferred from response labels.',
      '- GSE135222: pretreatment sampling and treatment-start PFS are verified in Jung et al. 2019 Methods. Original Supplementary Table 3 explicitly labels progression event=1 and censoring=0; all 27 GEO event indicators agree. GEO day values equal source PFS months times 30 for all 27 patients. This is PFS, not OS or ORR. The original durable-clinical-benefit category includes some stable disease and is not mapped to objective response.',
      '- All six programmes and 24 tests are disclosed. Eight Memory/Dysfunction tests have a separate BH correction; neither family controls the original manuscript analyses. Fixed scoring uses the cohort-wide expression distribution, including unknown-response patients, but no outcome labels. This is cohort-based exploration, not a deployable single-patient predictor.',
      '- Bulk-adapted scores reflect cell mixture as well as transcriptional state. Cross-cancer associations cannot establish a CD8-cell-intrinsic mechanism or treatment-specific predictive benefit. Cox estimates are unadjusted associations; PH rank-test P values are included.',
      '- Small cohorts and event counts do not justify high-dimensional model selection or a clinical nomogram. Median splits appear only in descriptive KM plots and were not optimised against survival.',
      '', '## Sources','',
      '- Hugo et al. Cell 2016;165:35–44. https://doi.org/10.1016/j.cell.2016.02.065; GEO https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220',
      '- Riaz et al. Cell 2017;171:934–949.e16. https://doi.org/10.1016/j.cell.2017.09.028; GEO https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061',
      '- Jung et al. Nature Communications 2019;10:4278. https://doi.org/10.1038/s41467-019-12159-9; GEO https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135222',
      '- Kim et al. Clinical Epigenetics 2020;12:119. https://doi.org/10.1186/s13148-020-00907-4',
      '- Gene identifiers: NCBI Homo_sapiens.gene_info (downloaded 2026-10-27) and the original single-cell Ensembl feature map.','']
    (OUT/'Clinical_Outcome_Results_EN.md').write_text('\n'.join(lines),encoding='utf-8')
    methods='''## Methods: exploratory external clinical endpoint analyses

We evaluated the original six gene programmes in three public bulk tumour transcriptome cohorts. The cohorts were selected post hoc to examine clinical endpoints beyond residual viable tumour fraction. Only baseline biopsies were included in response analyses. For GSE78220, log1p-transformed FPKM values from two biopsies of Pt27 were averaged before standardisation, yielding one observation per patient. We calculated each bulk-adapted programme score as the mean of gene-wise z-scores clipped to ±3. At least half of the predefined genes were required. Gene-wise standardisation was performed within each cohort without using outcome labels. Complete or partial response defined objective response; stable or progressive disease defined nonresponse. Unknown response labels were excluded from response tests. We calculated fixed-direction AUCs, stratified patient-bootstrap 95% confidence intervals (3,000 replicates), and two-sided Mann–Whitney tests. We fitted separate univariable Cox models to continuous programme scores standardised to one SD, excluding missing survival records. Proportional hazards were evaluated using a rank-transformed Schoenfeld residual test. We applied Benjamini–Hochberg correction across the 24 programme–endpoint tests and separately across the eight tests of Memory and Dysfunction. Kaplan–Meier plots used outcome-independent median cutoffs for visualisation. Analyses used patients as the independent units.

## Results: draft, pending endpoint-source completion

'''
    for name,row in flow.set_index('cohort').iterrows():
        methods+=f"{name} contributed {row.baseline_patients} baseline patients"
        if row.response_evaluable:methods+=f", including {row.responders}/{row.response_evaluable} objective responders"
        if row.survival_evaluable:methods+=f", with {row.survival_evaluable} survival-evaluable patients and {row.events} events"
        methods+='.'+'\n\n'
    for _,r in res[res.programme.isin(PRIMARY)].iterrows():
        if r.endpoint=='Response':methods+=f"The {r.programme} score in {r.cohort} had an AUC of {r.AUC:.3f} (95% CI {r.AUC_low:.3f}–{r.AUC_high:.3f}; P={r.p:.4g}; eight-test q={r.q_primary_8_tests:.4g}).\n\n"
        else:methods+=f"For {r.endpoint} in {r.cohort}, the unadjusted HR per SD of {r.programme} was {r.HR:.3f} (95% CI {r.CI_low:.3f}–{r.CI_high:.3f}; P={r.p:.4g}; eight-test q={r.q_primary_8_tests:.4g}; PH-test P={r.PH_rank_p:.4g}).\n\n"
    methods+='The eight Memory/Dysfunction tests did not meet the BH threshold of 0.05. None of the full 24-test family met this threshold. These results do not establish a clinical response or survival predictor. Bulk-adapted scores include tumour-composition effects and do not isolate CD8-cell-intrinsic state. GSE135222 baseline timing, PFS definition, day scale and event coding were verified against Jung et al. 2019 Methods and Supplementary Table 3. The exact origin of GSE78220 author-recorded OS is not specified in the downloaded GEO fields and remains an endpoint qualification.\n'
    (OUT/'Manuscript_Clinical_Extension_DRAFT_EN.md').write_text(methods,encoding='utf-8')
    manifest=[]
    for p in RAW.iterdir():
        if p.is_file():manifest.append({'file':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    (OUT/'raw_file_manifest.json').write_text(json.dumps(manifest,indent=2))
    import importlib.metadata
    versions={n:importlib.metadata.version(n) for n in ['numpy','pandas','scipy','matplotlib','scikit-learn','statsmodels','lifelines','openpyxl']}
    versions['python']=platform.python_version()
    (OUT/'software_versions.json').write_text(json.dumps(versions,indent=2))

if __name__=='__main__':
    if '--render-existing' in sys.argv:
        datasets={n:pd.read_csv(OUT/'tables'/f'{n}_patient_scores_and_outcomes.csv',index_col=0) for n in ['GSE78220','GSE91061','GSE135222']}
        genes={n:pd.read_csv(OUT/'source_data'/f'{n}_marker_zscores.csv',index_col=0) for n in datasets}
        results=pd.read_csv(OUT/'tables'/'all_programme_outcome_statistics.csv')
        flow=pd.read_csv(OUT/'tables'/'patient_flow.csv')
    else:
        datasets,genes=prepare()
        results,flow=analyse(datasets)
    print(flow.to_string(index=False),flush=True)
    print(results[results.programme.isin(PRIMARY)].to_string(index=False),flush=True)
    write_report(datasets,results,flow)
    make_figures(datasets,genes,results)
    shutil.copy2(__file__,OUT/'scripts'/'clinical_outcomes.py')
    plan=Path(__file__).with_name('outcome_analysis_plan.md')
    if plan.exists():shutil.copy2(plan,OUT/'scripts'/plan.name)
    print('Saved',OUT,flush=True)
