from pathlib import Path
import sys,os,json,io,re
ROOT=Path(r'D:\Downloads\Virtual_Biotech_Study_20261025\biology_rebuild_20261026');PREV=ROOT.parent
sys.path.insert(0,str(ROOT/'packages'));os.environ['MPLCONFIGDIR']=str(ROOT/'work')
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.patches import Rectangle,PathPatch
from matplotlib.path import Path as MPath
from matplotlib.lines import Line2D
from PIL import Image,ImageDraw
from pypdf import PdfWriter
R=ROOT/'results';F=ROOT/'figures';S=ROOT/'figure_source_data';F.mkdir(exist_ok=True);S.mkdir(exist_ok=True)
CO=['BCC','NSCLC','ORAL'];CN={'BCC':'BCC','NSCLC':'Lung','ORAL':'Oral','BREAST':'Breast'};CC={'BCC':'#3D6CA8','NSCLC':'#CF6A64','ORAL':'#31988C','BREAST':'#875DA6'}
P=['Memory','Cytotoxic','Dysfunction','Cycling','Interferon','Residency','Mixed'];PC=dict(zip(P,['#399C92','#D5A23F','#AC4F73','#906BB1','#5488B8','#D97549','#BFC6CE']));AB=['Mem','Cyto','Dys','Cyc','IFN','Res','Mix']
TARGETS=json.loads((PREV/'work/panel.json').read_text())['targets'];DEFS=json.loads((ROOT/'work/programme_definitions.json').read_text())['programmes'];CELLS={}
INK='#262B31';GREY='#B4BBC5';LIGHT='#E6E9ED';DIV=LinearSegmentedColormap.from_list('div',['#315F96','#B6CBE1','#FFFFFF','#EDB9B4','#BB4F63']);SEQ=LinearSegmentedColormap.from_list('seq',['#F4F2F6','#CFBDD7','#8B589F','#492F63'])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':8.5,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none','text.color':INK,'axes.labelcolor':INK,'axes.edgecolor':INK,'xtick.color':INK,'ytick.color':INK,'savefig.facecolor':'white'})
MAN=[]
def read(n):return pd.read_csv(R/(n+'.csv'))
def cells(c):
 if c not in CELLS:CELLS[c]=pd.read_csv(ROOT/'data'/(c+'_cells.csv.gz'))
 return CELLS[c]
def label(ax,ch,title):
 pos=ax.get_position();ax.figure.text(max(.012,pos.x0-.035),pos.y1+.012,ch,fontsize=11,fontweight='bold',va='bottom');ax.set_title(title,loc='left',pad=13)
def clean(ax):ax.spines[['top','right']].set_visible(False);ax.tick_params(width=.6,length=3)
def hm(ax,a,x=None,y=None,v=1,cmap=DIV):
 a=np.asarray(a,float);im=ax.pcolormesh(np.ma.masked_invalid(a),cmap=cmap,vmin=-v if cmap==DIV else 0,vmax=v,edgecolors='white',linewidth=.35);ax.set(ylim=(a.shape[0],0),xlim=(0,a.shape[1]));ax.spines[:].set_visible(False);ax.tick_params(length=0)
 ax.set_xticks(np.arange(a.shape[1])+.5,x if x is not None else ['']*a.shape[1]);ax.set_yticks(np.arange(a.shape[0])+.5,y if y is not None else ['']*a.shape[0]);return im
def cb(fig,im,rect,title,ticks=None):
 c=fig.colorbar(im,cax=fig.add_axes(rect),orientation='horizontal');c.set_label(title,fontsize=7.5);c.outline.set_visible(False);c.ax.tick_params(labelsize=7,length=2)
 if c.solids is not None:c.solids.set_rasterized(False)
 if ticks is not None:c.set_ticks(ticks)
def legend(fig,programmes=False,y=.015):
 co=PC if programmes else CC;names=P if programmes else CO
 fig.legend([Line2D([],[],marker='o',ls='',c=co[x],ms=4.5) for x in names],[x if programmes else CN[x] for x in names],loc='lower center',bbox_to_anchor=(.5,y),ncol=4 if programmes else 3,frameon=False,fontsize=7)
def umap(ax,m,color='programme',g=None):
 ax.spines[:].set_visible(False);ax.set_xticks([]);ax.set_yticks([])
 if g:
  z=m['log_'+g] if 'log_'+g in m else m[g];order=np.argsort(z.to_numpy());im=ax.scatter(m.UMAP1.iloc[order],m.UMAP2.iloc[order],c=z.iloc[order],s=1.5,cmap=SEQ,vmin=0,vmax=float(np.quantile(z,.99)) or 1,lw=0,rasterized=False);return im
 if color=='programme':
  for p in ['Mixed']+P[:-1]:
   u=m[m.programme.eq(p)];ax.scatter(u.UMAP1,u.UMAP2,c=PC[p],s=1.3,lw=0,alpha=.75,rasterized=False)
 elif color=='time':
  for t,col in [('pre','#BBC3CE'),('post','#805785')]:
   u=m[m.time.eq(t)];ax.scatter(u.UMAP1,u.UMAP2,c=col,s=1.4,lw=0,alpha=.65,rasterized=False)
 else:
  lev=sorted(m[color].unique());cm=plt.get_cmap('tab20')
  for j,k in enumerate(lev):
   u=m[m[color].eq(k)];ax.scatter(u.UMAP1,u.UMAP2,c=[cm(j%20)],s=1.3,lw=0,alpha=.7,rasterized=False)
   if color=='leiden':ax.text(u.UMAP1.median(),u.UMAP2.median(),str(k),ha='center',va='center',fontsize=6,fontweight='bold',bbox={'fc':'white','ec':'none','alpha':.65,'pad':.4})
 ax.set_xlim(m.UMAP1.min()-.6,m.UMAP1.max()+.6);ax.set_ylim(m.UMAP2.min()-.6,m.UMAP2.max()+.6)
def export(fig,name,title,panels,source):
 for e in ['png','svg','pdf']:fig.savefig(F/(name+'.'+e),dpi=300)
 b=io.BytesIO();fig.savefig(b,format='png',dpi=600);b.seek(0);Image.open(b).convert('RGB').save(F/(name+'.tif'),dpi=(600,600),compression='tiff_lzw');MAN.append({'figure':name,'title':title,'panels':panels,'source':source,'width_mm':183,'height_mm':float(fig.get_figheight()*25.4)});plt.close(fig);print('FIGURE',name,flush=True)
def save(d,n):d.to_csv(S/(n+'.csv'),index=False)
def forest(ax,u,names,col='mean',scale=100):
 for j,c in enumerate(CO):
  z=u[u.cohort.eq(c)].set_index('programme').reindex(names);y=np.arange(len(names))+(j-1)*.22
  ax.errorbar(z[col]*scale,y,xerr=[(z[col]-z.low)*scale,(z.high-z[col])*scale],fmt=['o','s','^'][j],ms=3.5,c=CC[c],lw=.8)
 ax.axvline(0,c=GREY,lw=.7);ax.set_yticks(range(len(names)),names);ax.set_ylim(len(names)-.5,-.5);clean(ax)
def f1():
 fig=plt.figure(figsize=(183/25.4,8.8))
 for j,c in enumerate(CO):
  m=cells(c);a=fig.add_axes([.065+j*.317,.68,.27,.235]);label(a,'abc'[j],f'{CN[c]}  |  {m.leiden.nunique()} clusters');umap(a,m,'leiden')
  b=fig.add_axes([.065+j*.317,.38,.27,.235]);label(b,'def'[j],CN[c]+' CXCL13');im=umap(b,m,g='CXCL13');cb(fig,im,[.09+j*.317,.345,.20,.008],'log(1 + CP10k)')
 ax=fig.add_axes([.09,.12,.86,.15]);label(ax,'g','Cell programmes across individual paired samples');pr=read('programme_patient');frames=[]
 for c in CO:
  u=pr[pr.cohort.eq(c)].pivot(index=['patient','time'],columns='programme',values='fraction').reindex(columns=P);u=u.reindex(sorted(u.index,key=lambda x:(x[0],0 if x[1]=='pre' else 1)));frames.append(u)
 w=pd.concat(frames);base=np.zeros(len(w))
 for p in P:ax.bar(range(len(w)),w[p],bottom=base,width=.92,color=PC[p],lw=0);base+=w[p].to_numpy()
 ax.set(xlim=(-1,len(w)),ylim=(0,1),yticks=[0,.5,1],ylabel='Cell fraction');labs=[p.replace('su','')+(' P' if t=='pre' else ' T') for p,t in w.index];ax.set_xticks(range(len(w)),labs,rotation=90,fontsize=5.8);clean(ax)
 legend(fig,True,-.005);export(fig,'Fig1','Cohort-resolved CD8 transcriptional atlas',7,['full_gene_counts','programme_patient.csv','atlas_coordinates'])
def f2():
 fig=plt.figure(figsize=(183/25.4,8.7));genes=sum([v[:4] for v in DEFS.values()],[]);rows=[];blocks=[];xlabs=[];bounds=[]
 for c in CO:
  m=cells(c);group=m.groupby(['patient','leiden']);means=group[['log_'+g for g in genes]].mean();means=means[group.size()>=5].groupby('leiden').mean();means=means.reindex(sorted(means.index,key=int));z=means.to_numpy().T;z=(z-np.nanmean(z,axis=1,keepdims=True))/np.maximum(np.nanstd(z,axis=1,keepdims=True),1e-6);blocks.append(z);xlabs.extend([str(k) for k in means.index]);bounds.append(z.shape[1]);rows.extend({'cohort':c,'cluster':p,'gene':g,'mean_logCP10k':means.loc[p,'log_'+g]} for p in means.index for g in genes)
 a=np.hstack(blocks);ax=fig.add_axes([.14,.51,.81,.40]);label(ax,'a','Marker profiles of within-cohort clusters');im=hm(ax,a,xlabs,genes,v=2)
 for x in np.cumsum(bounds)[:-1]:ax.axvline(x,c='white',lw=2.5)
 off=0
 for j,c in enumerate(CO):ax.text((off+bounds[j]/2)/sum(bounds),1.02,CN[c],transform=ax.transAxes,ha='center',c=CC[c],fontsize=8,weight='bold');off+=bounds[j]
 for y in [4,8,12,16,20]:ax.axhline(y,c='white',lw=2)
 ax.tick_params(axis='x',labelsize=6.5);ax.tick_params(axis='y',labelsize=7);plt.setp(ax.get_xticklabels(),rotation=55,ha='right')
 cb(fig,im,[.36,.425,.39,.011],'Within-cohort gene profile z score',[-2,0,2])
 eff=read('programme_effects');b=fig.add_axes([.15,.105,.31,.23]);label(b,'b','Paired programme redistribution');forest(b,eff[eff.endpoint.eq('fraction')],P[:-1]);b.set_xlabel('Post − pre fraction (pp)')
 cax=fig.add_axes([.65,.105,.30,.23]);label(cax,'c','Paired continuous programmes');forest(cax,eff[eff.endpoint.eq('score')],P[:-1],scale=1);cax.set_xlabel('Post − pre score')
 legend(fig,False,-.005);save(pd.DataFrame(rows),'Fig2_marker_profile');export(fig,'Fig2','Programme markers and patient-paired state changes',3,['programme_effects.csv','per-cell marker expression'])
def f3():
 fig=plt.figure(figsize=(183/25.4,8.9))
 for j,c in enumerate(CO):
  ax=fig.add_axes([.095+j*.31,.735,.25,.18]);label(ax,'abc'[j],CN[c]+' paired pseudobulk');u=pd.read_csv(R/(c+'_paired_DESeq2.csv'),index_col=0);q=u.padj.fillna(1).clip(lower=1e-200);sig=q<.05;color=np.where(sig&u.log2FoldChange.gt(0),'#C45C70',np.where(sig,'#4475A5','#CFD3DA'));ax.scatter(u.log2FoldChange,-np.log10(q),c=color,s=3,alpha=.7,lw=0);ax.axhline(-np.log10(.05),c=GREY,lw=.6,ls='--');ax.axvline(0,c=LIGHT,lw=.6);ax.set(xlabel='Post / pre log₂ fold change',ylabel='−log₁₀ adjusted P' if j==0 else '');clean(ax);ax.text(.02,.95,f'{int(sig.sum()):,} '+('gene' if sig.sum()==1 else 'genes')+' at FDR < 0.05',va='top',transform=ax.transAxes,fontsize=6.7)
 terms=['INTERFERON_GAMMA_RESPONSE','INTERFERON_ALPHA_RESPONSE','INFLAMMATORY_RESPONSE','TNFA_SIGNALING_VIA_NFKB','IL2_STAT5_SIGNALING','IL6_JAK_STAT3_SIGNALING','COMPLEMENT','APOPTOSIS','E2F_TARGETS','G2M_CHECKPOINT','MYC_TARGETS_V1','DNA_REPAIR','OXIDATIVE_PHOSPHORYLATION','GLYCOLYSIS','HYPOXIA','MTORC1_SIGNALING'];allg=[]
 for c in CO:
  u=read(c+'_Hallmark_GSEA');u['cohort']=c;allg.append(u)
 gs=pd.concat(allg);gs['short']=gs.Term.str.replace('HALLMARK_','',regex=False);w=gs.pivot(index='short',columns='cohort',values='NES').reindex(index=terms,columns=CO).astype(float);q=gs.pivot(index='short',columns='cohort',values='FDR q-val').reindex(index=terms,columns=CO).astype(float)
 ax=fig.add_axes([.28,.365,.24,.28]);label(ax,'d','Hallmark pathway context');im=hm(ax,w,[CN[c] for c in CO],[t.replace('_',' ').title() for t in terms],v=3);ax.tick_params(axis='y',labelsize=6.6)
 for y,x in np.argwhere(q.to_numpy()<.05):ax.plot(x+.5,y+.5,'o',ms=2.5,mfc='none',mec=INK,mew=.6)
 cb(fig,im,[.30,.315,.21,.009],'Normalized enrichment score',[-3,0,3])
 for j,c in enumerate(CO):
  ax=fig.add_axes([.66,.57-j*.125,.29,.075]);label(ax,'efg'[j],CN[c]+' enrichment');u=pd.read_csv(R/(c+'_GSEA_curves.csv.gz'))
  for term,col,name in [('HALLMARK_INTERFERON_GAMMA_RESPONSE','#8C5EA2','IFNγ'),('HALLMARK_E2F_TARGETS','#359A8B','E2F')]:
   z=u[u.term.eq(term)];ax.plot(z['rank']/max(z['rank'].max(),1),z.running_ES,c=col,lw=1,label=name)
  ax.axhline(0,c=LIGHT,lw=.6);ax.set(xlim=(0,1));clean(ax);ax.tick_params(labelsize=6)
  if j==0:ax.legend(frameon=False,fontsize=6,ncol=2,loc='lower right')
  if j==2:ax.set_xlabel('Fraction of ranked genes',fontsize=7)
  else:ax.tick_params(labelbottom=False)
 genes=['TCF7','CCR7','LTB','CXCL13','TOX','GZMB','NKG7','ISG15','IFIT3','MKI67','TOP2A','IL7R'];blocks=[]
 for c in CO:
  x=pd.read_csv(ROOT/'data'/(c+'_normalized_pseudobulk.csv.gz'),index_col=0);des=pd.read_csv(ROOT/'data'/(c+'_pseudobulk_design.csv')).set_index('sample');order=des.assign(time_order=des.condition.map({'pre':0,'post':1})).sort_values(['patient','time_order']).index;x=np.log1p(x.reindex(index=order,columns=genes));z=(x-x.mean())/x.std().replace(0,1);blocks.append(z.to_numpy().T)
 ax=fig.add_axes([.14,.095,.81,.13]);label(ax,'h','Paired patient transcript profiles');im=hm(ax,np.hstack(blocks),None,genes,v=2);ax.tick_params(axis='y',labelsize=6)
 ax.set_xticks(np.arange(50)+.5,['P','T']*25,fontsize=5.5)
 for x in [22,38]:ax.axvline(x,c='white',lw=2.5)
 cb(fig,im,[.39,.047,.31,.008],'Within-cohort log-expression z score',[-2,0,2]);save(gs,'Fig3_Hallmark_results');export(fig,'Fig3','Patient-paired transcriptional and pathway remodelling',8,['paired_DESeq2','Hallmark_GSEA','normalized_pseudobulk'])
def f4():
 fig=plt.figure(figsize=(183/25.4,8.8));sp=read('decomposition_support');de=read('decomposition_effects');dp=read('decomposition_patient')
 ax=fig.add_axes([.10,.745,.39,.17]);label(ax,'a','Programme support in each patient');clean(ax)
 ix=0
 for c in CO:
  u=sp[sp.cohort.eq(c)];x=np.arange(ix,ix+len(u));ax.vlines(x,u.pre_coverage*100,u.post_coverage*100,colors=CC[c],lw=1);ax.scatter(x,u.pre_coverage*100,fc='white',ec=CC[c],s=14,lw=.7);ax.scatter(x,u.post_coverage*100,c=CC[c],s=13,lw=0);ix+=len(u)+1
 ax.axhline(50,c=GREY,ls='--',lw=.7);ax.set(ylim=(0,103),ylabel='Cells in supported states (%)',xticks=[],xlabel='Individual paired patients')
 bx=fig.add_axes([.68,.745,.27,.17]);label(bx,'b','Supported versus full change');clean(bx)
 for c,m in zip(CO,['o','s','^']):
  u=dp[dp.cohort.eq(c)].groupby('target')[['full_total','supported_total']].mean()*100;bx.scatter(u.full_total,u.supported_total,c=CC[c],marker=m,s=17,ec='white',lw=.3)
 bx.plot([-20,20],[-20,20],c=GREY,ls='--',lw=.6);bx.set(xlabel='Full-population change (pp)',ylabel='Supported change (pp)')
 limit=max(15,float(de['mean'].abs().max()*100)*1.2)
 for j,c in enumerate(CO):
  ax=fig.add_axes([.14+j*.295,.15,.24,.48]);label(ax,'cde'[j],CN[c]+' change components');clean(ax)
  u=de[de.cohort.eq(c)]
  for term,col,dy in [('composition','#9274AD',-.17),('within_programme','#31988C',.17)]:
   z=u[u.term.eq(term)].set_index('target').loc[TARGETS];ax.barh(np.arange(24)+dy,z['mean']*100,height=.28,color=col,alpha=.9)
  z=u[u.term.eq('supported_total')].set_index('target').loc[TARGETS];ax.scatter(z['mean']*100,range(24),c=INK,s=11,marker='D',zorder=4);ax.axvline(0,c=GREY,lw=.6);ax.set(ylim=(23.7,-.7),xlim=(-limit,limit),yticks=range(24),yticklabels=TARGETS if j==0 else ['']*24,xlabel='Patient-equal contribution (pp)')
 fig.legend([Rectangle((0,0),1,1,fc='#9274AD'),Rectangle((0,0),1,1,fc='#31988C'),Line2D([],[],marker='D',c=INK,ls='',ms=4)],['Composition','Within programme','Supported total'],frameon=False,ncol=3,loc='lower center',bbox_to_anchor=(.5,.035),fontsize=7)
 export(fig,'Fig4','Separating programme redistribution from receptor change',5,['decomposition_support.csv','decomposition_patient.csv','decomposition_effects.csv'])
def ribbon_flow(ax,mat):
 # Exact patient-equal weights; node heights and ribbons share the same scale.
 a=mat.to_numpy();gap=.025;scale=1-gap*(len(P)-1);left=a.sum(1)*scale;right=a.sum(0)*scale;yl=np.cumsum(np.r_[0,left[:-1]])+np.arange(len(P))*gap;yr=np.cumsum(np.r_[0,right[:-1]])+np.arange(len(P))*gap;ol=yl.copy();orr=yr.copy()
 for i in range(7):
  for j in range(7):
   w=a[i,j]*scale
   if w<=0:continue
   y0,y1=ol[i],orr[j];verts=[(.08,y0),(.43,y0),(.57,y1),(.92,y1),(.92,y1+w),(.57,y1+w),(.43,y0+w),(.08,y0+w),(.08,y0)];codes=[MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.LINETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.CLOSEPOLY];ax.add_patch(PathPatch(MPath(verts,codes),fc=PC[P[i]],ec='none',alpha=.4));ol[i]+=w;orr[j]+=w
 for i,p in enumerate(P):
  ax.add_patch(Rectangle((.045,yl[i]),.035,left[i],fc=PC[p],lw=0));ax.add_patch(Rectangle((.92,yr[i]),.035,right[i],fc=PC[p],lw=0))
 ax.set(xlim=(0,1),ylim=(1,0));ax.axis('off');ax.text(.055,1.03,'Pre',transform=ax.transAxes,fontsize=7);ax.text(.84,1.03,'Post',transform=ax.transAxes,fontsize=7)
def f5():
 fig=plt.figure(figsize=(183/25.4,8.6));fl=read('shared_clone_links');ct=read('clone_history_programme');cl=read('clonality_patient')
 for j,c in enumerate(CO):
  ax=fig.add_axes([.075+j*.31,.64,.255,.265]);u=fl[fl.cohort.eq(c)];n=u.patient.nunique();mat=u.groupby(['pre_programme','post_programme']).patient_fraction.sum().unstack(fill_value=0).reindex(index=P,columns=P,fill_value=0)/n;assert np.isclose(mat.to_numpy().sum(),1);label(ax,'abc'[j],CN[c]+f'  |  {n} patients');ribbon_flow(ax,mat)
  save(mat.reset_index(),'Fig5_'+c+'_flow')
 ax=fig.add_axes([.13,.295,.31,.23]);label(ax,'d','Programme bias by clone history');data=[]
 for c in CO:
  u=ct[ct.cohort.eq(c)&ct.n_cells.ge(10)];w=u.pivot(index=['patient','programme'],columns='history',values='fraction').dropna();d=(w['Newly detected']-w['Previously detected']).groupby('programme').mean().reindex(P);data.append(d*100)
 im=hm(ax,np.asarray(data).T,[CN[c] for c in CO],P,v=20);cb(fig,im,[.14,.235,.29,.009],'New minus previous fraction (pp)',[-20,0,20])
 for j,(metric,title) in enumerate([('clonality','Observed TCR clonality'),('top10_fraction','Top-ten clone occupancy')]):
  ax=fig.add_axes([.63,.405-j*.235,.32,.13]);label(ax,'ef'[j],title);clean(ax)
  for k,c in enumerate(CO):
   w=cl[cl.cohort.eq(c)].pivot(index='patient',columns='time',values=metric)[['pre','post']].dropna()
   for _,r in w.iterrows():ax.plot([3*k,3*k+1],r,c=CC[c],lw=.6,alpha=.5,marker='o',ms=2)
   ax.plot([3*k,3*k+1],w.mean(),c=CC[c],lw=2,marker='o',ms=3.5)
  ax.set(ylim=(0,1),yticks=[0,.5,1],xticks=[0,1,3,4,6,7],xticklabels=['Pre','Post']*3)
 legend(fig,True,.035);export(fig,'Fig5','Patient-scoped clonal sharing across transcriptional programmes',6,['shared_clone_links.csv','clone_history_programme.csv','clonality_patient.csv'])
def f6():
 fig=plt.figure(figsize=(183/25.4,8.7));m=cells('BREAST')
 for j,(field,title) in enumerate([(None,'Breast CD8 programmes'),('IL7R_adt_corrected','CD127 corrected ADT'),('PDCD1_adt_corrected','PD-1 corrected ADT')]):
  ax=fig.add_axes([.065+j*.317,.68,.27,.23]);label(ax,'abc'[j],title)
  if field:
   z=m[field];ord=np.argsort(z);ax.scatter(m.UMAP1.iloc[ord],m.UMAP2.iloc[ord],c=z.iloc[ord],s=2,cmap=SEQ,lw=0,vmin=float(np.quantile(z,.01)),vmax=float(np.quantile(z,.99)));ax.axis('off');im=matplotlib.cm.ScalarMappable(norm=Normalize(np.quantile(z,.01),np.quantile(z,.99)),cmap=SEQ);cb(fig,im,[.085+j*.317,.63,.22,.009],'Corrected signal')
  else:umap(ax,m)
 cor=pd.read_csv(PREV/'results/protein_rna_correlations.csv');cor=cor[cor.grouping.eq('donor')];med=cor.groupby('target').rna_adt_rho.median().sort_values();genes=med.index.tolist()
 ax=fig.add_axes([.15,.19,.31,.32]);label(ax,'d','Within-donor RNA–ADT agreement');clean(ax)
 for i,g in enumerate(genes):
  z=cor[cor.target.eq(g)].rna_adt_rho;ax.scatter(z,[i]*len(z),s=15,c=GREY,ec='white',lw=.3);ax.plot([med[g]-.01,med[g]+.01],[i,i],c=CC['BREAST'],lw=3)
 ax.set(yticks=range(17),yticklabels=genes,xlabel='Spearman ρ',ylim=(-.6,16.6));ax.axvline(0,c=LIGHT,lw=.6)
 ax=fig.add_axes([.65,.19,.30,.32]);label(ax,'e','Independent donor tissue contrasts');clean(ax);dif=pd.read_csv(PREV/'results/protein_paired_donor_differences.csv');dif=dif[dif['mode'].eq('primary_cell_weighted')]
 for i,g in enumerate(genes):
  z=dif[dif.target.eq(g)].tumour_minus_peritumour;ax.scatter(z,i+np.linspace(-.12,.12,len(z)),c=PC['Residency'],s=13,ec='white',lw=.3);ax.scatter(z.mean(),i,c=INK,s=14,marker='D')
 ax.set(yticks=range(17),yticklabels=genes,xlabel='Tumor − peritumoral corrected ADT',ylim=(-.6,16.6));ax.axvline(0,c=LIGHT,lw=.6)
 legend(fig,True,.035);export(fig,'Fig6','Independent surface-protein measurements in breast CD8 cells',5,['BREAST_cells','protein_rna_correlations.csv','protein_paired_donor_differences.csv'])
def f7():
 fig=plt.figure(figsize=(183/25.4,8.8));rn=read('receptor_programme_rho_summary');ad=read('protein_programme_patient_rho');ad=ad.groupby(['target','programme']).rho.median().unstack().reindex(index=TARGETS,columns=P[:-1]);matrices=[]
 for j,c in enumerate(CO+['BREAST']):
  ax=fig.add_axes([.15+j*.21,.32,.17,.58]);label(ax,'abcd'[j],CN[c]+(' protein' if c=='BREAST' else ' RNA'))
  w=ad if c=='BREAST' else rn[rn.cohort.eq(c)].pivot(index='target',columns='programme',values='median').reindex(index=TARGETS,columns=P[:-1]);matrices.append(w);im=hm(ax,w,AB[:-1],TARGETS if j==0 else None,v=.5);ax.tick_params(axis='x',labelsize=6);plt.setp(ax.get_xticklabels(),rotation=55,ha='right');ax.tick_params(axis='y',labelsize=7)
  if c=='BREAST':
   for y,x in np.argwhere(w.isna().to_numpy()):ax.text(x+.5,y+.5,'×',ha='center',va='center',color=GREY,fontsize=6)
  save(w.reset_index(),'Fig7_'+c+'_programme_associations')
 cb(fig,im,[.34,.265,.43,.01],'Median within-patient Spearman ρ',[-.5,0,.5])
 selected=[('Memory','IL7R'),('Dysfunction','PDCD1'),('Dysfunction','LAG3')]
 ax=fig.add_axes([.15,.08,.80,.09]);label(ax,'e','Receptor-independent programmes and protein context');clean(ax)
 for k,c in enumerate(CO+['BREAST']):
  for i,(p,g) in enumerate(selected):
   vals=read('receptor_programme_patient_rho');z=ad # direct donor points use saved table
   v=read('protein_programme_patient_rho').query('programme==@p and target==@g').rho if c=='BREAST' else vals.query('cohort==@c and programme==@p and target==@g').rho
   xx=i*5+k;ax.scatter(xx+np.linspace(-.11,.11,len(v)),v,c=CC[c],s=13,ec='white',lw=.3);ax.plot([xx-.25,xx+.25],[v.median()]*2,c=INK,lw=1.2)
 ax.axhline(0,c=LIGHT,lw=.6);ax.set(xticks=[1.5,6.5,11.5],xticklabels=['Memory / IL7R','Dysfunction / PDCD1','Dysfunction / LAG3'],ylabel='Spearman ρ',ylim=(-.2,.8));fig.legend([Line2D([],[],marker='o',ls='',c=CC[c],ms=4) for c in CO+['BREAST']],[CN[c] for c in CO+['BREAST']],frameon=False,ncol=4,loc='lower center',bbox_to_anchor=(.5,.0),fontsize=7)
 export(fig,'Fig7','Conserved and assay-dependent receptor–programme relationships',5,['receptor_programme_patient_rho.csv','protein_programme_patient_rho.csv'])
def s1():
 fig=plt.figure(figsize=(183/25.4,8.8));axs=fig.subplots(3,3);fig.subplots_adjust(left=.08,right=.97,bottom=.08,top=.93,wspace=.35,hspace=.42)
 for j,c in enumerate(CO):
  m=cells(c)
  for i,col,title in [(0,'leiden','Leiden clusters'),(1,'patient','Patients'),(2,'time','Sampling time')]:
   ax=axs[i,j];label(ax,'abcdefghi'[i*3+j],CN[c]+' '+title);umap(ax,m,col)
  lev=sorted(m.patient.unique());pos=axs[1,j].get_position();fig.legend([Line2D([],[],marker='o',ls='',c=plt.get_cmap('tab20')(k%20),ms=2.5) for k in range(len(lev))],lev,loc='upper center',bbox_to_anchor=((pos.x0+pos.x1)/2,pos.y0-.002),ncol=3,frameon=False,fontsize=5.1,handletextpad=.25,columnspacing=.6,labelspacing=.2)
 fig.legend([Line2D([],[],marker='o',ls='',c=x,ms=3) for x in ['#BBC3CE','#805785']],['Pre','Post'],loc='lower center',bbox_to_anchor=(.5,.012),ncol=2,frameon=False,fontsize=7)
 export(fig,'FigS1','Unsupervised clusters and patient/time context',9,['atlas_coordinates'])
def s2():
 fig=plt.figure(figsize=(183/25.4,8.8));genes=sum(list(DEFS.values()),[]);parts=[]
 for c in CO:
  m=cells(c);u=m.groupby(['patient','programme'])[['log_'+g for g in genes]].mean().groupby('programme').mean().reindex(P);z=u.to_numpy().T;z=(z-np.nanmean(z,axis=1,keepdims=True))/np.maximum(np.nanstd(z,axis=1,keepdims=True),1e-6);parts.append(z)
 ax=fig.add_axes([.15,.20,.80,.69]);label(ax,'a','Complete receptor-independent marker programmes');im=hm(ax,np.hstack(parts),AB*3,genes,v=2);ax.tick_params(axis='y',labelsize=7);plt.setp(ax.get_xticklabels(),rotation=55,ha='right')
 for x in [7,14]:ax.axvline(x,c='white',lw=2)
 for j,c in enumerate(CO):fig.text(.15+.80*(j+.5)/3,.895,CN[c],ha='center',fontsize=7,color=CC[c],fontweight='bold')
 cb(fig,im,[.37,.11,.40,.011],'Within-cohort mean-profile z score',[-2,0,2]);export(fig,'FigS2','Complete marker programme profiles',1,['programme_definitions','full marker expression'])
def s3():
 fig=plt.figure(figsize=(183/25.4,8.8));gs=pd.concat([pd.read_csv(R/(c+'_Hallmark_GSEA.csv')).assign(cohort=c) for c in CO])
 def canonical(s):
  s=re.sub(r'[^A-Z0-9]+','_',s.upper()).strip('_')
  for a,b in [('TNF_ALPHA','TNFA'),('NF_KB','NFKB'),('IL_2','IL2'),('IL_6','IL6'),('G2_M','G2M')]:s=s.replace(a,b)
  return s if s.startswith('HALLMARK_') else 'HALLMARK_'+s
 terms=sorted(canonical(line.split('\t')[0]) for line in (ROOT/'references/MSigDB_Hallmark_2020.gmt').read_text().splitlines());w=gs.pivot(index='Term',columns='cohort',values='NES').reindex(index=terms,columns=CO).astype(float);q=gs.pivot(index='Term',columns='cohort',values='FDR q-val').reindex(index=terms,columns=CO).astype(float)
 ax=fig.add_axes([.53,.13,.36,.78]);label(ax,'a','Complete Hallmark enrichment');im=hm(ax,w,[CN[c] for c in CO],[s.replace('HALLMARK_','').replace('_',' ').title() for s in terms],v=3);ax.tick_params(axis='y',labelsize=6.4)
 for y,x in np.argwhere(q.to_numpy()<.05):ax.plot(x+.5,y+.5,'o',ms=2.5,mfc='none',mec=INK,mew=.6)
 for y,x in np.argwhere(w.isna().to_numpy()):ax.text(x+.5,y+.5,'×',ha='center',va='center',color=GREY,fontsize=6)
 cb(fig,im,[.53,.065,.36,.01],'Normalized enrichment score',[-3,0,3]);export(fig,'FigS3','Whole-family pathway enrichment results',1,['Hallmark_GSEA'])
def s4():
 fig=plt.figure(figsize=(183/25.4,8.8));support=read('decomposition_support');dp=read('decomposition_patient');axs=fig.subplots(2,3);fig.subplots_adjust(left=.12,right=.97,bottom=.12,top=.92,wspace=.4,hspace=.5)
 for j,c in enumerate(CO):
  u=dp[dp.cohort.eq(c)].groupby('patient')[['pre_coverage','post_coverage']].first();ax=axs[0,j];label(ax,'abc'[j],CN[c]+' support');clean(ax)
  for _,r in u.iterrows():ax.plot([0,1],r[['pre_coverage','post_coverage']]*100,c=CC[c],marker='o',lw=.8,ms=3,alpha=.6)
  ax.set(xticks=[0,1],xticklabels=['Pre','Post'],ylim=(0,100),ylabel='Cells retained (%)')
  ax=axs[1,j];label(ax,'def'[j],CN[c]+' partition sensitivity');clean(ax);z=read('decomposition_sensitivity_patient');z=z[z.cohort.eq(c)&z.minimum_state_cells.eq(5)].groupby(['state_definition','target']).within_programme.mean().unstack('state_definition');ax.scatter(z.programme*100,z.leiden*100,s=18,c=CC[c],ec='white',lw=.3);lo=min(z.min().min()*100,-1);hi=max(z.max().max()*100,1);ax.plot([lo,hi],[lo,hi],ls='--',c=GREY,lw=.7);ax.set(xlabel='Programme within-component (pp)',ylabel='Cluster within-component (pp)')
 export(fig,'FigS4','Supported-population coverage and state-granularity sensitivity',6,['decomposition_patient.csv','decomposition_sensitivity_patient.csv'])
def s5():
 fig=plt.figure(figsize=(183/25.4,8.8));d=read('external_antigen_clone_detection');genes=['IL7R','TCF7','PDCD1','LAG3','ENTPD1','CXCL13','GZMB','TOX'];donors=[]
 for p,u in d.groupby('patient'):
  n=u.groupby('antigen').size()
  if all(n.get(k,0)>=2 for k in ['MANA','viral']):donors.append(p)
 d=d[d.patient.isin(donors)].sort_values(['patient','antigen','clone']);ax=fig.add_axes([.15,.64,.79,.27]);label(ax,'a',f'Functionally labelled clones  |  {len(donors)} qualified patients');im=hm(ax,d[genes].to_numpy().T,['M' if k=='MANA' else 'V' for k in d.antigen],genes,v=1,cmap=SEQ);cb(fig,im,[.37,.575,.40,.01],'Clone-level transcript detection fraction',[0,.5,1])
 cursor=0
 for j,p in enumerate(donors):
  n=int(d.patient.eq(p).sum());fig.text(.15+.79*(cursor+n/2)/len(d),.917,p,ha='center',fontsize=7,color=['#5A86AC','#C77668'][j%2]);cursor+=n
  if cursor<len(d):ax.axvline(cursor,color='white',lw=2)
 for i,g in enumerate(genes):
  ax=fig.add_axes([.10+(i%4)*.23,.30-(i//4)*.20,.17,.13]);label(ax,'bcdefghi'[i],g);clean(ax)
  for j,p in enumerate(donors):
   u=d[d.patient.eq(p)].groupby('antigen')[g].mean();ax.plot([0,1],[u['viral'],u['MANA']],c=['#5A86AC','#C77668'][j%2],marker='o',lw=1.0,ms=3)
  ax.set(xticks=[0,1],xticklabels=['Viral','MANA'],ylim=(0,1),yticks=[0,.5,1])
 fig.legend([Line2D([],[],marker='o',ls='-',c=['#5A86AC','#C77668'][j%2],ms=3) for j in range(len(donors))],donors,loc='lower center',bbox_to_anchor=(.5,.025),ncol=2,frameon=False,fontsize=7)
 export(fig,'FigS5','External antigen-labelled clone context',9,['external_antigen_clone_detection.csv'])
def s6():
 fig=plt.figure(figsize=(183/25.4,8.8));m=cells('BREAST');axs=fig.subplots(2,3);fig.subplots_adjust(left=.08,right=.97,bottom=.08,top=.93,wspace=.3,hspace=.4)
 for j,(col,title) in enumerate([('leiden','Unsupervised clusters'),('patient','Independent donors'),('tissue','Tissue labels')]):
  ax=axs[0,j];label(ax,'abc'[j],title);umap(ax,m,col)
  if col in ['patient','tissue']:
   lev=sorted(m[col].unique());pos=ax.get_position();fig.legend([Line2D([],[],marker='o',ls='',c=plt.get_cmap('tab20')(k%20),ms=3) for k in range(len(lev))],lev,loc='upper center',bbox_to_anchor=((pos.x0+pos.x1)/2,pos.y0-.003),ncol=2,frameon=False,fontsize=6,handletextpad=.3,columnspacing=.7,labelspacing=.3)
 for j,g in enumerate(['TCF7','CXCL13','MKI67']):
  ax=axs[1,j];label(ax,'def'[j],g);im=umap(ax,m,g=g);cb(fig,im,[.10+j*.307,.055,.23,.008],'log(1 + CP10k)')
 export(fig,'FigS6','Independent breast atlas and marker context',6,['BREAST_cells'])
def s7():
 fig=plt.figure(figsize=(183/25.4,8.8));rna=read('RNA_partial_programme_rho');pro=read('protein_partial_programme_rho')
 for j,c in enumerate(CO+['BREAST']):
  ax=fig.add_axes([.15+j*.21,.30,.17,.59]);label(ax,'abcd'[j],CN[c]+(' ADT' if c=='BREAST' else ' RNA'));u=pro if c=='BREAST' else rna[rna.cohort.eq(c)];w=u.groupby(['target','programme']).rho.median().unstack().reindex(index=TARGETS,columns=P[:-1]);im=hm(ax,w,AB[:-1],TARGETS if j==0 else None,v=.5);plt.setp(ax.get_xticklabels(),rotation=55,ha='right',fontsize=6);ax.tick_params(axis='y',labelsize=7)
  for y,x in np.argwhere(w.isna().to_numpy()):ax.text(x+.5,y+.5,'×',ha='center',va='center',color=GREY,fontsize=6)
 cb(fig,im,[.34,.22,.43,.011],'Median partial rank correlation',[-.5,0,.5]);fig.text(.15,.10,'RNA: adjusted for library depth and sampling time.\nADT: adjusted for RNA depth, ADT depth, tissue and technical lane.',fontsize=8,linespacing=1.5)
 export(fig,'FigS7','Associations after measured technical and sampling adjustment',4,['RNA_partial_programme_rho.csv','protein_partial_programme_rho.csv'])
def s8():
 fig=plt.figure(figsize=(183/25.4,8.8))
 for j,c in enumerate(CO):
  curves=pd.read_csv(R/(c+'_GSEA_curves.csv.gz'));gs=read(c+'_Hallmark_GSEA').set_index('Term');rank=pd.read_csv(R/(c+'_Wald_ranking.csv'),index_col=0).iloc[:,0]
  for i,(term,name) in enumerate([('HALLMARK_INTERFERON_GAMMA_RESPONSE','Interferon-γ response'),('HALLMARK_E2F_TARGETS','E2F targets')]):
   x0=.095+j*.31;y0=.60-i*.43;width=.25;u=curves[curves.term.eq(term)];xx=u['rank'].to_numpy()/u['rank'].max();ax=fig.add_axes([x0,y0,width,.24]);label(ax,'abcdef'[i*3+j],CN[c]+' '+name);clean(ax);ax.plot(xx,u.running_ES,c='#36998B',lw=1.4);ax.axhline(0,c=GREY,lw=.6);ax.set(xlim=(0,1),ylabel='Running enrichment score' if j==0 else '');ax.tick_params(labelbottom=False);nes=float(gs.loc[term,'NES']);q=float(gs.loc[term,'FDR q-val']);ax.text(.02,.03,f'NES = {nes:.2f}; '+('q < 0.001' if q<.001 else f'q = {q:.3f}'),transform=ax.transAxes,fontsize=6.5)
   rug=fig.add_axes([x0,y0-.029,width,.021]);rug.vlines(xx[u.hit.to_numpy()],0,1,color=INK,lw=.3);rug.set(xlim=(0,1));rug.axis('off')
   metric=fig.add_axes([x0,y0-.09,width,.052]);metric.fill_between(np.linspace(0,1,len(rank)),rank.to_numpy(),0,color='#B9C4D0',lw=0);metric.axhline(0,c=GREY,lw=.5);metric.set(xlim=(0,1),xlabel='Fraction of ranked genes');metric.tick_params(labelsize=6,length=2);clean(metric)
 export(fig,'FigS8','Full running enrichment profiles with gene hits and rank metrics',6,['GSEA_curves','Wald_ranking','Hallmark_GSEA'])
def finish():
 (ROOT/'figure_manifest.json').write_text(json.dumps(MAN,indent=2));w=PdfWriter()
 for item in MAN:w.append(F/(item['figure']+'.pdf'),outline_item=item['figure']+' '+item['title'])
 w.write(ROOT/'Biology_Figure_Atlas.pdf');tiles=[]
 for it in MAN:
  im=Image.open(F/(it['figure']+'.png')).convert('RGB');im.thumbnail((580,740));tile=Image.new('RGB',(620,790),'white');tile.paste(im,((620-im.width)//2,30));ImageDraw.Draw(tile).text((20,10),it['figure'],fill='black');tiles.append(tile)
 sheet=Image.new('RGB',(2480,790*((len(tiles)+3)//4)),'#E8EBEF')
 for i,t in enumerate(tiles):sheet.paste(t,((i%4)*620,(i//4)*790))
 sheet.save(ROOT/'Biology_Figure_Overview.png')
if __name__=='__main__':
 jobs=[f1,f2,f3,f4,f5,f6,f7,s1,s2,s3,s4,s5,s6,s7,s8];chosen=set(sys.argv[1:])
 for f in jobs:
  if not chosen or f.__name__ in chosen:f()
 if not chosen:finish()
