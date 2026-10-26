from pathlib import Path
import os,sys,json,io,string
ROOT=Path(r'D:\Downloads\Virtual_Biotech_Study_20261025\biology_rebuild_20261026')
REV=ROOT/'reference_expansion_20261026';CLIN=ROOT/'clinical_validation_20261026';OUT=ROOT/'multipanel_rebuild_20261026'
sys.path.insert(0,str(ROOT/'packages'));os.environ['MPLCONFIGDIR']=str(OUT/'work')
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,to_rgba,Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle,PathPatch,Circle,FancyArrowPatch
from matplotlib.path import Path as MPath
from scipy.ndimage import gaussian_filter
from scipy.cluster.hierarchy import linkage,leaves_list,dendrogram
from scipy import stats
from PIL import Image,ImageDraw
CO=['BCC','NSCLC','ORAL'];ALL=CO+['BREAST'];CN={'BCC':'BCC','NSCLC':'Lung','ORAL':'Oral','BREAST':'Breast'}
P=['Memory','Cytotoxic','Dysfunction','Cycling','Interferon','Residency','Mixed'];AB=['Mem','Cyto','Dys','Cyc','IFN','Res','Mix']
PC=dict(zip(P,['#91D650','#42BBA6','#F49A83','#B55DCE','#5B9FD0','#E2B54A','#C7C7C7']))
CC=dict(zip(ALL,['#7EBA40','#36A9A2','#EC8A75','#AD6AC7']));PRE='#4DBBD5';POST='#E64B35';INK='#242424';GRAY='#B8BEC4';GREY=GRAY
DIV=LinearSegmentedColormap.from_list('reference_div',['#3B4CC0','#FAFAFA','#D73027']);DIV.set_bad('#E6E6E6')
ROSE=LinearSegmentedColormap.from_list('expression_density',['#FFFFFF','#FBD5DD','#F68095','#DB1F48']);ROSE.set_bad('white')
DEFS=json.loads((ROOT/'work/programme_definitions.json').read_text())['programmes'];TARGETS=json.loads((ROOT.parent/'work/panel.json').read_text())['targets'];CACHE={}
plt.rcParams.update({'font.family':'Arial','font.size':7.4,'axes.titlesize':8,'axes.labelsize':7.3,'xtick.labelsize':7,'ytick.labelsize':7,'axes.linewidth':.55,'pdf.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white','savefig.facecolor':'white','axes.spines.top':False,'axes.spines.right':False,'xtick.major.size':2.5,'ytick.major.size':2.5,'text.color':INK,'axes.labelcolor':INK,'axes.edgecolor':INK,'legend.frameon':False})
def read(name):return pd.read_csv(REV/'results'/(name+'.csv'))
def old(name):return pd.read_csv(ROOT/'results'/(name+'.csv'))
def source(d,name,index=False):
 p=OUT/'figure_source_data'/(name+'.csv');p.parent.mkdir(parents=True,exist_ok=True);d.to_csv(p,index=index);return str(p)
def cells(c):
 if c not in CACHE:
  d=pd.read_csv(ROOT/'data'/(c+'_cells.csv.gz'));p=REV/'data'/(c+'_additional_markers.csv.gz')
  if p.exists():
   e=pd.read_csv(p);keep=['cell_id']+[k for k in e if k!='cell_id' and k not in d];d=d.merge(e[keep],on='cell_id',validate='one_to_one')
  CACHE[c]=d
 return CACHE[c]
def clean(ax):ax.spines[['top','right']].set_visible(False);ax.tick_params(width=.55,length=2.5,pad=2)
def panel(ax,i,title):
 label=string.ascii_lowercase[i] if isinstance(i,int) else str(i)
 pos=ax.get_position();ax.figure.text(max(.012,pos.x0-.034),pos.y1+.011,label,weight='bold',fontsize=10,va='bottom');ax.set_title(title,loc='left',pad=8,fontsize=7.9);clean(ax)
def grid(title,rows,cols,height=7.8):
 fig,axs=plt.subplots(rows,cols,figsize=(183/25.4,height),squeeze=False);fig.subplots_adjust(left=.09,right=.965,bottom=.12,top=.94,wspace=.48,hspace=.65);return fig,axs
def note(fig,text):fig.text(.065,.018,text,fontsize=6.5,va='bottom',linespacing=1.3)
def legend(fig,kind='programme',y=.04):
 colors=PC if kind=='programme' else {'Pre':PRE,'Post':POST}
 fig.legend([Line2D([],[],marker='o',ls='',color=v,ms=4) for v in colors.values()],list(colors),loc='lower center',bbox_to_anchor=(.52,y),ncol=7 if kind=='programme' else 2,columnspacing=1.,handletextpad=.3,fontsize=7)
def export(fig,stem,title,panelids,caption,sources):
 ids=list(string.ascii_lowercase[:panelids]) if isinstance(panelids,int) else list(panelids)
 for p in ['figures','figure_source_data','scripts','work']:(OUT/p).mkdir(parents=True,exist_ok=True)
 for t in fig.findobj(match=matplotlib.text.Text):t.set_text(t.get_text().translate(str.maketrans('₀₁₂₃₄₅₆₇₈₉','0123456789')))
 for ext in ['png','svg','pdf']:fig.savefig(OUT/'figures'/(stem+'.'+ext),dpi=300)
 b=io.BytesIO();fig.savefig(b,format='png',dpi=600);b.seek(0);Image.open(b).convert('RGB').save(OUT/'figures'/(stem+'.tif'),compression='tiff_lzw',dpi=(600,600))
 rec={'figure':stem,'title':title,'panels':len(ids),'panel_ids':ids,'source':sources,'width_mm':183,'height_mm':fig.get_figheight()*25.4,'caption':caption}
 (OUT/'work'/(stem+'_metadata.json')).write_text(json.dumps(rec,indent=2),encoding='utf-8');plt.close(fig);print('FIGURE',stem,len(ids),'panels',flush=True);return rec
def spatial(ax,d,values):
 x=d.UMAP1.to_numpy();y=d.UMAP2.to_numpy();dx=np.ptp(x)*.06;dy=np.ptp(y)*.06;extent=(x.min()-dx,x.max()+dx,y.min()-dy,y.max()+dy)
 h,xe,ye=np.histogram2d(x,y,bins=110,range=[extent[:2],extent[2:]]);w,_,_=np.histogram2d(x,y,bins=[xe,ye],weights=np.maximum(np.asarray(values,float),0));h=gaussian_filter(h,1.6);w=gaussian_filter(w,1.6);z=w/max(float(w.max()),1e-12)
 ax.imshow(np.where(h>h.max()*.01,z,np.nan).T,extent=extent,origin='lower',cmap=ROSE,vmin=0,vmax=1,interpolation='bilinear',aspect='auto')
 ax.contour((xe[1:]+xe[:-1])/2,(ye[1:]+ye[:-1])/2,h.T,levels=[h.max()*.025],colors='#777777',linewidths=.4,linestyles='dashed');ax.set(xticks=[],yticks=[],xlim=extent[:2],ylim=extent[2:]);ax.spines[:].set_visible(False);return z
def umap(ax,d,color='programme',g=None):
 if g:return spatial(ax,d,d['log_'+g])
 if color=='programme':
  for p in ['Mixed']+P[:-1]:
   u=d[d.programme.eq(p)];ax.scatter(u.UMAP1,u.UMAP2,c=PC[p],s=1.3,lw=0,alpha=.85,rasterized=True)
 elif color=='time':
  for t,col in [('pre',PRE),('post',POST)]:
   u=d[d.time.eq(t)];ax.scatter(u.UMAP1,u.UMAP2,c=col,s=1.3,lw=0,alpha=.6,rasterized=True)
 else:
  for j,(k,u) in enumerate(d.groupby(color)):ax.scatter(u.UMAP1,u.UMAP2,c=[plt.get_cmap('tab20')(j%20)],s=1.3,lw=0,rasterized=True)
 ax.set(xticks=[],yticks=[]);ax.spines[:].set_visible(False)
def heat(ax,m,v=2,cmap=DIV,xlabels=True,ylabels=True):
 m=pd.DataFrame(m);im=ax.imshow(m.to_numpy(),cmap=cmap,vmin=-v if cmap is DIV else 0,vmax=v,aspect='auto',interpolation='nearest');ax.set_xticks(np.arange(m.shape[1]),m.columns,rotation=90) if xlabels else ax.set_xticks([]);ax.set_yticks(np.arange(m.shape[0]),m.index) if ylabels else ax.set_yticks([]);ax.tick_params(length=0,labelsize=7);ax.spines[:].set_visible(False);return im
def zscore(a):
 a=np.asarray(a,float);return (a-np.nanmean(a,axis=1,keepdims=True))/np.maximum(np.nanstd(a,axis=1,keepdims=True),1e-8)
def violin(ax,groups,colors,positions=None,width=.75):
 pos=list(range(len(groups))) if positions is None else positions
 for p,v,col in zip(pos,groups,colors):
  v=np.asarray(v,float);v=v[np.isfinite(v)]
  if len(v)>2 and np.ptp(v)>1e-8:
   body=ax.violinplot([v],positions=[p],widths=width,showextrema=False,points=70)['bodies'][0];body.set_facecolor(col);body.set_edgecolor('none');body.set_alpha(.8)
  if len(v):
   q1,med,q3=np.quantile(v,[.25,.5,.75]);ax.plot([p,p],[q1,q3],c=INK,lw=1);ax.plot([p-.12,p+.12],[med,med],c=INK,lw=.8)
 clean(ax)
def paired_lines(ax,u,value,scale=1,colors=None,labels=('Pre','Post')):
 times=['pre','post'] if set(u.time.unique())<=set(['pre','post']) else ['Peri_Healthy','Tumour'];w=u.pivot(index='patient',columns='time',values=value).reindex(columns=times).dropna()
 for p,row in w.iterrows():
  col=colors.get(p,GRAY) if colors else GRAY;ax.plot([0,1],row.to_numpy()*scale,c=col,lw=.6,alpha=.8);ax.scatter([0,1],row.to_numpy()*scale,s=12,c=[PRE,POST] if colors is None else [col,col],edgecolors='white',linewidths=.3,zorder=3)
 ax.set(xticks=[0,1],xticklabels=labels,xlim=(-.3,1.3));clean(ax);return w
def cbar(fig,im,rect,title='',ticks=None,orientation='vertical'):
 cb=fig.colorbar(im,cax=fig.add_axes(rect),orientation=orientation);cb.outline.set_visible(False);cb.ax.tick_params(labelsize=6.6,length=2,pad=2)
 if title:cb.ax.set_title(title,fontsize=6.6,pad=3)
 if ticks is not None:cb.set_ticks(ticks)
 return cb
