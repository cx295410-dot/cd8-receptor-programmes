from pathlib import Path
import os,sys,json,io,string
ROOT=Path(r'D:\Downloads\Virtual_Biotech_Study_20261025\biology_rebuild_20261026');OUT=ROOT/'reference_expansion_20261026'
sys.path.insert(0,str(ROOT/'packages'));os.environ['MPLCONFIGDIR']=str(OUT/'work')
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Circle,FancyArrowPatch,Rectangle
from scipy.ndimage import gaussian_filter
from scipy.cluster.hierarchy import linkage,leaves_list,dendrogram
from scipy import stats
from PIL import Image,ImageDraw
from pypdf import PdfWriter
CO=['BCC','NSCLC','ORAL'];ALL=CO+['BREAST'];CN={'BCC':'BCC','NSCLC':'Lung','ORAL':'Oral','BREAST':'Breast'}
P=['Memory','Cytotoxic','Dysfunction','Cycling','Interferon','Residency','Mixed']
PC=dict(zip(P,['#91D650','#42BBA6','#F49A83','#B55DCE','#5B9FD0','#E2B54A','#C7C7C7']))
CC=dict(zip(ALL,['#7EBA40','#36A9A2','#EC8A75','#AD6AC7']));PRE='#4DBBD5';POST='#E64B35';INK='#242424';GRAY='#B8BEC4'
DIV=LinearSegmentedColormap.from_list('reference_div',['#3B4CC0','#FAFAFA','#D73027']);DIV.set_bad('#E6E6E6')
ROSE=LinearSegmentedColormap.from_list('expression_density',['#FFFFFF','#FBD5DD','#F68095','#DB1F48'])
DEFS=json.loads((ROOT/'work/programme_definitions.json').read_text())['programmes'];TARGETS=json.loads((ROOT.parent/'work/panel.json').read_text())['targets'];CACHE={}
MAN=json.loads((OUT/'expanded_figure_manifest.json').read_text()) if (OUT/'expanded_figure_manifest.json').exists() else []
CAPTIONS=json.loads((OUT/'Expanded_Figure_Legends_EN.json').read_text()) if (OUT/'Expanded_Figure_Legends_EN.json').exists() else {}
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.linewidth':.6,'pdf.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white','savefig.facecolor':'white','axes.spines.top':False,'axes.spines.right':False,'xtick.major.size':2.5,'ytick.major.size':2.5,'text.color':INK,'axes.labelcolor':INK,'axes.edgecolor':INK,'legend.frameon':False})
def read(name):return pd.read_csv(OUT/'results'/(name+'.csv'))
def old(name):return pd.read_csv(ROOT/'results'/(name+'.csv'))
def source(d,name):d.to_csv(OUT/'figure_source_data'/(name+'.csv'),index=False)
def cells(c):
 if c not in CACHE:
  d=pd.read_csv(ROOT/'data'/(c+'_cells.csv.gz'));extra=OUT/'data'/(c+'_additional_markers.csv.gz')
  if extra.exists():e=pd.read_csv(extra)
  else:
   import anndata as ad
   a=ad.read_h5ad(ROOT/'data'/(c+'_counts.h5ad'),backed='r');gs=[g for g in ['GZMK','GZMA','CXCR4','CXCR3','CCL5','IFNG'] if g in a.var_names and 'log_'+g not in d]
   idx=a.obs_names.get_indexer(d.cell_id.astype(str));assert (idx>=0).all();raw=a[:,gs].X.toarray()[idx];a.file.close();e=pd.DataFrame({'cell_id':d.cell_id})
   for j,g in enumerate(gs):e['raw_'+g]=raw[:,j];e['log_'+g]=np.log1p(raw[:,j]*10000/d.library_UMI.to_numpy())
   e.to_csv(extra,index=False)
  d=d.merge(e,on='cell_id',validate='one_to_one');CACHE[c]=d
 return CACHE[c]
def clean(ax):
 ax.spines[['top','right']].set_visible(False);ax.tick_params(width=.6,length=2.5,pad=2)
def panel(ax,i,title):
 label=string.ascii_lowercase[i] if i<26 else 'a'+str(i-25)
 ax.annotate(label,xy=(0,1),xycoords='axes fraction',xytext=(-17,7),textcoords='offset points',weight='bold',fontsize=10,ha='left',va='bottom');ax.set_title(title,loc='left',pad=7,fontsize=7.8);clean(ax)
def grid(title,rows,cols,height=8.1):
 fig,axs=plt.subplots(rows,cols,figsize=(183/25.4,height),squeeze=False)
 fig.subplots_adjust(left=.085,right=.965,bottom=.12,top=.91,wspace=.48,hspace=.70)
 fig.text(.05,.974,title,fontsize=10,weight='bold',va='top');return fig,axs
def note(fig,text):fig.text(.065,.025,text,fontsize=6.5,va='bottom',linespacing=1.35)
def legend(fig,kind='programme',y=.065):
 colors=PC if kind=='programme' else {'Pre':PRE,'Post':POST}
 fig.legend([Line2D([],[],marker='o',ls='',color=v,ms=4) for v in colors.values()],list(colors),loc='lower center',bbox_to_anchor=(.52,y),ncol=7 if kind=='programme' else 2,columnspacing=1.2,handletextpad=.3,fontsize=6.2)
def sanitize(fig):
 for text in fig.findobj(match=matplotlib.text.Text):
  text.set_text(text.get_text().translate(str.maketrans('₀₁₂₃₄₅₆₇₈₉','0123456789')))
def export(fig,stem,title,panels,caption,sources):
 sanitize(fig)
 for ext in ['png','svg','pdf']:fig.savefig(OUT/'figures'/(stem+'.'+ext),dpi=300)
 buf=io.BytesIO();fig.savefig(buf,format='png',dpi=600);buf.seek(0);Image.open(buf).convert('RGB').save(OUT/'figures'/(stem+'.tif'),compression='tiff_lzw',dpi=(600,600))
 MAN[:]=[r for r in MAN if r['figure']!=stem];MAN.append({'figure':stem,'title':title,'panels':panels,'source':sources,'width_mm':183,'height_mm':fig.get_figheight()*25.4});MAN.sort(key=lambda r:int(r['figure'].replace('FigS','')));CAPTIONS[stem]=caption
 plt.close(fig);(OUT/'expanded_figure_manifest.json').write_text(json.dumps(MAN,indent=2),encoding='utf-8');(OUT/'Expanded_Figure_Legends_EN.json').write_text(json.dumps(CAPTIONS,indent=2),encoding='utf-8');print('FIGURE',stem,panels,flush=True)
def spatial(ax,d,values):
 x=d.UMAP1.to_numpy();y=d.UMAP2.to_numpy();padx=np.ptp(x)*.06;pady=np.ptp(y)*.06;extent=(x.min()-padx,x.max()+padx,y.min()-pady,y.max()+pady)
 h,xe,ye=np.histogram2d(x,y,bins=110,range=[extent[:2],extent[2:]]);w,_,_=np.histogram2d(x,y,bins=[xe,ye],weights=np.maximum(np.asarray(values,float),0))
 h=gaussian_filter(h,1.6);w=gaussian_filter(w,1.6);z=w/max(float(w.max()),1e-12);mask=h>h.max()*.01
 ax.imshow(np.where(mask,z,np.nan).T,extent=extent,origin='lower',cmap=ROSE,vmin=0,vmax=1,interpolation='bilinear',aspect='auto')
 if h.max()>0:ax.contour((xe[1:]+xe[:-1])/2,(ye[1:]+ye[:-1])/2,h.T,levels=[h.max()*.025],colors='#666666',linewidths=.45,linestyles='dashed')
 ax.set(xticks=[],yticks=[],xlim=extent[:2],ylim=extent[2:]);ax.spines[:].set_visible(False)
 return z
def violin(ax,groups,colors,positions=None,width=.75):
 pos=list(range(len(groups))) if positions is None else positions
 for p,v,col in zip(pos,groups,colors):
  v=np.asarray(v,float);v=v[np.isfinite(v)]
  if len(v)>2 and np.ptp(v)>1e-8:
   body=ax.violinplot([v],positions=[p],widths=width,showextrema=False,points=70)['bodies'][0];body.set_facecolor(col);body.set_edgecolor('none');body.set_alpha(.8)
  if len(v):
   q1,med,q3=np.quantile(v,[.25,.5,.75]);ax.plot([p,p],[q1,q3],c=INK,lw=1.1);ax.plot([p-.12,p+.12],[med,med],c=INK,lw=.8)
 clean(ax)
def paired_lines(ax,u,value,scale=1,colors=None,labels=('Pre','Post')):
 times=['pre','post'] if set(u.time.unique())<=set(['pre','post']) else ['Peri_Healthy','Tumour']
 w=u.pivot(index='patient',columns='time',values=value).reindex(columns=times).dropna()
 for j,(p,row) in enumerate(w.iterrows()):
  color=colors.get(p,GRAY) if colors else GRAY;ax.plot([0,1],row.to_numpy()*scale,c=color,lw=.6,alpha=.8,zorder=1)
  ax.scatter([0,1],row.to_numpy()*scale,s=12,c=[PRE,POST] if colors is None else [color,color],edgecolors='white',linewidths=.3,zorder=3)
 ax.set(xticks=[0,1],xticklabels=labels,xlim=(-.3,1.3));clean(ax);return w
def heat(ax,m,v=2,cmap=DIV,xlabels=True,ylabels=True):
 im=ax.imshow(m.to_numpy(),cmap=cmap,vmin=-v,vmax=v,aspect='auto',interpolation='nearest')
 ax.set_xticks(np.arange(m.shape[1]),m.columns if xlabels else [],rotation=90) if xlabels else ax.set_xticks([])
 ax.set_yticks(np.arange(m.shape[0]),m.index if ylabels else []) if ylabels else ax.set_yticks([])
 ax.tick_params(length=0,labelsize=6);ax.spines[:].set_visible(False);return im
def zscore(a):return (a-a.mean(axis=1,keepdims=True))/np.maximum(a.std(axis=1,keepdims=True),1e-8)
def atlas():
 manifests=[]
 for n in ['existing_figure_manifest.json','expanded_figure_manifest.json']:
  if (OUT/n).exists():manifests.extend(json.loads((OUT/n).read_text()))
 (OUT/'figure_manifest.json').write_text(json.dumps(manifests,indent=2));w=PdfWriter();tiles=[]
 for rec in manifests:
  stem=rec['figure'];w.append(OUT/'figures'/(stem+'.pdf'),outline_item=stem+' '+rec['title']);im=Image.open(OUT/'figures'/(stem+'.png')).convert('RGB');im.thumbnail((520,690));tile=Image.new('RGB',(550,730),'white');tile.paste(im,((550-im.width)//2,25));ImageDraw.Draw(tile).text((12,8),stem,fill='black');tiles.append(tile)
 w.write(OUT/'Complete_Reference_Style_Atlas.pdf')
 for start in range(0,len(tiles),12):
  sheet=Image.new('RGB',(2200,730*3),'#E5E5E5')
  for i,t in enumerate(tiles[start:start+12]):sheet.paste(t,((i%4)*550,(i//4)*730))
  sheet.save(OUT/'work'/f'Figure_Contact_{start//12+1}.png')
 print('ATLAS',len(manifests),'figures',sum(r['panels'] for r in manifests),'panels',flush=True)
