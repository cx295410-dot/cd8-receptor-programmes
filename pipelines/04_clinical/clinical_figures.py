"""Reference-aligned vector and high-resolution clinical figures from verified patient tables."""
from pathlib import Path
import sys,os,json,io,shutil,string
ROOT=Path('D:/Downloads/Virtual_Biotech_Study_20261025/biology_rebuild_20261026');O=ROOT/'clinical_validation_20261026'
sys.path.insert(0,str(ROOT/'packages'));os.environ['MPLCONFIGDIR']=str(O/'work')
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage,leaves_list,dendrogram
from PIL import Image
from pypdf import PdfWriter
SPEC=json.loads((O/'analysis_protocol_frozen.json').read_text());P=list(SPEC['programmes']);T=SPEC['targets']
COL=dict(zip(P,['#91D650','#42BBA6','#F49A83','#B55DCE','#5B9FD0','#E2B54A']))
BLUE='#4DBBD5';RED='#E64B35';INK='#242424';GRAY='#A8ADB3'
DIV=LinearSegmentedColormap.from_list('reference_signed',['#3B4CC0','#FAFAFA','#D73027'])
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6,'ytick.labelsize':6,'axes.linewidth':.6,'pdf.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white','axes.spines.top':False,'axes.spines.right':False,'xtick.major.size':2.5,'ytick.major.size':2.5,'legend.frameon':False,'text.color':INK})
M=pd.read_csv(O/'results/patient_clinical_features.csv');A=pd.read_csv(O/'results/clinical_associations_prespecified.csv');S=json.loads((O/'results/clinical_summary.json').read_text());MAN=[];CAP={}
def panel(ax,letter,title):
 ax.annotate(letter,(0,1),xycoords='axes fraction',xytext=(-18,22),textcoords='offset points',fontsize=10,weight='bold',va='bottom');ax.set_title(title,loc='left',pad=22)
def grid(title,rows,cols,height=9):
 f,a=plt.subplots(rows,cols,figsize=(183/25.4,height),squeeze=False);f.subplots_adjust(left=.10,right=.96,top=.90,bottom=.095,wspace=.55,hspace=.72);f.text(.045,.98,title,weight='bold',fontsize=10,va='top');return f,a
def get(feature,subset='primary'):return A[(A.feature==feature)&(A.subset==subset)].iloc[0]
def scat(ax,u,col,feature='Memory',stats=None,label=False):
 for primary in [True,False]:
  v=u[u.primary_analysis.eq(primary)]
  ax.scatter(v[col],v.residual_tumor_pct,s=26,c=np.where(v.MPR,BLUE,RED) if primary else 'none',edgecolor=np.where(v.MPR,BLUE,RED),linewidth=.8,zorder=3)
  if label:
   for _,r in v.iterrows():ax.annotate(r.patient,(r[col],r.residual_tumor_pct),xytext=(3,4),textcoords='offset points',fontsize=5.5)
 ax.set(ylim=(-5,111),yticks=[0,25,50,75,100],ylabel='Residual viable tumor (%)',xlabel=feature+' score')
 if stats is not None:
  line=f"n = {int(stats['n'])}; r = {stats['rho']:.2f}; P = {stats['p']:.3f}"
  if 'q_secondary' in stats and pd.notna(stats['q_secondary']):line+=f"; q = {stats['q_secondary']:.3f}"
  ax.text(0,1.015,line,transform=ax.transAxes,va='bottom',fontsize=6.2,clip_on=False)
def legend(f,context=False,y=.025):
 handles=[Line2D([],[],marker='o',ls='',c=BLUE,ms=4),Line2D([],[],marker='o',ls='',c=RED,ms=4)];labels=['MPR (residual tumor <=10%)','Non-MPR (>10%)']
 if context:handles.append(Line2D([],[],marker='o',ls='',mfc='white',mec=GRAY,ms=4));labels.append('Previously used functional donor')
 f.legend(handles,labels,loc='lower center',bbox_to_anchor=(.53,y),ncol=3 if context else 2,fontsize=6,columnspacing=1)
def export(f,stem,title,n,caption,files):
 for ext in ['png','svg','pdf']:f.savefig(O/'figures'/(stem+'.'+ext),dpi=300)
 b=io.BytesIO();f.savefig(b,format='png',dpi=600);b.seek(0);Image.open(b).convert('RGB').save(O/'figures'/(stem+'.tif'),compression='tiff_lzw',dpi=(600,600))
 MAN.append(dict(figure=stem,title=title,panels=n,source=files,width_mm=183,height_mm=f.get_figheight()*25.4));CAP[stem]=caption;plt.close(f);print('FIGURE',stem,flush=True)
def main():
 f,axs=grid('Pathology-linked CD8 programmes in resected lung tumors',3,2,height=9.5)
 ax=axs[0,0];panel(ax,'a','Original pathological response');u=M.sort_values(['residual_tumor_pct','patient']).reset_index(drop=True)
 ax.bar(np.arange(len(u)),u.residual_tumor_pct,color=np.where(u.MPR,BLUE,RED),width=.72,edgecolor=np.where(u.primary_analysis,'none',INK),linewidth=.8)
 zeros=u[u.residual_tumor_pct.eq(0)];ax.scatter(zeros.index,np.zeros(len(zeros)),c=BLUE,edgecolors=np.where(zeros.primary_analysis,BLUE,INK),linewidths=.8,s=18,zorder=3)
 ax.axhline(10,c=GRAY,lw=.6,ls='--');ax.set(xticks=np.arange(len(u)),xticklabels=u.patient,ylim=(-3,109),yticks=[0,25,50,75,100],ylabel='Residual viable tumor (%)');ax.tick_params(axis='x',rotation=90,labelsize=5.8)
 ax=axs[0,1];panel(ax,'b','Memory programme and residual tumor');scat(ax,M[M.primary_analysis],'score_Memory',stats=get('Memory'))
 ax=axs[1,0];panel(ax,'c','Prespecified programme associations');z=A[(A.subset=='primary')&(A.kind=='programme')].set_index('feature').loc[P]
 for i,(name,r) in enumerate(z.iterrows()):ax.plot([r.low,r.high],[i,i],c=COL[name],lw=1.4);ax.scatter(r.rho,i,c=COL[name],s=24,zorder=3)
 ax.axvline(0,c=GRAY,lw=.6,ls='--');ax.set(yticks=np.arange(6),yticklabels=P,xlim=(-1.05,1.05),xticks=[-1,-.5,0,.5,1],xlabel='Patient-level Spearman r (95% bootstrap CI)');ax.invert_yaxis()
 ax=axs[1,1];panel(ax,'d','IL7R and Memory within each patient');r=pd.read_csv(O/'results/IL7R_Memory_patient_correlations.csv').sort_values('IL7R_Memory_rho');ax.scatter(r.IL7R_Memory_rho,np.arange(len(r)),c='#42BBA6',s=18);ax.axvline(0,c=GRAY,lw=.6,ls='--');ax.set(yticks=np.arange(len(r)),yticklabels=r.patient,xlim=(-.05,1),xlabel='Within-patient Spearman r');ax.tick_params(axis='y',labelsize=5.8)
 ax=axs[2,0];panel(ax,'e','Fixed receptor panel and residual tumor');q=A[(A.subset=='primary')&(A.kind=='target')].set_index('feature').loc[T];y=np.arange(len(T));ax.hlines(y,0,q.rho,color=np.where(q.rho>0,RED,BLUE),lw=.8);ax.scatter(q.rho,y,c=np.where(q.q_secondary<.05,INK,np.where(q.rho>0,RED,BLUE)),s=11,zorder=3);ax.axvline(0,c=GRAY,lw=.6);ax.set(yticks=y,yticklabels=T,xlim=(-1,1),xticks=[-1,-.5,0,.5,1],xlabel='Patient-level Spearman r');ax.invert_yaxis();ax.tick_params(axis='y',labelsize=5.5)
 ax=axs[2,1];panel(ax,'f','Primary association sensitivity');ss=pd.read_csv(O/'results/primary_sensitivities.csv');rows=[('Primary, 13 patients',get('Memory')),('All 15 patients',get('Memory','all_patients')),('Equal library weight',ss[ss.analysis=='equal_library'].iloc[0]),('Rank-based score',ss[ss.analysis=='rank_scoring'].iloc[0])]
 for i,(name,r) in enumerate(rows):ax.plot([r['low'],r['high']],[i,i],c='#36A9A2',lw=1.4);ax.scatter(r['rho'],i,c='#36A9A2',s=22)
 ax.axvline(0,c=GRAY,lw=.6,ls='--');ax.set(yticks=np.arange(4),yticklabels=[x[0] for x in rows],xlim=(-1.05,1.05),xticks=[-1,-.5,0,.5,1],xlabel='Spearman r (95% bootstrap CI)');ax.invert_yaxis()
 legend(f,y=.012)
 caption=f"a, Original residual viable tumor percentages for 15 patients with primary resection tumor GEX and author CD8 annotations. Black bar or point outlines identify the two previously used functional-clone donors; zero-residual tumors are shown by dots. Dashed line marks the original MPR threshold of 10%. b, Prespecified Memory score association in the 13 other patients. Each point is one patient; r is Spearman correlation and P is a two-sided 99,999-permutation P value. c, All six fixed programmes with marginal 95% patient-bootstrap intervals (10,000 resamples); no simultaneous-coverage claim. d, IL7R–Memory correlations within each of 15 patients, displayed without cell-level P values. e, All 24 prespecified target-receptor correlations in the primary set. Black dots, if present, denote BH q<0.05 in the 29-test secondary family (five nonprimary programmes and 24 targets). f, Primary Memory association and prespecified alternative analyses with the same bootstrap scheme. Expression was measured after neoadjuvant anti-PD-1 treatment, at surgery; these are contemporaneous pathology associations. Cell totals were {S['CD8_cells']:,} across 51 contributing libraries (52 audited) and {S['primary_CD8_cells']:,} in the primary patient set. MPR classification supplies color only; all association tests use the continuous endpoint."
 export(f,'Fig8','Pathology-linked CD8 programmes in resected lung tumors',6,caption,['patient_clinical_features.csv','clinical_associations_prespecified.csv','IL7R_Memory_patient_correlations.csv','primary_sensitivities.csv'])
def programmes():
 f,a=grid('Programme associations with continuous pathological outcome',3,2,9.3)
 for i,p in enumerate(P):
  ax=a.flat[i];r=get(p);panel(ax,string.ascii_lowercase[i],p);scat(ax,M,'score_'+p,p,stats=r)
 legend(f,True,y=.012)
 export(f,'FigS37','Programme associations with continuous pathological outcome',6,'a–f, Patient-mean scores for all six frozen receptor-independent programmes against original residual viable tumor percentage. Filled points are the 13 primary patients; open points are MD01-004 and MD01-005, previously shown in the functional-clone comparison and excluded from primary statistics. Statistics use only the 13 filled-point patients and the continuous outcome. The Memory comparison is the prespecified primary endpoint. Other programme P values belong to the 29-test BH secondary family. Color denotes original MPR classification; no expression threshold or direction was optimized. All samples are primary resection tumors; metastatic recurrence, normal lung, lymph node and stimulated cultures are excluded.',['patient_clinical_features.csv','clinical_associations_prespecified.csv'])
def clustered():
 f=plt.figure(figsize=(183/25.4,8.3));f.text(.045,.98,'Patient-resolved marker and receptor profiles',weight='bold',fontsize=10,va='top')
 for i,genes in enumerate([sum(SPEC['programmes'].values(),[]),T]):
  x0=.105+i*.455;w=.305;y0=.21;h=.60
  v=M.set_index('patient')[['log_'+g for g in genes]].T;v.index=genes;z=(v.to_numpy()-v.to_numpy().mean(1,keepdims=True))/np.maximum(v.to_numpy().std(1,keepdims=True),1e-8);rl=linkage(z,method='average',metric='euclidean');cl=linkage(z.T,method='average',metric='euclidean');ri=leaves_list(rl);ci=leaves_list(cl)
  da=f.add_axes([x0,.862,w,.075]);dendrogram(cl,ax=da,no_labels=True,color_threshold=0,above_threshold_color=INK);da.axis('off');f.text(x0-.025,.95,string.ascii_lowercase[i],weight='bold',fontsize=10,va='bottom');f.text(x0,.95,'Programme markers' if i==0 else 'Fixed receptor targets',fontsize=8,va='bottom')
  for collection in da.collections:collection.set_linewidth(.6)
  da=f.add_axes([x0-.063,y0,.055,h]);dendrogram(rl,ax=da,orientation='left',no_labels=True,color_threshold=0,above_threshold_color=INK);da.axis('off');da.invert_yaxis()
  for collection in da.collections:collection.set_linewidth(.6)
  ax=f.add_axes([x0,y0,w,h]);im=ax.imshow(np.clip(z[np.ix_(ri,ci)],-2,2),cmap=DIV,vmin=-2,vmax=2,aspect='auto',interpolation='nearest');ax.set(xticks=np.arange(len(ci)),xticklabels=v.columns[ci],yticks=np.arange(len(ri)),yticklabels=v.index[ri]);ax.tick_params(axis='x',rotation=90,length=0,labelsize=6);ax.tick_params(axis='y',length=0,labelsize=5.5);ax.yaxis.tick_right();ax.spines[:].set_visible(False)
  strip=f.add_axes([x0,.826,w,.022]);c=M.set_index('patient').loc[v.columns[ci]];strip.imshow(c.residual_tumor_pct.to_numpy()[None,:],cmap='Blues',vmin=0,vmax=100,aspect='auto');strip.set(xticks=[],yticks=[]);strip.spines[:].set_visible(False)
  for j,primary in enumerate(c.primary_analysis):
   if not primary:strip.plot(j,0,'o',mfc='white',mec=INK,ms=3,mew=.5)
  cb=f.add_axes([x0,.135,w,.012]);f.colorbar(im,cax=cb,orientation='horizontal',ticks=[-2,0,2]);cb.set_xlabel('Patient mean expression, row z score',fontsize=6.5,labelpad=2);cb.tick_params(labelsize=6)
  pd.DataFrame(z[np.ix_(ri,ci)],index=v.index[ri],columns=v.columns[ci]).to_csv(O/'results'/f'FigS38_{i}_heatmap_matrix.csv')
 f.text(.11,.04,'Top strip: residual viable tumor, 0–100% (white to blue).\nOpen strip markers: the two previously used functional donors.',fontsize=6.8)
 export(f,'FigS38','Patient-resolved marker and receptor profiles',2,'a, All 38 fixed programme markers. b, All 24 receptor targets excluded from those programme definitions. Values are patient-mean log1p(CP10k), standardized within gene across the 15 pathology-matched patients and clipped to ±2 for display. Both axes use average-linkage Euclidean hierarchical clustering calculated on the unclipped standardized matrix, with independently displayed row and column dendrograms. Top strips show original residual viable tumor percentage on a fixed 0–100% scale. Open strip circles identify the two previously used functional donors. Both primary and sensitivity patients are displayed; heatmap ordering is descriptive and does not determine any statistical test. No differential-expression gene selection was performed.',['patient_clinical_features.csv','FigS38_0_heatmap_matrix.csv','FigS38_1_heatmap_matrix.csv'])
def sensitivity():
 f,a=grid('Patient and sampling sensitivity of the pathology association',3,2,9.4);primary=M[M.primary_analysis];ss=pd.read_csv(O/'results/primary_sensitivities.csv')
 ax=a[0,0];panel(ax,'a','All pathology-matched patients');scat(ax,M,'score_Memory',stats=get('Memory','all_patients'))
 ax=a[0,1];panel(ax,'b','Equal weighting of tumor libraries');scat(ax,primary,'equal_library_score_Memory',stats=ss[ss.analysis=='equal_library'].iloc[0])
 ax=a[1,0];panel(ax,'c','Rank-based Memory programme');scat(ax,primary,'rank_Memory',stats=ss[ss.analysis=='rank_scoring'].iloc[0])
 ax=a[1,1];panel(ax,'d','Leave-one-patient-out analysis');u=pd.read_csv(O/'results/primary_leave_one_out.csv');ax.scatter(u.rho,np.arange(len(u)),c='#36A9A2',s=17);ax.axvline(get('Memory').rho,c=GRAY,ls='--',lw=.7);ax.axvline(0,c=GRAY,lw=.5);ax.set(yticks=np.arange(len(u)),yticklabels=u.omitted_patient,xlim=(-1,1),xlabel='Memory–residual tumor Spearman r');ax.tick_params(axis='y',labelsize=5.8)
 ax=a[2,0];panel(ax,'e','Sampled CD8 depth and pathology');scat(ax,M,'n_cells','CD8 cells sampled');ax.set_xscale('log');ax.set_xlabel('Author-annotated tumor CD8 cells (log scale)')
 ax=a[2,1];panel(ax,'f','IL7R expression and pathology');scat(ax,primary,'log_IL7R','Mean IL7R log(1 + CP10k)',stats=get('IL7R'));ax.set_xlabel('Patient-mean IL7R log(1 + CP10k)')
 legend(f,True,y=.012)
 export(f,'FigS39','Patient and sampling sensitivity of the pathology association',6,'a, Memory programme and continuous residual tumor percentage in all 15 patients. b, Primary-set analysis after giving each GEX library equal weight within a patient. c, Primary-set analysis using mean within-cohort percentile ranks of the six frozen Memory markers. d, All 13 leave-one-patient-out Spearman coefficients; the dashed line is the full primary estimate. e, Sampled CD8 cell number and pathology, with log-scaled cell counts; no depth filtering was optimized. f, Prespecified target IL7R and residual tumor percentage in the 13 primary patients. Statistics in a–c and f use two-sided patient-level permutation tests. Filled/open symbols retain the primary/context split and color denotes MPR. Sampling and pathology association do not establish pretreatment prediction or intervention effects.',['patient_clinical_features.csv','primary_sensitivities.csv','primary_leave_one_out.csv','clinical_associations_prespecified.csv'])
if __name__=='__main__':
 main();programmes();clustered();sensitivity();(O/'figure_manifest.json').write_text(json.dumps(MAN,indent=2));(O/'Figure_Legends_EN.json').write_text(json.dumps(CAP,indent=2));w=PdfWriter()
 for r in MAN:w.append(O/'figures'/(r['figure']+'.pdf'),outline_item=r['figure']+' '+r['title'])
 w.write(O/'Clinical_Validation_Figure_Atlas.pdf');shutil.copy2(__file__,O/'scripts'/Path(__file__).name)
