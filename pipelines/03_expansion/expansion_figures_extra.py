from expansion_plot_core import *
from scipy.ndimage import gaussian_filter1d

def networks():
 sh=read('programme_clone_sharing');cl=read('clone_patient_visit');rep=read('repertoire_patient_visit');fig,axs=grid('Patient-scoped clonal sharing and visit abundance',3,3,8.2)
 theta=np.linspace(0,2*np.pi,7,endpoint=False)+np.pi/2;xy=np.c_[np.cos(theta),np.sin(theta)]
 for i,c in enumerate(CO):
  u=sh[sh.cohort.eq(c)&sh.time.eq('post')];m=u.groupby(['programme1','programme2']).jaccard.median().unstack().reindex(index=P,columns=P)
  ax=axs[i,0];panel(ax,i*3,CN[c]+' observed sharing');ax.axis('off');ax.set(xlim=(-1.5,1.5),ylim=(-1.45,1.45))
  for j in range(7):
   for k in range(j):
    w=m.iloc[j,k]
    if np.isfinite(w) and w>0:ax.plot(xy[[j,k],0],xy[[j,k],1],lw=.3+min(4,w*12),c='#8B8B8B',alpha=.65,zorder=1)
  for j,p in enumerate(P):ax.scatter(*xy[j],s=95,c=PC[p],ec='white',lw=.8,zorder=3);ax.text(xy[j,0]*1.28,xy[j,1]*1.23,p,fontsize=5.8,ha='center',va='center')
  ax=axs[i,1];panel(ax,i*3+1,CN[c]+' median Jaccard');im=ax.imshow(m,cmap='Blues',vmin=0,vmax=1);ax.set_xticks(range(7),['Mem','Cyto','Dys','Cyc','IFN','Res','Mix'],rotation=90);ax.set_yticks(range(7),['Mem','Cyto','Dys','Cyc','IFN','Res','Mix']);ax.tick_params(labelsize=5.5,length=0);ax.spines[:].set_visible(False)
  ax=axs[i,2];panel(ax,i*3+2,CN[c]+' clone abundance');v=cl[cl.cohort.eq(c)];w=v.pivot(index=['patient','scoped_clone'],columns='time',values='frequency').fillna(0);tot=rep[rep.cohort.eq(c)].pivot(index='patient',columns='time',values='assigned_cells');patients=w.index.get_level_values('patient');x=np.maximum(w.pre.to_numpy(),.5/tot.loc[patients,'pre'].to_numpy());y=np.maximum(w.post.to_numpy(),.5/tot.loc[patients,'post'].to_numpy());ax.scatter(x,y,c=np.where(w.pre.eq(0),POST,PRE),s=3,alpha=.42,lw=0,rasterized=True);ax.set(xscale='log',yscale='log',xlabel='Pre clone fraction',ylabel='Post clone fraction');lims=(min(x.min(),y.min())*.7,max(x.max(),y.max())*1.4);ax.plot(lims,lims,c=GRAY,lw=.7,ls='--');ax.set(xlim=lims,ylim=lims)
  source(m.reset_index(),f'FigS22_{c}_median_sharing')
 cb=fig.colorbar(im,cax=fig.add_axes([.973,.42,.008,.17]));cb.set_ticks([0,.5,1]);cb.ax.tick_params(labelsize=5,length=2);cb.outline.set_visible(False)
 note(fig,'Post-treatment networks: median patient Jaccard; ≥2 cells per clone/programme. Abundance axes use a 0.5-cell display floor.')
 export(fig,'FigS22','Observed clonal sharing networks',9,'a–i, Within each patient and post-treatment visit, a clone enters a programme set if at least two of its cells occupy that programme. Edges and heatmaps show the median patient Jaccard index (intersection/union) for each pair of programme sets. Empty unions remain missing. The networks are undirected observed sharing graphs. Clone-abundance panels plot pre/post assigned-cell fractions for patient-scoped clonotypes, with a half-cell floor used only to display absent observations on logarithmic axes. Red points were not detected at baseline; blue points were detected. No antigen assignment or differentiation direction is inferred.',['programme_clone_sharing','clone_patient_visit','repertoire_patient_visit'])

def clone_heatmaps():
 cl=read('clone_patient_visit');genes=['IL7R','TCF7','CXCL13','PDCD1','LAG3','ENTPD1','TIGIT','CTLA4','GZMB','PRF1','NKG7','TOX','MKI67','STMN1','ZNF683','HSPA1A'];fig,axs=grid('Patient-balanced dominant-clone transcript profiles',3,2,9);fig.subplots_adjust(left=.12,right=.91,wspace=.40,hspace=.65,bottom=.13)
 for i,c in enumerate(CO):
  u=cl[cl.cohort.eq(c)];rank=u.groupby(['patient','scoped_clone']).frequency.max().reset_index();chosen=rank.sort_values(['patient','frequency'],ascending=[True,False]).groupby('patient').head(3).scoped_clone.tolist();lookup={r.scoped_clone:r.patient for r in u.itertuples()};mats=[u[u.time.eq(t)].set_index('scoped_clone').reindex(chosen)[genes].T for t in ['pre','post']];stack=np.hstack([m.to_numpy() for m in mats]);mu=np.nanmean(stack,axis=1,keepdims=True);sd=np.maximum(np.nanstd(stack,axis=1,keepdims=True),1e-8);z=[(m.to_numpy()-mu)/sd for m in mats];ro=leaves_list(linkage(np.nan_to_num(np.hstack(z)),method='average',optimal_ordering=True))
  for j,t in enumerate(['pre','post']):
   ax=axs[i,j];panel(ax,i*2+j,CN[c]+' · '+t);w=pd.DataFrame(z[j][ro],index=np.array(genes)[ro],columns=chosen);im=heat(ax,w,xlabels=False);ax.tick_params(labelsize=5.5);previous=None
   for k,cl_id in enumerate(chosen):
    patient=lookup[cl_id]
    if previous and patient!=previous:ax.axvline(k-.5,c='white',lw=1)
    previous=patient
   ax.set_xlabel(f'{len(chosen)} clones; top three within each patient',fontsize=6);w.to_csv(OUT/'figure_source_data'/(f'FigS23_{c}_{t}_clone_z.csv'))
 fig.colorbar(im,cax=fig.add_axes([.94,.43,.012,.18])).outline.set_visible(False);note(fig,'Clone selection uses the maximum pre/post within-patient frequency. Gray columns: clone not observed at that visit.')
 export(fig,'FigS23','Dominant-clone expression heatmaps',6,'a–f, The three clonotypes with the largest maximum pre/post assigned-cell frequency are retained from each eligible patient. Gene expression is averaged over the observed cells of each clone/visit; row standardization uses the combined observed pre/post clone profiles. Rows share an average-linkage Euclidean ordering across the two visits. White separators identify patient blocks. Missing clone/visit measurements are gray rather than zero. Selected clone identifiers and all expression values accompany the figures.',['clone_patient_visit','patient-scoped selected clone IDs'])

def clinical_context():
 d=cells('ORAL');pr=old('programme_patient');pr=pr[pr.cohort.eq('ORAL')];meta=d[['patient','Path_response','regimen']].drop_duplicates();assert meta.patient.is_unique;colors={'Low':'#56B4E9','Medium':'#D36ABB','High':'#C02B3D'};mapping={p:colors[r] for p,r in zip(meta.patient,meta.Path_response)};fig,axs=grid('Original oral clinical categories and patient-level changes',3,3,8.2)
 for i,p in enumerate(P[:-1]):
  ax=axs.flat[i];panel(ax,i,p);u=pr[pr.programme.eq(p)];paired_lines(ax,u,'fraction',100,mapping);ax.set_ylabel('Programme fraction (%)');ax.set_ylim(bottom=0)
 ax=axs.flat[6];panel(ax,6,'Patient programme shifts');mem=pr[pr.programme.eq('Memory')].pivot(index='patient',columns='time',values='fraction');dys=pr[pr.programme.eq('Dysfunction')].pivot(index='patient',columns='time',values='fraction');x=(mem.post-mem.pre)*100;y=(dys.post-dys.pre)*100
 for row in meta.itertuples():ax.scatter(x[row.patient],y[row.patient],c=mapping[row.patient],marker='o' if row.regimen=='Mono' else '^',s=28);ax.annotate(row.patient,(x[row.patient],y[row.patient]),fontsize=5.5,xytext=(3,3),textcoords='offset points')
 ax.axhline(0,c=GRAY,lw=.6);ax.axvline(0,c=GRAY,lw=.6);ax.set(xlabel='Memory change (pp)',ylabel='Dysfunction change (pp)')
 ax=axs.flat[7];panel(ax,7,'Independent patients by category');counts=meta.Path_response.value_counts().reindex(['Low','Medium','High'],fill_value=0);ax.bar(range(3),counts,color=[colors[x] for x in counts.index],width=.55);ax.set(xticks=range(3),xticklabels=counts.index,ylabel='Patients',yticks=range(5));ax.set_ylim(0,max(counts)+1)
 ax=axs.flat[8];panel(ax,8,'Regimen and pathological response');tab=pd.crosstab(meta.regimen,meta.Path_response).reindex(index=['Mono','Combo'],columns=['Low','Medium','High'],fill_value=0);base=np.zeros(2)
 for p in tab.columns:ax.bar(range(2),tab[p],bottom=base,c=colors[p],label=p) if False else ax.bar(range(2),tab[p],bottom=base,color=colors[p],width=.55);base+=tab[p].to_numpy()
 ax.set(xticks=range(2),xticklabels=tab.index,ylabel='Patients',yticks=range(5));fig.legend([Line2D([],[],marker='o',ls='',c=v,ms=4) for v in colors.values()],list(colors),loc='lower center',bbox_to_anchor=(.5,.065),ncol=3,fontsize=6.5);note(fig,'Original pathology groups are retained; High contains one patient. These panels are descriptive clinical context.')
 source(meta,'FigS24_original_clinical_categories');export(fig,'FigS24','Oral pathology and regimen context',9,'a–f, Individual patient programme fractions before and after treatment, colored by the original Low/Medium/High pathological response labels. g, Paired memory and dysfunction-programme changes; point shape denotes Mono or Combo source regimen. h–i, Independent patient counts, not cell counts, for the original clinical groups. High contains one patient, so these displays retain source categories and do not fit a response classifier or estimate an independent response effect. All eight patients in the separate primary lung paired cohort were originally marked Response=Yes; no lung responder/nonresponder comparison is created.',['programme_patient.csv','clinical_patient_visit','original source clinical metadata'])

def functional_validation():
 d=read('external_functional_clone_logexpr');auc=read('external_within_donor_auc');r=read('external_within_donor_roc');genes=['CXCL13','IL7R','TCF7','PDCD1','LAG3','ENTPD1','GZMB','TOX'];gc=dict(zip(genes,['#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#B09C85']));donors=sorted(d.patient.unique());fig,axs=grid('Independent functionally labelled clonotypes',3,4,8);fig.subplots_adjust(wspace=.65,hspace=.80,bottom=.16)
 for i,p in enumerate(donors):
  ax=axs.flat[i];panel(ax,i,p+' · marker ROC')
  for g in genes:
   z=r[r.patient.eq(p)&r.gene.eq(g)];ax.step(z.fpr,z.tpr,where='post',c=gc[g],lw=.9)
  ax.plot([0,1],[0,1],c=GRAY,ls='--',lw=.6);u=d[d.patient.eq(p)];ax.set(xlim=(0,1),ylim=(0,1),xlabel='False-positive rate',ylabel='True-positive rate');ax.text(.03,.05,f'{u.antigen.eq("MANA").sum()} MANA / {u.antigen.eq("viral").sum()} viral',transform=ax.transAxes,fontsize=5.8)
 ax=axs.flat[2];panel(ax,2,'Within-donor AUC');m=auc.pivot(index='gene',columns='patient',values='AUC').reindex(genes);im=ax.imshow(m,cmap=DIV,vmin=0,vmax=1,aspect='auto');ax.set_xticks(range(2),['004','005']);ax.set_yticks(range(8),genes,fontsize=5.5);ax.tick_params(length=0)
 for i in range(8):
  for j in range(2):ax.text(j,i,f'{m.iloc[i,j]:.2f}',ha='center',va='center',fontsize=5.7)
 ax.spines[:].set_visible(False)
 ax=axs.flat[3];panel(ax,3,'Clone-level score context')
 for ant,col in [('viral',PRE),('MANA',POST)]:
  for pi,p in enumerate(donors):
   u=d[d.antigen.eq(ant)&d.patient.eq(p)];ax.scatter(u.memory_score,u.exhaustion_score,s=24,c=col,marker=['o','^'][pi],ec='white',lw=.4)
 ax.set(xlabel='Memory rank score',ylabel='Exhaustion rank score')
 for j,g in enumerate(genes):
  ax=axs.flat[j+4];panel(ax,j+4,g)
  for pi,p in enumerate(donors):
   for ai,ant in enumerate(['viral','MANA']):
    u=d[d.patient.eq(p)&d.antigen.eq(ant)];x=ai+(pi-.5)*.15;ax.scatter(x+np.linspace(-.045,.045,len(u)),u[g],c=[PRE,POST][ai],s=14,marker=['o','^'][pi],ec='white',lw=.3);ax.plot([x-.06,x+.06],[u[g].mean()]*2,c=INK,lw=.9)
  ax.set(xticks=[0,1],xticklabels=['Viral','MANA'],xlim=(-.35,1.35),ylabel='Clone mean log expression')
 fig.legend([Line2D([],[],c=gc[g],lw=1) for g in genes],genes,loc='lower center',bbox_to_anchor=(.5,.075),ncol=4,fontsize=6,columnspacing=1.3);note(fig,'Source functional labels; 15 clones, 404 cells, two donors. ROC is descriptive within donor; positive class is MANA.')
 export(fig,'FigS25','Independent functional-clone expression and ranking',12,'a–b, Within-donor univariate ROC curves for mean log(1 + CP10k) transcript expression in independently functionally labelled MANA versus viral clones. Higher expression always ranks MANA as positive; curve direction is never reversed to maximize AUC. Each retained clone has at least three observed cells, and each qualified donor has at least two clones of each antigen class. c, Exact AUC values; all are independently verified by positive–negative clone pair concordance, including ties. d, Pre-existing source memory/exhaustion rank scores for the same clones, with circle/triangle identifying donors 004/005 and blue/red identifying viral/MANA labels. e–l, Individual clone means with within-donor/class mean segments. These 15 clones from two donors give descriptive functional context, without model training, a generalizable diagnostic-performance claim, or new antigen assignments.',['external_functional_clone_logexpr','external_within_donor_auc','external_within_donor_roc','original experimental antigen labels'])

def protein_density():
 d=cells('BREAST');genes=['IL7R','PDCD1','LAG3'];fig,axs=grid('Independent breast RNA and surface-protein context',3,3,8.3);rows=[]
 for j,g in enumerate(genes):
  for i,(col,title) in enumerate([('log_'+g,'RNA'),(g+'_adt_corrected','Corrected ADT')]):
   ax=axs[i,j];panel(ax,i*3+j,g+' · '+title);spatial(ax,d,d[col])
  ax=axs[2,j];panel(ax,6+j,g+' RNA / ADT');ax.scatter(d['log_'+g],d[g+'_adt_corrected'],s=2,c=GRAY,alpha=.25,lw=0,rasterized=True);rhos=[]
  for pi,(p,u) in enumerate(d.groupby('patient')):
   rho=stats.spearmanr(u['log_'+g],u[g+'_adt_corrected']).statistic;rhos.append(rho);rows.append({'donor':p,'gene':g,'rho':rho,'n_cells':len(u)});ax.scatter(u['log_'+g].mean(),u[g+'_adt_corrected'].mean(),s=22,c=[plt.get_cmap('tab10')(pi)],ec='white',lw=.4)
  ax.set(xlabel='RNA log(1 + CP10k)',ylabel='Corrected ADT');ax.text(.04,.96,f'Median donor ρ = {np.median(rhos):.3f}',transform=ax.transAxes,fontsize=6,va='top')
 cb=fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,1),cmap=ROSE),cax=fig.add_axes([.82,.073,.12,.009]),orientation='horizontal');cb.set_ticks([0,1],labels=['Min','Max']);cb.ax.tick_params(labelsize=5,length=2);cb.outline.set_visible(False)
 note(fig,'Densities are scaled within panel; corrected-ADT density weights use positive signal. Scatter plots retain negative corrected values.')
 source(pd.DataFrame(rows),'FigS26_donor_RNA_ADT_correlations');export(fig,'FigS26','Spatial RNA and protein measurement context',9,'a–f, RNA-expression and positive corrected-ADT weighted spatial densities for IL7R, PDCD1 and LAG3 on the independent breast CD8 embedding. Only the density weights truncate negative corrected ADT at zero. g–i, Individual-cell RNA/ADT scatter plots retain the full corrected-ADT range, with larger donor-mean points. Correlations are computed within each of the five donors and summarized by their median; cell counts do not substitute for donor replication. Density magnitudes are scaled within each panel and cannot rank proteins across antibodies.',['BREAST_cells','corrected ADT','within-donor correlations'])

def thresholds():
 d=read('detection_threshold_patient');fig,axs=grid('RNA-detection threshold sensitivity',3,3,7.5)
 for i,c in enumerate(CO):
  for j,g in enumerate(['CXCL13','CD8A','IL7R']):
   ax=axs[i,j];panel(ax,i*3+j,CN[c]+' · '+g);u=d[d.cohort.eq(c)&d.gene.eq(g)]
   for t,col in [('pre',PRE),('post',POST)]:
    for p,v in u[u.time.eq(t)].groupby('patient'):ax.plot(v.min_UMI,v.fraction*100,c=col,lw=.5,alpha=.25)
    means=u[u.time.eq(t)].groupby('min_UMI').fraction.mean();ax.plot(means.index,means*100,c=col,lw=1.3,marker='o',ms=3)
   ax.set(xticks=[1,2,3],xlabel='Minimum recorded UMI',ylabel='RNA-detected cells (%)',ylim=(0,100))
 legend(fig,'time',.055);note(fig,'Thin curves: individual patients. Thick curves: equal-patient means; these are technical detection thresholds.')
 export(fig,'FigS27','RNA-detection sensitivity',9,'a–i, Fractions of retained CD8 cells meeting one-, two- or three-UMI detection thresholds for CXCL13, CD8A and IL7R in each treatment cohort and visit. Thin lines are individual patients; thick lines are patient-equal means. Every patient contributes at every threshold, so fractions are monotonically non-increasing as the threshold increases. These technical sensitivity displays do not define tumour specificity or surface-protein positivity.',['detection_threshold_patient'])

def qc():
 d=read('qc_patient_visit');fig,axs=grid('Patient-level library and sampling context',4,3,8.8)
 for i,c in enumerate(ALL):
  for j,(a,b,title) in enumerate([('library_UMI','library_UMI','Library UMI'),('detected_genes','rna_n_genes','Detected genes'),('mitochondrial_fraction','rna_mito_fraction','Mitochondrial fraction')]):
   field=a if c!='BREAST' else b;u=d[d.cohort.eq(c)&d.metric.eq(field)];ax=axs[i,j];panel(ax,i*3+j,CN[c]+' · '+title);paired_lines(ax,u,'median',100 if j==2 else 1,labels=('Peri','Tumour') if c=='BREAST' else ('Pre','Post'));ax.set_ylabel('Median (%)' if j==2 else 'Patient/visit median');ax.set_ylim(bottom=0)
 note(fig,'One point per independent patient/visit. These summaries describe the frozen source-defined CD8 populations.')
 export(fig,'FigS28','Patient-level quality and sampling diagnostics',12,'a–l, Patient/visit median total gene-expression UMI, number of detected genes, and mitochondrial-count fraction for each source-defined CD8 population. Breast comparisons use peri-tumour versus tumour tissue, not treatment visits. All original analysis cells are retained; these diagnostics neither introduce an outcome-dependent cell filter nor assert that cohort-specific technical differences have been removed.',['qc_patient_visit','frozen source quality annotations'])

def differential_context():
 import itertools
 de={c:pd.read_csv(ROOT/'results'/(c+'_paired_DESeq2.csv'),index_col=0) for c in CO};gs={c:old(c+'_Hallmark_GSEA').set_index('Term') for c in CO};fig,axs=grid('Cohort-resolved transcriptional and pathway concordance',3,3,8)
 for j,(a,b) in enumerate(itertools.combinations(CO,2)):
  m=de[a][['log2FoldChange','padj']].join(de[b][['log2FoldChange','padj']],lsuffix='_a',rsuffix='_b',how='inner').replace([np.inf,-np.inf],np.nan).dropna(subset=['log2FoldChange_a','log2FoldChange_b']);ax=axs[0,j];panel(ax,j,CN[a]+' / '+CN[b]);sig=(m.padj_a<.05)|(m.padj_b<.05);ax.scatter(m.log2FoldChange_a,m.log2FoldChange_b,c=np.where(sig,POST,GRAY),s=2,alpha=.5,lw=0,rasterized=True);rho=stats.spearmanr(m.log2FoldChange_a,m.log2FoldChange_b).statistic;ax.axhline(0,c=GRAY,lw=.5);ax.axvline(0,c=GRAY,lw=.5);ax.set(xlabel=CN[a]+' log₂ fold change',ylabel=CN[b]+' log₂ fold change');ax.text(.02,.96,f'{len(m):,} genes; ρ={rho:.2f}',transform=ax.transAxes,fontsize=6,va='top');source(m.reset_index(),f'FigS29_{a}_{b}_gene_concordance')
  ax=axs[1,j];panel(ax,3+j,'Hallmark NES · '+CN[a]+' / '+CN[b]);m=gs[a][['NES','FDR q-val']].join(gs[b][['NES','FDR q-val']],lsuffix='_a',rsuffix='_b').dropna();ax.scatter(m.NES_a,m.NES_b,c=np.where((m['FDR q-val_a']<.05)|(m['FDR q-val_b']<.05),POST,GRAY),s=12,lw=0);ax.axhline(0,c=GRAY,lw=.5);ax.axvline(0,c=GRAY,lw=.5);ax.set(xlabel=CN[a]+' NES',ylabel=CN[b]+' NES');ax.text(.02,.96,f'{len(m)} matched sets',transform=ax.transAxes,fontsize=6,va='top')
 for j,c in enumerate(CO):
  ax=axs[2,j];panel(ax,6+j,CN[c]+' gene-level discoveries');m=de[c];counts=[int(((m.padj<.05)&(m.log2FoldChange<0)).sum()),int(((m.padj<.05)&(m.log2FoldChange>0)).sum())];ax.bar([0,1],counts,color=[PRE,POST],width=.55);ax.set(xticks=[0,1],xticklabels=['Down','Up'],ylabel='Genes at FDR < 0.05');ax.set_ylim(0,max(counts)*1.25+1)
  for x,n in enumerate(counts):ax.text(x,n+.03*(max(counts)+1),str(n),ha='center',fontsize=7)
 note(fig,'Frozen patient-blocked DESeq2 and Hallmark GSEA results. Red points have a discovery in either contributing cohort.')
 export(fig,'FigS29','Cross-cohort gene and pathway context',9,'a–c, Gene-level post/pre log2 fold changes for genes tested in both contributing patient-blocked models. Descriptive Spearman correlations summarize concordance across genes and do not create patient replication. d–f, NES for pathways eligible in both studies; red identifies GSEA q < 0.05 in either cohort. g–i, Numbers of up- and downregulated genes at within-cohort gene FDR < 0.05. These displays reuse the frozen differential-expression and enrichment results without refitting or selecting pathways for significance.',['frozen paired_DESeq2','frozen Hallmark_GSEA'])

def heterogeneity():
 d=read('programme_heterogeneity');fig,axs=grid('Programme composition and patient heterogeneity',3,3,8)
 for i,c in enumerate(CO):
  u=d[d.cohort.eq(c)];w=u.pivot(index=['patient','time'],columns='programme',values='fraction').reindex(columns=P);w=w.reindex(sorted(w.index,key=lambda x:(x[0],x[1]=='post')));ax=axs[i,0];panel(ax,i*3,CN[c]+' paired compositions');base=np.zeros(len(w))
  for p in P:ax.bar(range(len(w)),w[p],bottom=base,color=PC[p],width=.9,lw=0);base+=w[p].to_numpy()
  ax.set(ylim=(0,1),ylabel='Cell fraction',xticks=np.arange(len(w)),xticklabels=[p.replace('su','')+('P' if t=='pre' else 'T') for p,t in w.index]);ax.tick_params(axis='x',rotation=90,labelsize=4.8)
  ax=axs[i,1];panel(ax,i*3+1,CN[c]+' programme changes');m=u[u.programme.eq('Memory')].pivot(index='patient',columns='time',values='fraction');z=u[u.programme.eq('Dysfunction')].pivot(index='patient',columns='time',values='fraction');ax.scatter((m.post-m.pre)*100,(z.post-z.pre)*100,c=CC[c],s=20,ec='white',lw=.4);ax.axhline(0,c=GRAY,lw=.6);ax.axvline(0,c=GRAY,lw=.6);ax.set(xlabel='Memory change (pp)',ylabel='Dysfunction change (pp)')
  ax=axs[i,2];panel(ax,i*3+2,CN[c]+' composition entropy');entropy=(-(w*np.log(w.where(w>0,1))).sum(axis=1)).rename('entropy').reset_index();paired_lines(ax,entropy,'entropy');ax.set(ylabel='Shannon programme entropy',ylim=(0,np.log(7)+.1))
 legend(fig,y=.055);note(fig,'Programmes are operational transcript features; Mixed is retained. Programme entropy concerns observed composition.')
 export(fig,'FigS30','Patient-level programme heterogeneity',9,'a–i, Patient/visit programme compositions, paired memory-versus-dysfunction fraction changes, and Shannon entropy over the seven programme fractions. Each point in the change and entropy panels represents one patient or paired sample. The same programme definitions and complete retained-cell populations are used throughout; no shared cell-type taxonomy is inferred from these operational programmes.',['programme_heterogeneity'])

def dotplots():
 genes=sum(DEFS.values(),[]);fig,axs=grid('Complete programme-marker expression and detection',2,2,9.1);fig.subplots_adjust(left=.14,right=.92,wspace=.45,hspace=.30,bottom=.13,top=.92);rows=[]
 for k,c in enumerate(ALL):
  d=cells(c);parts=[]
  for (p,prog),u in d.groupby(['patient','programme']):
   if len(u)>=5:parts.append(pd.DataFrame({'gene':genes,'mean_log':u[['log_'+g for g in genes]].mean().to_numpy(),'fraction':u[['raw_'+g for g in genes]].gt(0).mean().to_numpy(),'patient':p,'programme':prog}))
  ag=pd.concat(parts).groupby(['gene','programme'])[['mean_log','fraction']].mean();m=ag.mean_log.unstack().reindex(index=genes,columns=P);z=zscore(m.to_numpy());ax=axs.flat[k];panel(ax,k,CN[c]);X,Y=np.meshgrid(np.arange(7),np.arange(len(genes)));f=ag.fraction.unstack().reindex(index=genes,columns=P).to_numpy();im=ax.scatter(X.ravel(),Y.ravel(),s=(f.ravel()*31+1),c=z.ravel(),cmap=DIV,vmin=-2,vmax=2,lw=.15,edgecolors='#EEEEEE');ax.set(xticks=range(7),xticklabels=P,yticks=range(len(genes)),yticklabels=genes,ylim=(len(genes)-.5,-.5));ax.tick_params(axis='x',rotation=90,labelsize=5.5);ax.tick_params(axis='y',labelsize=5.2,length=0);ax.spines[:].set_visible(False);rows.extend(ag.reset_index().assign(cohort=c).to_dict('records'))
 cb=fig.colorbar(im,cax=fig.add_axes([.95,.55,.01,.20]));cb.outline.set_visible(False);cb.ax.set_title('z',fontsize=6)
 handles=[plt.scatter([],[],s=v*31+1,c='#777777',lw=0) for v in [.25,.5,1]];fig.legend(handles,['25%','50%','100%'],loc='lower center',bbox_to_anchor=(.53,.048),ncol=3,fontsize=6,title='Detected fraction',title_fontsize=6)
 note(fig,'Color: within-gene programme z score. Dot area: patient-equal detected fraction; each contributing group contains ≥5 cells.')
 source(pd.DataFrame(rows),'FigS31_patient_equal_marker_dotplot');export(fig,'FigS31','Complete marker dotplots',4,'a–d, All 38 designated programme markers across the seven operational programmes in each cohort. Mean log expression and RNA-detection fractions are first computed within patient/programme groups with at least five cells, then averaged equally across contributing patients. Dot area represents detection fraction; color represents the within-gene z score of patient-equal programme means. These are marker-definition and expression summaries, with no cell-level hypothesis test.',['per-cell markers','patient-equal programme summaries'])

def coexpression():
 d=read('target_coexpression_patient');mats={c:d[d.cohort.eq(c)].groupby(['gene1','gene2']).rho.median().unstack().reindex(index=TARGETS,columns=TARGETS) for c in ALL};order=leaves_list(linkage(np.hstack([mats[c].fillna(0).to_numpy() for c in CO]),method='average',optimal_ordering=True));fig,axs=grid('Patient-level receptor coexpression structure',2,2,8.6);fig.subplots_adjust(left=.14,right=.91,bottom=.15,wspace=.38,hspace=.40)
 for k,c in enumerate(ALL):
  ax=axs.flat[k];panel(ax,k,CN[c]+' RNA');w=mats[c].iloc[order,order];im=heat(ax,w,v=1);ax.tick_params(labelsize=4.7);w.to_csv(OUT/'figure_source_data'/(f'FigS32_{c}_median_target_rho.csv'))
 fig.colorbar(im,cax=fig.add_axes([.945,.52,.012,.18])).outline.set_visible(False);note(fig,'Each matrix summarizes within-patient Spearman correlations. The RNA-derived gene order is shared across all four cohorts.')
 export(fig,'FigS32','Receptor RNA coexpression matrices',4,'a–d, Median within-patient Spearman correlations among all 24 receptor transcripts in each cohort. The three treatment-cohort median matrices jointly determine one average-linkage Euclidean row/column order, which is shared with the independent breast RNA matrix. Undefined correlations remain missing. Diagonal entries are self-correlations; these descriptive association matrices do not establish receptor regulation or interaction.',['target_coexpression_patient'])

def crosswalks():
 d=read('author_state_crosswalk');fig,axs=grid('Source annotations and operational programme crosswalks',2,2,6.6);fig.subplots_adjust(left=.17,right=.92,wspace=.50,hspace=.65,bottom=.17)
 for k,c in enumerate(ALL):
  ax=axs.flat[k];panel(ax,k,CN[c]);u=d[d.cohort.eq(c)];m=u.pivot(index='author_state',columns='programme',values='fraction').reindex(columns=P);im=ax.imshow(m,cmap=ROSE,vmin=0,vmax=1,aspect='auto');ax.set_xticks(range(7),P,rotation=90,fontsize=6);ax.set_yticks(range(len(m)),m.index,fontsize=6);ax.tick_params(length=0);ax.spines[:].set_visible(False)
  for i in range(len(m)):
   for j in range(7):ax.text(j,i,f'{m.iloc[i,j]*100:.0f}',ha='center',va='center',fontsize=5.5)
 fig.colorbar(im,cax=fig.add_axes([.95,.52,.012,.20])).outline.set_visible(False);note(fig,'Cells are grouped by the original source annotation; entries are row percentages. Operational programmes retain mixed features.')
 export(fig,'FigS33','Author annotation and programme crosswalks',4,'a–d, Within each original author-defined CD8 annotation, the percentage of retained cells assigned to each operational transcript programme. Each row sums to 100% before rounding. The analysis preserves source annotations and does not treat different original taxonomies as a shared ground truth. Numeric oral source states remain numeric rather than being renamed using unvalidated functional labels.',['author_state_crosswalk'])

def expression_histograms():
 fig,axs=grid('CD8, CXCL13 and IL7R expression distributions',4,3,8.8)
 for i,c in enumerate(ALL):
  d=cells(c)
  for j,g in enumerate(['CD8A','CXCL13','IL7R']):
   ax=axs[i,j];panel(ax,i*3+j,CN[c]+' · '+g);x=d['log_'+g].to_numpy();hist,edges=np.histogram(x,bins=100,range=(0,max(x.max(),1)),density=True);sm=gaussian_filter1d(hist,1,mode='nearest');ax.plot((edges[:-1]+edges[1:])/2,sm,c='#C95785' if g=='CD8A' else '#3A9E87',lw=1.2);ax.axvline(0,c=GRAY,ls='--',lw=.6);ax.set(xlabel='log(1 + CP10k)',ylabel='Smoothed density',xlim=(0,max(x.max(),1)));ax.text(.97,.94,f'Zero: {(x==0).mean()*100:.1f}%',transform=ax.transAxes,ha='right',va='top',fontsize=6)
 note(fig,'Full source-defined CD8 populations; no expression-based reclassification. Zero fractions are shown explicitly.')
 export(fig,'FigS34','Expression sparsity and density diagnostics',12,'a–l, Smoothed 100-bin expression histograms (Gaussian sigma one bin) for CD8A, CXCL13 and IL7R, using the full within-cohort expression ranges. The fraction of recorded zero values is shown in every panel. These displays describe transcript sparsity and do not redefine the original CD8 populations or impose a CXCL13-derived tumour-specificity label.',['frozen per-cell transcript expression'])

def paga():
 d=read('paga_programme_connectivity');fig,axs=grid('Undirected transcriptional neighbourhood connectivity',2,4,6.1);fig.subplots_adjust(left=.07,right=.96,wspace=.32,hspace=.35,bottom=.18);theta=np.linspace(0,2*np.pi,7,endpoint=False)+np.pi/2;xy=np.c_[np.cos(theta),np.sin(theta)]
 for k,c in enumerate(ALL):
  m=d[d.cohort.eq(c)].pivot(index='programme1',columns='programme2',values='connectivity').reindex(index=P,columns=P);ax=axs[0,k];panel(ax,k,CN[c]+' · PAGA');ax.axis('off');ax.set(xlim=(-1.5,1.5),ylim=(-1.5,1.5))
  for i in range(7):
   for j in range(i):
    w=m.iloc[i,j]
    if w>0:ax.plot(xy[[i,j],0],xy[[i,j],1],c='#888888',lw=.2+3*w,alpha=.65)
  for i,p in enumerate(P):ax.scatter(*xy[i],s=55,c=PC[p],ec='white',lw=.5,zorder=3);ax.text(xy[i,0]*1.25,xy[i,1]*1.25,['Mem','Cyto','Dys','Cyc','IFN','Res','Mix'][i],ha='center',va='center',fontsize=5.5)
  ax=axs[1,k];panel(ax,k+4,CN[c]+' · programme map');u=cells(c)
  for p in ['Mixed']+P[:-1]:
   v=u[u.programme.eq(p)];ax.scatter(v.UMAP1,v.UMAP2,c=PC[p],s=1,lw=0,rasterized=True)
  ax.set(xticks=[],yticks=[]);ax.spines[:].set_visible(False)
 legend(fig,y=.065);note(fig,'PAGA uses 15-neighbour graphs from frozen 30-PC coordinates; edge width represents undirected graph connectivity.')
 export(fig,'FigS35','Transcriptional graph connectivity',8,'a–d, Scanpy PAGA connectivities between the unchanged operational programme groups, computed on 15-neighbour graphs from each cohort’s original 30 principal components (random seed 20261026). All positive edges are displayed, with width proportional to connectivity; node positions are fixed for comparison. e–h, The same programme identities on the unchanged source UMAP coordinates. These are undirected neighbourhood summaries, without a root, pseudotime, lineage inference or observed differentiation trajectory.',['paga_programme_connectivity','frozen PCA and UMAP coordinates'])

def replacement_context():
 d=read('clone_replacement');rep=read('repertoire_patient_visit');fig,axs=grid('Observed clone replacement and repertoire sampling',3,3,7.7);rows=[]
 for i,c in enumerate(CO):
  u=d[d.cohort.eq(c)].merge(rep[rep.cohort.eq(c)&rep.time.eq('pre')][['patient','hill_D1']],on='patient',validate='one_to_one')
  for j,(x,title) in enumerate([('pre_cells','Baseline assigned CD8 cells'),('post_cells','Post-treatment assigned CD8 cells'),('hill_D1','Baseline effective clone number')]):
   ax=axs[i,j];panel(ax,i*3+j,CN[c]+' · sampling context');ax.scatter(u[x],u.new_post_clone_fraction*100,s=25,c=CC[c],ec='white',lw=.4);rho=stats.spearmanr(u[x],u.new_post_clone_fraction).statistic;ax.set(xlabel=title,ylabel='Newly detected post clones (%)',ylim=(0,105));ax.text(.04,.96,f'n={len(u)}; ρ={rho:.2f}',transform=ax.transAxes,fontsize=6,va='top');rows.append({'cohort':c,'x_measure':x,'rho':rho,'n_patients':len(u)})
 note(fig,'Baseline absence is observed non-detection. Sampling depth and repertoire complexity remain visible; correlations are descriptive.')
 source(pd.DataFrame(rows),'FigS36_sampling_correlations');export(fig,'FigS36','Clone replacement and sampling context',9,'a–i, Patient-level associations of the fraction of post-treatment clonotypes not detected at baseline with baseline assigned CD8 cell count, post-treatment assigned CD8 cell count, and baseline Shannon effective clone number. Every point is one independent patient; Spearman rho is descriptive and no patient is multiplied by its clone count. These displays quantify sampling context. Newly detected clones are not automatically newly generated clones, and baseline assigned CD8 count is not a whole-tumour T-cell infiltration measurement.',['clone_replacement','repertoire_patient_visit'])

if __name__=='__main__':
 jobs={'networks':networks,'clones':clone_heatmaps,'clinical':clinical_context,'functional':functional_validation,'protein':protein_density,'thresholds':thresholds,'qc':qc,'DE':differential_context,'heterogeneity':heterogeneity,'dotplots':dotplots,'coexpression':coexpression,'crosswalks':crosswalks,'histograms':expression_histograms,'paga':paga,'replacement':replacement_context}
 for key in sys.argv[1:] or jobs:jobs[key]()
 atlas()
