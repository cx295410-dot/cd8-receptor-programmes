from expansion_plot_core import *

def densities():
 genes=['IL7R','GZMK','CXCL13','PDCD1','LAG3','ENTPD1','MKI67','HSPA1A']
 for k,c in enumerate(ALL):
  d=cells(c);fig,axs=grid(CN[c]+' CD8 transcript-expression densities',2,4,5.4);fig.subplots_adjust(hspace=.36,wspace=.18,bottom=.17)
  for j,g in enumerate(genes):
   ax=axs.flat[j];panel(ax,j,g);z=spatial(ax,d,d['log_'+g]);np.savez_compressed(OUT/'figure_source_data'/(f'FigS{k+9}_{g}_density.npz'),relative_density=z)
  cax=fig.add_axes([.79,.09,.15,.013]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,1),cmap=ROSE),cax=cax,orientation='horizontal');cb.set_ticks([0,1],labels=['Min','Max']);cb.outline.set_visible(False)
  note(fig,f'{len(d):,} cells; {d.patient.nunique()} independent patients/donors. Dashed contours: total-cell embedding density.')
  export(fig,f'FigS{k+9}',CN[c]+' gene-expression density atlas',8,'a–h, Expression-weighted spatial densities for eight transcript features on the unchanged within-cohort UMAP. Non-negative log(1 + CP10k) weights were binned on a 110 × 110 grid, smoothed with a Gaussian kernel (sigma 1.6 bins), and scaled to each panel maximum. Dashed contours show the total-cell density at 2.5% of its maximum. These displays locate transcript-enriched regions; panel intensities are not cross-gene abundance estimates. The full source cell population is used, without changing programme or antigen labels.',['frozen per-cell expression','within-cohort UMAP','additional full-count markers'])

def cell_heatmaps():
 genes=sum(DEFS.values(),[])
 for k,c in enumerate(ALL):
  d=cells(c);chosen=[]
  for _,u in d.groupby(['programme','patient','time'],sort=True):chosen.append(u.sort_values('cell_id').sample(n=min(30,len(u)),random_state=20261026))
  s=pd.concat(chosen);s['prog_order']=s.programme.map({p:i for i,p in enumerate(P)});s=s.sort_values(['prog_order','patient','time','cell_id']);raw=s[['log_'+g for g in genes]].to_numpy().T;z=zscore(raw);order=leaves_list(linkage(z,method='average',optimal_ordering=True))
  fig=plt.figure(figsize=(183/25.4,8.3));fig.text(.055,.97,CN[c]+' cell-resolved programme and receptor profiles',weight='bold',fontsize=10)
  ax=fig.add_axes([.17,.43,.73,.43]);panel(ax,0,f'Balanced display: {len(s):,} of {len(d):,} CD8 cells');ax.set_title(f'Balanced display: {len(s):,} of {len(d):,} CD8 cells',loc='left',pad=23)
  mat=pd.DataFrame(z[order],index=np.asarray(genes)[order],columns=s.cell_id);im=heat(ax,mat,xlabels=False);ax.tick_params(axis='y',labelsize=5.8)
  tax=fig.add_axes([.055,.43,.040,.43]);dendrogram(linkage(z,method='average',optimal_ordering=True),orientation='left',no_labels=True,ax=tax,color_threshold=0,above_threshold_color=INK);tax.invert_yaxis();tax.axis('off')
  for collection in tax.collections:collection.set_linewidth(.4)
  strips=fig.add_axes([.17,.867,.73,.013]);strips.imshow(np.array([[to_rgba(PC[p]) for p in s.programme]]),aspect='auto');strips.axis('off')
  mids=[];names=[];offset=0
  for p in P:
   n=int(s.programme.eq(p).sum())
   if n:ax.axvline(offset-.5,color='white',lw=.8);mids.append(offset+n/2-.5);names.append(p)
   offset+=n
  ax.set_xticks(mids,names,rotation=90,fontsize=6);cbar=fig.colorbar(im,cax=fig.add_axes([.925,.66,.012,.17]));cbar.set_ticks([-2,0,2]);cbar.ax.set_title('z',fontsize=6);cbar.outline.set_visible(False)
  s[['cell_id','patient','time','programme']].to_csv(OUT/'figure_source_data'/(f'FigS{k+13}_display_cells.csv'),index=False);mat.to_csv(OUT/'figure_source_data'/(f'FigS{k+13}_cell_z.csv.gz'))
  # Patient-equal within-programme target profile, requiring five cells per contributing patient/programme.
  parts=[]
  for (p,prog),u in d.groupby(['patient','programme']):
   if len(u)>=5:parts.append(pd.Series(u[['log_'+g for g in TARGETS]].mean().to_numpy(),index=TARGETS,name=(p,prog)))
  profiles=pd.DataFrame(parts);profiles.index=pd.MultiIndex.from_tuples(profiles.index,names=['patient','programme']);means=profiles.groupby('programme').mean().reindex(P).T;arr=zscore(means.to_numpy());ro=leaves_list(linkage(arr,method='average',optimal_ordering=True));w=pd.DataFrame(arr[ro],index=np.asarray(TARGETS)[ro],columns=P)
  ax2=fig.add_axes([.20,.09,.65,.23]);panel(ax2,1,'Patient-equal receptor profiles');im=heat(ax2,w);ax2.tick_params(labelsize=5.7);fig.colorbar(im,cax=fig.add_axes([.88,.15,.012,.12])).outline.set_visible(False);w.to_csv(OUT/'figure_source_data'/(f'FigS{k+13}_target_programme_z.csv'))
  note(fig,'Top: maximum 30 cells per patient/visit/programme; programme blocks retained. Bottom: equal patient weighting.')
  export(fig,f'FigS{k+13}',CN[c]+' cell-level expression heatmaps',2,'a, A deterministic balanced display subset contains at most 30 cells from each patient/visit/programme stratum; cells are grouped by programme, patient and visit. All 38 programme markers are shown as row z scores of log(1 + CP10k), capped at ±2. Average-linkage Euclidean clustering determines marker order. No new cell labels were learned. b, Mean target expression is first calculated within each patient/programme having at least five cells, then averaged equally over contributing patients. These target profiles are standardized across programmes for each target. Both panels are descriptive expression summaries; full population counts remain the denominators of the inferential analyses.',['full per-cell expression','programme definitions','balanced display cell IDs'])

def marker_distributions():
 genes=['CXCL13','IL7R','GZMK','PDCD1','LAG3','MKI67'];fig,axs=grid('Transcript distributions across receptor-independent programmes',4,6,9);fig.subplots_adjust(left=.075,right=.985,wspace=.55,hspace=.8,bottom=.11)
 rows=[]
 for i,c in enumerate(ALL):
  d=cells(c)
  for j,g in enumerate(genes):
   ax=axs[i,j];panel(ax,i*6+j,CN[c]+' · '+g);violin(ax,[d.loc[d.programme.eq(p),'log_'+g] for p in P],[PC[p] for p in P])
   means=d.groupby(['patient','programme'])['log_'+g].agg(['mean','size']).reset_index();means=means[means['size']>=5]
   for pi,p in enumerate(P):
    v=means[means.programme.eq(p)];ax.scatter(np.full(len(v),pi),v['mean'],s=3,c=INK,alpha=.45,zorder=4)
   rows.extend(means.assign(cohort=c,gene=g).to_dict('records'));ax.set_xticks(range(7),['Mem','Cyto','Dys','Cyc','IFN','Res','Mix'],rotation=90,fontsize=5.2)
   if j==0:ax.set_ylabel('log(1 + CP10k)')
 legend(fig,y=.047);note(fig,'Violin: individual cells. Black points: patient/programme means (at least five cells); no cell-level P values.')
 source(pd.DataFrame(rows),'FigS17_patient_gene_programme_means');export(fig,'FigS17','Programme-stratified transcript distributions',24,'a–x, Single-cell distributions of six transcript features across seven operational programmes in BCC, lung, oral and independent breast cohorts. Violin widths display density; black segments show median and interquartile range. Black points are patient/programme means with at least five contributing cells. Cell distributions are descriptive and do not supply independent patient replication. The programme definitions exclude all 24 target receptors; distributions of programme-defining genes illustrate the definitions rather than an independent validation.',['per-cell expression','programme labels','patient-level means'])

def paired_scores():
 d=old('programme_patient');ef=old('programme_effects');fig,axs=grid('Patient-paired transcriptional programme scores',3,6,7.5);fig.subplots_adjust(wspace=.75,hspace=.75,bottom=.10)
 for i,c in enumerate(CO):
  for j,p in enumerate(P[:-1]):
   ax=axs[i,j];panel(ax,i*6+j,CN[c]+' · '+p);u=d[d.cohort.eq(c)&d.programme.eq(p)];w=paired_lines(ax,u,'mean_score');ymin,ymax=ax.get_ylim();ax.set_ylim(ymin,ymax+(ymax-ymin)*.18);q=ef.query('cohort==@c and programme==@p and endpoint=="score"').iloc[0].q;ax.text(.02,.96,f'n={len(w)}; q={q:.2f}',transform=ax.transAxes,fontsize=5.8,va='top')
   if j==0:ax.set_ylabel('Mean programme score')
 legend(fig,'time',.045);note(fig,'Lines connect the same patient. BH q values refer to six paired score endpoints within each cohort.')
 export(fig,'FigS18','Patient-paired programme scores',18,'a–r, Patient-level mean programme scores before and after treatment. Lines connect paired samples, with each patient receiving equal inferential weight. Scores and exact two-sided sign-flip tests are the frozen original estimates. BH adjustment is over the six continuous-score endpoints within each cohort. No score endpoint meets within-cohort q < 0.05.',['programme_patient.csv','programme_effects.csv'])

def target_pairs():
 d=read('gene_patient_visit');ef=read('target_paired_effects')
 for part in range(2):
  genes=TARGETS[part*12:(part+1)*12];fig,axs=grid('Patient-paired receptor detection · panel '+str(part+1),4,3,9.0);fig.subplots_adjust(hspace=.80,wspace=.48,bottom=.12)
  for j,g in enumerate(genes):
   ax=axs.flat[j];panel(ax,j,g)
   for k,c in enumerate(CO):
    u=d[d.cohort.eq(c)&d.gene.eq(g)];w=u.pivot(index='patient',columns='time',values='detected_fraction').dropna()
    for p,row in w.iterrows():ax.plot([k-.15,k+.15],row[['pre','post']].to_numpy()*100,c=CC[c],lw=.45,alpha=.5);ax.scatter([k-.15,k+.15],row[['pre','post']].to_numpy()*100,c=[PRE,POST],s=6,zorder=3,edgecolors='none')
   ax.set(xticks=range(3),xticklabels=['BCC','Lung','Oral'],xlim=(-.45,2.45));ax.set_ylim(bottom=0)
   if j%3==0:ax.set_ylabel('RNA-detected cells (%)')
  legend(fig,'time',.055);note(fig,'Each segment is one patient. Target families include all 24 receptors; all paired BH q values exceed 0.05.')
  source(d[d.gene.isin(genes)],f'FigS{part+19}_patient_target_detection');export(fig,f'FigS{part+19}','Complete paired receptor detection family '+str(part+1),12,'a–l, Fractions of cells with at least one recorded UMI for each target, displayed separately for paired patient visits in the three treatment cohorts. Colors identify visit, and lines connect the same patient. The numerator is RNA detection and is not a protein-positive threshold. Exact paired sign-flip tests and patient bootstrap intervals are supplied in target_paired_effects.csv, with 24 targets comprising each cohort/endpoint BH family. No target meets q < 0.05. These additional target analyses were specified for the present exploratory expansion after inspection of the earlier study results.',['gene_patient_visit','target_paired_effects'])

def repertoire():
 d=read('repertoire_patient_visit');cl=read('clone_patient_visit');fig,axs=grid('Clonal expansion and sampled repertoire structure',3,4,8.2);fig.subplots_adjust(hspace=.7,wspace=.52)
 for i,c in enumerate(CO):
  u=d[d.cohort.eq(c)];v=cl[cl.cohort.eq(c)]
  ax=axs[i,0];panel(ax,i*4,CN[c]+' clone-size spectrum');cats=['1','2–4','5–19','≥20'];rows=[]
  for (p,t),z in v.groupby(['patient','time']):
   cat=pd.cut(z.n_cells,[0,1,4,19,np.inf],labels=cats);counts=z.groupby(cat,observed=False).n_cells.sum()/z.n_cells.sum();rows.extend({'patient':p,'time':t,'category':b,'fraction':counts[b]} for b in cats)
  z=pd.DataFrame(rows)
  for ti,t in enumerate(['pre','post']):
   means=z[z.time.eq(t)].groupby('category').fraction.mean().reindex(cats);ax.bar(np.arange(4)+(ti-.5)*.32,means*100,width=.30,color=[PRE,POST][ti])
  ax.set(xticks=range(4),xticklabels=cats,ylabel='Assigned cells (%)',xlabel='Cells per clone',ylim=(0,100))
  for j,(value,title) in enumerate([('expanded_cell_fraction','Cells in expanded clones'),('hill_D1','Shannon effective clones'),('tcr_coverage','TCR assignment coverage')],1):
   ax=axs[i,j];panel(ax,i*4+j,CN[c]+' · '+title);paired_lines(ax,u,value,scale=1 if value=='hill_D1' else 100);ax.set_ylabel('Effective clone number' if value=='hill_D1' else '%');ax.set_ylim(bottom=0)
 legend(fig,'time',.055);note(fig,'Clones are scoped to patient; oral P23 is excluded from paired TCR analyses. Fractions use assigned-cell denominators.')
 export(fig,'FigS21','Sampled clonal repertoire diversity and expansion',12,'a–l, Clone-size spectra and patient-paired sampled-repertoire summaries. Spectra average patient-level fractions equally. Expanded clones contain at least two observed cells. Shannon effective clone number is exp(−sum p log p), without a claim of equal repertoire sampling depth. TCR assignment coverage uses all retained CD8 cells as denominator; clone fractions use assigned CD8 cells. Oral P23 has no baseline TCR library and is excluded. Six paired repertoire endpoints per cohort, bootstrap intervals and BH q values are supplied in the source tables.',['repertoire_patient_visit','clone_patient_visit','repertoire_paired_effects'])

if __name__=='__main__':
 jobs={'densities':densities,'cell_heatmaps':cell_heatmaps,'distributions':marker_distributions,'scores':paired_scores,'targets':target_pairs,'repertoire':repertoire}
 for key in sys.argv[1:] or jobs:
  jobs[key]()
 atlas()
