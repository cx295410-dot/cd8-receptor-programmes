from largepanel_core import *
MAN=[]

def panel_details(stem):
 def item(identity,src,meaning):return {'identity':identity,'sources':[str(x) for x in src],'quantitative_meaning':meaning}
 out={}
 if stem=='Fig1':
  for j,c in enumerate(CO):
   out[string.ascii_lowercase[j]]=item(CN[c]+' CD8 programme UMAP',[ROOT/'data'/(c+'_cells.csv.gz')],'Within-cohort full-gene embedding; programme color, all retained cells')
   out[string.ascii_lowercase[3+j]]=item(CN[c]+' CXCL13 weighted density',[ROOT/'data'/(c+'_cells.csv.gz'),OUT/'figure_source_data'/(f'Fig1_{c}_CXCL13_density.npz')],'Gaussian-smoothed log-expression-weighted density, normalized to its panel maximum')
  out['g']=item('Patient programme composition',[ROOT/'results/programme_patient.csv',OUT/'figure_source_data/Fig1_patient_composition.csv'],'Within-patient/visit cell fractions, seven programmes')
 if stem=='Fig2':
  for j,c in enumerate(CO):out[string.ascii_lowercase[j]]=item(CN[c]+' all-marker programme profiles',[ROOT/'data'/(c+'_cells.csv.gz'),OUT/'figure_source_data'/('Fig2_'+c+'_marker_z.csv')],'All 38 markers; patient-equal means within programme; within-cohort row z score; shared gene clustering')
  out['d']=item('BCC source-state crosswalk',[REV/'results/author_state_crosswalk.csv',OUT/'figure_source_data/Fig2_BCC_source_crosswalk.csv'],'Within-source-state programme percentages')
  out['e']=item('BCC cell-resolved marker matrix',[OUT/'figure_source_data/Fig2_BCC_display_cells.csv',OUT/'figure_source_data/Fig2_BCC_display_cell_z.csv'],'Twelve markers, gene-wise z score over fixed balanced display cells')
  out['f']=item('Paired programme fraction changes',[ROOT/'results/programme_patient.csv'],'One point per patient; post minus pre Memory and Dysfunction fractions in percentage points')
  out['g']=item('BCC GZMK by programme',[ROOT/'data/BCC_cells.csv.gz',REV/'data/BCC_additional_markers.csv.gz',OUT/'figure_source_data/Fig2_BCC_GZMK_patient_programme.csv'],'Single-cell log-expression distribution and patient/programme mean points')
 if stem=='Fig3':
  for j,c in enumerate(CO):out[string.ascii_lowercase[j]]=item(CN[c]+' paired DESeq2 volcano',[ROOT/'results'/(c+'_paired_DESeq2.csv')],'Frozen patient-blocked post/pre log2 effect and gene BH-adjusted P value')
  out['d']=item('Complete Hallmark family',[ROOT/'references/MSigDB_Hallmark_2020.gmt']+[ROOT/'results'/(c+'_Hallmark_GSEA.csv') for c in CO],'All 50 GMT terms; frozen NES and GSEA q < 0.05 circles; missing eligible sets gray')
  for k,term in enumerate(['IFNG','E2F']):
   out[string.ascii_lowercase[4+k]]=item(term+' running enrichment',[ROOT/'results'/(c+'_GSEA_curves.csv.gz') for c in CO],'Frozen running ES over fraction of full Wald-ranked transcriptome')
   out[string.ascii_lowercase[6+k]]=item(term+' leading-edge pseudobulk expression',[OUT/'figure_source_data'/('Fig3_'+term+'_leadingedge_z.csv')]+[ROOT/'data'/(c+'_normalized_pseudobulk.csv.gz') for c in CO]+[ROOT/'data'/(c+'_pseudobulk_design.csv') for c in CO],'Fourteen recorded leading-edge genes; within-cohort log-normalized-count row z score; patient-adjacent visits')
 if stem=='Fig4':
  for j,c in enumerate(CO):
   out[string.ascii_lowercase[j]]=item(CN[c]+' all-target programme profiles',[ROOT/'data'/(c+'_cells.csv.gz'),OUT/'figure_source_data'/('Fig4_'+c+'_receptor_profile_z.csv')],'24 transcripts; patient-equal programme means; within-cohort row z score; shared gene clustering')
   out[string.ascii_lowercase[3+j]]=item(CN[c]+' measured target coexpression network',[REV/'results/target_coexpression_patient.csv',OUT/'figure_source_data'/('Fig4_'+c+'_coexpression_edges.csv')],'12 fixed display targets; median within-patient Spearman correlations, edges abs(rho) >= 0.15')
  out['g']=item('BCC PDCD1 weighted density',[ROOT/'data/BCC_cells.csv.gz',OUT/'figure_source_data/Fig4_BCC_PDCD1_density.npz'],'Panel-normalized log-expression-weighted spatial density')
  out['h']=item('Additional state transcript profiles',[REV/'data'/(c+'_additional_markers.csv.gz') for c in CO]+[OUT/'figure_source_data/Fig4_additional_state_gene_z.csv'],'Six non-defining transcript profiles; patient-equal programme mean, within-cohort row z score')
 return out

def reg(fig,stem,title,n,caption,sources):
 rec=export(fig,stem,title,list(string.ascii_lowercase[:n]),caption,sources);rec['panel_details']=panel_details(stem);(OUT/'work'/(stem+'_metadata.json')).write_text(json.dumps(rec,indent=2),encoding='utf-8');MAN.append(rec)
def newfig(height=8):return plt.figure(figsize=(183/25.4,height))
def balanced(d,genes):
 sizes=d.groupby(['patient','programme']).size();m=d.groupby(['patient','programme'])[['log_'+g for g in genes]].mean();m=m.loc[sizes.ge(5)].groupby('programme').mean().reindex(P).T;m.index=genes;return m

def clustered(ax,m,labels=True,cluster=True,v=2):
 z=zscore(m);order=leaves_list(linkage(np.nan_to_num(z),method='average',optimal_ordering=True)) if cluster else np.arange(len(z));w=pd.DataFrame(z[order],index=m.index[order],columns=m.columns);im=heat(ax,w,v=v,ylabels=labels);ax.tick_params(axis='x',labelsize=6.5);ax.tick_params(axis='y',labelsize=6.6)
 if cluster:
  pos=ax.get_position();da=ax.figure.add_axes([pos.x0-.016,pos.y0,.013,pos.height]);dendrogram(linkage(np.nan_to_num(z),method='average',optimal_ordering=True),orientation='left',no_labels=True,ax=da,color_threshold=0,above_threshold_color=INK);da.invert_yaxis();da.axis('off')
  for coll in da.collections:coll.set_linewidth(.35)
 return im,w

def f1():
 fig=newfig(7.35);src=[]
 for j,c in enumerate(CO):
  d=cells(c);ax=fig.add_axes([.065+j*.317,.735,.27,.215]);panel(ax,j,CN[c]+f' CD8 cells ({len(d):,})');umap(ax,d)
  ax.text(.01,.01,f'{d.patient.nunique()} patients',transform=ax.transAxes,fontsize=6.9)
  ax=fig.add_axes([.065+j*.317,.445,.27,.215]);panel(ax,3+j,CN[c]+' CXCL13 density');z=spatial(ax,d,d.log_CXCL13);np.savez_compressed(OUT/'figure_source_data'/(f'Fig1_{c}_CXCL13_density.npz'),relative_density=z);src.append(str(ROOT/'data'/(c+'_cells.csv.gz')))
 cbar(fig,plt.cm.ScalarMappable(norm=Normalize(0,1),cmap=ROSE),[.81,.401,.13,.009],'Relative density',[0,1],'horizontal')
 pr=old('programme_patient');pieces=[];cohort_ids=[]
 for c in CO:
  w=pr[pr.cohort.eq(c)].pivot(index=['patient','time'],columns='programme',values='fraction').reindex(columns=P);w=w.reindex(sorted(w.index,key=lambda x:(x[0],x[1]=='post')));pieces.append(w);cohort_ids.extend([c]*len(w))
 w=pd.concat(pieces);ax=fig.add_axes([.09,.145,.86,.165]);panel(ax,6,'Patient-resolved programme composition');base=np.zeros(len(w))
 for p in P:ax.bar(range(len(w)),w[p],bottom=base,color=PC[p],width=.94,lw=0);base+=w[p].to_numpy()
 off=0
 for c,piece in zip(CO,pieces):
  if off:ax.axvline(off-.5,c='white',lw=2.5)
  ax.text((off+len(piece)/2)/len(w),1.03,CN[c],transform=ax.transAxes,ha='center',fontsize=7.2);off+=len(piece)
 ax.set(ylim=(0,1),xlim=(-.6,len(w)-.4),yticks=[0,.5,1],ylabel='Cell fraction');ax.set_xticks(range(len(w)),[p.replace('su','')+('P' if t=='pre' else 'T') for p,t in w.index],rotation=90,fontsize=6.1);legend(fig,y=.032)
 src.append(source(w.reset_index().assign(cohort=cohort_ids),'Fig1_patient_composition'))
 reg(fig,'Fig1','Cohort-resolved CD8 transcriptional atlas',7,'a–c, Within-cohort UMAPs of retained author-defined CD8 cells from BCC, lung and oral treatment cohorts, colored by the seven receptor-independent operational programmes. Embeddings were calculated from the full gene-expression matrices and are shown separately for each cohort. d–f, CXCL13 expression-weighted densities on the same embeddings. Log(1 + CP10k) weights were binned on 110 × 110 grids, Gaussian-smoothed (sigma 1.6 bins) and scaled to each panel maximum; dashed contours indicate total-cell density at 2.5% of its maximum. Density intensity is relative within panel. g, Programme fractions among all retained cells in each paired patient/visit. Adjacent P and T bars indicate pretreatment and post-treatment samples; cohort boundaries are separated in white. Programme colors are shared across panels. Source cell tables supply cell IDs, original annotations, patient identity, programme labels and UMAP coordinates.',src)

def f2():
 fig=newfig(8.6);genes=sum(DEFS.values(),[]);src=[];mats={c:balanced(cells(c),genes) for c in CO};stack=np.hstack([zscore(mats[c]) for c in CO]);order=leaves_list(linkage(stack,method='average',optimal_ordering=True));ordered=np.asarray(genes)[order]
 for j,c in enumerate(CO):
  ax=fig.add_axes([.145+j*.287,.565,.238,.385]);panel(ax,j,CN[c]+' marker programmes');m=mats[c];w=pd.DataFrame(zscore(m)[order],index=ordered,columns=AB);im=heat(ax,w,ylabels=j==0);ax.tick_params(axis='y',labelsize=6.4);ax.tick_params(axis='x',labelsize=6.5);src.append(source(w.reset_index(names='gene'),'Fig2_'+c+'_marker_z'))
 da=fig.add_axes([.045,.565,.027,.385]);dendrogram(linkage(stack,method='average',optimal_ordering=True),orientation='left',no_labels=True,ax=da,color_threshold=0,above_threshold_color=INK);da.invert_yaxis();da.axis('off')
 for coll in da.collections:coll.set_linewidth(.35)
 cbar(fig,im,[.971,.69,.009,.17],'z',[-2,0,2])
 cw=read('author_state_crosswalk');u=cw[cw.cohort.eq('BCC')];m=u.pivot(index='author_state',columns='programme',values='fraction').reindex(columns=P);m.columns=AB;ax=fig.add_axes([.145,.305,.335,.16]);panel(ax,3,'Original BCC states to programmes');im=heat(ax,m,v=1,cmap=ROSE);ax.tick_params(axis='y',labelsize=6.7);ax.tick_params(axis='x',labelsize=6.4)
 for i in range(len(m)):
  for j in range(7):ax.text(j,i,f'{100*m.iloc[i,j]:.0f}',fontsize=6.2,ha='center',va='center')
 cbar(fig,im,[.493,.34,.008,.08],'Fraction',[0,1]);src.append(source(m.reset_index(),'Fig2_BCC_source_crosswalk'))
 d=cells('BCC');sample=pd.concat([u.sort_values('cell_id').sample(n=min(12,len(u)),random_state=20261026) for _,u in d.groupby(['patient','time','programme'],sort=True)]).copy();sample['po']=sample.programme.map({p:i for i,p in enumerate(P)});sample=sample.sort_values(['po','patient','time','cell_id']);sel=['TCF7','CCR7','GZMB','PRF1','CXCL13','TOX','MKI67','TOP2A','ISG15','IFIT3','ZNF683','CD69'];mat=pd.DataFrame(zscore(sample[['log_'+g for g in sel]].to_numpy().T),index=sel,columns=sample.cell_id);ax=fig.add_axes([.66,.305,.29,.16]);panel(ax,4,'BCC cell-resolved marker identity');heat(ax,mat,xlabels=False);ax.tick_params(axis='y',labelsize=6.4);strip=fig.add_axes([.66,.296,.29,.006]);strip.imshow(np.array([[to_rgba(PC[p]) for p in sample.programme]]),aspect='auto');strip.axis('off');src.append(source(mat.reset_index(names='gene'),'Fig2_BCC_display_cell_z'));src.append(source(sample[['cell_id','patient','time','programme']],'Fig2_BCC_display_cells'))
 ax=fig.add_axes([.145,.092,.335,.13]);panel(ax,5,'Patient-paired programme shifts');pr=old('programme_patient')
 for c in CO:
  u=pr[pr.cohort.eq(c)];a=u[u.programme.eq('Memory')].pivot(index='patient',columns='time',values='fraction');b=u[u.programme.eq('Dysfunction')].pivot(index='patient',columns='time',values='fraction');ax.scatter((a.post-a.pre)*100,(b.post-b.pre)*100,s=19,c=CC[c],ec='white',lw=.4,label=CN[c])
 ax.axhline(0,c=GRAY,lw=.5);ax.axvline(0,c=GRAY,lw=.5);ax.set(xlabel='Memory change (pp)',ylabel='Dysfunction change (pp)');ax.legend(fontsize=6.2,ncol=1,loc='center left',bbox_to_anchor=(1.01,.5),handletextpad=.2,columnspacing=.6)
 ax=fig.add_axes([.66,.092,.29,.13]);panel(ax,6,'BCC GZMK distribution by programme');violin(ax,[d.loc[d.programme.eq(p),'log_GZMK'] for p in P],[PC[p] for p in P]);u=d.groupby(['patient','programme']).log_GZMK.agg(['mean','size']);u=u[u['size']>=5].reset_index()
 for j,p in enumerate(P):ax.scatter(np.full(sum(u.programme.eq(p)),j),u.loc[u.programme.eq(p),'mean'],s=4,c=INK,alpha=.45,zorder=4)
 ax.set(xticks=range(7),xticklabels=AB,ylabel='log(1 + CP10k)');ax.tick_params(axis='x',rotation=90,labelsize=6.5);src.append(source(u,'Fig2_BCC_GZMK_patient_programme'));legend(fig,y=.013)
 reg(fig,'Fig2','Receptor-independent programme identity and patient heterogeneity',7,'a–c, All 38 designated programme markers in BCC, lung and oral CD8 cells. Mean log(1 + CP10k) expression is calculated within patient/programme groups with at least five cells, then averaged equally across patients. Rows are standardized within each cohort and share an average-linkage Euclidean gene ordering learned from the three cohort profiles; color is capped at ±2. d, Original BCC source annotations crossed with operational programmes. Entries are within-source-state percentages before rounding, and each row sums to 100%. e, Representative BCC single-cell marker matrix, with at most 12 cells per patient/visit/programme stratum (seed 20261026), grouped by programme and shown as gene-wise z scores. The color strip gives programme identity; the supplied cell-ID table records the display subset. f, Patient-paired changes in Memory and Dysfunction cell fractions; each point is one patient and values are percentage points. g, Single-cell GZMK distributions across BCC programmes. Black intervals mark medians and interquartile ranges; black points are patient/programme means with at least five cells. GZMK is an additional transcript feature outside the 38 defining markers. The programme marker matrix describes programme identity, while f and g show patient heterogeneity and an independent transcript feature.',src+[str(ROOT/'results/programme_patient.csv'),str(REV/'data/BCC_additional_markers.csv.gz')])

def canonical(s):
 import re
 s=re.sub(r'[^A-Z0-9]+','_',s.upper()).strip('_')
 for a,b in [('TNF_ALPHA','TNFA'),('NF_KB','NFKB'),('IL_2','IL2'),('IL_6','IL6'),('G2_M','G2M')]:s=s.replace(a,b)
 return s if s.startswith('HALLMARK_') else 'HALLMARK_'+s

def leading_matrix(term,limit=14):
 lists={};score={};blocks=[];names=[];bounds=[]
 for c in CO:
  g=old(c+'_Hallmark_GSEA').set_index('Term');lists[c]=str(g.loc[term,'Lead_genes']).split(';');rank=pd.read_csv(ROOT/'results'/(c+'_Wald_ranking.csv'),index_col=0).iloc[:,0]
  for k in lists[c]:score[k]=max(score.get(k,0),abs(rank.get(k,0)))
 genes=sorted(score,key=score.get,reverse=True)[:limit]
 for c in CO:
  x=pd.read_csv(ROOT/'data'/(c+'_normalized_pseudobulk.csv.gz'),index_col=0);de=pd.read_csv(ROOT/'data'/(c+'_pseudobulk_design.csv')).set_index('sample');de=de.assign(torder=de.condition.map({'pre':0,'post':1})).sort_values(['patient','torder']);x=np.log1p(x.reindex(index=de.index,columns=genes));blocks.append(zscore(x.to_numpy().T));names.extend([c+':'+str(i) for i in de.index]);bounds.append(len(de))
 return pd.DataFrame(np.hstack(blocks),index=genes,columns=names),bounds,lists

def f3():
 fig=newfig(8.25);src=[]
 for j,c in enumerate(CO):
  ax=fig.add_axes([.09+j*.312,.785,.255,.155]);panel(ax,j,CN[c]+' paired transcriptome');u=pd.read_csv(ROOT/'results'/(c+'_paired_DESeq2.csv'),index_col=0);q=u.padj.fillna(1).clip(lower=1e-200);sig=q<.05;ax.scatter(u.log2FoldChange,-np.log10(q),c=np.where(sig&u.log2FoldChange.gt(0),POST,np.where(sig,PRE,'#CED1D4')),s=2.5,alpha=.75,lw=0,rasterized=True);ax.axhline(-np.log10(.05),c=GRAY,lw=.55,ls='--');ax.axvline(0,c='#DDDDDD',lw=.5);ax.set(xlabel='Post / pre log₂ fold change',ylabel='−log₁₀ adjusted P' if j==0 else '');ax.text(.02,.96,f'{sig.sum():,} FDR < 0.05',fontsize=6.5,transform=ax.transAxes,va='top');src.append(str(ROOT/'results'/(c+'_paired_DESeq2.csv')))
 gs=pd.concat([old(c+'_Hallmark_GSEA').assign(cohort=c) for c in CO]);terms=sorted(canonical(s.split('\t')[0]) for s in (ROOT/'references/MSigDB_Hallmark_2020.gmt').read_text().splitlines());w=gs.pivot(index='Term',columns='cohort',values='NES').reindex(index=terms,columns=CO).astype(float);q=gs.pivot(index='Term',columns='cohort',values='FDR q-val').reindex(index=terms,columns=CO).astype(float);w.columns=[CN[c] for c in CO];w.index=[t.replace('HALLMARK_','').replace('_',' ').capitalize().replace('Il2','IL2').replace('Il6','IL6').replace('stat3','STAT3').replace('stat5','STAT5').replace('G2m','G2M').replace('Myc','MYC').replace('Mtorc1','MTORC1').replace('Tnfa','TNFA').replace('nfkb','NFKB').replace('E2f','E2F').replace('Kras','KRAS').replace('P53','p53').replace('Pi3k','PI3K').replace('Tgf','TGF').replace('Uv','UV').replace('Wnt','WNT').replace('Notch','NOTCH') for t in terms];ax=fig.add_axes([.325,.105,.20,.595]);panel(ax,3,'Full Hallmark family');im=heat(ax,w,v=3);ax.tick_params(axis='y',labelsize=6.3);ax.tick_params(axis='x',labelsize=7)
 for y,x in np.argwhere(q.to_numpy()<.05):ax.plot(x,y,'o',ms=2.8,mfc='none',mec=INK,mew=.55)
 cb=cbar(fig,im,[.33,.051,.19,.009],'',[-3,0,3],'horizontal');fig.text(.31,.0555,'NES',fontsize=6.6,ha='right',va='center');src.append(source(gs,'Fig3_full_Hallmark_results'))
 for k,(term,ttl) in enumerate([('HALLMARK_INTERFERON_GAMMA_RESPONSE','Interferon-γ response'),('HALLMARK_E2F_TARGETS','E2F targets')]):
  ax=fig.add_axes([.705,.595-k*.157,.25,.085]);panel(ax,4+k,ttl)
  for c in CO:
   u=pd.read_csv(ROOT/'results'/(c+'_GSEA_curves.csv.gz'));u=u[u.term.eq(term)];ax.plot(u['rank']/max(u['rank'].max(),1),u.running_ES,c=CC[c],lw=1,label=CN[c])
  ax.axhline(0,c=GRAY,lw=.5);ax.set(xlim=(0,1),ylabel='Running ES');ax.tick_params(labelsize=6.5)
  if k==0:fig.legend([Line2D([],[],c=CC[c],lw=1.2) for c in CO],[CN[c] for c in CO],fontsize=6.6,ncol=3,loc='center',bbox_to_anchor=(.826,.725),columnspacing=.8,handlelength=1.1)
  else:ax.set_xlabel('')
  mat,bounds,lists=leading_matrix(term);ax=fig.add_axes([.705,.255-k*.185,.25,.135]);panel(ax,6+k,('IFN-γ' if k==0 else 'E2F')+' leading-edge genes');heat(ax,mat,xlabels=False);ax.tick_params(axis='y',labelsize=6.3);off=0
  for c,n in zip(CO,bounds):
   if off:ax.axvline(off-.5,c='white',lw=1.6)
   ax.text((off+n/2)/sum(bounds),-.055,CN[c],transform=ax.transAxes,ha='center',va='top',fontsize=6.5);off+=n
  src.append(source(mat.reset_index(names='gene'),'Fig3_'+('IFNG' if k==0 else 'E2F')+'_leadingedge_z'));(OUT/'figure_source_data'/('Fig3_'+str(k)+'_recorded_leadingedge.json')).write_text(json.dumps(lists,indent=2))
 note(fig,'Heatmap circles: GSEA q < 0.05. Gray: pathway ineligible. Expression heatmaps: within-cohort gene z scores (±2).')
 reg(fig,'Fig3','Patient-paired transcriptional and pathway remodelling',8,'a–c, Frozen patient-blocked DESeq2 post/pre effects in BCC, lung and oral pseudobulks. Red and blue indicate positive and negative effects at gene-level BH FDR < 0.05; gray indicates other tested genes. Dashed horizontal lines mark adjusted P = 0.05. d, All 50 Hallmark sets from the recorded GMT, with normalized enrichment scores from the original full Wald-statistic rankings. Circles indicate within-cohort GSEA q < 0.05; gray entries were ineligible under the original gene-set size filters. e,f, Recorded running enrichment scores for interferon-γ response and E2F targets, plotted over the full ranked gene lists in each cohort. g,h, Patient-paired normalized pseudobulk profiles of 14 recorded leading-edge genes per pathway. Display genes are the union of the three recorded leading-edge lists, ordered by their maximum absolute Wald statistic and truncated to 14; this display rule does not change enrichment results. Columns retain patient-adjacent pretreatment/post-treatment ordering within each cohort; white dividers separate cohorts. Expression is log(1 + normalized count), standardized across samples within each cohort and gene. All differential-expression and enrichment estimates are unchanged.',src+[str(ROOT/'results'/(c+'_GSEA_curves.csv.gz')) for c in CO]+[str(ROOT/'data'/(c+'_normalized_pseudobulk.csv.gz')) for c in CO])

def f4():
 fig=newfig(8.1);src=[];mats={c:balanced(cells(c),TARGETS) for c in CO};stack=np.hstack([zscore(mats[c]) for c in CO]);order=leaves_list(linkage(stack,method='average',optimal_ordering=True))
 for j,c in enumerate(CO):
  m=mats[c];w=pd.DataFrame(zscore(m)[order],index=np.asarray(TARGETS)[order],columns=AB);ax=fig.add_axes([.15+j*.285,.61,.232,.33]);panel(ax,j,CN[c]+' receptor by programme');im=heat(ax,w,ylabels=j==0);ax.tick_params(axis='y',labelsize=6.6);ax.set_xticklabels(AB,rotation=90,fontsize=6.6);src.append(source(w.reset_index(names='target'),'Fig4_'+c+'_receptor_profile_z'))
 da=fig.add_axes([.035,.61,.027,.33]);dendrogram(linkage(stack,method='average',optimal_ordering=True),orientation='left',no_labels=True,ax=da,color_threshold=0,above_threshold_color=INK);da.invert_yaxis();da.axis('off')
 for coll in da.collections:coll.set_linewidth(.35)
 cbar(fig,im,[.968,.69,.01,.17],'z',[-2,0,2]);co=read('target_coexpression_patient');genes=['IL7R','PDCD1','LAG3','HAVCR2','TIGIT','CTLA4','ENTPD1','TNFRSF9','ICOS','CD27','CD28','CXCR6'];theta=np.linspace(0,2*np.pi,len(genes),endpoint=False)+np.pi/2;xy=np.c_[np.cos(theta),np.sin(theta)]
 for j,c in enumerate(CO):
  ax=fig.add_axes([.075+j*.31,.315,.265,.205]);panel(ax,3+j,CN[c]+' target coexpression');u=co[co.cohort.eq(c)].groupby(['gene1','gene2']).rho.median().unstack().reindex(index=genes,columns=genes);edges=[]
  for a in range(len(genes)):
   for b in range(a):
    rho=u.iloc[a,b]
    if np.isfinite(rho) and abs(rho)>=.15:ax.plot(xy[[a,b],0],xy[[a,b],1],c=POST if rho>0 else PRE,alpha=.4,lw=.4+abs(rho)*4);edges.append({'cohort':c,'gene1':genes[a],'gene2':genes[b],'median_within_patient_rho':rho})
  for a,g in enumerate(genes):
   ax.scatter(*xy[a],s=31,c='#777777',ec='white',lw=.5,zorder=3)
   xpos,ypos=xy[a];ha='center'
   if xpos<-.1:xpos-=.10;ha='right'
   elif xpos>.1:xpos+=.10;ha='left'
   else:ypos+=.14*np.sign(ypos)
   ax.text(xpos,ypos,g,ha=ha,va='center',fontsize=6.5)
  ax.set(xlim=(-1.55,1.55),ylim=(-1.34,1.38));ax.axis('off');src.append(source(pd.DataFrame(edges,columns=['cohort','gene1','gene2','median_within_patient_rho']),'Fig4_'+c+'_coexpression_edges'))
 ax=fig.add_axes([.075,.067,.265,.16]);panel(ax,6,'BCC PDCD1 expression density');z=spatial(ax,cells('BCC'),cells('BCC').log_PDCD1);np.savez_compressed(OUT/'figure_source_data/Fig4_BCC_PDCD1_density.npz',relative_density=z)
 genes2=['GZMK','GZMA','CCL5','IFNG','CXCR3','CXCR4'];blocks=[]
 for c in CO:
  m=balanced(cells(c),genes2);blocks.append(zscore(m))
 w=pd.DataFrame(np.hstack(blocks),index=genes2,columns=AB*3);ax=fig.add_axes([.50,.067,.45,.16]);panel(ax,7,'Additional state-associated transcript profiles');im=heat(ax,w);ax.tick_params(axis='x',labelsize=6.3);ax.tick_params(axis='y',labelsize=7)
 for x in [6.5,13.5]:ax.axvline(x,c='white',lw=2)
 for j,c in enumerate(CO):ax.text((j+.5)/3,1.02,CN[c],transform=ax.transAxes,ha='center',fontsize=7)
 src.append(source(w.reset_index(names='gene'),'Fig4_additional_state_gene_z'))
 reg(fig,'Fig4','Patient-balanced receptor and transcriptional state context',8,'a–c, All 24 receptor transcripts across the seven operational programmes in BCC, lung and oral cells. Within-patient/programme means require at least five cells and receive equal patient weighting. Gene profiles are standardized across programmes. Joint average-linkage Euclidean clustering of the three cohort profiles gives the shared row order; the thin dendrogram shows this ordering and colors are capped at ±2. d–f, Observed coexpression networks for a fixed 12-target display set: IL7R, PDCD1, LAG3, HAVCR2, TIGIT, CTLA4, ENTPD1, TNFRSF9, ICOS, CD27, CD28 and CXCR6. An edge is drawn when the absolute median within-patient single-cell Spearman correlation is at least 0.15. Edge width is proportional to its absolute correlation; red and blue mark positive and negative associations. All display nodes are retained, including isolated targets. These are transcript coexpression associations rather than regulatory or physical interactions. g, PDCD1 expression-weighted density on the unchanged BCC embedding, computed and scaled as in Fig. 1d–f. h, Six additional state-associated transcripts outside the 38 defining markers, summarized with the same patient-equal procedure as a–c and standardized within each cohort. The three seven-column blocks follow the shared programme order. Source data record every displayed profile and network edge.',src+[str(REV/'results/target_coexpression_patient.csv')]+[str(ROOT/'data'/(c+'_cells.csv.gz')) for c in CO]+[str(REV/'data'/(c+'_additional_markers.csv.gz')) for c in CO])

if __name__=='__main__':
 for p in ['figures','figure_source_data','scripts','work']:(OUT/p).mkdir(parents=True,exist_ok=True)
 jobs={'f1':f1,'f2':f2,'f3':f3,'f4':f4}
 for key in sys.argv[1:] or jobs:jobs[key]()
 # Independent manifests never alter the other worker's figures.
 existing=json.loads((OUT/'manifest_main_1_4.json').read_text()) if (OUT/'manifest_main_1_4.json').exists() else []
 existing=[r for r in existing if r['figure'] not in {r['figure'] for r in MAN}]+MAN;existing=sorted(existing,key=lambda r:int(r['figure'][3:]));(OUT/'manifest_main_1_4.json').write_text(json.dumps(existing,indent=2),encoding='utf-8');(OUT/'legends_main_1_4.json').write_text(json.dumps({r['figure']:r['caption'] for r in existing},indent=2),encoding='utf-8')

