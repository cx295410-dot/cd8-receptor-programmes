"""Source-derived main figures 5–8; no result or antigen label is inferred for display."""
from largepanel_core import *
import shutil

MAN=[];CAP={}
def axat(f,r,i,title):
 a=f.add_axes(r);panel(a,i,title);return a
def record(f,stem,title,n,caption,sources):
 r=export(f,stem,title,n,caption,sources);MAN.append(r);CAP[stem]=caption
def rowheat(ax,m,v=2,cluster=True,xlabels=True):
 m=pd.DataFrame(m);z=pd.DataFrame(zscore(m),index=m.index,columns=m.columns)
 if cluster and len(z)>1:z=z.iloc[leaves_list(linkage(np.nan_to_num(z),method='average',optimal_ordering=True))]
 return heat(ax,z,v=v,xlabels=xlabels),z
def densitykey(f,rect):
 cbar(f,plt.cm.ScalarMappable(norm=Normalize(0,1),cmap=ROSE),rect,'Relative expression density',[0,1],orientation='horizontal')
def flow(ax,mat):
 a=mat.to_numpy();gap=.025;s=1-gap*6;l=a.sum(1)*s;r=a.sum(0)*s;yl=np.cumsum(np.r_[0,l[:-1]])+np.arange(7)*gap;yr=np.cumsum(np.r_[0,r[:-1]])+np.arange(7)*gap;ol=yl.copy();orr=yr.copy()
 for i in range(7):
  for j in range(7):
   w=a[i,j]*s
   if w<=0:continue
   y0,y1=ol[i],orr[j];verts=[(.07,y0),(.43,y0),(.57,y1),(.93,y1),(.93,y1+w),(.57,y1+w),(.43,y0+w),(.07,y0+w),(.07,y0)];codes=[MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.LINETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.CLOSEPOLY];ax.add_patch(PathPatch(MPath(verts,codes),fc=PC[P[i]],ec='none',alpha=.42));ol[i]+=w;orr[j]+=w
 for i,p in enumerate(P):
  ax.add_patch(Rectangle((.035,yl[i]),.035,l[i],fc=PC[p],lw=0));ax.add_patch(Rectangle((.93,yr[i]),.035,r[i],fc=PC[p],lw=0))
 ax.set(xlim=(0,1),ylim=(1,0));ax.axis('off');ax.text(.035,-.09,'Pre',transform=ax.transAxes,fontsize=7);ax.text(.85,-.09,'Post',transform=ax.transAxes,fontsize=7)

def f5():
 f=plt.figure(figsize=(183/25.4,9.0));fl=old('shared_clone_links');cl=read('clone_patient_visit');sh=read('programme_clone_sharing');ss=[]
 for j,c in enumerate(CO):
  u=fl[fl.cohort.eq(c)];m=u.groupby(['pre_programme','post_programme']).patient_fraction.sum().unstack(fill_value=0).reindex(index=P,columns=P,fill_value=0)/u.patient.nunique();assert np.isclose(m.to_numpy().sum(),1)
  a=axat(f,[.08+j*.31,.724,.255,.208],j,CN[c]+f' shared clones · n={u.patient.nunique()}');flow(a,m);ss.append(source(m.reset_index(),'Fig5_'+c+'_flow'))
  u=sh[sh.cohort.eq(c)&sh.time.eq('post')];m=u.groupby(['programme1','programme2']).jaccard.median().unstack().reindex(index=P,columns=P)
  a=axat(f,[.10+j*.31,.421,.215,.195],j+3,CN[c]+' programme sharing');im=heat(a,m,v=1,cmap=plt.get_cmap('Blues'));a.set_xticklabels(AB,rotation=60,ha='right');a.set_yticklabels(AB);ss.append(source(m.reset_index(),'Fig5_'+c+'_sharing'))
 cbar(f,im,[.95,.445,.008,.13],'Jaccard',[0,.5,1])
 # A fixed number per patient prevents large repertoires dominating expression display.
 genes=['IL7R','TCF7','CXCL13','PDCD1','LAG3','ENTPD1','GZMB','NKG7','TOX','MKI67','ZNF683','HSPA1A'];blocks=[];meta=[]
 for c in CO:
  u=cl[cl.cohort.eq(c)];sel=u.groupby(['patient','scoped_clone']).frequency.max().reset_index().sort_values(['patient','frequency'],ascending=[True,False]).groupby('patient').head(1)
  for r in sel.itertuples():
   v=u[u.scoped_clone.eq(r.scoped_clone)].set_index('time').reindex(['pre','post']);blocks.extend([v.loc[t,genes].to_numpy(float) for t in ['pre','post']]);meta.extend([{'cohort':c,'patient':r.patient,'scoped_clone':r.scoped_clone,'time':t} for t in ['pre','post']])
 m=pd.DataFrame(np.asarray(blocks).T,index=genes,columns=[f'{r["cohort"]}:{r["patient"]}:{r["time"]}' for r in meta]);a=axat(f,[.105,.112,.48,.22],6,'Dominant-clone expression');im,z=rowheat(a,m,xlabels=False);a.set_xlabel('One dominant clone per patient; pre / post');
 for k in range(2,m.shape[1],2):a.axvline(k-.5,c='white',lw=.4)
 ss += [source(z.reset_index(),'Fig5_dominant_clone_z'),source(pd.DataFrame(meta),'Fig5_dominant_clone_selection')];cbar(f,im,[.60,.145,.009,.14],'Row z',[-2,0,2])
 a=axat(f,[.74,.112,.215,.22],7,'Clone-size distribution');rows=[];cols=['#DCEED0','#91D650','#42BBA6','#556C8D']
 for ci,c in enumerate(CO):
  for ti,t in enumerate(['pre','post']):
   vals=[]
   for p,u in cl[cl.cohort.eq(c)&cl.time.eq(t)].groupby('patient'):
    cats=pd.cut(u.n_cells,[0,1,4,19,np.inf],labels=False);v=np.array([u.loc[cats.eq(k),'n_cells'].sum()/u.n_cells.sum() for k in range(4)]);vals.append(v);rows.append({'cohort':c,'patient':p,'time':t,**dict(zip(['one','two_four','five_nineteen','twenty_plus'],v))})
   mean=np.mean(vals,axis=0);bottom=0
   for k,v in enumerate(mean):a.bar(ci*2.5+ti,v*100,bottom=bottom,color=cols[k],width=.65);bottom+=v*100
 a.set(xticks=[.5,3,5.5],xticklabels=['BCC','Lung','Oral'],ylim=(0,100),ylabel='Assigned CD8 cells (%)');a.text(.5,-.20,'Within each cohort: pre, post',ha='center',transform=a.transAxes,fontsize=7);f.legend([Rectangle((0,0),1,1,fc=v) for v in cols],['1','2–4','5–19','≥20 cells / clone'],ncol=4,loc='lower center',bbox_to_anchor=(.75,.018),fontsize=6.8);ss.append(source(pd.DataFrame(rows),'Fig5_repertoire_size_patient'))
 legend(f,y=.031)
 caption='a–c, Pre/post programme occupancy of patient-scoped clonotypes detected at both visits. Dominant programme labels and ribbon weights retain the original shared-clone analysis; each patient contributes equally. Ribbon width sums to one within cohort and indicates observed sharing, not a differentiation trajectory. d–f, Median post-treatment within-patient Jaccard overlap between programme clonotype sets. A clonotype enters a programme set with at least two observed cells; empty unions remain missing. g, Mean log(1 + CP10k) expression of the dominant clonotype from each patient, selected by its largest within-patient frequency across either visit. Each patient supplies one clone and adjacent pre/post columns. Gene rows are standardized across observed profiles and clustered by average-linkage Euclidean distance; gray denotes an unobserved clone/visit. h, Clone-size composition, averaged equally over patients, with pre and post bars within each cohort. Oral P23 has no baseline TCR library and is excluded from paired repertoire analysis. All clone identities are scoped to patient.'
 record(f,'Fig5','Clonal sharing connects transcriptional programmes',8,caption,ss+[str(ROOT/'results/shared_clone_links.csv'),str(REV/'results/clone_patient_visit.csv'),str(REV/'results/programme_clone_sharing.csv')])

def f6():
 d=cells('BREAST');f=plt.figure(figsize=(183/25.4,9));ss=[]
 panels=[(None,'Breast CD8 programmes'),('log_IL7R','IL7R RNA density'),('IL7R_adt_corrected','CD127 protein density'),('log_PDCD1','PDCD1 RNA density'),('PDCD1_adt_corrected','PD-1 protein density')]
 for i,(g,t) in enumerate(panels):
  a=axat(f,[.075+(i%3)*.31,.734-(i//3)*.29,.26,.21],i,t)
  if g:spatial(a,d,d[g])
  else:umap(a,d)
 cor=pd.read_csv(ROOT.parent/'results/protein_rna_correlations.csv');cor=cor[cor.grouping.eq('donor')];genes=cor.groupby('target').rna_adt_rho.median().sort_values().index.tolist();m=cor.pivot(index='target',columns='donor',values='rna_adt_rho').reindex(genes)
 a=axat(f,[.748,.427,.195,.224],5,'RNA / ADT concordance');im=heat(a,m,v=.6);a.tick_params(labelsize=6.8);a.set_xticklabels(m.columns,rotation=60,ha='right');cbar(f,im,[.956,.455,.008,.15],'Spearman ρ',[-.6,0,.6]);ss.append(source(m.reset_index(),'Fig6_donor_RNA_ADT_rho'))
 dif=pd.read_csv(ROOT.parent/'results/protein_paired_donor_differences.csv');dif=dif[dif['mode'].eq('primary_cell_weighted')];m=dif.pivot(index='target',columns='donor',values='tumour_minus_peritumour').reindex(genes)
 a=axat(f,[.135,.101,.29,.235],6,'Tumor − peritumor protein');v=max(.5,np.nanmax(np.abs(m)));im=heat(a,m,v=v);a.set_xticklabels(m.columns,rotation=60,ha='right');cbar(f,im,[.445,.14,.009,.15],'ADT difference',[-round(v,1),0,round(v,1)]);ss.append(source(m.reset_index(),'Fig6_paired_tissue_ADT'))
 par=old('protein_partial_programme_rho');m=par.groupby(['target','programme']).rho.median().unstack().reindex(index=genes,columns=P[:-1]);a=axat(f,[.63,.101,.29,.235],7,'Adjusted protein / programme');im=heat(a,m,v=.5);a.set_xticklabels(AB[:-1],rotation=60,ha='right');cbar(f,im,[.94,.14,.009,.15],'Partial ρ',[-.5,0,.5]);ss.append(source(m.reset_index(),'Fig6_adjusted_ADT_programme'))
 legend(f,y=.025);densitykey(f,[.375,.704,.18,.008])
 caption='a, Independent breast CD8 transcript programmes on the full-transcriptome embedding. b–e, Binned expression-weighted spatial densities for IL7R RNA/CD127 protein and PDCD1 RNA/PD-1 protein. Densities are scaled within each panel; corrected-ADT weights are truncated at zero for density display alone. f, Within-donor RNA/ADT Spearman correlation for every matched receptor antibody in five independent donors. g, Corrected-ADT tumor minus peritumor mean within each donor; donor means use the original primary cell-weighted tissue summaries. h, Median within-donor partial rank correlation of corrected ADT with each frozen transcript programme, adjusted for RNA depth, ADT depth, tissue and technical lane. All measured matched antibodies are retained and unmeasured antibodies are not imputed. Density magnitudes cannot compare abundance across antibodies; patient matrices preserve donor replication.'
 record(f,'Fig6','RNA and surface protein resolve complementary receptor context',8,caption,ss+[str(ROOT/'data/BREAST_cells.csv.gz'),str(ROOT.parent/'results/protein_rna_correlations.csv'),str(ROOT.parent/'results/protein_paired_donor_differences.csv'),str(ROOT/'results/protein_partial_programme_rho.csv')])

def f7():
 ex=read('external_functional_clone_logexpr');assert len(ex)==15 and ex.n_cells.sum()==404 and set(ex.patient)=={'MD01-004','MD01-005'}
 path=Path('D:/CD8_Immunotherapy_20261023/analysis_extensions/GSE176021_antigen_validation/tables/tumor_CD8_expression_TCR_functional_labels.csv.gz');e=pd.read_csv(path,low_memory=False).rename(columns={'UMAP_1':'UMAP1','UMAP_2':'UMAP2'});e=e[e.patient.isin(ex.patient.unique())].copy();keys=ex[['patient','clone','antigen']].rename(columns={'clone':'clone_beta'});q=e.merge(keys,on=['patient','clone_beta','antigen'],how='inner',validate='many_to_one');assert len(q)==404
 f=plt.figure(figsize=(183/25.4,9.1));ss=[];a=axat(f,[.075,.744,.26,.20],0,'Functional labels on source UMAP');a.scatter(e.UMAP1,e.UMAP2,s=.8,c='#D9D9D9',lw=0,rasterized=True)
 for ant,col in [('viral',PRE),('MANA',POST)]:
  u=q[q.antigen.eq(ant)];a.scatter(u.UMAP1,u.UMAP2,s=7,c=col,lw=0,rasterized=True,label=ant)
 a.axis('off');a.legend(loc='lower left',fontsize=7,markerscale=2,handletextpad=.2)
 for i,g in enumerate(['IL7R','CXCL13'],1):a=axat(f,[.075+i*.31,.744,.26,.20],i,g+' RNA density');spatial(a,e,e['log_'+g])
 densitykey(f,[.38,.694,.2,.008])
 genes=['IL7R','TCF7','PDCD1','LAG3','ENTPD1','CXCL13','GZMB','TOX'];ex=ex.sort_values(['patient','antigen','clone']).reset_index(drop=True);ids=[f'{r.patient[-3:]} {"M" if r.antigen=="MANA" else "V"}{i+1}' for i,r in ex.iterrows()];m=ex[genes].T;m.columns=ids
 a=axat(f,[.12,.398,.42,.222],3,'Functionally measured clonotype profiles');im,z=rowheat(a,m);a.tick_params(labelsize=6.8);ss.append(source(z.reset_index(),'Fig7_functional_clone_z'));ss.append(source(ex.assign(display_id=ids),'Fig7_functional_clone_key'));cbar(f,im,[.557,.43,.009,.15],'Row z',[-2,0,2])
 tab=pd.crosstab(q.author_state,q.antigen).reindex(columns=['viral','MANA'],fill_value=0);tab=tab.div(tab.sum(axis=0),axis=1)*100;a=axat(f,[.775,.398,.17,.222],4,'Original CD8 states');im=heat(a,tab,v=100,cmap=plt.get_cmap('Blues'));a.set_xticklabels(['Viral','MANA'],rotation=0);a.tick_params(labelsize=6.8);cbar(f,im,[.957,.43,.008,.15],'% cells',[0,50,100]);ss.append(source(tab.reset_index(),'Fig7_functional_author_state_percent'))
 for i,g in enumerate(['IL7R','CXCL13']):
  a=axat(f,[.11+i*.295,.105,.20,.20],i+5,g+' · clone means')
  for pi,p in enumerate(['MD01-004','MD01-005']):
   for ai,ant in enumerate(['viral','MANA']):
    u=ex[ex.patient.eq(p)&ex.antigen.eq(ant)];x=ai+(pi-.5)*.17;a.scatter(x+np.linspace(-.045,.045,len(u)),u[g],c=[PRE,POST][ai],marker=['o','^'][pi],s=23,ec='white',lw=.4);a.plot([x-.075,x+.075],[u[g].mean()]*2,c=INK,lw=.8)
  a.set(xticks=[0,1],xticklabels=['Viral','MANA'],xlim=(-.32,1.32),ylabel='Mean log(1 + CP10k)')
 auc=read('external_within_donor_auc');m=auc.pivot(index='gene',columns='patient',values='AUC').reindex(genes);a=axat(f,[.78,.105,.15,.20],7,'Fixed-direction AUC');im=heat(a,m-.5,v=.5);a.set_yticklabels(genes);a.set_xticklabels(['004','005'],rotation=0)
 for i in range(len(genes)):
  for j in range(2):a.text(j,i,f'{m.iloc[i,j]:.2f}',ha='center',va='center',fontsize=6.8)
 ss.append(source(m.reset_index(),'Fig7_fixed_direction_AUC'));ss.append(source(q[['cell_id','patient','clone_beta','antigen','author_state','UMAP1','UMAP2']+[f'log_{g}' for g in genes]],'Fig7_qualified_functional_cells'))
 f.legend([Line2D([],[],marker=m,ls='',color=INK,ms=4) for m in ['o','^']],['MD01-004','MD01-005'],ncol=2,loc='lower center',bbox_to_anchor=(.34,.022),fontsize=7)
 caption='a, Source author UMAP for the two qualified functional donors, with independently measured viral and mutation-associated neoantigen (MANA) clonotypes overlaid on gray donor CD8 cells. Gray cells carry no inferred antigen specificity. b–c, IL7R and CXCL13 expression-weighted density in all available donor CD8 cells, scaled separately per panel. d, Clone-mean log(1 + CP10k) expression, standardized within gene and hierarchically ordered by average-linkage Euclidean distance. Columns identify 15 source-labelled clonotypes, comprising 404 observed cells from MD01-004 and MD01-005. Retained clones have at least three cells matched by patient-scoped TRB nucleotide sequence, and each donor has at least two clones of each antigen class. e, Original author-state composition among qualified labelled cells, normalized within antigen class across both donors. f–g, Individual clone means; circle and triangle distinguish donors and short bars show each donor/class mean. h, Within-donor univariate ROC AUC, with higher expression always ranking MANA positive. Directions are never reversed to maximize AUC. These fixed-direction rankings describe this measured clone set; no classifier is trained and no expression-derived antigen assignment is made.'
 record(f,'Fig7','Measured antigen labels anchor functional clonotype context',8,caption,ss+[str(path),str(REV/'results/external_functional_clone_logexpr.csv'),str(REV/'results/external_within_donor_auc.csv')])

def f8():
 d=pd.read_csv(CLIN/'data/CAUSHI_CD8_cells_frozen_scores.csv.gz').rename(columns={'UMAP_1':'UMAP1','UMAP_2':'UMAP2'});M=pd.read_csv(CLIN/'results/patient_clinical_features.csv');pr=M[M.primary_analysis].copy();assert len(d)==151601 and len(pr)==13 and len(d[d.patient.isin(pr.patient)])==125887
 f=plt.figure(figsize=(183/25.4,9.15));ss=[]
 for i,title in enumerate(['Original author CD8 states','Frozen CD8 programmes','IL7R RNA density','CXCL13 RNA density']):
  a=axat(f,[.06+i*.24,.775,.215,.165],i,title)
  if i==0:
   states=sorted(d.author_state.unique());colors=plt.get_cmap('tab20')
   for j,state in enumerate(states):
    u=d[d.author_state.eq(state)];a.scatter(u.UMAP1,u.UMAP2,c=[colors(j)],s=.7,lw=0,rasterized=True);a.text(u.UMAP1.median(),u.UMAP2.median(),state.replace('Effector','Eff.'),fontsize=5.8,ha='center',va='center',bbox=dict(fc='white',ec='none',alpha=.65,pad=.2))
   a.axis('off')
  elif i==1:umap(a,d)
  else:spatial(a,d,d['log_'+['IL7R','CXCL13'][i-2]])
 densitykey(f,[.705,.725,.19,.008])
 de=pd.read_csv(OUT/'results/clinical_continuous_DESeq2.csv');a=axat(f,[.10,.423,.35,.225],4,'Continuous residual-tumor transcriptome');u=de[de.pvalue.notna()].copy();u['minuslogp']=-np.log10(np.maximum(u.pvalue,1e-300));sig=u.padj.lt(.05);a.scatter(u.log2FoldChange,u.minuslogp,s=4,c=np.where(sig,np.where(u.log2FoldChange>0,POST,PRE),'#D2D2D2'),lw=0,rasterized=True);a.axvline(0,c=GRAY,lw=.55);a.set(xlabel='log2 fold change per residual fraction 0 → 1',ylabel='−log10 nominal P');a.text(.02,.90,f'{len(de):,} tested genes; {sig.sum()} BH q<0.05',transform=a.transAxes,va='top',fontsize=7)
 for k,r in enumerate(u.sort_values('padj').head(5).itertuples()):a.annotate(r.gene,(r.log2FoldChange,r.minuslogp),xytext=[(5,5),(-36,15),(7,-13),(7,7),(8,-12)][k],textcoords='offset points',fontsize=6.5,arrowprops=dict(arrowstyle='-',lw=.4,color=GRAY))
 ss.append(source(de,'Fig8_continuous_DESeq2'))
 gs=pd.read_csv(OUT/'results/clinical_Hallmark_GSEA.csv');terms=['INTERFERON_GAMMA_RESPONSE','INTERFERON_ALPHA_RESPONSE','IL2_STAT5_SIGNALING','E2F_TARGETS','INFLAMMATORY_RESPONSE','TNFA_SIGNALING_VIA_NFKB','MYC_TARGETS_V1','OXIDATIVE_PHOSPHORYLATION'];g=gs.set_index('Term').reindex(['HALLMARK_'+x for x in terms]);a=axat(f,[.765,.423,.18,.225],5,'Fixed Hallmark context');y=np.arange(len(g));a.barh(y,g.NES,color=np.where(g.NES>0,POST,PRE),height=.58);a.axvline(0,c=GRAY,lw=.5);a.set(yticks=y,yticklabels=['IFN-γ response','IFN-α response','IL2–STAT5','E2F targets','Inflammatory','TNF–NF-κB','MYC targets V1','Ox. phosph.'],xlabel='Normalized enrichment');a.invert_yaxis();a.set_xlim(-2.5,2.5)
 for i,q in enumerate(g['FDR q-val']):a.text(2.6,i,f'{q:.3f}' if q>=.001 else '<0.001',va='center',fontsize=6.7,clip_on=False)
 a.text(2.6,-1.02,'FDR q',fontsize=6.7);ss.append(source(g.reset_index(),'Fig8_fixed_Hallmark_context'))
 norm=pd.read_csv(OUT/'data/clinical_normalized_pseudobulk.csv.gz',index_col=0);genes=de.sort_values(['padj','pvalue']).gene.head(20).tolist();order=pr.sort_values(['residual_tumor_pct','patient']).patient.tolist();mat=np.log1p(norm.loc[order,genes]).T;a=axat(f,[.12,.107,.46,.20],6,'Patient pseudobulk expression');im,z=rowheat(a,mat);a.set_xticklabels([p.replace('MD01-','M').replace('NY','N') for p in order],rotation=60,ha='right');a.tick_params(labelsize=6.8);ss.append(source(z.reset_index(),'Fig8_discovery_gene_patient_z'));cbar(f,im,[.597,.135,.009,.14],'Row z',[-2,0,2]);
 strip=f.add_axes([.12,.352,.46,.014]);strip.imshow(pr.set_index('patient').loc[order].residual_tumor_pct.to_numpy()[None,:],cmap='Blues',vmin=0,vmax=100,aspect='auto');strip.axis('off');f.text(.12,.373,'Residual viable tumor: 0 → 100%',fontsize=6.8)
 a=axat(f,[.765,.118,.18,.215],7,'Frozen Memory primary test');a.scatter(pr.score_Memory,pr.residual_tumor_pct,c=np.where(pr.MPR,PRE,POST),s=27,ec='white',lw=.4);a.set(xlabel='Patient Memory score',ylabel='Residual viable tumor (%)',ylim=(-5,105));a.text(.01,.33,'n=13; Spearman r=0.670\nPermutation P=0.0152',transform=a.transAxes,va='top',fontsize=6.8);a.legend([Line2D([],[],marker='o',ls='',color=PRE,ms=3),Line2D([],[],marker='o',ls='',color=POST,ms=3)],['MPR','Non-MPR'],loc='lower right',fontsize=6.5,handletextpad=.1);ss.append(source(pr[['patient','score_Memory','residual_tumor_pct','MPR']],'Fig8_frozen_primary_patient'))
 legend(f,y=.026)
 caption='a–d, Original author UMAP coordinates for 151,601 CD8 cells from 15 pathology-matched patients, shown by author state, frozen operational programme and IL7R/CXCL13 expression-weighted spatial density. Density magnitudes are scaled within panel. The 2,602 source CD8 records lacking patient identity remain unassigned and are excluded. e, Exploratory genome-wide association of patient pseudobulk UMI counts with continuous residual viable tumor fraction in 13 primary patients, adjusted for source patient-prefix stratum. Positive log2 fold change denotes greater residual viable tumor at resection; effect size corresponds to fraction 0 to 1. BH correction spans all 14,129 tested genes; colored points have q<0.05 and labels identify the five smallest q values. f, A fixed eight-pathway Hallmark set spanning interferon, IL2–STAT5, proliferation, inflammation and metabolism, including nonsignificant results. Printed values are permutation-based GSEA FDR q values; complete Hallmark results accompany the figure. g, Normalized pseudobulk expression for the 20 genes with the smallest adjusted-model q values, row-standardized and average-linkage clustered; primary patients are ordered by residual viable tumor and the strip shows its continuous percentage. This discovery-selected display is exploratory. h, Unchanged prespecified Memory association in the 13 primary patients (125,887 cells), excluding MD01-004 and MD01-005: Spearman r=0.670383, two-sided 99,999-permutation P=0.0152. Color denotes original MPR (≤10% residual viable tumor). None of the 29 prespecified secondary tests passes BH q<0.05; all-15-patient sensitivity remains separate. Expression is measured after neoadjuvant treatment at surgery, and the endpoint is residual viable tumor percentage, not tumor diameter or baseline prediction. The new genome-wide analyses are separate exploratory follow-up to the frozen primary test.'
 record(f,'Fig8','Resection CD8 biology tracks continuous residual viable tumor',8,caption,ss+[str(CLIN/'data/CAUSHI_CD8_cells_frozen_scores.csv.gz'),str(CLIN/'results/clinical_associations_prespecified.csv'),str(OUT/'results/clinical_Hallmark_GSEA.csv'),str(OUT/'data/clinical_pseudobulk_design.csv')])

if __name__=='__main__':
 jobs={'5':f5,'6':f6,'7':f7,'8':f8}
 for k in (sys.argv[1:] or list(jobs)):jobs[k]()
 MAN=[json.loads((OUT/'work'/f'Fig{k}_metadata.json').read_text(encoding='utf-8')) for k in jobs if (OUT/'work'/f'Fig{k}_metadata.json').exists()];CAP={r['figure']:r['caption'] for r in MAN}
 (OUT/'manifest_main_5_8.json').write_text(json.dumps(MAN,indent=2),encoding='utf-8');(OUT/'legends_main_5_8.json').write_text(json.dumps(CAP,indent=2),encoding='utf-8');shutil.copy2(__file__,OUT/'scripts'/Path(__file__).name)
