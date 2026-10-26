"""Whole source-derived supplementary units; numbered to the multipanel study."""
import expansion_figures as E
import expansion_figures_extra as X
from largepanel_core import *
import textwrap,shutil
MAN=[];CAP={};CUR='';GENERATED=[];USED=[]
def src(d,name,index=False):
 p=source(d,CUR+'_'+name,index);GENERATED.append(p);return p
def rd(name):
 p=REV/'results'/(name+'.csv');USED.append(str(p));return pd.read_csv(p)
def oldrd(name):
 p=ROOT/'results'/(name+'.csv');USED.append(str(p));return pd.read_csv(p)
def ce(c):USED.append(str(ROOT/'data'/(c+'_cells.csv.gz')));return cells(c)
def pa(ax,i,title):
 width=ax.get_position().width;panel(ax,i,'\n'.join(textwrap.wrap(title,24 if width<.19 else 32)))
def gr(title,rows,cols,height=9):
 fig,axs=grid(title,rows,cols,max(height,9.2 if rows==3 else 10 if rows==4 else 7.4));fig.subplots_adjust(left=.105,right=.95,top=.90,bottom=.14,wspace=.65,hspace=.82);return fig,axs
def ex(fig,stem,title,n,caption,sources):
 if CUR=='FigS17':
  fig.delaxes(fig.axes[-1]);cbar(fig,fig.axes[1].images[0],[.44,.065,.18,.008],'Median Jaccard',[0,.5,1],orientation='horizontal')
 if CUR=='FigS22':
  a=fig.axes[8];a.clear();pa(a,8,'Oral · retained CD8 cells');u=rd('qc_patient_visit');u=u[u.cohort.eq('ORAL')&u.metric.eq('library_UMI')];paired_lines(a,u,'n_cells');a.set_ylabel('Retained CD8 cells');a.set_ylim(bottom=0);src(u,'oral_sampling_counts');caption=caption.replace('a–l, Patient/visit','a–h and j–l, Patient/visit');caption+=' Panel i shows retained oral CD8 sampling counts because mitochondrial fractions are unavailable in the original oral source; missing mitochondrial measurements are not imputed.'
 rec=export(fig,CUR,title,n,caption,list(dict.fromkeys(USED+GENERATED)));MAN.append(rec);CAP[CUR]=caption
for module in [E,X]:
 for name,obj in {'OUT':OUT,'read':rd,'old':oldrd,'source':src,'cells':ce,'panel':pa,'grid':gr,'export':ex,'heat':heat,'zscore':zscore,'paired_lines':paired_lines,'legend':legend,'note':note}.items():setattr(module,name,obj)

def breast():
 d=ce('BREAST');f=plt.figure(figsize=(183/25.4,10));genes=sum(DEFS.values(),[]);parts=[]
 for _,u in d.groupby(['programme','patient','time']):parts.append(u.sort_values('cell_id').sample(n=min(30,len(u)),random_state=20261026))
 s=pd.concat(parts).copy();s=s.assign(po=s.programme.map({p:i for i,p in enumerate(P)})).sort_values(['po','patient','time','cell_id']);mat=pd.DataFrame(zscore(s[['log_'+g for g in genes]].to_numpy().T),index=genes,columns=s.cell_id);mat=mat.iloc[leaves_list(linkage(mat,method='average',optimal_ordering=True))]
 a=f.add_axes([.135,.61,.77,.285]);panel(a,0,f'Breast cell expression · {len(s):,} balanced cells');a.set_title(f'Breast cell expression · {len(s):,} balanced cells',loc='left',pad=24);im=heat(a,mat,xlabels=False);a.tick_params(axis='y',labelsize=6.5);src(mat,'cell_marker_z',True);src(s[['cell_id','patient','time','programme']],'balanced_cells');cbar(f,im,[.927,.70,.009,.12],'Row z',[-2,0,2]);a.set_xlabel('Programme blocks; patient / tissue strata within each block')
 strip=f.add_axes([.135,.901,.77,.012]);strip.imshow(np.array([[to_rgba(PC[p]) for p in s.programme]]),aspect='auto');strip.axis('off')
 parts=[]
 for (p,prog),u in d.groupby(['patient','programme']):
  if len(u)>=5:parts.append(pd.Series(u[['log_'+g for g in TARGETS]].mean().to_numpy(),index=TARGETS,name=(p,prog)))
 profiles=pd.DataFrame(parts);profiles.index=pd.MultiIndex.from_tuples(profiles.index,names=['patient','programme']);m=profiles.groupby('programme').mean().reindex(P).T;z=pd.DataFrame(zscore(m),index=m.index,columns=AB);z=z.iloc[leaves_list(linkage(z,method='average',optimal_ordering=True))];a=f.add_axes([.12,.265,.285,.26]);panel(a,1,'Patient-equal receptor RNA');im=heat(a,z);a.set_xticklabels(AB,rotation=60,ha='right');src(z,'patient_equal_receptor_z',True);cbar(f,im,[.42,.33,.009,.13],'Row z',[-2,0,2])
 cor=pd.read_csv(ROOT.parent/'results/protein_rna_correlations.csv');USED.append(str(ROOT.parent/'results/protein_rna_correlations.csv'));cor=cor[cor.grouping.eq('donor')];m=cor.pivot(index='target',columns='donor',values='rna_adt_rho');m=m.loc[m.median(axis=1).sort_values().index];a=f.add_axes([.65,.265,.255,.26]);panel(a,2,'All matched RNA / ADT assays');im=heat(a,m,v=.6);a.set_xticklabels(m.columns,rotation=60,ha='right');cbar(f,im,[.927,.33,.009,.13],'Spearman ρ',[-.6,0,.6]);src(m,'all_assay_concordance',True)
 for j,g in enumerate(['IL7R','PDCD1','LAG3']):
  a=f.add_axes([.105+j*.31,.065,.235,.115]);panel(a,j+3,g+' RNA / protein');a.scatter(d['log_'+g],d[g+'_adt_corrected'],s=1,c=GRAY,alpha=.18,lw=0,rasterized=True);rows=[]
  for pi,(p,u) in enumerate(d.groupby('patient')):
   rho=stats.spearmanr(u['log_'+g],u[g+'_adt_corrected']).statistic;rows.append({'patient':p,'gene':g,'n':len(u),'rho':rho});a.scatter(u['log_'+g].mean(),u[g+'_adt_corrected'].mean(),c=[plt.get_cmap('tab10')(pi)],s=20,ec='white',lw=.4)
  a.set(xlabel='RNA log(1 + CP10k)',ylabel='Corrected ADT');a.text(.03,.96,f'Donor median ρ={np.median([r["rho"] for r in rows]):.2f}',transform=a.transAxes,fontsize=6.8,va='top');src(pd.DataFrame(rows),g+'_donor_agreement')
 caption='a, All 38 frozen programme markers in a deterministic display subset of at most 30 cells per patient/tissue/programme stratum. Columns retain programme blocks and patient/tissue order; rows show log(1 + CP10k) z scores and average-linkage Euclidean clustering. b, All 24 receptor transcripts averaged within patient/programme groups containing at least five cells, then averaged equally over patients and row-standardized. d–f, Full-range individual-cell RNA versus corrected ADT for IL7R/CD127, PDCD1/PD-1 and LAG3, including negative corrected protein values. Larger colored points are donor means; printed correlations summarize the five independent donor-specific correlations. c, All 17 matched receptor RNA/ADT correlations within each donor. Cell subsampling is for display alone; no antibody or gene is added by imputation.'
 ex(f,CUR,'Breast cell profiles and matched RNA/protein context',6,caption,[])

def robustness():
 M=pd.read_csv(CLIN/'results/patient_clinical_features.csv');A=pd.read_csv(CLIN/'results/clinical_associations_prespecified.csv');S=pd.read_csv(CLIN/'results/primary_sensitivities.csv');L=pd.read_csv(CLIN/'results/primary_leave_one_out.csv');USED.extend([str(CLIN/'results'/n) for n in ['patient_clinical_features.csv','clinical_associations_prespecified.csv','primary_sensitivities.csv','primary_leave_one_out.csv']]);f,axs=gr('',4,2,10.5)
 sets=[(M,'score_Memory','All 15 pathology-matched patients',A[A.feature.eq('Memory')&A.subset.eq('all_patients')].iloc[0]),(M[M.primary_analysis],'equal_library_score_Memory','Equal weighting of tumor libraries',S[S.analysis.eq('equal_library')].iloc[0]),(M[M.primary_analysis],'rank_Memory','Frozen rank-based Memory score',S[S.analysis.eq('rank_scoring')].iloc[0]),(M,'n_cells','Sampled CD8 depth and pathology',None),(M[M.primary_analysis],'log_IL7R','IL7R and residual viable tumor',A[A.feature.eq('IL7R')&A.subset.eq('primary')].iloc[0])]
 for i,(u,x,title,r) in zip([0,1,2,4,5],sets):
  a=axs.flat[i];pa(a,i,title)
  for primary in [True,False]:
   v=u[u.primary_analysis.eq(primary)];a.scatter(v[x],v.residual_tumor_pct,c=np.where(v.MPR,PRE,POST) if primary else 'none',edgecolors=np.where(v.MPR,PRE,POST),s=20,lw=.7)
  a.set(ylabel='Residual viable tumor (%)',xlabel={'score_Memory':'Memory score','equal_library_score_Memory':'Equal-library Memory score','rank_Memory':'Memory rank score','n_cells':'Sampled CD8 cells','log_IL7R':'Mean IL7R log(1 + CP10k)'}[x],ylim=(-5,108))
  if x=='n_cells':a.set_xscale('log')
  if r is not None:a.text(.02,.40,f'n={int(r["n"])}; r={r["rho"]:.3f}; P={r["p"]:.4f}',va='top',transform=a.transAxes,fontsize=6.7,bbox=dict(fc='white',ec='none',alpha=.8,pad=.2))
 a=axs.flat[3];pa(a,3,'Leave-one-patient-out Memory');a.scatter(L.rho,range(len(L)),c=CC['NSCLC'],s=15);a.axvline(.670383,c=GRAY,lw=.6,ls='--');a.set(yticks=range(len(L)),yticklabels=L.omitted_patient,xlim=(-1,1),xlabel='Spearman r after omission');a.tick_params(axis='y',labelsize=6.4)
 de=pd.read_csv(OUT/'results/clinical_continuous_DESeq2.csv').set_index('gene');un=pd.read_csv(OUT/'results/clinical_continuous_unadjusted_DESeq2.csv').set_index('gene');m=de[['log2FoldChange','padj']].join(un[['log2FoldChange','padj']],lsuffix='_adjusted',rsuffix='_unadjusted');a=axs.flat[6];pa(a,6,'Full-gene adjustment sensitivity');a.scatter(m.log2FoldChange_adjusted,m.log2FoldChange_unadjusted,c=np.where(m.padj_adjusted.lt(.05)|m.padj_unadjusted.lt(.05),POST,GRAY),s=2,lw=0,alpha=.5,rasterized=True);a.axhline(0,c=GRAY,lw=.5);a.axvline(0,c=GRAY,lw=.5);a.set(xlabel='Prefix-adjusted log2 fold change',ylabel='Unadjusted log2 fold change');src(m,'gene_model_comparison',True)
 aa=set(de.index[de.padj.lt(.05)]);bb=set(un.index[un.padj.lt(.05)]);counts=[len(aa-bb),len(aa&bb),len(bb-aa)];a=axs.flat[7];pa(a,7,'Gene-level BH discovery overlap');a.bar(range(3),counts,color=[PRE,'#9769A8',POST],width=.6);a.set(xticks=range(3),xticklabels=['Adjusted\nonly','Both','Unadjusted\nonly'],ylabel='Genes at BH q<0.05',ylim=(0,max(counts)*1.2+1));src(pd.DataFrame({'category':['adjusted_only','both','unadjusted_only'],'n':counts}),'BH_discovery_overlap')
 USED.extend([str(OUT/'results'/n) for n in ['clinical_continuous_DESeq2.csv','clinical_continuous_unadjusted_DESeq2.csv']]);note(f,'Filled points: primary patients. Open points: prior functional donors. Endpoint: resection residual viable tumor, not baseline prediction.')
 ex(f,CUR,'Clinical association and genome-wide adjustment robustness',8,'a–c, Memory association in all 15 pathology-matched patients, with equal weighting of libraries in the 13 primary patients, and using the frozen rank-based score. d, All 13 leave-one-patient-out Memory coefficients; the dashed line is the unchanged primary estimate. e, Sampled CD8 cell count against continuous residual viable tumor percentage on a log depth axis. f, Prespecified IL7R expression association in the primary set. Filled symbols are primary patients and open symbols are MD01-004/005; color denotes original MPR. Statistics in a–c and f are patient-level Spearman correlations and two-sided permutation P values. g, Full-gene residual-fraction coefficients from the exploratory source-prefix-adjusted and unadjusted patient pseudobulk models. Color identifies BH q<0.05 in either model. h, Exact intersection and model-specific discovery counts, with BH correction applied to the complete 14,129-gene family in each model. This comparison does not change the primary model or the frozen Memory test.',[])

def genome():
 f=plt.figure(figsize=(183/25.4,9.8));models=[('clinical_continuous_DESeq2','Primary prefix adjusted'),('clinical_continuous_unadjusted_DESeq2','Unadjusted sensitivity'),('clinical_continuous_sexadjusted_DESeq2','Post hoc sex adjusted')];gfiles=['clinical_Hallmark_GSEA','clinical_Hallmark_GSEA_unadjusted','clinical_sexadjusted_Hallmark_GSEA'];gs=[]
 for j,(name,title) in enumerate(models):
  p=OUT/'results'/(name+'.csv');u=pd.read_csv(p);USED.append(str(p));a=f.add_axes([.09+j*.305,.79,.245,.145]);pa(a,j,title);a.scatter(u.log2FoldChange,-np.log10(u.pvalue.clip(lower=1e-300)),s=2,c=np.where(u.padj.lt(.05),np.where(u.log2FoldChange>0,POST,PRE),GRAY),lw=0,rasterized=True);a.axvline(0,c=GRAY,lw=.5);a.set(xlabel='log2 FC / residual fraction',ylabel='−log10 nominal P');a.text(.02,.96,f'{u.padj.lt(.05).sum()} / {len(u):,} BH q<0.05',transform=a.transAxes,va='top',fontsize=6.7)
  p=OUT/'results'/(gfiles[j]+'.csv');v=pd.read_csv(p).set_index('Term');gs.append(v);USED.append(str(p))
 order=sorted(set.intersection(*[set(g.index) for g in gs]));m=pd.DataFrame({name:g.NES.reindex(order) for name,g in zip(['Primary','Unadj.','Sex post hoc'],gs)});a=f.add_axes([.315,.10,.20,.58]);panel(a,3,'Complete eligible Hallmark pathways');im=heat(a,m,v=2.5);a.set_yticklabels([t.replace('HALLMARK_','').replace('_',' ').capitalize().replace('Dna ','DNA ').replace('E2f ','E2F ').replace('G2m ','G2M ').replace('Il2 stat5','IL2 STAT5').replace('Il6 jak stat3','IL6 JAK STAT3').replace('Kras ','KRAS ').replace('Mtorc1 ','mTORC1 ').replace('Myc ','MYC ').replace('Pi3k akt mtor','PI3K AKT mTOR').replace('Tgf ','TGF ').replace('Tnfa signaling via nfkb','TNFα signaling via NFκB').replace('Uv ','UV ').replace('Wnt ','WNT ') for t in order],fontsize=6.7);a.set_xticklabels(m.columns,rotation=60,ha='right',fontsize=6.8)
 for j,g in enumerate(gs):
  for i,t in enumerate(order):
   if g.loc[t,'FDR q-val']<.05:a.text(j,i,'•',ha='center',va='center',fontsize=7,color=INK)
 cbar(f,im,[.535,.125,.008,.15],'NES',[-2.5,0,2.5]);src(m,'complete_Hallmark_NES',True)
 curves=pd.read_csv(OUT/'results/clinical_GSEA_curves.csv.gz');USED.append(str(OUT/'results/clinical_GSEA_curves.csv.gz'));terms=['HALLMARK_MYC_TARGETS_V1','HALLMARK_PROTEIN_SECRETION'];norm=pd.read_csv(OUT/'data/clinical_normalized_pseudobulk.csv.gz',index_col=0);des=pd.read_csv(OUT/'data/clinical_pseudobulk_design.csv');orderp=des[des.primary_analysis].sort_values(['residual_tumor_pct','patient']).patient;USED.extend([str(OUT/'data/clinical_normalized_pseudobulk.csv.gz'),str(OUT/'data/clinical_pseudobulk_design.csv')])
 for k,t in enumerate(terms):
  u=curves[curves.term.eq(t)];a=f.add_axes([.685,.565-k*.145,.265,.09]);panel(a,k+4,['MYC targets V1 enrichment','Protein secretion enrichment'][k]);a.plot(u['rank']/len(u),u.running_ES,c=CC['NSCLC'],lw=1);a.axhline(0,c=GRAY,lw=.5);a.set(xlabel='Fraction of Wald-ranked genes',ylabel='Running ES');q=gs[0].loc[t,'FDR q-val'];a.text(.03,.18,f'NES={gs[0].loc[t,"NES"]:.2f}; FDR q={q:.3f}',transform=a.transAxes,va='top',fontsize=6.7);src(u,t+'_running_ES')
  genes=gs[0].loc[t,'Lead_genes'].split(';');de=pd.read_csv(OUT/'results/clinical_continuous_DESeq2.csv').set_index('gene');genes=de.reindex(genes).sort_values('stat').index[:10];mat=np.log1p(norm.loc[orderp,genes]).T;z=pd.DataFrame(zscore(mat),index=genes,columns=orderp);z=z.iloc[leaves_list(linkage(z,method='average',optimal_ordering=True))];a=f.add_axes([.69,.245-k*.16,.25,.105]);panel(a,k+6,['MYC leading-edge expression','Protein leading-edge expression'][k]);im=heat(a,z);a.set_xticklabels([p.replace('MD01-','M').replace('MD043-','D').replace('NY016-','N') for p in orderp],rotation=60,ha='right',fontsize=6.3);a.tick_params(axis='y',labelsize=6.5);src(z,t+'_leading_edge_z',True)
 cbar(f,im,[.955,.13,.008,.15],'Row z',[-2,0,2]);note(f,'Exploratory full-gene analyses. Sex adjustment is explicitly post hoc; the primary prefix-adjusted model is unchanged.')
 ex(f,CUR,'Complete clinical transcriptome and post hoc sex sensitivity',8,'a–c, Genome-wide residual viable tumor fraction coefficients from the primary source-prefix-adjusted model, unadjusted sensitivity, and explicitly post hoc source-prefix-plus-sex model. All models use 13 primary patient pseudobulks and the same 14,129-gene BH family; discoveries number 31, 18 and 27, respectively. Positive coefficients denote greater residual viable tumor at resection. Sex sensitivity was added after sex-linked genes appeared among the primary discoveries; it does not replace the primary model. d, All 48 size-eligible Hallmark pathways shared by the three models in alphabetical order. Color is normalized enrichment score and dots mark permutation-based GSEA FDR q<0.05. No pathway is selected for significance in this matrix. e–f, Full Wald-ranked running enrichment for MYC targets V1 and protein secretion in the primary model, with original GSEA statistics. g–h, Primary-model leading-edge expression for these two pathways. Ten genes per pathway are selected by most negative primary Wald statistic, with row-standardized log1p normalized pseudobulk counts and average-linkage clustering. Patient columns follow residual viable tumor percentage. These selected expression displays are exploratory; complete leading edges, models and source-prefix/sex distributions accompany the figure. Fifteen primary-model/sex-adjusted BH discoveries overlap. The frozen clinical Memory test remains unchanged.',[])

def run(k):
 global CUR,GENERATED,USED
 CUR='FigS'+str(k);GENERATED=[];USED=[]
 fun={16:E.repertoire,17:X.networks,18:X.clone_heatmaps,19:X.clinical_context,20:breast,21:X.thresholds,22:X.qc,23:robustness,24:genome,25:X.paga}[k];fun()
 # The imported whole-figure clone routine writes three matrices directly.
 if k==18:
  rec=MAN[-1]
  for c in CO:
   for t in ['pre','post']:
    p=OUT/'figure_source_data'/f'FigS23_{c}_{t}_clone_z.csv'
    if p.exists():rec['source'].append(str(p))
  (OUT/'work'/(CUR+'_metadata.json')).write_text(json.dumps(rec,indent=2),encoding='utf-8')
if __name__=='__main__':
 for k in map(int,sys.argv[1:] or range(16,26)):run(k)
 MAN=[json.loads((OUT/'work'/f'FigS{k}_metadata.json').read_text()) for k in range(16,26) if (OUT/'work'/f'FigS{k}_metadata.json').exists()];CAP={r['figure']:r['caption'] for r in MAN}
 (OUT/'manifest_supp_16_25.json').write_text(json.dumps(MAN,indent=2));(OUT/'legends_supp_16_25.json').write_text(json.dumps(CAP,indent=2));shutil.copy2(__file__,OUT/'scripts'/Path(__file__).name)
 for name in ['expansion_figures.py','expansion_figures_extra.py','expansion_plot_core.py']:shutil.copy2(Path(__file__).parent/name,OUT/'scripts'/name)
