from largepanel_core import *
from pypdf import PdfReader
import xml.etree.ElementTree as ET

def run():
 rn=old('receptor_programme_rho_summary');rp=old('receptor_programme_patient_rho');ad=old('protein_programme_patient_rho');adj=old('protein_partial_programme_rho');progs=P[:-1];abbr=AB[:-1];mats={};fullsources={'RNA':source(rp,'FigS26_all_RNA_patient_rho'),'ADT':source(ad,'FigS26_all_ADT_donor_rho'),'Adjusted':source(adj,'FigS26_all_adjusted_ADT_donor_rho')}
 for c in CO:mats[c]=rn[rn.cohort.eq(c)].pivot(index='target',columns='programme',values='median').reindex(index=TARGETS,columns=progs)
 mats['ADT']=ad.groupby(['target','programme']).rho.median().unstack().reindex(index=TARGETS,columns=progs);mats['Adjusted']=adj.groupby(['target','programme']).rho.median().unstack().reindex(index=TARGETS,columns=progs)
 fig=plt.figure(figsize=(183/25.4,7.7));dets={};paths={};contexts=['BCC','NSCLC','ORAL','ADT','Adjusted'];titles=['BCC RNA','Lung RNA','Oral RNA','Breast ADT','Adjusted breast ADT']
 for j,(c,title) in enumerate(zip(contexts,titles)):
  col=j if j<3 else j-3;row=0 if j<3 else 1;y=.55 if row==0 else .11;h=.39 if row==0 else .33;ax=fig.add_axes([.15+col*.285,y,.232,h]);panel(ax,j,title);w=mats[c].copy();w.columns=abbr;im=heat(ax,w,v=.5,ylabels=col==0);ax.tick_params(axis='y',labelsize=6.9);ax.tick_params(axis='x',labelsize=6.8)
  if row:
   for a,b in np.argwhere(w.isna().to_numpy()):ax.text(b,a,'×',ha='center',va='center',color='#999999',fontsize=6)
  paths[c]=source(w.reset_index(names='target'),'FigS26_'+c+'_median_rho');src=[paths[c],fullsources['RNA' if j<3 else c],ROOT/'results'/('receptor_programme_rho_summary.csv' if j<3 else 'protein_programme_patient_rho.csv' if j==3 else 'protein_partial_programme_rho.csv')]
  dets[string.ascii_lowercase[j]]={'identity':title+' receptor-programme correlations','sources':[str(p) for p in src],'quantitative_meaning':'Median within-patient single-cell Spearman correlation between target log RNA and receptor-independent programme score' if j<3 else 'Median within-donor Spearman correlation between corrected ADT and programme score' if j==3 else 'Median within-donor partial rank correlation adjusted for RNA depth, ADT depth, tissue and technical lane'}
 cbar(fig,im,[.957,.65,.006,.18],'ρ',[-.5,0,.5]);ax=fig.add_axes([.785,.11,.17,.33]);panel(ax,5,'Patient/donor correlations');pairs=[('Memory','IL7R'),('Dysfunction','PDCD1'),('Dysfunction','LAG3')];colors={**CC,'ADT':'#AD6AC7','Adjusted':'#543F7E'};rows=[]
 for i,(p,g) in enumerate(pairs):
  for j,c in enumerate(contexts):
   if c in CO:u=rp[rp.cohort.eq(c)&rp.programme.eq(p)&rp.target.eq(g)][['patient','rho']]
   elif c=='ADT':u=ad[ad.programme.eq(p)&ad.target.eq(g)][['donor','rho']].rename(columns={'donor':'patient'})
   else:u=adj[adj.programme.eq(p)&adj.target.eq(g)][['patient','rho']]
   v=u.rho.dropna();y=i*1.4+(j-2)*.16;ax.scatter(v,y+np.linspace(-.025,.025,len(v)),s=11,c=colors[c],ec='white',lw=.2,zorder=3);ax.plot([v.median()-.035,v.median()+.035],[y,y],c=INK,lw=1,zorder=4);rows.extend(u.assign(context=c,programme=p,target=g).to_dict('records'))
 lo=min(-.05,min(r['rho'] for r in rows if np.isfinite(r['rho']))-.04);hi=max(.65,max(r['rho'] for r in rows if np.isfinite(r['rho']))+.04);ax.axvline(0,c=GRAY,lw=.55);ax.set(xlim=(lo,hi),ylim=(3.35,-.55),yticks=[0,1.4,2.8],yticklabels=['IL7R / Mem','PDCD1 / Dys','LAG3 / Dys'],xlabel='Spearman ρ',xticks=[0,.3,.6]);ax.tick_params(axis='y',labelsize=6.6);ax.tick_params(axis='x',labelsize=6.8);rep=source(pd.DataFrame(rows),'FigS26_selected_patient_donor_correlations');dets['f']={'identity':'Independent patient and donor distributions for three receptor-programme pairs','sources':[rep]+list(fullsources.values())+[str(ROOT/'results/receptor_programme_patient_rho.csv'),str(ROOT/'results/protein_programme_patient_rho.csv'),str(ROOT/'results/protein_partial_programme_rho.csv')],'quantitative_meaning':'Every point is one measured independent patient/donor correlation; black segments are context-specific medians. Context order within each pair: BCC RNA, lung RNA, oral RNA, breast ADT, adjusted breast ADT.'}
 fig.legend([Line2D([],[],marker='o',c=colors[c],ls='',ms=3.8) for c in contexts],titles,loc='lower center',bbox_to_anchor=(.52,.015),ncol=5,fontsize=6.7,columnspacing=1.05,handletextpad=.35)
 caption='a–c, Median within-patient Spearman correlations between each of the 24 fixed target transcripts and the six continuous receptor-independent programme scores in BCC, lung and oral CD8 cells. d, Median within-donor Spearman correlations between the same programme scores and independently measured corrected breast ADT. e, Median within-donor partial rank correlations after adjustment for RNA library depth, ADT library depth, tissue and technical lane, using the frozen original adjustment model. All five maps retain the same 24-target row order and six-programme column order. Color is capped at ±0.5; gray cells marked × indicate unavailable or undefined protein associations. Seventeen of the 24 targets were measured by ADT, and the remaining targets are retained as missing. f, Independent patient/donor correlations for IL7R–Memory, PDCD1–Dysfunction and LAG3–Dysfunction in the five assay/model contexts. Each point represents one patient or donor; black segments mark medians. Within each receptor/programme pair, rows follow the context order shown in the legend. All correlations, missing measurements and adjusted estimates are reused unchanged from the original analysis. These correlations are distinct quantities from the programme-mean expression profiles in Fig. 4.'
 rec=export(fig,'FigS26','Complete measured receptor-programme correlation evidence',list('abcdef'),caption,sorted({p for v in dets.values() for p in v['sources']}));rec['panel_details']=dets;rec['previous_figures']=['original main Fig7','original FigS7 (adjusted breast ADT)'];(OUT/'work/FigS26_metadata.json').write_text(json.dumps(rec,indent=2),encoding='utf-8');(OUT/'manifest_supp_26.json').write_text(json.dumps([rec],indent=2),encoding='utf-8');(OUT/'legends_supp_26.json').write_text(json.dumps({'FigS26':caption},indent=2),encoding='utf-8')
 # Independent equality and Table3 rounding checks.
 expected={('Memory','IL7R'):[.338,.293,.271,.278,.245],('Dysfunction','PDCD1'):[.168,.290,.275,.135,.165],('Dysfunction','LAG3'):[.152,.229,.211,.032,.015]};checks={}
 for (p,g),numbers in expected.items():
  actual=[float(mats[c].loc[g,p]) for c in contexts];assert np.allclose(np.round(actual,3),numbers,atol=1e-12);checks[g+' / '+p]={'full_precision_medians':actual,'Table3_rounded':numbers,'matches':True}
 for c in CO:
  m=rp[rp.cohort.eq(c)].groupby(['target','programme']).rho.median().unstack().reindex(index=TARGETS,columns=progs);assert np.allclose(m,mats[c],equal_nan=True)
 assert set(ad.target)==set(adj.target) and ad.target.nunique()==17
 assert mats['ADT'].shape==mats['Adjusted'].shape==(24,6);assert int(mats['ADT'].isna().all(axis=1).sum())==7 and int(mats['Adjusted'].isna().all(axis=1).sum())==7
 for ext,dpi in [('png',300),('tif',600)]:
  im=Image.open(OUT/'figures'/('FigS26.'+ext));assert abs(float(im.info['dpi'][0])-dpi)<1
  if ext=='tif':assert im.mode=='RGB'
 assert abs(float(PdfReader(OUT/'figures/FigS26.pdf').pages[0].mediabox.width)*25.4/72-183)<.01
 ntext=len(ET.parse(OUT/'figures/FigS26.svg').findall('.//{http://www.w3.org/2000/svg}text'));assert ntext>100
 report={'panels':6,'all_source_paths_exist':all(Path(p).exists() for p in rec['source']),'Table3_checks':checks,'RNA_medians_recomputed_from_patient_tables':True,'measured_ADT_targets':17,'missing_ADT_targets':7,'unmeasured_protein_values_preserved_missing':True,'svg_editable_text_nodes':ntext,'PDF_183mm_PNG_300dpi_RGB_TIFF_600dpi':True};(OUT/'work/verification_supp_26.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':run()
