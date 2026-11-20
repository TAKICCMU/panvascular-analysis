"""Display-only merge of current CoVarNet and spatial figures.

No model, statistical test, normalization, row selection, or confidence interval
calculation. Original source values are read and plotted. Native panel letters
are excluded from standalone assets; included only in the review assembly.
"""
from pathlib import Path
import copy, hashlib, json, re, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
import pymupdf as fitz
from lxml import etree as ET

P=Path('/Users/taki/Desktop/投稿包_最新版_20260919')
O=Path(__file__).resolve().parents[1]
S=P/'Source_Data_Figures'
MM=72/25.4
INK='#26313A'; MUTED='#64717B'; GUIDE='#E5E9EC'
MOD={'CM01':'#526D8C','CM03':'#8D7298','CM09':'#248C85','Other':'#A6ADB4'}
GENE={'FN1':'#C1824A','SCARB1':'#2A7F78'}
EVID={'Positive':'#2A7F78','Inverse':'#C1824A','Not significant':'#747C86'}
plt.rcParams.update({'font.family':'Arial','font.size':6.5,'text.color':INK,
 'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':INK,
 'axes.edgecolor':'#87929B','axes.linewidth':.5,'xtick.major.width':.5,'ytick.major.width':.5,
 'xtick.major.size':2,'ytick.major.size':2,'svg.fonttype':'none','pdf.fonttype':42,
 'ps.fonttype':42,'svg.hashsalt':'PanvascularMergedFigure4_20261116','savefig.facecolor':'white'})
sources={}; extents=[]; geom={}; manifests=[]
def read(rel):
 f=S/rel;sources[str(f)]=hashlib.sha256(f.read_bytes()).hexdigest()
 shutil.copy2(f,O/'source_data'/f.name)
 return pd.read_csv(f,sep='\t',float_precision='round_trip')

N=read('Figure4/05_top_nodes_reference.tsv')
I=read('Figure4/FINAL_STATE_IDENTITY_MAP.tsv')
M=read('Figure4/14_final_module_evidence_and_tiers.tsv')
U=read('Figure5_final_companion/Figure5b_60_gene_unit_estimates.tsv')
Y=read('Figure5_final_companion/Figure5b_point_mapping.tsv')
POOL=read('Figure5_final_companion/Figure5b_pooled_FDR_I2_original.tsv')
F=read('Figure5_final_companion/Fig5c_all29.tsv')
ALL=read('Figure5_final_companion/all29_original_meta.tsv')
COH=read('Figure5_final_companion/module_coherence_meta.tsv')
DS=read('Figure5_final_companion/module_dataset_stats.tsv')
PROG=read('Figure5_final_companion/module_program_direction_meta.tsv')
ORDER=read('Figure3/Fig3d_e_all29.tsv').gene_symbol.tolist()
SEC=read('Figure5/section_input_audit.tsv')
workflow_source=P/'图件编辑资料/BioRender_当前素材/Figure4整图排版/Figure4a_source_and_scope.json'
sources[str(workflow_source)]=hashlib.sha256(workflow_source.read_bytes()).hexdigest()
shutil.copy2(workflow_source,O/'source_data'/workflow_source.name)
F=F.set_index('gene').loc[ORDER].reset_index()
assert len(F)==29 and len(U)==60 and U.biological_unit_id.nunique()==34
assert U.groupby('gene').size().to_dict()=={'FN1':34,'SCARB1':26}
assert M.loc[M.stable_structure & M.broad_context & M.coherent_within_study,'module'].tolist()==['CM09']
J=U.merge(Y,on=['gene','GSE','biological_unit_id'],validate='one_to_one',suffixes=('_raw','_display'))
assert np.allclose(J.beta_raw,J.beta_display,rtol=0,atol=1e-15)

def fig(w,h):
 f=plt.figure(figsize=(w/25.4,h/25.4)); f.mm=(w,h); return f
def ax(f,x,y,w,h):
 fw,fh=f.mm; return f.add_axes([x/fw,1-(y+h)/fh,w/fw,h/fh])
def tx(f,x,y,t,**kw):
 fw,fh=f.mm; return f.text(x/fw,1-y/fh,t,va='center',**kw)
def clean(a,left=True):
 a.spines[['top','right']].set_visible(False)
 if not left:a.spines['left'].set_visible(False)
 a.set_axisbelow(True)
def save(f,name):
 f.canvas.draw(); rr=f.canvas.get_renderer();fw,fh=f.canvas.get_width_height()
 for t in f.findobj(match=matplotlib.text.Text):
  if not t.get_visible() or not t.get_text():continue
  bb=t.get_window_extent(rr)
  extents.append({'asset':name,'text':t.get_text(),'font_pt':t.get_fontsize(),'bbox':[bb.x0,bb.y0,bb.x1,bb.y1],
   'outside':bool(bb.x0<-.5 or bb.y0<-.5 or bb.x1>fw+.5 or bb.y1>fh+.5)})
 for ext in ['pdf','svg','png']:f.savefig(O/'panels'/f'{name}.{ext}',dpi=450)
 geom[name]={'width_mm':f.mm[0],'height_mm':f.mm[1],'panel_letter_embedded':False}
 plt.close(f)

# a: parallel analysis streams. No implied selection of genes by CoVarNet.
f=fig(183,25); a=ax(f,0,0,183,25);a.set(xlim=(0,183),ylim=(25,0));a.axis('off')
def box(x,w,label,sub,fc):
 a.add_patch(FancyBboxPatch((x,5),w,16,boxstyle='round,pad=0.0,rounding_size=1.6',fc=fc,ec='none'))
 a.text(x+w/2,10,label,ha='center',va='center',fontsize=7.3,fontweight='bold')
 a.text(x+w/2,16,sub,ha='center',va='center',fontsize=6.5)
def arrow(x0,x1):a.add_patch(FancyArrowPatch((x0,13),(x1,13),arrowstyle='-|>',mutation_scale=7,lw=.7,color=MUTED))
box(3,39,'Cell fractions','223 profiles · 89 states','#F1F5F7')
arrow(43,48)
box(49,34,'CoVarNet','14 candidate modules','#EAF4F2')
a.plot([89,89],[5,21],c=GUIDE,lw=.6)
box(95,40,'Spatial data','49 sections · 29 candidates','#F1F5F7')
arrow(136,142)
box(143,37,'Association models','Expression and cell abundance','#F4F1F5')
save(f,'Figure4a')

# b: complete original top-five weights of the stable module and comparator.
f=fig(83,59)
display={
 'CM09':['GZMK cytotoxic memory T','ZNF683 cytotoxic T','GIMAP7 memory T','MKI67 proliferating T','IL7R memory T'],
 'CM03':['CXCL12 chemokine Fib','PI16 adventitial-like Fib','LAMA2 tissue-associated Fib','FAP activated ECM Fib','PROX1 lymphatic EC']}
rows=[]
for module,ys in [('CM09',7),('CM03',34)]:
 d=N[N.module.eq(module)].sort_values('rank').head(5).merge(I,left_on='cell_state',right_on='state_id',sort=False,validate='one_to_one')
 a=ax(f,49,ys,31,17.5);y=np.arange(5)
 a.barh(y,d.weight.to_numpy(),height=.63,fc=MOD[module],ec='none',zorder=3)
 a.set(xlim=(0,.3),ylim=(4.5,-.5));a.set_yticks([])
 a.set_xticks([0,.15,.3],['0','0.15','0.30'],fontsize=6.1)
 a.grid(axis='x',color=GUIDE,lw=.4);clean(a,left=False)
 a.tick_params(axis='x',pad=1.5)
 tx(f,49,ys-3.4,module,fontsize=7.2,fontweight='bold')
 tx(f,80,ys-3.4,'Stable' if module=='CM09' else 'Exploratory',ha='right',fontsize=6.2,color=MUTED)
 for j,(_,r) in enumerate(d.iterrows()):
  tx(f,47.5,ys+(j+.5)/5*17.5,display[module][j],fontsize=6.5,ha='right')
  rows.append({'module':module,'rank':r['rank'],'state_id':r.cell_state,'full_name':r.level2_annotation,'display_name':display[module][j],'weight':r.weight})
tx(f,64.5,57,'NMF weight',ha='center',fontsize=6.5)
pd.DataFrame(rows).to_csv(O/'source_data/Figure4b_display_and_original_identity.tsv',sep='\t',index=False)
save(f,'Figure4b')

# c: original beta values and deterministic display positions for all 60 estimates.
f=fig(83,80); a=ax(f,21,8,58,54)
gses=U.GSE.drop_duplicates().tolist();assert len(gses)==7
for i in range(7):
 if i%2==0:a.axhspan(i-.5,i+.5,fc='#F7F8FA',ec='none',zorder=0)
a.axvline(0,c='#94A0A9',lw=.6,zorder=0)
for gene,mark in [('FN1','o'),('SCARB1','^')]:
 d=Y[Y.gene.eq(gene)];a.scatter(d.beta,d.display_y,s=12,marker=mark,color=GENE[gene],edgecolors='none',zorder=3)
a.axhline(6.8,color=GUIDE,lw=.55)
for y,gene in [(7.5,'FN1'),(8.55,'SCARB1')]:
 r=POOL[POOL.gene.eq(gene)].iloc[0]
 a.plot([r.CI_low,r.CI_high],[y,y],c=GENE[gene],lw=1.05,zorder=3)
 a.scatter([r.beta],[y],marker='D',s=24,c=GENE[gene],edgecolors='none',zorder=4)
a.set(xlim=(-1.02,.55),ylim=(9.2,-.65));a.set_yticks(list(range(7))+[7.5,8.55],gses+['FN1 pooled','SCARB1 pooled'],fontsize=6.1)
a.set_xticks([-1,-.5,0,.5],['−1.0','−0.5','0','0.5'],fontsize=6.2)
a.tick_params(axis='y',length=0,pad=2);clean(a,left=False)
a.set_xlabel('Conditional Fib abundance coefficient',fontsize=6.4,labelpad=2.3)
for x,gene,mark in [(26,'FN1','o'),(49,'SCARB1','^')]:
 f.add_artist(Line2D([x/83],[1-3/80],marker=mark,ms=3.3,mfc=GENE[gene],mec='none',ls='',transform=f.transFigure))
 tx(f,x+2.7,3,gene,fontsize=6.5)
for y,gene in [(71,'FN1'),(75,'SCARB1')]:
 r=POOL[POOL.gene.eq(gene)].iloc[0]
 q=float(r.FDR_BH_all_29_pairs_x_3_terms)
 if q < 0.001:
  exp=int(np.floor(np.log10(q))); value=fr'${q/10**exp:.2f}\times10^{{{exp}}}$'
 else: value=f'{q:.4f}'
 tx(f,21,y,f'{gene}: FDR {value}; I² {r.I2:.1f}%',fontsize=6.3)
save(f,'Figure4c')

# d: all 29 candidates retain original order, coefficients, CI endpoints and FDR.
f=fig(96,143); a=ax(f,31,10,49,124)
a.axvline(0,c='#9AA4AC',lw=.65,zorder=0)
for i,r in F.iterrows():
 if i%2==0:a.axhspan(i-.5,i+.5,fc='#F7F8FA',ec='none',zorder=-1)
 cat='Not significant' if not r.FDR_BH_all_29_pairs_x_3_terms<.05 else ('Positive' if r.beta>0 else 'Inverse')
 a.plot([r.CI_low,r.CI_high],[i,i],c=EVID[cat],lw=.75,zorder=2)
 a.scatter([r.beta],[i],s=9.5,c=EVID[cat],edgecolors='none',zorder=3)
a.set(xlim=(-.46,.70),ylim=(28.5,-.5));a.set_yticks([])
a.set_xticks([-.4,-.2,0,.2,.4,.6],['−0.4','−0.2','0','0.2','0.4','0.6'],fontsize=6.1)
clean(a,left=False);a.set_xlabel('Conditional coefficient (95% CI)',fontsize=6.5,labelpad=2.4)
tx(f,1,7,'Gene',fontsize=6.4);tx(f,22,7,'Cell',fontsize=6.4);tx(f,94,7,'FDR',ha='right',fontsize=6.4)
short={'Fibroblast':'Fib','Endothelial':'EC','Endothelial cell':'EC','Macrophage':'Mac'}
for i,r in F.iterrows():
 y=10+(i+.5)/29*124
 tx(f,1,y,r.gene,fontsize=6.5);tx(f,22,y,short[r.cell_context],fontsize=6.2)
 tx(f,94,y,'<0.001' if r.FDR_BH_all_29_pairs_x_3_terms<.001 else f'{r.FDR_BH_all_29_pairs_x_3_terms:.3f}',fontsize=6.2,ha='right')
for x,label in [(6,'Positive'),(29,'Inverse'),(50,'Not significant')]:
 f.add_artist(Line2D([x/96],[1-2/143],marker='o',ms=2.7,mfc=EVID[label],mec='none',ls='',transform=f.transFigure))
 tx(f,x+2.5,2,label,fontsize=6.1)
save(f,'Figure4d')

# Extended Data 10a: original strict-parent null summary for all 14 modules.
def coh_panel(w,h,name):
 d=COH[(COH.posterior=='q05')&(COH.null_method=='strict_parent')]
 v=d.pivot(index='module',columns='metric',values='observed_median').sort_index()
 q=d.pivot(index='module',columns='metric',values='FDR_BH_28').sort_index()
 assert len(v)==14
 f=fig(w,h);a=ax(f,12,9,w-16,h-22)
 offsets={'CM01':(4,4),'CM02':(-7,-6),'CM03':(5,7),'CM04':(-9,3),'CM05':(-9,-10),'CM06':(-13,5),'CM07':(7,-4),
  'CM08':(6,-3),'CM09':(6,5),'CM10':(7,6),'CM11':(-7,5),'CM12':(-8,6),'CM13':(5,-6),'CM14':(5,-6)}
 for module,r in v.iterrows():
  sig=bool((q.loc[module]<.05).any());c=MOD.get(module,MOD['Other'])
  a.scatter(r.within_spot,r.spatial_lag,s=18,marker='s' if sig else 'o',color=c,edgecolors='none',zorder=3)
  a.annotate(module[-2:],(r.within_spot,r.spatial_lag),xytext=offsets[module],textcoords='offset points',fontsize=6.1,ha='center',va='center')
 a.set(xlim=(.25,.93),ylim=(.05,.5));a.set_xticks([.3,.5,.7,.9]);a.set_yticks([.1,.2,.3,.4,.5]);a.tick_params(labelsize=6.1)
 a.set_xlabel('Within-spot coherence',fontsize=6.4,labelpad=2)
 a.set_ylabel('Neighbourhood coherence',fontsize=6.4,labelpad=2);clean(a)
 f.add_artist(Line2D([12/w],[1-3/h],marker='s',ms=3.5,mfc=INK,mec='none',ls='',transform=f.transFigure))
 tx(f,15,3,'FDR < 0.05',fontsize=6.1)
 f.add_artist(Line2D([40/w],[1-3/h],marker='o',ms=3.3,mfc=INK,mec='none',ls='',transform=f.transFigure))
 tx(f,43,3,'FDR ≥ 0.05',fontsize=6.1)
 save(f,name)
coh_panel(65,69,'ExtendedData10a')

# ED10b/c/d/e preserve original source-count, candidate sensitivity and I2 data.
# Section-source summary is taken from the existing metadata, not inferred donors.
f=fig(65,45);a=ax(f,23,5,39,31)
gcol='GSE' if 'GSE' in SEC else 'study'
if gcol not in SEC:
 raise ValueError(SEC.columns.tolist())
ucol=next(c for c in ['biological_unit_id','analysis_unit_id'] if c in SEC)
cnt=SEC.groupby(gcol).agg(sections=(gcol,'size'),units=(ucol,'nunique')).reindex(gses)
a.barh(np.arange(7)-.18,cnt.sections,height=.32,color='#3B6F91',label='Sections')
a.barh(np.arange(7)+.18,cnt.units,height=.32,color='#B4C3CF',label='Analysis units')
a.set(ylim=(6.6,-.6),xlim=(0,20));a.set_yticks(range(7),gses,fontsize=6.0);a.set_xticks([0,10,20]);a.tick_params(axis='y',length=0,pad=2)
clean(a,left=False);a.set_xlabel('Count',fontsize=6.2,labelpad=1.6)
for x,t,c in [(5,'Sections (49)','#3B6F91'),(31,'Analysis units (35)','#B4C3CF')]:
 f.add_artist(Line2D([x/65],[1-2/45],marker='s',ms=3.1,mfc=c,mec='none',ls='',transform=f.transFigure));tx(f,x+2.4,2,t,fontsize=6)
save(f,'ExtendedData10b')

f=fig(65,39);a=ax(f,12,6,50,24)
# Histogram bins exactly match the retained figure's 0,20,40,60,80,100 grid.
a.hist(F.I2,bins=[0,20,40,60,80,100],color='#7F9AA3',edgecolor='white',lw=.5)
a.set(xlim=(0,100));a.set_xticks([0,50,100]);a.set_yticks([0,5,10]);a.set_xlabel('I² (%)',fontsize=6.4,labelpad=2)
a.set_ylabel('Candidates',fontsize=6.4,labelpad=2);clean(a)
save(f,'ExtendedData10d')

f=fig(68,170);a=ax(f,29,25,37,132)
for i in range(29):
 if i%2==0:a.axhspan(i-.5,i+.5,fc='#F7F8FA',ec='none')
a.axvline(0,c='#A2ADB5',lw=.6)
for sensitivity,dy,col,mark in [('all_sections',-.13,'#3B6F91','D'),('exclude_QC_fallback',.13,'#D89232','s')]:
 z=ALL[(ALL.posterior=='q05')&(ALL.term=='cell_abundance')&(ALL.sensitivity==sensitivity)].set_index('gene').loc[ORDER]
 for j,(_,r) in enumerate(z.iterrows()):
  a.plot([r.CI_low,r.CI_high],[j+dy,j+dy],c=col,lw=.55)
  a.scatter(r.beta,j+dy,s=8,marker=mark,c=col,edgecolors='none')
a.set(xlim=(-.46,.73),ylim=(28.5,-.5));a.set_yticks([]);a.set_xticks([-.25,0,.25,.5],['−0.25','0','0.25','0.50'],fontsize=6.0)
a.set_xlabel('Coefficient (95% CI)',fontsize=6.2,labelpad=2);clean(a,left=False)
for i,r in F.iterrows():tx(f,28,25+(i+.5)/29*132,f'{r.gene} | {short[r.cell_context]}',ha='right',fontsize=6.0)
for y,l,c,m in [(8,'All sections','#3B6F91','D'),(12,'Exclude QC-fallback','#D89232','s')]:
 f.add_artist(Line2D([9/68],[1-y/170],marker=m,ms=3.1,mfc=c,mec='none',ls='',transform=f.transFigure));tx(f,12,y,l,fontsize=6.1)
save(f,'ExtendedData10c')

f=fig(43,170);a=ax(f,1,25,40,132)
sets=[('mean','all_sections'),('mean','exclude_QC_fallback'),('q05','all_sections'),('q05','exclude_QC_fallback')]
matrix=[]
for post,sens in sets:
 z=ALL[(ALL.posterior==post)&(ALL.term=='cell_abundance')&(ALL.sensitivity==sens)].set_index('gene').loc[ORDER]
 matrix.append(z.FDR_BH_all_29_pairs_x_3_terms.to_numpy())
matrix=np.array(matrix).T
im=a.imshow(-np.log10(np.clip(matrix,np.finfo(float).tiny,1)),cmap='Blues',vmin=0,vmax=10,aspect='auto',interpolation='none',extent=(-.5,3.5,28.5,-.5))
a.set(ylim=(28.5,-.5));a.set_yticks([]);a.set_xticks([]);a.spines[:].set_visible(False)
for x in [.5,1.5,2.5]:a.axvline(x,c='white',lw=.4)
for j in np.arange(.5,28.5):a.axhline(j,c='white',lw=.4)
tx(f,11,11,'Mean',ha='center',fontsize=6.2);tx(f,31,11,'q05',ha='center',fontsize=6.2)
for j,t in enumerate(['All','Exclude QC-\nfallback','All','Exclude QC-\nfallback']):tx(f,1+(j+.5)*10,21,t,ha='center',fontsize=6.0,linespacing=1.0)
cb=f.colorbar(im,cax=ax(f,1,161,40,1.5),orientation='horizontal');cb.set_ticks([0,5,10],labels=['0','5','≥10']);cb.ax.tick_params(labelsize=6.0,pad=1,length=1)
tx(f,21,169,'−log10 FDR',ha='center',fontsize=6.0)
save(f,'ExtendedData10e')

# Supplementary panels retain the original four-setting count summary and all
# 35 dataset correlations / 5 saved medians, without a new ED number.
f=fig(180,71)
a=ax(f,53,9,39,42)
counts=[]
for j,(post,sens) in enumerate(sets):
 z=ALL[(ALL.posterior==post)&(ALL.term=='cell_abundance')&(ALL.sensitivity==sens)]
 pos=((z.FDR_BH_all_29_pairs_x_3_terms<.05)&(z.beta>0)).sum();neg=((z.FDR_BH_all_29_pairs_x_3_terms<.05)&(z.beta<0)).sum();ns=29-pos-neg
 counts.append([post,sens,int(pos),int(neg),int(ns)])
 left=0
 for val,col in [(pos,EVID['Positive']),(neg,EVID['Inverse']),(ns,'#D9DDE1')]:
  a.barh(j,val,left=left,height=.58,color=col,edgecolor='none');a.text(left+val/2,j,str(val),ha='center',va='center',fontsize=6.1);left+=val
 a.text(-1,j,('Mean' if post=='mean' else 'q05')+' | '+('All sections' if sens=='all_sections' else 'Exclude QC-fallback'),ha='right',va='center',fontsize=6.1)
a.set(ylim=(3.7,-.7),xlim=(0,29));a.set_yticks([]);a.set_xticks([0,10,20,29]);a.set_xlabel('Candidates',fontsize=6.3);clean(a,left=False)
for x,l,c in [(4,'Positive',EVID['Positive']),(31,'Inverse',EVID['Inverse']),(54,'Not significant','#D9DDE1')]:
 f.add_artist(Line2D([x/180],[1-61/71],marker='s',ms=3.8,mfc=c,mec='none',ls='',transform=f.transFigure));tx(f,x+2.8,61,l,fontsize=6.3)
pa=ax(f,130,9,38,42)
pairs=[('CM01','fibrous_cap_core','CM01 | fibrous cap'),('CM01','calcification','CM01 | calcification'),('CM02','inflammation','CM02 | inflammation'),('CM03','endothelial_activation','CM03 | EC activation'),('CM09','inflammation','CM09 | inflammation')]
out=[]
for j,(m,prog,label) in enumerate(pairs):
 vals=DS[(DS.posterior=='q05')&(DS.module==m)][['GSE','rho__'+prog]].sort_values('GSE')
 r=PROG[(PROG.posterior=='q05')&(PROG.module==m)&(PROG.program==prog)].iloc[0]
 ys=j+np.linspace(-.2,.2,len(vals));pa.scatter(vals['rho__'+prog],ys,s=10,c='#A3ABB2',edgecolors='none')
 pa.scatter([r.median_rho],[j],s=22,marker='D',c=MOD.get(m,MOD['Other']),edgecolors='none')
 tx(f,128,9+(j+.6)/5.2*42,label,ha='right',fontsize=6.1)
 tx(f,179,9+(j+.6)/5.2*42,f'{r.FDR_BH_84:.2f}',ha='right',fontsize=6.1)
 for (_,v),y in zip(vals.iterrows(),ys):out.append({'module':m,'program':prog,'GSE':v.GSE,'rho':v['rho__'+prog],'display_y':y,'saved_median':r.median_rho,'FDR':r.FDR_BH_84})
pa.axvline(0,c='#B4BCC3',lw=.6);pa.set(xlim=(-.6,.95),ylim=(4.6,-.6));pa.set_yticks([]);pa.set_xticks([-.5,0,.5],['−0.5','0','0.5']);pa.set_xlabel('Spearman ρ',fontsize=6.3);clean(pa,left=False)
tx(f,179,3,'FDR',ha='right',fontsize=6.2)
for x,t,mark,c in [(131,'Dataset','o','#A3ABB2'),(155,'Median','D',INK)]:
 f.add_artist(Line2D([x/180],[1-61/71],marker=mark,ms=3.2,mfc=c,mec='none',ls='',transform=f.transFigure));tx(f,x+2.8,61,t,fontsize=6.3)
pd.DataFrame(out).to_csv(O/'source_data/Supplementary_module_program_display.tsv',sep='\t',index=False)
pd.DataFrame(counts,columns=['posterior','section_set','positive','inverse','not_significant']).to_csv(O/'source_data/Supplementary_four_settings_counts.tsv',sep='\t',index=False)
save(f,'Supplementary_spatial_sensitivity_and_programmes')

def assemble(name,canvas,placements):
 ns='http://www.w3.org/2000/svg';W,H=canvas
 root=ET.Element('{'+ns+'}svg',nsmap={None:ns},width=f'{W}mm',height=f'{H}mm',viewBox=f'0 0 {W*MM} {H*MM}')
 ET.SubElement(root,'{'+ns+'}rect',width='100%',height='100%',fill='white')
 doc=fitz.open();pg=doc.new_page(width=W*MM,height=H*MM)
 for letter,(asset,x,y,w,h) in placements.items():
  src=fitz.open(O/'panels'/f'{asset}.pdf');pg.show_pdf_page(fitz.Rect(x*MM,y*MM,(x+w)*MM,(y+h)*MM),src,0)
  sr=ET.parse(str(O/'panels'/f'{asset}.svg')).getroot();vb=list(map(float,sr.attrib['viewBox'].split()))
  sx=w*MM/vb[2];sy=h*MM/vb[3];assert abs(sx-sy)<1e-5,(asset,sx,sy)
  group=ET.SubElement(root,'{'+ns+'}g',id=f'panel_{letter}',transform=f'translate({x*MM} {y*MM}) scale({sx})')
  # Scope IDs and references to avoid silent collisions in merged SVG.
  prefix=f'{letter}_'; ids={el.get('id'):prefix+el.get('id') for el in sr.iter() if el.get('id')}
  for el in sr.iter():
   if el.get('id'):el.set('id',ids[el.get('id')])
   for k,v in list(el.attrib.items()):
    for old,new in ids.items():
     v=v.replace(f'url(#{old})',f'url(#{new})')
     if v==f'#{old}':v=f'#{new}'
    el.set(k,v)
  for child in sr:group.append(copy.deepcopy(child))
  ly=y+2.4;lx=x+.4
  label=ET.SubElement(root,'{'+ns+'}text',x=str(lx*MM),y=str(ly*MM),style="font-family:Arial;font-weight:bold;font-size:8px;fill:#222831")
  label.text=letter
  pg.insert_text((lx*MM,ly*MM),letter,fontsize=8,fontname='hebo',color=(.13,.16,.19))
 ET.ElementTree(root).write(str(O/'review'/f'{name}.svg'),encoding='utf-8',xml_declaration=True)
 doc.save(O/'review'/f'{name}.pdf',garbage=4,deflate=True)
 pg.get_pixmap(matrix=fitz.Matrix(180/72,180/72),alpha=False).save(O/'review'/f'{name}.png')
 (O/'review'/f'{name}_placements.json').write_text(json.dumps({'canvas_mm':canvas,'panels':placements},indent=2))

placements={'a':['Figure4a',0,0,183,25],'b':['Figure4b',0,27,83,59],'c':['Figure4c',0,90,83,80],'d':['Figure4d',87,27,96,143]}
assemble('Figure4_merged',[183,170],placements)
assemble('ExtendedData_Figure10',[180,170],{'a':['ExtendedData10a',0,0,65,69],'b':['ExtendedData10b',0,73,65,45],'d':['ExtendedData10d',0,124,65,39],'c':['ExtendedData10c',67,0,68,170],'e':['ExtendedData10e',137,0,43,170]})
assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==v for p,v in sources.items())
pd.DataFrame(extents).to_csv(O/'review/text_extent_check.tsv',sep='\t',index=False)
(O/'panel_manifest.json').write_text(json.dumps(geom,indent=2))
(O/'source_hashes.json').write_text(json.dumps(sources,ensure_ascii=False,indent=2))
checks={'original_sources_unchanged':True,'candidate_count':len(F),'candidate_order':ORDER,
 'FN1_estimates':int((U.gene=='FN1').sum()),'SCARB1_estimates':int((U.gene=='SCARB1').sum()),'recorded_unit_IDs':int(U.biological_unit_id.nunique()),
 'unit_effect_values_preserved':True,'unit_SE_values_preserved_in_source':True,'unit_CI_newly_computed':False,
 'pooled_CI_endpoints_read_from_source':True,'module_weights_original':True,'scientific_analyses_rerun':False,
 'ED10c_e_shared_y_centres_mm':[25+(i+.5)/29*132 for i in range(29)],
 'ED10c_e_shared_ylim':[28.5,-.5],'all14_module_matchednull_retained':True,
 'supplementary_module_program_points':len(out),'supplementary_saved_medians':5,
 'text_outside':[x for x in extents if x['outside']]}
(O/'review/integrity_and_layout_checks.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2))
print(json.dumps({'panels':geom,'text_outside':checks['text_outside']},ensure_ascii=False))
