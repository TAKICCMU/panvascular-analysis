"""Panvascular main panels: render saved results, never fit/recompute science.

Bioinformatics_Visual_Library (private, cc1ccc9) supplies matrix primitives.
The grouped-marker and interval recipes inform the adapters below. Their
default Z scaling, heavy dot edges, headings and footnotes are not adopted.
"""
import sys, shutil, json, hashlib, textwrap, re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'reference'))
from matrix_figures import numeric_matrix, linear_dot_area, dot_matrix, matrix_image
import redraw_data_panels as h
from redraw_data_panels import np,pd,plt,E,fitz,Normalize,Line2D,Patch,Rectangle,LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle
P=h.P; S=h.S; C=h.C; L1=h.L1
h.O=ROOT
plt.rcParams.update({'font.size':8,'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,
 'axes.linewidth':.65,'grid.color':'#E7EAED','grid.linewidth':.45})
read=h.read; axmm=h.axmm; fig=h.fig; clean=h.clean; save=h.save
NAVY=C['method']['MAGMA']; GREY='#747C86'; LIGHT='#E7EAED'
RECORDS=h.records
def labels_wrap(xs,w=40):return ['\n'.join(textwrap.wrap(str(x),w,break_long_words=False)) for x in xs]
def colors_modules(ms):return [C['module'].get(x,C['module']['Other']) for x in ms]
def box(a,x,y,w,height,text,fc='#F2F5F6',fs=9):
 a.add_patch(FancyBboxPatch((x,y),w,height,boxstyle='round,pad=0.008,rounding_size=0.014',facecolor=fc,edgecolor='#B9C4C9',linewidth=.7))
 a.text(x+w/2,y+height/2,text,ha='center',va='center',fontsize=fs,linespacing=1.5)
def arrow(a,p,q,color=GREY):a.add_patch(FancyArrowPatch(p,q,arrowstyle='-|>',mutation_scale=9,linewidth=.8,color=color))
def schematic(n):
 f=fig(183,48);a=axmm(f,4,6,175,38);a.set(xlim=(0,1),ylim=(0,1));a.axis('off')
 if n==1:
  for x in [.03,.055,.08]:
   a.add_patch(Rectangle((x,.65),.13,.22,fc='white',ec='#A7B8C3',lw=.7))
   a.plot([x+.02,x+.10],[.81,.81],color=NAVY,lw=1)
   a.plot([x+.02,x+.09],[.75,.75],color='#BDCAD1',lw=.8)
  a.text(.14,.53,'Public human vascular datasets',ha='center',fontsize=8)
  a.text(.14,.28,'26 studies\n247 sample identifiers',ha='center',fontsize=9,linespacing=1.6)
  arrow(a,(.29,.60),(.34,.60));box(a,.36,.27,.28,.62,'Harmony integration\nL1 / L2 annotation',fs=10)
  a.text(.50,.08,'1,175,115 cells · 116 named states',ha='center',fontsize=9)
  arrow(a,(.66,.60),(.71,.60));box(a,.73,.48,.25,.40,'Atlas resource',fc='#EAF2F2',fs=10)
  a.text(.855,.25,'Genetics · expression\nModules · spatial evidence',ha='center',fontsize=8,linespacing=1.6)
  a.text(.855,.08,'Downstream analyses',ha='center',fontsize=8)
 elif n==2:
  a.set_xlim(-.015,1.015)  # Include rounded end-box edges in the visible frame.
  box(a,.005,.35,.15,.48,'8 vascular\nGWAS',fs=10);arrow(a,(.16,.59),(.20,.59))
  for x,y,t in [(.29,.80,'G'),(.225,.37,'A'),(.355,.37,'S')]:
   a.add_patch(Circle((x,y),.055,fc='white',ec=NAVY,lw=.85));a.text(x,y,t,ha='center',va='center',fontsize=9)
  # The orthogonal factors are separate; do not imply G -> A/S paths.
  a.text(.29,.12,'Genomic SEM',ha='center',fontsize=8)
  arrow(a,(.39,.59),(.43,.59));box(a,.45,.35,.18,.48,'15,639 evaluable\ngene–state pairs',fs=9)
  arrow(a,(.635,.59),(.67,.59));box(a,.685,.35,.13,.48,'Five joint\ncriteria',fs=9)
  arrow(a,(.82,.59),(.855,.59));box(a,.87,.35,.125,.48,'29\nprioritised pairs',fc='#EAF2F2',fs=9)
 else:
  # Tiny generic matrices are schematic elements, not synthetic measured data.
  for i in range(5):
   for j in range(7):a.add_patch(Rectangle((.01+j*.020,.44+i*.065),.017,.052,fc='#D7E3E9',ec='none'))
  a.text(.08,.25,'223 profiles\n89 retained states',ha='center',fontsize=9,linespacing=1.5)
  arrow(a,(.17,.60),(.225,.60));box(a,.24,.38,.19,.47,'NMF\n14 modules',fs=10)
  arrow(a,(.445,.60),(.49,.60));box(a,.51,.38,.23,.47,'Structure\nCoverage · coherence',fs=9)
  arrow(a,(.755,.60),(.80,.60));box(a,.82,.38,.17,.47,'Perturbations\nEvidence tiers',fc='#EAF2F2',fs=9)
 p=P/'图件编辑资料/editable_svg'/f'Figure{n}.svg'
 save(f,f'Figure{n}a',[p],'Compact editable process schematic; existing numbers and analytical relations retained')

def marker():
 b=S/'Current_Display_20260920/marker_audit_20260921';z=b/'per_gene_scaled_mean.tsv';d=b/'original_detection_fraction.tsv'
 Z=numeric_matrix(read(z).set_index('level1_annotation'));D=numeric_matrix(read(d).set_index('level1_annotation'),row_order=Z.index,column_order=Z.columns)
 assert Z.shape==D.shape==(14,28) and np.nanmin(D)>=0 and np.nanmax(D)<=1
 f=fig(173,85);a=axmm(f,41,8,103.6,53.2);rh=3.8;cw=3.7;diameter=.72*min(rh,cw)*72/25.4;smax=diameter**2
 vals=[]
 for i,st in enumerate(Z.index):
  for j,g in enumerate(Z.columns):vals.append((st,g,Z.iloc[i,j],D.iloc[i,j]))
 tab=pd.DataFrame(vals,columns=['state','gene','mean','fraction']);valid=tab['mean'].notna()&tab.fraction.notna()
 artist=dot_matrix(a,tab[valid],x='gene',y='state',x_order=Z.columns,y_order=Z.index,color='mean',areas=linear_dot_area(tab.loc[valid,'fraction'],domain=(0,1),min_area=0,max_area=smax),cmap=h.BL,vmin=0,vmax=1,color_space='Stored per-gene min–max mean, 14 categories',size_space='Original raw counts >0 fraction, within type',edgecolors='none',linewidth=0)
 a.set(xlim=(-.5,27.5),ylim=(13.5,-.5));a.set_xticks(range(28),Z.columns,rotation=90,fontstyle='italic');a.set_yticks(range(14),[s.replace(' cell','') for s in Z.index]);a.tick_params(length=0,pad=4);a.tick_params(axis='y',pad=7);a.spines[:].set_visible(False)
 for i,st in enumerate(Z.index):a.add_patch(Rectangle((-.92,i-.39),.16,.78,fc=L1[st],ec='none',clip_on=False))
 missing=~valid
 if missing.any():
  a.scatter([Z.columns.get_loc(x) for x in tab.loc[missing,'gene']],[Z.index.get_loc(x) for x in tab.loc[missing,'state']],marker='x',s=9,c=GREY)
 cb=f.colorbar(artist,cax=axmm(f,151,13,2.4,27));cb.set_ticks([0,.5,1]);cb.set_label('Per-gene scaled mean expression',labelpad=5);cb.outline.set_visible(False)
 la=axmm(f,43,76,100,7);la.axis('off');la.set(xlim=(0,100),ylim=(0,1));la.text(0,.5,'Detected cells (%)',va='center',fontsize=8)
 for x,v in zip([50,64,78,92],[.25,.5,.75,1]):la.scatter(x,.65,s=smax*v,c=GREY,ec='none');la.text(x,-.03,str(round(v*100)),ha='center',fontsize=7)
 save(f,'Figure1g',[z,d],'Library dot_matrix + linear_dot_area; preserved B values; borderless dots, independent identity stripe, explicit keys',{'matrix_cells':392,'area_pt2_max':smax,'diameter_fraction_of_pitch':.72,'zero_area_for_zero':True,'scaled_values_recomputed':False})

def summaries():
 p=S/'Figure1/L1_palette_bound.tsv';states=read(p).sort_values('order').L1.tolist();orig=P/'图件编辑资料/editable_svg/Figure1.svg'
 r=E.parse(str(orig));texts=lambda id:[''.join(t.itertext()) for t in r.xpath('//*[@id=$v]',v=id)[0].iter('{http://www.w3.org/2000/svg}text')]
 count=[int(x.replace(',','')) for x in texts('axes_5')[-16:]];studies=[int(x) for x in texts('axes_6')[-16:]]
 pd.DataFrame({'L1':states,'cells':count,'studies':studies}).to_csv(ROOT/'data_used/Figure1c_f_exact_display_values.tsv',sep='\t',index=False)
 for letter,values,lab in [('c',count,'Cells per lineage'),('f',studies,'Studies per lineage')]:
  f=fig(101,88);a=axmm(f,36,7,55,67);a.barh(range(16),values,height=.56,color=[L1[x] for x in states]);a.set_yticks(range(16),[x.replace(' cell','') for x in states]);a.set(ylim=(15.5,-.5),xlim=(0,max(values)*1.25),xlabel=lab);clean(a);a.grid(axis='x',lw=.4);a.set_axisbelow(True)
  for i,v in enumerate(values):a.text(v+max(values)*.025,i,f'{v:,}',va='center',fontsize=7.5)
  if letter=='c':a.set_xticks([0,100000,200000],['0','100k','200k'])
  else:a.set_xticks([0,10,20])
  save(f,'Figure1'+letter,[p,orig],'Shared lineage order; compact names, full counts and separate physical axes')
 import redraw_evidence_panels as ep
 ep.save=save;ep.tissue()

def genetic_matrix():
 p=S/'Figure2/Figure2_Fig3a_LDSC_genetic_correlations.tsv';d=read(p)
 print('GENETIC_COLUMNS',d.columns.tolist())
 # Adapter uses the same explicit disease aliases as the adopted panel.
 alias={'AA':'AD','TAAD':'TAA'};order=['AD','AAA','TAA','CAD','Stroke','PAD','CA','CAVS']
 mat=d.set_index('trait').rename(index=alias,columns=alias).loc[order,order]
 q=P/'Public_Release/Figure_Code/native_dependencies/Panvascular_FigureYa_20260919/code/genetics_assets/rg_with_P.tsv'
 pvals=read(q).pivot(index='trait1',columns='trait2',values='rg_p');pvals=pvals.combine_first(pvals.T)
 pv=pvals.rename(index=alias,columns=alias).loc[order,order]
 assert mat.shape==(8,8) and pv.notna().all().all()
 f=fig(105,91);a=axmm(f,14,8,66,66);im=matrix_image(a,mat,cmap='RdBu_r',vmin=-1,vmax=1,center=0,value_space='LDSC genetic correlation')
 for i in range(8):
  for j in range(i+1):
   a.text(j,i,f'{mat.iloc[i,j]:.2f}',ha='center',va='center',color='white' if abs(mat.iloc[i,j])>.7 else '#20252A',fontsize=8)
   if i!=j and pv.iloc[i,j]<.05:a.plot(j+.30,i-.28,'o',ms=1.8,c='#20252A')
 a.set_xticks(range(8),order,rotation=45,ha='right');a.set_yticks(range(8),order);a.tick_params(length=0);a.spines[:].set_visible(False)
 cb=f.colorbar(im,cax=axmm(f,85,13,2.5,37));cb.set_label('Genetic correlation (r$_g$)',labelpad=5);cb.outline.set_visible(False);f.text(.82,.20,'● Nominal\n   P < 0.05',fontsize=7)
 save(f,'Figure2b',[p,q],'Library matrix_image; numeric values centred; significance dots offset from numbers; all 8 traits retained')

def sem():
 p=S/'Figure2/Figure2_Fig3d_bifactor_loadings.tsv';d=read(p);order=['AA','AAA','TAAD','CAD','Stroke','PAD','CA','CAVS'];names=['AD','AAA','TAA','CAD','Stroke','PAD','CA','CAVS']
 f=fig(183,82);a=axmm(f,5,7,173,68);a.set(xlim=(-.02,1.02),ylim=(0,1));a.axis('off');xs=np.linspace(.055,.945,8)
 for x,name in zip(xs,names):box(a,x-.048,.39,.096,.16,name,fc='white',fs=9)
 centers={'G':(.50,.89),'AORTIC':(.19,.10),'ATHERO':(.69,.10)};co={'G':NAVY,'AORTIC':'#9D6979','ATHERO':'#2A8B86'}
 for group,(x,y) in centers.items():
  a.add_patch(FancyBboxPatch((x-.105,y-.062),.21,.124,boxstyle='round,pad=0.01,rounding_size=.055',fc='white',ec=co[group],lw=1))
  a.text(x,y,{'G':'Shared G','AORTIC':'Aortic','ATHERO':'Atherosclerotic'}[group],ha='center',va='center',fontsize=10)
 for _,r in d.iterrows():
  sx,sy=centers[r.lhs];tx=xs[order.index(r.rhs)];top=r.lhs=='G';start=(sx,sy-.075 if top else sy+.075);end=(tx,.56 if top else .38)
  a.add_patch(FancyArrowPatch(start,end,arrowstyle='-|>',mutation_scale=6,lw=.8,color=co[r.lhs],linestyle='--' if r.p_value>=.05 else '-'))
  # Place exact rounded loadings near the trait end, away from the dense fan.
  tt=.76 if top else .70;px=start[0]+(end[0]-start[0])*tt;py=start[1]+(end[1]-start[1])*tt
  a.text(px,py,f'{r.STD_All:.2f}',ha='center',va='center',fontsize=8,bbox=dict(fc='white',ec='none',pad=1.3))
 save(f,'Figure2c',[p],'Redrawn orthogonal factor graph; same standardised loadings, aliases and dashed P≥.05 paths')

def venn():
 p=Path('/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Content_Alignment_20260918/source_data/genetic_sets/Venn_gene_membership.tsv');d=read(p);assert d.gene_id.is_unique
 cols=['positional_set','eQTL_set','MAGMA_Bonferroni_set'];bits=d[cols].astype(int).astype(str).agg(''.join,axis=1);counts=bits.value_counts().to_dict()
 expected={'100':65,'010':1113,'001':31,'110':207,'101':10,'011':8,'111':224};assert counts==expected and len(d)==1658
 shutil.copy2(p,ROOT/'data_used/Figure2f_gene_membership.tsv')
 f=fig(96,93);a=axmm(f,5,8,85,80);a.set(xlim=(-1.9,1.9),ylim=(-1.9,1.65));a.set_aspect('equal');a.axis('off')
 for x,y,co in [(-.65,.40,C['method']['MAGMA']),(.65,.40,C['evidence']['support']),(0,-.65,'#B78D45')]:
  a.add_patch(Circle((x,y),1.12,facecolor=co,edgecolor=co,alpha=.10,lw=0));a.add_patch(Circle((x,y),1.12,fill=False,edgecolor=co,lw=1.35))
 for x,y,t in [(-.90,1.13,'Positional'),(.90,1.13,'eQTL'),(0,-1.38,'MAGMA\n(Bonferroni)')]:a.text(x,y,t,ha='center',va='center',fontsize=9)
 places={'100':(-1.15,.28),'010':(1.15,.28),'001':(0,-1.02),'110':(0,.85),'101':(-.64,-.54),'011':(.64,-.54),'111':(0,-.05)}
 for region,(x,y) in places.items():a.text(x,y,str(counts[region]),ha='center',va='center',fontsize=10)
 save(f,'Figure2f',[p],'Unique-gene membership validated; schematic circles with separate label zones and exact seven region counts',{'union_genes':1658,'triple_overlap':224,'region_counts':counts,'area_proportional':False})

def upset():
 base=Path('/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Nature_Reconstruction_20260918/source_data/genetics');p=base/'intersection_counts.tsv';q=base/'interval_membership.tsv';d=read(p);members=read(q)
 assert d['count'].tolist()==[128,41,25,13,2,2,1,1]
 for _,r in d.iterrows():assert ((members[['set0','set1','set2','set3']]==r[['set0','set1','set2','set3']]).all(axis=1)).sum()==r['count']
 f=fig(132,83);a=axmm(f,49,9,74,39);b=axmm(f,49,52,74,23,);xs=np.arange(8);a.bar(xs,d['count'],width=.57,color=NAVY);a.set(xlim=(-.6,7.6),ylim=(0,145),ylabel='Interval components');a.set_xticks([]);a.spines[['top','right','bottom']].set_visible(False)
 for x,v in enumerate(d['count']):a.text(x,v+3,str(v),ha='center',va='bottom',fontsize=8)
 for j,r in d.iterrows():
  active=np.flatnonzero(r[['set0','set1','set2','set3']].to_numpy(dtype=bool));b.scatter([j]*4,range(4),s=20,c='#E7EAED',ec='none');b.plot([j]*len(active),active,color=NAVY,lw=1);b.scatter([j]*len(active),active,s=25,c=NAVY,ec='none',zorder=3)
 b.set(xlim=(-.6,7.6),ylim=(3.5,-.5));b.set_yticks(range(4),['Primary r² ≥ 0.10','Primary r² ≥ 0.60','Leave-CAD r² ≥ 0.10','Leave-CAD r² ≥ 0.60']);b.set_xticks([]);b.tick_params(length=0);b.spines[:].set_visible(False)
 save(f,'Figure2e',[p,q],'Counts and original interval membership validated; aligned UpSet bars and dot intersections; enlarged category area',{'intersections':8,'interval_components':len(members)})

def genetics():
 # Accepted coordinates, identities, selected 19 states and original FDR bins.
 h.genetics()
 # Replace the generated b with wider separated method columns; leave a intact.
 RECORDS[:]=[r for r in RECORDS if not (r['figure']=='Figure3' and r['panel']=='b')]
 p=S/'Figure3/Fig3b_L2_selected19.tsv';d=read(p);b=h.svgaxes('Figure3','figure3_b_content');tx=[''.join(t.itertext()) for t in b.iter('{http://www.w3.org/2000/svg}text')];labs=tx[2:];assert len(labs)==19
 f=fig(140,110);a=axmm(f,85,13,34,80);a.set(xlim=(-.5,1.5),ylim=(18.5,-.5));a.set_xticks([0,1],['MAGMA','scDRS']);a.xaxis.tick_top();a.set_yticks(range(19),labs);a.tick_params(length=0,pad=5);a.spines[:].set_visible(False)
 sizes=[10,23,36,49]
 for i,r in d.iterrows():
  for j,key in enumerate(['FDR_BH','primary_r0.1__assoc_mcp_BH']):
   v=r[key];cat=0 if v>=.05 else 1 if v>=.01 else 2 if v>=.001 else 3;a.scatter(j,i,s=sizes[cat],c=L1[r['group'].split(' | ')[0]],ec='none')
 handles=[a.scatter([],[],s=sz,c=GREY,ec='none',label=la) for sz,la in zip(sizes,['FDR ≥ 0.05','0.01–<0.05','0.001–<0.01','<0.001'])];f.legend(handles=handles,loc='lower center',ncol=2,frameon=False,fontsize=7,labelspacing=.8)
 save(f,'Figure3b',[p,P/'图件编辑资料/editable_svg/Figure3.svg'],'Separate method columns with 17 mm centres; same selected 19 states and four original FDR area bins',{'rows':19,'points':38,'max_diameter_pt':7,'row_pitch_pt':80/19*72/25.4})

def pathways():
 b=S/'Current_Display_20260920/Figure3c_current';p=b/'MAGMA_display9.tsv';q=b/'pseudobulk_display9_all7_lineages.tsv';lp=b/'pathway_display_labels.tsv';m=read(p);e=read(q);lab=read(lp).set_index('pathway_id').display_label.to_dict()
 paths=e.pathway.drop_duplicates().tolist();states=e.cell_state.drop_duplicates().tolist();assert len(paths)==9 and len(states)==7
 f=fig(164,153);a=axmm(f,80,8,59,53.1);pitch=5.9
 for i,r in m.iterrows():
  a.plot([0,r.NES],[i,i],c=NAVY,lw=1.15);a.scatter(r.NES,i,s=27,c=NAVY,ec='none',zorder=3)
  f.text(145/164,1-(8+(i+.5)*pitch)/153,f'{r.padj:.3f}',va='center',fontsize=8)
 a.set_yticks(range(9),labels_wrap([lab[x] for x in m.pathway],48));a.set(ylim=(8.5,-.5),xlim=(0,2.15),xlabel='NES');clean(a);a.grid(axis='x',lw=.4);a.set_axisbelow(True)
 f.text(145/164,1-4/153,'Adjusted P',fontsize=8);f.text(80/164,1-2/153,'MAGMA GSEA',fontsize=8)
 v=e.pivot(index='pathway',columns='cell_state',values='NES').reindex(index=paths,columns=states);sig=e.pivot(index='pathway',columns='cell_state',values='padj').reindex(index=paths,columns=states);vm=float(e.NES.abs().max())
 a=axmm(f,80,86,59,49.5);im=matrix_image(a,v,cmap='RdBu_r',vmin=-vm,vmax=vm,center=0,value_space='Saved pseudobulk DL GSEA NES')
 yy,xx=np.where(sig.to_numpy()<.05);a.scatter(xx,yy,s=7,c='#20252A',ec='none');a.set_yticks(range(9),labels_wrap([lab[x] for x in paths],48));a.set_xticks(range(7),[{'Endothelial cell':'EC','Fibroblast':'Fib','Smooth muscle cell':'SMC','Macrophage':'Mac','Monocyte':'Mo','T cell':'T','Natural killer cell':'NK'}[s] for s in states]);a.xaxis.tick_top();a.tick_params(length=0,pad=4);a.spines[:].set_visible(False)
 f.text(80/164,1-76/153,'Pseudobulk GSEA',fontsize=8)
 for top,rh in [(8,pitch),(86,5.5)]:
  for j,label in enumerate(['GO BP','KEGG','Reactome']):
   f.text(15/164,1-(top+(j*3+1.5)*rh)/153,label,rotation=90,ha='center',va='center',fontsize=7.5)
   f.add_artist(Line2D([18/164]*2,[1-(top+j*3*rh+.6)/153,1-(top+(j*3+3)*rh-.6)/153],transform=f.transFigure,c=GREY,lw=.9))
 cb=f.colorbar(im,cax=axmm(f,145,91,2.6,25));cb.set_label('NES');cb.outline.set_visible(False);f.text(143/164,1-128/153,'● Adjusted\n   P < 0.05',fontsize=7)
 save(f,'Figure3c',[p,q,lp],'Horizontal interval/point enrichment recipe + intact NES matrix; fixed 9+9 pathways; explicit adjusted P column',{'MAGMA_rows':9,'PB_cells':63,'PB_color_limits':[-vm,vm]})

def forests_and_evidence():
 h.forests()
 # A shared physical row grid is declared for the forest and its evidence.
 # Standalone panels retain full names. Review composition uses the combined
 # view below, with only one copy of row labels and identical page y centers.
 p=P/'图件编辑资料/editable_svg/Figure3.svg';r=E.parse(str(p));g=r.xpath('//*[@id="FigureYa_vector_matrix_Figure3e"]')[0];els=[e for e in g if e.get('data-row')];assert len(els)==145
 d=read(S/'Figure3/Fig3d_e_all29.tsv');labs=[r.gene_symbol+'  '+{'Fibroblast':'Fib','Macrophage':'Mac','Endothelial cell':'EC'}[r.cell_state] for _,r in d.iterrows()]
 cells=pd.DataFrame([{'row':int(el.get('data-row'))-1,'col':int(el.get('data-col'))-1,'source_value':el.get('data-source-value'),'fill':el.get('fill')} for el in els]);cells.to_csv(ROOT/'data_used/Figure3e_145_frozen_cells.tsv',sep='\t',index=False)
 def draw_e(a,with_names):
  for _,v in cells.iterrows():a.add_patch(Rectangle((v.col-.41,v.row-.38),.82,.76,fc=v.fill,ec='none'))
  a.set(xlim=(-.5,4.5),ylim=(28.5,-.5));a.set_xticks(range(5),['Primary\nGWAS','Primary\nPB','Leave-\nCAD','Regional\nABF','Spatial']);a.xaxis.tick_top();a.set_yticks(range(29),labs if with_names else []);a.tick_params(length=0,pad=4);a.spines[:].set_visible(False)
 handles=[Patch(fc=C['evidence'][key],label=lab) for key,lab in [('support','Support'),('inverse','Inverse'),('no_support','No support')]]
 f=fig(106,137);a=axmm(f,35,22,60,102);draw_e(a,True);f.legend(handles=handles,ncol=3,loc='lower center',frameon=False,fontsize=7)
 save(f,'Figure3e',[p,S/'Figure3/Fig3d_e_all29.tsv'],'Full 145-cell evidence matrix; same row pitch as the forest, separate category key',{'rows':29,'cells':145,'row_pitch_mm':102/29})
 f=fig(184,137);a=axmm(f,34,16,77,102);b=axmm(f,115,16,60,102)
 for i,v in d.iterrows():a.plot([v.meta_CI_low,v.meta_CI_high],[i,i],c=L1[v.cell_state],lw=.85);a.scatter(v.meta_logFC,i,s=17,c=L1[v.cell_state],ec='none',zorder=3)
 a.set_yticks(range(29),labs);a.set(ylim=(28.5,-.5),xlabel='Meta log2 fold change (95% CI)');a.axvline(0,c=GREY,lw=.55);clean(a);draw_e(b,False)
 f.text(.015,.98,'d',va='top',weight='bold',fontsize=10);f.text(115/184,.98,'e',va='top',weight='bold',fontsize=10)
 f.legend(handles=handles,ncol=3,loc='lower right',frameon=False,fontsize=7)
 f.canvas.draw();ya=a.transData.transform(np.column_stack([np.zeros(29),np.arange(29)]))[:,1];yb=b.transData.transform(np.column_stack([np.zeros(29),np.arange(29)]))[:,1];assert np.array_equal(ya,yb)
 for ext in ['pdf','svg']:f.savefig(ROOT/f'review/Figure3de_shared_rows.{ext}')
 plt.close(f);doc=fitz.open(ROOT/'review/Figure3de_shared_rows.pdf');doc[0].get_pixmap(matrix=fitz.Matrix(2,2)).save(ROOT/'review/Figure3de_shared_rows.png')
 (ROOT/'review/shared_rows.json').write_text(json.dumps({'rows':29,'pitch_mm':102/29,'data_top_mm':16,'data_height_mm':102,'data_gap_mm':4,'max_y_delta_pixels':float(np.abs(ya-yb).max()),'checked_genes':[d.gene_symbol.iloc[0],d.gene_symbol.iloc[17],d.gene_symbol.iloc[-1]]},indent=2))
 import redraw_evidence_panels as ep
 ep.save=save;ep.coloc()

def modules():
 h.modules()
 p=S/'Figure4/14_final_module_evidence_and_tiers.tsv';d=read(p).set_index('module');keys=['adjacent_rank_median','leave_one_study','leave_one_tissue','condition_only','hierarchical_study_bootstrap','transformation'];labs=['Adjacent rank','Leave one study','Leave one tissue','Condition only','Study / donor bootstrap','Transformation']
 f=fig(122,71);a=axmm(f,53,10,59,42)
 for i,(key,lab) in enumerate(zip(keys,labs)):
  vals=[d.loc[x,key] for x in ['CM09','CM03']];a.plot(vals,[i-.1,i+.1],color='#CAD1D5',lw=1.2,zorder=1)
  for module,offset in [('CM09',-.1),('CM03',.1)]:a.scatter(d.loc[module,key],i+offset,c=C['module'][module],s=30,ec='none',label=module if i==0 else None,zorder=3)
 a.set_yticks(range(6),labs);a.set(ylim=(5.5,-.5),xlim=(0,1.04),xlabel='Recorded matching similarity');a.set_xticks([0,.25,.5,.75,1]);clean(a);a.grid(axis='x',lw=.4);f.legend(*a.get_legend_handles_labels(),loc='lower center',ncol=2,frameon=False)
 save(f,'Figure4e',[p],'Two-module perturbation comparison using existing six metrics; same 0–1 scale; all values retained')

def spatial_quant():
 p=S/'Figure5_final_companion/Figure5b_60_gene_unit_estimates.tsv';q=S/'Figure5_final_companion/Figure5b_pooled_FDR_I2_original.tsv';d=read(p);meta=read(q);assert len(d)==60 and d.biological_unit_id.nunique()==34
 studies=d.GSE.drop_duplicates().tolist();f=fig(143,99);a=axmm(f,23,12,80,68);cols={'FN1':C['evidence']['inverse'],'SCARB1':C['evidence']['support']};markers={'FN1':'o','SCARB1':'^'}
 for i,st in enumerate(studies):
  for gene,offset in [('FN1',-.14),('SCARB1',.14)]:
   z=d[d.GSE.eq(st)&d.gene.eq(gene)].sort_values('biological_unit_id');ys=i+offset+np.linspace(-.07,.07,len(z));a.scatter(z.beta,ys,s=18,c=cols[gene],marker=markers[gene],ec='white',lw=.25,zorder=3)
 a.axhline(6.6,color=LIGHT,lw=.65)
 for j,gene in enumerate(['FN1','SCARB1']):
  r=meta[meta.gene.eq(gene)].iloc[0];y=7.25+j;a.plot([r.CI_low,r.CI_high],[y,y],color=cols[gene],lw=1.25);a.scatter(r.beta,y,s=40,c=cols[gene],marker='D',ec='none',zorder=3)
  yy=1-(12+(y+.5)/9.1*68)/99
  f.text(110/143,yy,f'{r.FDR_BH_all_29_pairs_x_3_terms:.3g}',va='center',fontsize=8);f.text(128/143,yy,f'{r.I2:.1f}%',va='center',fontsize=8)
 a.axvline(0,c=GREY,lw=.7);a.set(ylim=(8.6,-.5),xlim=(-1.02,.55),xlabel='Conditional Fib-abundance coefficient');a.set_yticks(list(range(7))+[7.25,8.25],studies+['FN1 pooled','SCARB1 pooled']);a.set_xticks([-1,-.5,0,.5]);clean(a)
 hy=1-(12+(6.85+.5)/9.1*68)/99;f.text(110/143,hy,'FDR',fontsize=8);f.text(128/143,hy,'I²',fontsize=8)
 handles=[Line2D([],[],marker=markers[g],color=cols[g],ls='',label=g) for g in ['FN1','SCARB1']]
 f.legend(handles=handles,loc='upper center',ncol=2,frameon=False,bbox_to_anchor=(.48,.98));f.legend(handles=[Line2D([],[],marker='D',c=GREY,ms=5,label='Pooled estimate and 95% CI')],loc='lower center',frameon=False)
 save(f,'Figure5b',[p,q],'Source-stratified 60 gene–unit estimates; pooled FDR/I2 in independent columns; full opposite-direction estimates retained',{'estimates':60,'recorded_unit_IDs':34,'unit_SE_not_drawn_as_CI':True})

def spatial_modules():
 p=S/'Figure5_final_companion/all29_original_meta.tsv';m=read(p)
 f=fig(122,62);a=axmm(f,57,11,60,34);rows=[('mean','all_sections'),('mean','exclude_QC_fallback'),('q05','all_sections'),('q05','exclude_QC_fallback')];counts=[]
 for i,(post,qc) in enumerate(rows):
  z=m[m.posterior.eq(post)&m.sensitivity.eq(qc)&m.term.eq('cell_abundance')];assert len(z)==29
  sig=z.FDR_BH_all_29_pairs_x_3_terms.lt(.05);v=[int((sig&z.beta.gt(0)).sum()),int((sig&z.beta.lt(0)).sum()),int((~sig).sum())];counts.append(v);left=0
  for value,key in zip(v,['support','inverse','no_support']):
   a.barh(i,value,left=left,height=.56,color=C['evidence'][key],ec='none');a.text(left+value/2,i,str(value),ha='center',va='center',fontsize=8);left+=value
 a.set(xlim=(0,29),ylim=(3.5,-.5),xlabel='Candidates');a.set_xticks([0,10,20,29]);a.set_yticks([]);clean(a)
 for i,(post,qc) in enumerate(rows):
  y=1-(11+34*(i+.5)/4)/62;f.text(4/122,y,'Mean' if post=='mean' else 'q05',va='center',fontsize=8);f.text(23/122,y,'All sections' if qc=='all_sections' else 'Exclude QC-\nfallback',va='center',fontsize=8)
 f.text(4/122,1-5/62,'Posterior',fontsize=8);f.text(23/122,1-5/62,'Section set',fontsize=8)
 f.legend(handles=[Patch(fc=C['evidence'][k],label=l) for k,l in [('support','Positive'),('inverse','Inverse'),('no_support','No support')]],loc='lower center',ncol=3,frameon=False)
 assert counts==[[6,3,20],[7,4,18],[6,4,19],[6,4,19]]
 save(f,'Figure5d',[p],'Read existing FDR classifications, retain four 29-candidate totals; row fields and count scale aligned',{'counts':counts,'new_correction':False})
 p=S/'Figure5/module_spatial/module_coherence_meta.tsv';d=read(p);z=d[d.posterior.eq('q05')&d.null_method.eq('strict_parent')];x=z[z.metric.eq('within_spot')].set_index('module');y=z[z.metric.eq('spatial_lag')].set_index('module');assert len(x)==len(y)==14
 f=fig(104,86);a=axmm(f,17,12,76,56)
 offsets={'CM01':(-7,6),'CM02':(-9,-7),'CM03':(5,6),'CM04':(7,4),'CM05':(-8,-7),'CM06':(-12,7),'CM07':(7,-6),'CM08':(7,5),'CM09':(6,5),'CM10':(-10,7),'CM11':(-10,5),'CM12':(-10,7),'CM13':(5,-7),'CM14':(7,-7)}
 for module,r in x.iterrows():
  xx,yy=r.observed_median,y.loc[module,'observed_median'];sig=(z[z.module.eq(module)].FDR_BH_28<.05).any()
  a.scatter(xx,yy,c=C['module'].get(module,C['module']['Other']),marker='s' if sig else 'o',s=35,ec='none',zorder=3)
  a.annotate(module[-2:],(xx,yy),xytext=offsets[module],textcoords='offset points',fontsize=8,va='center',arrowprops=dict(arrowstyle='-',color='#B7BFC5',lw=.45))
 a.set(xlabel='Within-spot coherence',ylabel='Spatial-lag coherence',xlim=(.25,.92),ylim=(.08,.50));a.spines[['top','right']].set_visible(False)
 f.legend(handles=[Line2D([],[],marker='s',ls='',c=GREY,label='Any metric FDR < 0.05'),Line2D([],[],marker='o',ls='',c=GREY,label='Other')],loc='upper center',ncol=2,frameon=False,fontsize=7)
 save(f,'Figure5e',[p],'14 original module coordinates; leadered labels; original any-metric FDR shape and module colors',{'module_count':14,'coordinate_fields':'observed_median of within_spot and spatial_lag, q05 strict_parent','shape_source':'reconstruct_spatial.py:119, original any-metric FDR<.05'})
 p=S/'Figure5/module_spatial/module_dataset_stats.tsv';q=S/'Figure5/module_spatial/module_program_direction_meta.tsv';d=read(p);m=read(q);pairs=[('CM01','fibrous_cap_core','fibrous cap'),('CM01','calcification','calcification'),('CM02','inflammation','inflammation'),('CM03','endothelial_activation','EC activation'),('CM09','inflammation','inflammation')]
 f=fig(139,70);a=axmm(f,53,14,66,40)
 for i,(mod,prog,label) in enumerate(pairs):
  z=d[d.posterior.eq('q05')&d.module.eq(mod)].sort_values('GSE');assert len(z)==7
  a.scatter(z['rho__'+prog],i+np.linspace(-.12,.12,7),s=20,c='#99A1A6',ec='white',lw=.25,zorder=2)
  r=m[m.posterior.eq('q05')&m.module.eq(mod)&m.program.eq(prog)].iloc[0];a.scatter(r.median_rho,i,s=45,c=C['module'].get(mod,C['module']['Other']),marker='D',ec='none',zorder=3)
  f.text(124/139,1-(14+40*(i+.5)/5)/70,f'{r.FDR_BH_84:.2f}',va='center',fontsize=8)
 a.axvline(0,color=GREY,lw=.6);a.set(ylim=(4.5,-.5),xlim=(-.58,.58),xlabel='Spearman ρ');a.set_xticks([-.5,0,.5]);a.set_yticks(range(5),[f'{mod} | {label}' for mod,pr,label in pairs]);clean(a)
 f.text(124/139,1-10/70,'FDR',fontsize=8)
 f.legend(handles=[Line2D([],[],marker='o',ls='',c='#99A1A6',label='Dataset'),Line2D([],[],marker='D',ls='',c=GREY,label='Median')],loc='upper center',ncol=2,frameon=False,bbox_to_anchor=(.60,.98),fontsize=8)
 save(f,'Figure5f',[p,q],'Complete 35 dataset estimates and 5 stored medians; labels, keys and FDR in independent columns',{'dataset_points':35,'stored_medians':5,'new_meta':False})

def inherited_panels():
 # Retain reliable native point/image layers, coordinates and annotations. The
 # independent source SVGs already contain complete context and tight margins.
 for n,letters in [(1,'be'),(2,'d'),(5,'a')]:
  for letter in letters:
   name=f'Figure{n}{letter}';src=P/f'独立小图_最新版/Figure{n}'
   for ext in ['pdf','svg','png']:shutil.copy2(src/f'{name}.{ext}',ROOT/f'Figure{n}'/f'{name}.{ext}')
   doc=fitz.open(ROOT/f'Figure{n}/{name}.pdf')
   RECORDS.append({'figure':f'Figure{n}','panel':letter,'mode':'preserved_native_geometry','source_files':[str(src/f'{name}.svg')],'change':'Existing frozen geometry, native image/point layer and all statistical encodings retained; full context preserved','width_mm':doc[0].rect.width/72*25.4,'height_mm':doc[0].rect.height/72*25.4,'details':{}})

def main():
 try:
  schematic(1);schematic(2);schematic(4);marker();summaries();genetic_matrix();sem();venn();upset();genetics();pathways();forests_and_evidence();modules();spatial_quant();spatial_modules();inherited_panels()
 finally:
  (ROOT/'review/panel_manifest.json').write_text(json.dumps(RECORDS,ensure_ascii=False,indent=2));(ROOT/'review/source_hashes.json').write_text(json.dumps(h.sources,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
