"""Draw standalone panels from existing, adopted summary statistics only.
No expression processing, tests, fitting, ranking or marker selection is run.
"""
from pathlib import Path
import json,hashlib,textwrap,re,copy
import numpy as np,pandas as pd
from lxml import etree as E
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle,Patch
import pymupdf as fitz
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');O=P/'独立小图_最新版';S=P/'Source_Data_Figures'
C=json.loads((P/'图件编辑资料/semantic_colors_current.json').read_text());L1=C['L1'];mm=1/25.4
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':7,'xtick.labelsize':7,'ytick.labelsize':7,'svg.fonttype':'none','pdf.fonttype':42,'text.color':'#20252A','axes.labelcolor':'#20252A','axes.edgecolor':'#647078','xtick.color':'#20252A','ytick.color':'#20252A','axes.linewidth':.55,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
BL=LinearSegmentedColormap.from_list('adopted_marker',plt.colormaps['Blues'](np.linspace(.15,.95,256)),N=256)
records=[];sources={}
def read(f):
 f=Path(f);sources[str(f)]=hashlib.sha256(f.read_bytes()).hexdigest();return pd.read_csv(f,sep='\t',float_precision='round_trip')
def axmm(f,x,y,w,h):
 W,H=f.get_size_inches()*25.4;return f.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def fig(w,h):return plt.figure(figsize=(w*mm,h*mm))
def clean(a):a.spines[['top','right','left']].set_visible(False);a.tick_params(axis='y',length=0)
def save(f,name,inputs,changes,details=None):
 root=re.sub('[a-z]$','',name);letter=name[-1];folder=O/root;folder.mkdir(exist_ok=True)
 # Anchor the letter to the rendered content, not an unused canvas corner.
 # Only the outer frame changes; shared-axis geometry remains untouched.
 f.canvas.draw();bb=f.get_tightbbox(f.canvas.get_renderer());W,H=f.get_size_inches()
 f.text(bb.x0/W,(bb.y1+2.2/25.4)/H,letter,ha='left',va='bottom',weight='bold',fontsize=10)
 for ext in ['pdf','svg']:f.savefig(folder/f'{name}.{ext}',bbox_inches='tight',pad_inches=1.5/25.4)
 plt.close(f);d=fitz.open(folder/f'{name}.pdf');d[0].get_pixmap(matrix=fitz.Matrix(2,2)).save(folder/f'{name}.png')
 records.append({'figure':root,'panel':letter,'mode':'redrawn_from_existing_summary','source_files':[str(x) for x in inputs],'change':changes,'width_mm':d[0].rect.width/72*25.4,'height_mm':d[0].rect.height/72*25.4,'details':details or {}})
 print('redrawn',name,flush=True)

def svgaxes(n,id):
 r=E.parse(str(P/'图件编辑资料/editable_svg'/f'{n}.svg'));return r.xpath('//*[@id=$v]',v=id)[0]
def axis_labels(n,id):
 a=svgaxes(n,id);axes=[g for g in a if 'matplotlib.axis' in g.get('id','')]
 return [[''.join(t.itertext()) for t in g.iter('{http://www.w3.org/2000/svg}text')] for g in axes]

def marker(name,d,states,genes,labels,vmax,scaled=False):
 nr,nc=len(states),len(genes);rh=3.25;cw=3.15;left=67 if nr>20 else 58;right=28;top=10;bottom=30
 w=left+nc*cw+right;h=top+nr*rh+bottom;f=fig(w,h);a=axmm(f,left,top,nc*cw,nr*rh)
 mean=d.pivot(index='state_id',columns='gene_id',values='color_value').reindex(index=states,columns=genes)
 det=d.pivot(index='state_id',columns='gene_id',values='detection_fraction').reindex(index=states,columns=genes)
 assert mean.shape==det.shape==(nr,nc);assert ((det.dropna()>=0)&(det.dropna()<=1)).all().all()
 X,Y=np.meshgrid(np.arange(nc),np.arange(nr));valid=np.isfinite(mean)&np.isfinite(det)&(det>0)
 diameter=.72*min(rh,cw)*72/25.4;smax=diameter**2
 a.scatter(X[valid],Y[valid],s=(smax*det.to_numpy())[valid],c=mean.to_numpy()[valid],cmap=BL,norm=Normalize(0,vmax),linewidths=0,edgecolors='none')
 missing=~np.isfinite(mean)|~np.isfinite(det)
 if missing.to_numpy().any():a.scatter(X[missing],Y[missing],marker='x',s=9,c='#747C86',linewidths=.6)
 a.set(xlim=(-.5,nc-.5),ylim=(nr-.5,-.5));a.set_xticks(range(nc),genes,rotation=90,fontstyle='italic');a.set_yticks(range(nr),labels);a.tick_params(length=0,pad=4);a.tick_params(axis='y',pad=9)
 for sp in a.spines.values():sp.set_visible(False)
 parents=[]
 index=read(S/'ExtendedData3_4_5/ED345_all_named_state_index.tsv').set_index('level2_annotation')
 for i,st in enumerate(states):
  parent=st if st in L1 else index.loc[st,'level1_annotation'];parents.append(parent)
  co=L1[st] if st in L1 else C['L2'][C['L2_name_to_id'][st]]
  a.add_patch(Rectangle((-.77,i-.43),.12,.86,fc=co,ec='none',clip_on=False))
 for i in range(1,nr):
  if nr>20 and parents[i]!=parents[i-1]:a.axhline(i-.5,color='#E2E6E8',lw=.55,zorder=0)
 # Group bars remain distinct from individual state colours.
 if nr>20:
  for pa in dict.fromkeys(parents):
   ix=np.where(np.array(parents)==pa)[0];y0=top+rh*ix.min();height=rh*len(ix)
   bar=axmm(f,8,y0,1,height-.4);bar.set_facecolor(L1[pa]);bar.set_xticks([]);bar.set_yticks([]);[sp.set_visible(False) for sp in bar.spines.values()]
   f.text(5/w,1-(y0+height/2)/h,pa.replace(' cell',''),rotation=90,ha='center',va='center',fontsize=7)
 cax=axmm(f,left+nc*cw+7,top+5,2.3,min(nr*rh-9,25));cb=f.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,vmax),cmap=BL),cax=cax);cb.outline.set_visible(False);cb.ax.tick_params(length=2)
 cb.set_label('Per-gene scaled mean expression' if scaled else 'Mean log1p(CP10k)',labelpad=4)
 la=axmm(f,(w-85)/2,top+nr*rh+20,85,8);la.set_xlim(0,100);la.set_ylim(0,10);la.axis('off');la.text(0,9,'Detected cells (%)',va='center',fontsize=7)
 for j,v in enumerate([.25,.5,.75,1]):
  x=48+j*14;la.scatter(x,9,s=smax*v,c='#747C86',edgecolors='none');la.text(x,2,str(int(v*100)),ha='center',va='center',fontsize=7)
 if missing.to_numpy().any():la.text(0,1,'× Not measured',fontsize=7)
 return f,{'rows':nr,'columns':nc,'cells':nr*nc,'vmin':0,'vmax':vmax,'maximum_diameter_pt':diameter,'row_spacing_pt':rh*72/25.4,'column_spacing_pt':cw*72/25.4,'zero_detection_count':int((det==0).sum().sum()),'NA_cells':int(missing.sum().sum()),'area_mapping':'s = smax * original_detection_fraction; zero = 0','state_order':states,'gene_order':genes}

def markers():
 mf=P/'Public_Release/Figure_Code/native_dependencies/Panvascular_Marker_Communication_Enrichment_20260920/source_data/display_marker_values.tsv';M=read(mf)
 for n,letter,ai in [(2,'b',1),(3,'b',1),(4,'b',1),(4,'d',4),(5,'b',1)]:
  genes,labs=axis_labels(f'ExtendedData_Figure{n}',f'newED{n}_axes_{ai}')
  d=M[M.display_label.isin(labs)&M.gene_id.isin(genes)].copy();assert len(d)==len(labs)*len(genes)
  look=d.drop_duplicates('display_label').set_index('display_label').state_id.to_dict();states=[look[l] for l in labs];caps=d.vmax.unique();assert len(caps)==1
  d['color_value']=d.mean_log1p
  f,checks=marker(f'ExtendedData_Figure{n}{letter}',d,states,genes,labs,float(caps[0]));save(f,f'ExtendedData_Figure{n}{letter}',[mf],'fixed complete matrix; 3.25 mm row pitch; separate names and state stripe; borderless proportional circles; compact independent keys',checks)
 zf=S/'Current_Display_20260920/marker_audit_20260921/per_gene_scaled_mean.tsv';df=S/'Current_Display_20260920/marker_audit_20260921/original_detection_fraction.tsv'
 Z=read(zf).set_index('level1_annotation');D=read(df).set_index('level1_annotation');rows=[]
 for st in Z.index:
  for g in Z.columns:rows.append({'state_id':st,'gene_id':g,'color_value':Z.loc[st,g],'detection_fraction':D.loc[st,g]})
 f,ch=marker('Figure1g',pd.DataFrame(rows),list(Z.index),list(Z.columns),[s.replace(' cell','') for s in Z.index],1,True);save(f,'Figure1g',[zf,df],'adopted B values read directly; no expression or min-max recalculation; row/column pitch optimized',ch)

def lineage():
 path=S/'Figure1/L1_palette_bound.tsv';p=read(path).sort_values('order');states=list(p.L1)
 # Exact counts and study totals printed in adopted axes; no million-cell recomputation.
 r=E.parse(str(P/'图件编辑资料/editable_svg/Figure1.svg'));texts=lambda id:[''.join(t.itertext()) for t in r.xpath('//*[@id=$v]',v=id)[0].iter('{http://www.w3.org/2000/svg}text')]
 counts=np.array([int(t.replace(',','')) for t in texts('axes_5')[-16:]]);studies=np.array([int(t) for t in texts('axes_6')[-16:]])
 for letter,values,xlabel in [('c',counts,'Cells per lineage'),('f',studies,'Studies per lineage')]:
  f=fig(110,85);a=axmm(f,35,8,65,65);a.barh(range(16),values,height=.66,color=[L1[s] for s in states]);a.set_yticks(range(16),[s.replace(' cell','') for s in states]);a.invert_yaxis();a.set_xlabel(xlabel);a.set_xlim(0,max(values)*1.19);a.set_axisbelow(True);a.grid(axis='x',color='#E5E8EA',lw=.4);clean(a)
  for i,v in enumerate(values):a.text(v+max(values)*.015,i,f'{v:,}',va='center',fontsize=7)
  if letter=='c':a.set_xticks([0,100000,200000],['0','100k','200k'])
  else:a.set_xticks([0,13,26])
  save(f,'Figure1'+letter,[P/'图件编辑资料/editable_svg/Figure1.svg',path],'complete names restored on standalone panel; aligned values; fixed lineage order',{'values':values.tolist(),'state_order':states})
 path=S/'Figure1/L1_tissue_sample_mean_fraction.tsv';d=read(path);print('tissuecols',d.columns.tolist())

def genetics():
 p=S/'Figure3/Fig3a_L1.tsv';d=read(p);f=fig(95,84);a=axmm(f,18,10,65,53)
 for _,v in d.iterrows():
  sig=v.FDR_BH<.05 and v['primary_r0.1__assoc_mcp_BH']<.05;co=L1[v.group]
  a.scatter(v.BETA_STD,v['primary_r0.1__assoc_mcz'],s=23 if not sig else 32,facecolor=co if sig else 'white',edgecolor=co,lw=.8,zorder=3)
  if sig or v.group=='Smooth muscle cell':
   lab={'Endothelial cell':'Endothelial','Smooth muscle cell':'Smooth muscle'}.get(v.group,v.group)
   off={'Pericyte':(-6,2),'Fibromyocyte':(-8,5),'Smooth muscle':(7,7),'Endothelial':(-5,-9)}[lab]
   a.annotate(lab,(v.BETA_STD,v['primary_r0.1__assoc_mcz']),xytext=off,textcoords='offset points',ha='left' if lab=='Smooth muscle' else 'right',fontsize=7,arrowprops=dict(arrowstyle='-',lw=.4,color='#747C86'))
 a.axhline(0,c='#D9DDE1',lw=.5);a.axvline(0,c='#D9DDE1',lw=.5);a.set(xlim=(-.0215,.04),ylim=(-3.3,6.85),xlabel='MAGMA standardised beta',ylabel='scDRS association Z');a.set_xticks([-.02,0,.02,.04]);a.set_yticks([-2,0,2,4,6]);f.legend(handles=[Line2D([],[],ls='',marker='o',mfc='#747C86',mec='#747C86',label='Both FDR < 0.05'),Line2D([],[],ls='',marker='o',mfc='white',mec='#747C86',label='Not both FDR < 0.05')],loc='lower center',ncol=1,frameon=False,fontsize=7)
 save(f,'Figure3a',[p],'larger data area; explicit hollow/filled key; 15 unchanged coordinates',{'points':15})
 p=S/'Figure3/Fig3b_L2_selected19.tsv';d=read(p);f=fig(126,104);a=axmm(f,85,12,24,76);a.set(xlim=(-.65,1.65),ylim=(18.5,-.5));a.set_xticks([0,1],['MAGMA','scDRS']);a.xaxis.tick_top();a.tick_params(length=0,pad=5)
 labs=axis_labels('Figure3','figure3_b_content') if False else None
 # Labels from approved SVG are a display dictionary, not new state names.
 b=svgaxes('Figure3','figure3_b_content');tx=[''.join(t.itertext()) for t in b.iter('{http://www.w3.org/2000/svg}text')];labs=tx[2:];assert len(labs)==19
 a.set_yticks(range(19),labs)
 sizes=[10,23,36,49]
 for i,r in d.iterrows():
  for j,key in enumerate(['FDR_BH','primary_r0.1__assoc_mcp_BH']):
   v=r[key];cat=0 if v>=.05 else 1 if v>=.01 else 2 if v>=.001 else 3
   a.scatter(j,i,s=sizes[cat],color=L1[r['group'].split(' | ')[0]],edgecolors='none')
 for s in a.spines.values():s.set_visible(False)
 # Explicit area categories kept; diameter7pt < row11.3pt.
 handles=[plt.scatter([],[],s=sz,c='#747C86',edgecolors='none',label=la) for sz,la in zip(sizes,['FDR ≥ 0.05','0.01–<0.05','0.001–<0.01','<0.001'])];f.legend(handles=handles,loc='lower center',ncol=2,frameon=False,labelspacing=.8)
 save(f,'Figure3b',[p,P/'图件编辑资料/editable_svg/Figure3.svg'],'19 rows with adequate pitch; names close to method columns; categorical FDR area key',{'rows':19,'point_diameter_pt_max':7,'row_pitch_pt':76/19*72/25.4})

def forests():
 p=S/'Figure3/Fig3d_e_all29.tsv';d=read(p);f=fig(112,127);a=axmm(f,40,10,60,102)
 for i,r in d.iterrows():
  a.plot([r.meta_CI_low,r.meta_CI_high],[i,i],c=L1[r.cell_state],lw=.8);a.scatter(r.meta_logFC,i,s=16,c=L1[r.cell_state],zorder=3)
 a.set_yticks(range(29),[r.gene_symbol+'  '+{'Fibroblast':'Fib','Macrophage':'Mac','Endothelial cell':'EC'}[r.cell_state] for _,r in d.iterrows()]);a.set(ylim=(28.5,-.5),xlim=(float(d.meta_CI_low.min())-.12,float(d.meta_CI_high.max())+.12),xlabel='Meta log2 fold change (95% CI)');a.axvline(0,c='#AAB1B7',lw=.55);clean(a);save(f,'Figure3d',[p],'29 aligned rows; labels adjacent to forest; full stored CI endpoints',{'rows':29})
 p2=S/'Figure5_final_companion/Fig5c_all29.tsv';s=read(p2).set_index('gene').loc[d.gene_symbol].reset_index();f=fig(118,128);a=axmm(f,35,10,61,102)
 for i,r in s.iterrows():
  q=r.FDR_BH_all_29_pairs_x_3_terms;co=C['evidence']['support'] if q<.05 and r.beta>0 else C['evidence']['inverse'] if q<.05 else C['NO_SUPPORT_POINT'];a.plot([r.CI_low,r.CI_high],[i,i],c=co,lw=.85);a.scatter(r.beta,i,s=16,c=co,zorder=3);f.text(100/118,1-(10+102*(i+.5)/29)/128,('<0.001' if q<.001 else f'{q:.3f}'),va='center',fontsize=7)
 f.text(100/118,1-7/128,'FDR',fontsize=7);a.set_yticks(range(29),[r.gene+'  '+{'Fibroblast':'Fib','Macrophage':'Mac','Endothelial cell':'EC'}[r.cell_context] for _,r in s.iterrows()]);a.set(ylim=(28.5,-.5),xlim=(float(s.CI_low.min())-.05,float(s.CI_high.max())+.05),xlabel='Conditional coefficient (95% CI)');a.axvline(0,c='#AAB1B7',lw=.55);clean(a)
 f.legend(handles=[Line2D([],[],marker='o',ls='-',ms=3,c=c,label=l) for c,l in [(C['evidence']['support'],'Positive'),(C['evidence']['inverse'],'Inverse'),(C['NO_SUPPORT_POINT'],'No support')]],loc='lower center',ncol=3,frameon=False,fontsize=7)
 save(f,'Figure5c',[p,p2],'full29 independent forest with row labels and FDR; opaque gray data; compact support key',{'rows':29,'CI':'stored endpoints, not recalculated'})

def modules():
 p=S/'Figure4/14_final_module_evidence_and_tiers.tsv';d=read(p);co=lambda m:C['module'].get(m,C['module']['Other'])
 f=fig(160,40);a=axmm(f,39,5,112,20);mat=d[['stable_structure','broad_context','coherent_within_study']].to_numpy().T.astype(int);cm=LinearSegmentedColormap.from_list('evidence',[C['NO_SUPPORT_FILL'],C['evidence']['support']],N=2);a.imshow(mat,cmap=cm,vmin=0,vmax=1,aspect='auto',interpolation='nearest');a.set_xticks(np.arange(-.5,14,1),minor=True);a.set_yticks(np.arange(-.5,3,1),minor=True);a.grid(which='minor',color='white',linewidth=.8);a.tick_params(which='minor',length=0);a.set_xticks(range(14),d.module,fontsize=7);a.set_yticks(range(3),['Structure','Context breadth','Within-study coherence']);a.tick_params(length=0);[sp.set_visible(False) for sp in a.spines.values()];f.legend(handles=[Patch(fc=C['evidence']['support'],label='Criterion met'),Patch(fc=C['NO_SUPPORT_FILL'],label='Not met')],loc='lower center',bbox_to_anchor=(.595,0),ncol=2,frameon=False);save(f,'Figure4b',[p],'matrix width allocated per module; key centered below data; original Boolean gates')
 p2=S/'Figure4/05_top_nodes_reference.tsv';n=read(p2);p3=S/'Figure4/FINAL_STATE_IDENTITY_MAP.tsv';idx=read(p3);print('identitycols',idx.columns.tolist())
 # Original displayed top-five state strings are retained from adopted vector text.
 for letter,module,axis in [('c','CM09','axes_6'),('d','CM03','axes_7')]:
  v=n[n.module.eq(module)].sort_values('rank').head(5)
  # identity map field names are explicit below after inspection.
  names=[]
  for st in v.cell_state:
   q=idx[idx['cell_state'].eq(st)] if 'cell_state' in idx else idx[idx['state_id'].eq(st)]
   names.append(str(q.iloc[0]['level2_annotation']))
  labs=['\n'.join(textwrap.wrap(x,38)) for x in names];f=fig(136,65);a=axmm(f,70,9,58,40);a.barh(range(5),v.weight,height=.60,color=co(module));a.set_yticks(range(5),labs);a.set(xlim=(0,.3),ylim=(4.5,-.5),xlabel='NMF weight');a.set_xticks([0,.15,.3]);a.grid(axis='x',color='#E4E8EA',lw=.45);a.set_axisbelow(True);clean(a);f.text(99/136,1-3/65,module,ha='center',va='top',fontsize=8,weight='bold');save(f,'Figure4'+letter,[p2,p3],'shared 0–0.30 scale retained; enlarged 58mm data width and 8mm rows; module centred above data',{'module':module,'rows':5})
 f=fig(154,89);keys=['effective_n_studies','effective_n_tissues','median_pooled_rho'];labels=['Effective studies','Effective tissues','Median pooled\nSpearman ρ']
 for j,(key,lab) in enumerate(zip(keys,labels)):
  a=axmm(f,15+45*j,6,38,59);a.barh(range(14),d[key],height=.6,color=[co(m) for m in d.module]);a.set(ylim=(13.5,-.5),xlabel=lab);a.set_yticks(range(14),d.module if j==0 else []);a.grid(axis='x',color='#E4E8EA',lw=.4);a.set_axisbelow(True);clean(a)
  if j==2:a.set_xlim(0,1)
 f.legend(handles=[Patch(fc=co(m),label=m if m!='Other' else 'Other modules') for m in ['CM01','CM03','CM09','Other']],loc='lower center',ncol=4,frameon=False);save(f,'Figure4f',[p],'shared row centres across three evidence columns; full module key retained')

def qc():
 p=S/'ExtendedData1/QC247_sample_medians.tsv';d=read(p);f=fig(166,56)
 for j,(key,lab) in enumerate(zip(['n_genes_by_counts','total_counts','pct_counts_mt','doublet_score'],['Detected genes','Total UMIs','Mitochondrial counts (%)','Doublet score'])):
  a=axmm(f,15+j*40,8,27,33);v=d[key].dropna();a.boxplot(v,positions=[0],widths=.48,showfliers=False,patch_artist=True,boxprops=dict(fc='#D9E5EC',ec='#3B6F91',lw=.65),medianprops=dict(c='#20252A',lw=.85),whiskerprops=dict(c='#747C86'),capprops=dict(c='#747C86'));a.scatter(np.linspace(-.2,.2,len(v)),v,s=3,c='#3B6F91',alpha=.4,edgecolors='none');a.set_xticks([]);a.set_xlabel(lab);a.set_xlim(-.45,.45);a.spines[['top','right','bottom']].set_visible(False)
 save(f,'ExtendedData_Figure1a',[p],'same247 sample medians shown individually with summary distribution; no QC recomputation',{'sample_count':247})
 # ranks in existing sorted counts are descriptive display, no change to observations.
 f=fig(89,65);a=axmm(f,17,8,64,42);vals=d.n_cells.sort_values(ascending=False).to_numpy();a.plot(np.arange(1,248),vals,c='#3B6F91',lw=.75);a.scatter(np.arange(1,248),vals,s=4,c='#3B6F91',edgecolors='none');a.set(yscale='log',xlabel='Sample rank by retained cell count',ylabel='Retained cells',xlim=(1,247));save(f,'ExtendedData_Figure1b',[p],'247 retained counts; readable rank axes; original log scale')
 p2=S/'ExtendedData1/benchmark_publication_CPU_only.tsv';b=read(p2);f=fig(110,68);a=axmm(f,24,8,78,45)
 colours=['#5D859B','#A5B8BC','#2A7F78']
 for j,key in enumerate(['biology_score','batch_score','overall_score']):a.scatter(b[key],np.arange(8)+(j-1)*.22,s=19 if j<2 else 24,color=colours[j],marker=['o','s','D'][j],label=['Biology','Batch','Overall'][j],edgecolors='none')
 names={'scpoli':'scPoli','scanvi':'scANVI','harmony':'Harmony','scvi':'scVI','scalex':'SCALEX','pca':'PCA','liger':'LIGER','scanorama':'Scanorama'};a.set_yticks(range(8),[names.get(x.lower(),x) for x in b.method]);a.set(ylim=(7.5,-.5),xlim=(0,1),xlabel='scIB score');a.grid(axis='x',c='#E4E8EA',lw=.5);clean(a);f.legend(*a.get_legend_handles_labels(),loc='lower center',ncol=3,frameon=False);save(f,'ExtendedData_Figure1c',[p2],'three score types with redundant shapes; full8 methods; original ranking')
 v=b[b.method.str.lower().str.contains('harmony')].iloc[0];f=fig(89,49);a=axmm(f,22,7,56,28);vals=[v.biology_score,v.batch_score,v.overall_score];a.barh(range(3),vals,color=colours,height=.58);a.set_yticks(range(3),['Biology','Batch','Overall']);a.set(ylim=(2.5,-.5),xlim=(0,1),xlabel='scIB score');clean(a)
 for i,val in enumerate(vals):a.text(val+.02,i,f'{val:.3f}',va='center')
 save(f,'ExtendedData_Figure1d',[p2],'individual Harmony scores with explicit0–1 scale and direct values')

if __name__=='__main__':
 try:
  markers();lineage();genetics();forests();modules();qc()
 finally:
  (O/'review/redraw_manifest.json').write_text(json.dumps(records,ensure_ascii=False,indent=2));(O/'review/source_hashes.json').write_text(json.dumps(sources,ensure_ascii=False,indent=2))
