"""Figure 1 layout only: frozen coordinates and stored summary statistics.

No inference, embedding, annotation, normalization or expression aggregation.
Native panel sizes are also the BioRender placement sizes (mm).
"""
from pathlib import Path
import os, io, json, copy, base64, hashlib
os.environ['DYLD_FALLBACK_LIBRARY_PATH']='/Users/taki/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/lib'
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[k]='4'
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,FancyBboxPatch,FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap,Normalize
from PIL import Image
from lxml import etree as E
import fitz
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

P=Path('/Users/taki/Desktop/投稿包_最新版_20260919')
O=Path(__file__).parent; S=P/'Source_Data_Figures'
(O/'panels').mkdir(exist_ok=True); (O/'review').mkdir(exist_ok=True)
MM=72/25.4; NS='http://www.w3.org/2000/svg'; XL='http://www.w3.org/1999/xlink'; Q=lambda x:'{'+NS+'}'+x
mpl.rcParams.update({'font.family':'Arial','font.size':6.2,'axes.labelsize':6.2,'xtick.labelsize':6,'ytick.labelsize':6,
 'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.45,'axes.spines.top':False,'axes.spines.right':False,
 'text.color':'#222831','axes.labelcolor':'#222831','xtick.color':'#222831','ytick.color':'#222831','savefig.facecolor':'white'})

palette=pd.read_csv(S/'Figure1/L1_palette_bound.tsv',sep='\t').sort_values('order')
palette['hex']=['#356D8A','#C09554','#718C65','#8D739E','#AD6879','#2A8C85','#AE7658','#7D91B2','#658A98','#56799D','#B77455','#9C8CAD','#B78D8D','#8A995B','#A4B4B8','#C5C8CC']
L1=palette.L1.tolist(); pal=palette.set_index('L1').hex.to_dict();palette.to_csv(O/'L1_display_palette.tsv',sep='\t',index=False)
short={x:x.replace(' cell','') for x in L1}
TP={'Aorta':'#708F9C','Coronary artery':'#B98261','Carotid artery':'#8A9C76','Pulmonary artery':'#9486A6','Aortic valve':'#B6A47D','Other systemic arteries':'#6C9C9A'}
pd.DataFrame(TP.items(),columns=['Tissue','hex']).to_csv(O/'Tissue_display_palette.tsv',sep='\t',index=False)
D=pd.read_csv(S/'Figure1/full_cells_frozen_umap.tsv.gz',sep='\t');assert len(D)==1175115 and D._index.is_unique
T=pd.read_csv(S/'Current_Display_20260920/Figure1b_tissue_sample_mapping.tsv',sep='\t');assert T.sample_uid.is_unique
D['Tissue']=D.sample_uid.map(T.set_index('sample_uid').Tissue);assert D.Tissue.notna().all()
extent=[D.UMAP1.min(),D.UMAP1.max(),D.UMAP2.min(),D.UMAP2.max()]
dx=extent[1]-extent[0];dy=extent[3]-extent[2]
LIM=[extent[0]-.04*dx,extent[1]+.04*dx,extent[2]-.04*dy,extent[3]+.04*dy]
ST=pd.read_csv('/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Figure_Redesign_20260917_v3/config/state_id_to_hex.tsv',sep='\t')
qual=['#377E9A','#B98D4F','#7C6EAA','#4D9B8D','#B56870','#768B51','#6486B6','#B67D58','#9B7FA1','#4C7779','#AD9D67','#718C9C','#B8808F','#889869','#526C9A','#8C7155','#6E9EAA']
sub=ST[ST.level1_annotation.isin(['Macrophage','Fibroblast','Endothelial cell'])].copy()
for par,z in sub.groupby('level1_annotation'):
 for j,ix in enumerate(z.sort_values('within_parent_display_order').index):sub.loc[ix,'hex']=qual[j]
sub.to_csv(O/'L2_display_palette.tsv',sep='\t',index=False)
lp=sub.set_index(['level1_annotation','level2_annotation']).hex.to_dict()
WINDOWS={'Macrophage':[-1.5,10.1,8.1,20.3],'Fibroblast':[-11.,1.,.7,12.5],'Endothelial cell':[-11.,-4.9,-3.9,2.2]}
manifest={};stats={};geom=[]

def new(w,h):
 f=plt.figure(figsize=(w/25.4,h/25.4));f._mm=(w,h);return f
def ax(f,x,y,w,h):
 W,H=f._mm;a=f.add_axes([x/W,1-(y+h)/H,w/W,h/H]);return a
def tx(f,x,y,s,size=6.2,**kw):
 W,H=f._mm;return f.text(x/W,1-y/H,s,fontsize=size,va='center',**kw)
def letter(f,s):
 # Panel letters belong to the assembly, never to individual image files.
 pass
def arrow(f,x1,y1,x2,y2):
 W,H=f._mm;f.add_artist(FancyArrowPatch((x1/W,1-y1/H),(x2/W,1-y2/H),transform=f.transFigure,arrowstyle='-|>',mutation_scale=5,lw=.65,color='#64717C'))
def save(f,name):
 w,h=f._mm
 f.savefig(O/'panels'/f'Figure1{name}.pdf',dpi=900)
 f.savefig(O/'panels'/f'Figure1{name}.svg',dpi=900)
 f.savefig(O/'panels'/f'Figure1{name}.png',dpi=900)
 f.canvas.draw();rr=f.canvas.get_renderer();boxes=[]
 for a in f.axes:
  p=a.get_position();b=a.get_tightbbox(rr).transformed(f.transFigure.inverted());boxes.append({'data_mm':[p.x0*w,(1-p.y1)*h,p.width*w,p.height*h],'content_mm':[b.x0*w,(1-b.y1)*h,b.width*w,b.height*h]})
 manifest[name]={'width_mm':w,'height_mm':h,'axes':boxes,'raster_dpi':900}
 plt.close(f)

# a: reuse licensed native BioRender assets, repositioned without stretching.
native=E.parse(str(O/'source_assets/Figure1a_BioRender_native.svg'))
ims={e.get('id'):Image.open(io.BytesIO(base64.b64decode(e.get('{'+XL+'}href').split(',')[1]))).convert('RGBA') for e in native.iter(Q('image'))}
def icon(f,key,x,y,w,h):
 a=ax(f,x,y,w,h);a.set_zorder(2);a.imshow(ims[key]);a.set_aspect('equal');a.axis('off')
f=new(183,29);letter(f,'a')
tx(f,27,3.5,'Public human vascular datasets',6.6,weight='bold',ha='center')
icon(f,'image_59',5,7,5.8,13);icon(f,'image_63',6.1,7,4.2,13)
icon(f,'image_11',14,8,5,5);icon(f,'image_15',14,14,5,5)
tx(f,21,10.5,'26 studies',6.6,weight='bold');tx(f,21,16.5,'247 sample identifiers',6.3)
arrow(f,48,12,54,12)
tx(f,83,3.5,'Integration and annotation',6.6,weight='bold',ha='center')
for key,x,w in [('image_55',57,9),('image_51',70,6),('image_47',80,9),('image_43',96,11)]:icon(f,key,x,7,w,8)
for x1,x2 in [(66,69),(76.5,79),(90,95)]:arrow(f,x1,11,x2,11)
for x,s in [(61,'Counts and\nmetadata'),(73,'QC'),(84.5,'Harmony'),(101.5,'L1 / L2\nannotation')]:tx(f,x,18,s,6,ha='center',linespacing=1.0)
arrow(f,110,12,116,12)
tx(f,151,3.5,'Panvascular atlas',6.6,weight='bold',ha='center')
icon(f,'image_31',120,8,12,7);tx(f,137,10,'1,175,115 cells',8,weight='bold')
tx(f,150,17,'16 L1 categories · 116 named L2 states',6.1,ha='center')
W,H=f._mm;f.add_artist(FancyBboxPatch((3/W,1-28.7/H),177/W,6.3/H,boxstyle='round,pad=0.0015',transform=f.transFigure,fc='#EDF2F5',ec='none'))
tx(f,6,25.5,'Downstream analyses',6.3,weight='bold')
for key,x,label in [('image_19',44,'Genetics'),('image_7',76,'Expression'),('image_23',111,'Cell communication'),('image_67',153,'Spatial')]:icon(f,key,x,23,5,5);tx(f,x+6.5,25.5,label,6.1)
save(f,'a')

# b: reviewed compact callouts and the six-class tissue legend.
from figure1b_layout import draw_panel
stats['b']=draw_panel(D,pal,TP,LIM,S,O)
manifest['b']={'width_mm':107,'height_mm':51,'raster_dpi':900,
 'axes':[], 'layout_source':'figure1b_layout.py'}

# c and f: identical 16 row centres, saved values only.
C=pd.read_csv(P/'图件编辑资料/可视化库重绘_当前/data_used/Figure1c_f_exact_display_values.tsv',sep='\t').set_index('L1').loc[L1]
assert C.cells.sum()==1175115
for panel,w,datax,dataw,value,caption in [('c',60,18.5,29,'cells','Cells per lineage'),('f',13,1,8.5,'studies','Studies per lineage')]:
 f=new(w,51);letter(f,panel);a=ax(f,datax,5.5,dataw,40.5);y=np.arange(16)
 a.barh(y,C[value],color=[pal[s] for s in L1],height=.67);a.set_ylim(15.5,-.5);a.set_yticks(y)
 a.set_yticklabels([short[s] for s in L1] if panel=='c' else [],fontsize=6)
 a.tick_params(axis='y',length=0,pad=2);a.tick_params(axis='x',length=2,pad=1.5,labelsize=5.8);a.spines[['left','right','top']].set_visible(False)
 if panel=='c':a.set_xlim(0,270000);a.set_xticks([0,100000,200000],['0','100k','200k'])
 else:a.set_xlim(0,28);a.set_xticks([])
 for j,v in enumerate(C[value]):a.text(1.025,j,f'{v:,}' if panel=='c' else str(v),transform=a.get_yaxis_transform(),va='center',fontsize=6)
 tx(f,datax+dataw/2 if panel=='c' else 6.5,50 if panel=='c' else 48.3,caption if panel=='c' else 'Studies per\nlineage',6,ha='center',linespacing=.9)
 save(f,panel)

# e: identical old coordinate windows and full locator; new display colours only.
f=new(122,35);letter(f,'e');target=D[D.level1_annotation.isin(WINDOWS)].copy()
target['display_hex']=[lp.get((a,b),'#A4AAAF') for a,b in zip(target.level1_annotation,target.level2_annotation)]
loc=ax(f,1,7,21,26);loc.scatter(target.UMAP1,target.UMAP2,c=target.display_hex,s=.025,lw=0,alpha=.75,rasterized=True);loc.set_aspect('equal');loc.axis('off');loc.margins(.045)
for j,(parent,lim) in enumerate(WINDOWS.items()):
 x=26+32*j;z=target[target.level1_annotation.eq(parent)];a=ax(f,x,6.4,30,28)
 ix=np.random.default_rng(20260917).permutation(len(z));a.scatter(z.UMAP1.to_numpy()[ix],z.UMAP2.to_numpy()[ix],c=z.display_hex.to_numpy()[ix],s=.07,lw=0,alpha=.77,rasterized=True);a.set(xlim=lim[:2],ylim=lim[2:]);a.set_aspect('equal');a.axis('off')
 tx(f,x+15,3.2,'ABC'[j]+'  '+short[parent],6.6,weight='bold',ha='center')
 loc.add_patch(Rectangle((lim[0],lim[2]),lim[1]-lim[0],lim[3]-lim[2],fill=False,ec='#424C55',lw=.5))
 loc.annotate('ABC'[j],xy=((lim[0]+lim[1])/2,(lim[2]+lim[3])/2),xytext=[(.82,.96),(.05,.60),(.04,.13)][j],textcoords='axes fraction',ha='center',va='center',fontsize=6,weight='bold',bbox=dict(fc='white',ec='none',pad=.15),arrowprops=dict(arrowstyle='-',lw=.3,color='#424C55'))
 outside=~z.UMAP1.between(lim[0],lim[1])|~z.UMAP2.between(lim[2],lim[3]);stats[parent]={'cells':len(z),'window':lim,'inside':int((~outside).sum()),'outside':int(outside.sum()),'all_in_locator':True}
save(f,'e')

# d: original tissue mean fractions, no recomputation.
f=new(58,35);letter(f,'d');comp=pd.read_csv(S/'Figure1/Fig1d.tsv',sep='\t').set_index('tissue').reindex(columns=L1)
a=ax(f,29.5,4,27,25.3);left=np.zeros(len(comp))
for k in L1:a.barh(range(len(comp)),comp[k],left=left,color=pal[k],height=.8);left+=comp[k]
a.set_ylim(len(comp)-.5,-.5);a.set_yticks(range(len(comp)),comp.index,fontsize=6);a.set_xlim(0,1);a.set_xticks([0,.5,1]);a.tick_params(axis='y',length=0,pad=2);a.tick_params(axis='x',length=2,pad=1,labelsize=6);a.spines[['top','right','left']].set_visible(False);tx(f,43,33,'Mean sample fraction',6,ha='center');save(f,'d')

# g: B display adopted by author: existing per-gene values; original fractions.
path=S/'Current_Display_20260920/marker_audit_20260921'
V=pd.read_csv(path/'per_gene_scaled_mean.tsv',sep='\t',index_col=0);F=pd.read_csv(path/'original_detection_fraction.tsv',sep='\t',index_col=0)
assert V.shape==(14,28) and V.index.equals(F.index) and V.columns.equals(F.columns)
assert F.max().max()<=1 and F.min().min()>=0
cm=LinearSegmentedColormap.from_list('blue_teal_sequential',['#EFF4F5','#BAD4D8','#76A8B4','#367789','#173F5C'])
f=new(183,49);letter(f,'g');a=ax(f,23,3,140,34)
x,y=np.meshgrid(range(28),range(14));smax=(.70*min(140/28,34/14)*MM)**2;valid=V.notna().to_numpy()&F.notna().to_numpy()&(F.to_numpy()>0)
a.scatter(x[valid],y[valid],s=smax*F.to_numpy()[valid],c=V.to_numpy()[valid],cmap=cm,norm=Normalize(0,1),edgecolors='none',linewidths=0)
a.set(xlim=(-.5,27.5),ylim=(13.5,-.5));a.set_xticks(range(28),V.columns,rotation=90,fontstyle='italic',fontsize=6.5);a.set_yticks(range(14),[short[k] for k in V.index],fontsize=6.4);a.tick_params(length=0,pad=2);a.tick_params(axis='y',pad=6.5)
for sp in a.spines.values():sp.set_visible(False)
for j,s in enumerate(V.index):a.add_patch(Rectangle((-.86,j-.37),.12,.74,fc=pal[s],ec='none',clip_on=False))
tx(f,173.2,5,'Per-gene scaled\nmean expression',6.1,ha='center',linespacing=1.1)
cbax=ax(f,167,10,2.2,17);cb=f.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap=cm),cax=cbax);cb.set_ticks([0,.5,1]);cb.ax.tick_params(labelsize=6,length=2,pad=2);cb.outline.set_linewidth(.35)
tx(f,173,31,'Detected\ncells (%)',6.1,ha='center',linespacing=1)
la=ax(f,165,34.5,17,13);la.set(xlim=(0,17),ylim=(13,0));la.axis('off')
for j,z in enumerate([.25,.5,.75,1]):la.scatter(3,1.5+j*3,s=smax*z,c='#617E89',edgecolors='none');la.text(7,1.5+j*3,str(int(z*100)),fontsize=6.1,va='center')
stats['g']={'shape':list(V.shape),'size_area_max_pt2':smax,'mapping':'area = max area * original fraction; zero area at zero','colour':'stored per-gene min-max mean, 0–1','no_new_statistics':True}
save(f,'g')

# Native composition: only translations, no anisotropic scaling or hidden old panels.
positions={'a':(0,0),'b':(0,31),'c':(109,31),'f':(170,31),'e':(0,84),'d':(124,84),'g':(0,121)}
root=E.Element(Q('svg'),nsmap={None:NS,'xlink':XL},width='183mm',height='170mm',viewBox=f'0 0 {183*MM} {170*MM}')
doc=fitz.open();page=doc.new_page(width=183*MM,height=170*MM)
for name,(xx,yy) in positions.items():
 r=E.parse(str(O/'panels'/f'Figure1{name}.svg')).getroot();ids={e.get('id'):f'{name}_{e.get("id")}' for e in r.iter() if e.get('id')}
 for e in r.iter():
  for k,v in list(e.attrib.items()):
   if k=='id':e.set(k,ids[v])
   else:
    for old,newid in ids.items():v=v.replace('url(#'+old+')','url(#'+newid+')');v=('#'+newid) if v=='#'+old else v
    e.set(k,v)
 g=E.SubElement(root,Q('g'),id='panel_'+name,transform=f'translate({xx*MM} {yy*MM})')
 for e in list(r):g.append(e)
 d=fitz.open(O/'panels'/f'Figure1{name}.pdf');w,h=manifest[name]['width_mm'],manifest[name]['height_mm'];page.show_pdf_page(fitz.Rect(xx*MM,yy*MM,(xx+w)*MM,(yy+h)*MM),d,0);d.close()
 manifest[name]['origin_mm']=[xx,yy]
# Independent labels in the local composite mirror BioRender-native text objects.
fontfile='/System/Library/Fonts/Supplemental/Arial Bold.ttf'
page.insert_font(fontname='PanelArialBold',fontfile=fontfile)
for name,(xx,yy) in positions.items():
 px=(xx+1)*MM;py=(yy+2.2)*MM+2.8
 page.insert_text((px,py),name,fontname='PanelArialBold',fontsize=8,color=(.133,.157,.192))
 t=E.SubElement(root,Q('text'),id='panel_letter_'+name,x=str(px),y=str(py),fill='#222831',**{'font-family':'Arial','font-weight':'bold','font-size':'8'})
 t.text=name
 manifest[name]['panel_letter']={'included_in_panel_file':False,'assembly_anchor_mm':[xx+1,yy+2.2],'font_pt':8,'biorender_role':'independent text object'}
doc.save(O/'Figure1.pdf',garbage=4,deflate=True);page.get_pixmap(dpi=200).save(O/'review/Figure1.png');doc.close()
E.ElementTree(root).write(str(O/'Figure1.svg'),xml_declaration=True,encoding='utf-8')
(O/'panel_placement.json').write_text(json.dumps(manifest,indent=2))
(O/'preserved_data_checks.json').write_text(json.dumps(stats,indent=2,default=lambda x:x.item()))
print('DONE',json.dumps(stats,default=lambda x:x.item()))
