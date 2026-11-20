from pathlib import Path
import os,sys,copy,json,io,re
os.environ['DYLD_FALLBACK_LIBRARY_PATH']='/Users/taki/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/lib'
sys.path.insert(0,'/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Local_Touchup_20260918/tools')
import matplotlib as mpl;mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Circle,Polygon
import numpy as np,pandas as pd,fitz,cairosvg,xml.etree.ElementTree as E
from PIL import Image,ImageDraw
import cairosvg.image
from cairosvg.surface import cairo
cairosvg.image.IMAGE_RENDERING.update({"pixelated":cairo.FILTER_NEAREST,"crisp-edges":cairo.FILTER_NEAREST})
B=Path('/Users/taki/Desktop/Panvascular_FigureYa_20260919');O=Path(__file__).resolve().parents[1]
S=Path('/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Nature_Reconstruction_20260918/SourceData/spatial');A=Path('/Users/taki/Desktop/Codex任务1/spatial_visium_geo_local/extracted')
NS='http://www.w3.org/2000/svg';N={'s':NS};E.register_namespace('',NS);E.register_namespace('xlink','http://www.w3.org/1999/xlink');MM=72/25.4
mpl.rcParams.update({'font.family':'Arial','font.size':6,'text.color':'#222222','svg.fonttype':'none','pdf.fonttype':42})
params=json.load(open(S/'Figure5a_display_parameters.json'));rois=json.load(open(O/'audit/spatial_roi_definitions.json'))['rois']
cm=LinearSegmentedColormap.from_list('original_positive_spatial',['#F2B447','#D46A22','#A32C36','#541B43'])
f=plt.figure(figsize=(183/25.4,170/25.4));display=[]
for i,p in enumerate(params):
 gene=p['gene'];gse='GSE248608' if gene=='LIPA' else 'GSE314851';stem=p['section'];image=Image.open(A/gse/(stem+'_tissue_hires_image.png'))
 z=pd.read_csv(O/'audit'/f'spatial_counts_{gene}_plotted_spots.tsv',sep='\t');v=z.plotted_normalised_value.to_numpy();x=z.plot_hires_pixel_x.to_numpy();y=z.plot_hires_pixel_y.to_numpy();flag=np.zeros(len(z),bool)
 if gene=='LIPA':
  r=pd.read_csv(O/'audit/spatial_counts_roi_LIPA_central_clear_space_center_included.tsv',sep='\t');ids=set(r.loc[r.entire_physical_spot_disk_inside_roi,'barcode']);flag=z.barcode.isin(ids).to_numpy();assert flag.sum()==7
 norm=Normalize(p['positive_cap02'],p['positive_cap98'],clip=True)
 a=f.add_axes([(5+i*47)/183,1-49/170,39/183,39/170]);a.imshow(image,alpha=1,interpolation='nearest')
 zero=np.isfinite(v)&(v==0);positive=np.isfinite(v)&(v>0);missing=~np.isfinite(v)
 a.scatter(x[zero&~flag],y[zero&~flag],s=.65,c='#CED3D6',lw=0,rasterized=True)
 a.scatter(x[positive&~flag],y[positive&~flag],s=1.15,c=v[positive&~flag],cmap=cm,norm=norm,lw=0,rasterized=True)
 # Hollow retains the exact original expression colour, and is not recoded as zero or missing.
 colors=[cm(norm(val)) if val>0 else '#CED3D6' for val in v[flag]]
 a.scatter(x[flag],y[flag],s=3.0,facecolors='none',edgecolors=colors,lw=.5,rasterized=True)
 a.scatter(x[missing],y[missing],s=7,c='#333',marker='x',lw=.6,rasterized=True)
 a.set_xlim(p['xlim']);a.set_ylim(p['ylim']);a.set_aspect('equal');a.axis('off')
 display.append({'gene':gene,'n_spots':len(z),'n_positive':int(positive.sum()),'n_zero':int(zero.sum()),'n_flagged_hollow':int(flag.sum()),'flagged_barcodes':z.loc[flag,'barcode'].tolist(),'vmin':norm.vmin,'vmax':norm.vmax,'xlim':p['xlim'],'ylim':p['ylim'],'image_alpha':1,'positive_size_pt2':1.15,'zero_size_pt2':.65,'hollow_size_pt2':3,'expression_interpolation':False,'smoothing':False,'image_interpolation':'nearest','all_coordinates_and_values_unchanged':True})
 # Expanded review: unmodified image, original expression display, new display; physical spot footprint in audit only.
 for roi in [q for q in rois if q['gene']==gene and q['role']=='suspected_acellular_space']:
  rr=pd.read_csv(O/'audit'/f"spatial_counts_roi_{roi['id']}_center_included.tsv",sep='\t');poly=np.array(roi['polygon']);rad=float(rr.physical_spot_radius_hires_pixels.iloc[0]);pad=30
  fig,axs=plt.subplots(1,3,figsize=(10.5,3.7));limits=[poly[:,0].min()-pad,poly[:,0].max()+pad,poly[:,1].min()-pad,poly[:,1].max()+pad]
  for j,ax in enumerate(axs):
   ax.imshow(image,alpha=.68 if j==1 else 1,interpolation='nearest');ax.set_xlim(limits[:2]);ax.set_ylim(limits[3],limits[2]);ax.add_patch(Polygon(poly,fill=False,edgecolor='#247AA5',lw=.8))
   if j>0:
    for _,r in rr.iterrows():
     val=r.plotted_normalised_value;col=cm(norm(val)) if val>0 else '#CED3D6';hollow=j==2 and bool(r.entire_physical_spot_disk_inside_roi)
     ax.add_patch(Circle((r.source_hires_pixel_x,r.source_hires_pixel_y),rad*.65,facecolor='none' if hollow else col,edgecolor=col,lw=1))
   if j==0:
    for k,(_,r) in enumerate(rr.iterrows(),1):
     ax.add_patch(Circle((r.source_hires_pixel_x,r.source_hires_pixel_y),rad,fill=False,edgecolor='#222',lw=.6));ax.text(r.source_hires_pixel_x,r.source_hires_pixel_y,str(k),ha='center',va='center',fontsize=7,color='black')
   ax.set_xlabel('Hires image x (pixels)');ax.set_ylabel('Hires image y (pixels)');ax.set_title(['Original H&E + physical spot footprint','Previous expression encoding','Current encoding; hollow = H&E flag'][j],fontsize=8)
  fig.tight_layout();fig.savefig(O/'review'/f"spatial_ROI_{roi['id']}.png",dpi=180);fig.savefig(O/'review'/f"spatial_ROI_{roi['id']}.pdf");plt.close(fig)
leg=f.legend(handles=[Line2D([],[],marker='o',ls='',mfc='#CED3D6',mec='none',label='Zero',ms=3),Line2D([],[],marker='o',ls='',mfc='none',mec='#555',label='H&E flag',ms=3)],loc='upper left',bbox_to_anchor=(31/183,1-56/170),ncol=2,frameon=False,borderaxespad=0,fontsize=6,columnspacing=.9,handletextpad=.3)
f.savefig(O/'review/Figure5a_replacement.svg',dpi=400);plt.close(f)
new=E.parse(O/'review/Figure5a_replacement.svg').getroot();ids={el.get('id'):'st_audit_'+el.get('id') for el in new.iter() if el.get('id')}
for el in new.iter():
 for key,val in list(el.attrib.items()):
  if key=='id':el.set(key,ids[val])
  else:
   for old,n in ids.items():val=val.replace('url(#'+old+')','url(#'+n+')');val=n if False else val
   if val.startswith('#') and val[1:] in ids:val='#'+ids[val[1:]]
   el.set(key,val)
tree=E.parse(B/'editable_svg/Figure5.svg');root=tree.getroot();fg=root.find("s:g[@id='figure_1']",N)
for target,source in [('axes_1','axes_1'),('axes_3','axes_2'),('legend_1','legend_1')]:
 old=fg.find(f"s:g[@id='{target}']",N);replacement=new.find(f".//s:g[@id='st_audit_{source}']",N);replacement.set('id',target);idx=list(fg).index(old);fg.remove(old);fg.insert(idx,copy.deepcopy(replacement))
# Collect definitions for new images/clips; originals remain usable by untouched panels.
for d in new.findall('.//s:defs',N):root.append(copy.deepcopy(d))
def bbox(el):
 rr=E.Element(root.tag,root.attrib);defs=E.SubElement(rr,'{'+NS+'}defs')
 for d in root.findall('.//s:defs',N):
  for ch in d:defs.append(copy.deepcopy(ch))
 rr.append(copy.deepcopy(el));ar=np.array(Image.open(io.BytesIO(cairosvg.svg2png(bytestring=E.tostring(rr),scale=2,background_color='white'))).convert('RGB'));yy,xx=np.where(ar.min(2)<245)
 return np.array([xx.min(),yy.min(),xx.max()+1,yy.max()+1])*.375 if len(xx) else None
def shift(el,dx,dy):el.set('transform',f'translate({dx:.5f} {dy:.5f}) '+el.get('transform',''))
def union(bs):return np.array([min(b[0] for b in bs),min(b[1] for b in bs),max(b[2] for b in bs),max(b[3] for b in bs)])
def get(id):return fg.find(f"s:g[@id='{id}']",N)
# Separate f legend from e axis text; move e as a coherent block 1.5 mm up.
shift(get('axes_9'),0,-1.5*MM);shift(get('text_156'),0,-1.5*MM)
# Align legend left edge with its visible labelled data block.
legend_records=[]
for lg,target,y in [('legend_2','FigureYa6_native_forest',None),('legend_3','text_138',None),('legend_4','axes_10',None)]:
 el=get(lg);bb=bbox(el);tb=bbox(get(target));dx=tb[0]-bb[0];dy=0
 if lg=='legend_4':dy=396-bb[1]
 shift(el,dx,dy);legend_records.append({'legend':lg,'dx_pt':dx,'dy_pt':dy,'target':target})
# Explicit content membership: no automated reassignment of nearby labels.
panels={'a':(['axes_1','axes_2','axes_3','axes_4','text_133','text_134','text_135','legend_1'],(2,100)), 'b':(['axes_5','axes_6'],(101,181)), 'c':(['FigureYa6_native_forest','legend_2'],(2,100)), 'd':(['axes_8','legend_3']+[f'text_{i}' for i in range(138,148)],(102,181)), 'e':(['axes_9','text_156'],(102,181)), 'f':(['axes_10','legend_4']+[f'text_{i}' for i in range(150,156)],(102,181))}
records=[]
for panel,(members,(l,r)) in panels.items():
 bs=[bbox(get(q)) for q in members];bb=union(bs);dx=(l+r)*MM/2-(bb[0]+bb[2])/2
 for q in members:shift(get(q),dx,0)
 records.append({'panel':panel,'before_bbox_pt':bb.tolist(),'dx_pt':dx,'after_bbox_pt':(bb+np.array([dx,0,dx,0])).tolist(),'centering':'horizontal visible-content bbox; vertical bands retained to preserve gutters'})
for lg in root.findall('.//s:g',N):
 if lg.get('id','').startswith('legend_'):
  for t in lg.findall('.//s:text',N):
   assert 'text-anchor: middle' not in t.get('style','') and 'text-anchor: end' not in t.get('style','')
tree.write(O/'editable_svg/Figure5.svg',encoding='utf-8',xml_declaration=True);cairosvg.svg2pdf(url=str(O/'editable_svg/Figure5.svg'),write_to=str(O/'main/Figure5.pdf'))
(O/'audit/Figure5_display_and_alignment.json').write_text(json.dumps({'spatial_display':display,'panels':records,'legends':legend_records},indent=2))
for side,path in [('before',B/'main/Figure5.pdf'),('after',O/'main/Figure5.pdf')]:
 d=fitz.open(path);d[0].get_pixmap(matrix=fitz.Matrix(1.8,1.8)).save(O/'review'/f'Figure5_{side}.png')
ims=[Image.open(O/'review'/f'Figure5_{s}.png') for s in ['before','after']];out=Image.new('RGB',(ims[0].width*2+20,ims[0].height+30),'white');dr=ImageDraw.Draw(out)
for i,(im,label) in enumerate(zip(ims,['BEFORE','AFTER'])):out.paste(im,(i*(im.width+20),30));dr.text((i*(im.width+20)+5,5),label,fill='black')
out.save(O/'review/Figure5_comparison.png')
print('Figure5 complete')
