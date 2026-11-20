from pathlib import Path
import os,sys,json,copy,hashlib
os.environ['DYLD_FALLBACK_LIBRARY_PATH']='/Users/taki/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/lib'
sys.path.insert(0,'/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Local_Touchup_20260918/tools')
import matplotlib as mpl;mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,LinearSegmentedColormap
from matplotlib.lines import Line2D
import numpy as np,pandas as pd,fitz,cairosvg,xml.etree.ElementTree as E
import cairosvg.image
from cairosvg.surface import cairo
from PIL import Image,ImageDraw
cairosvg.image.IMAGE_RENDERING.update({'pixelated':cairo.FILTER_NEAREST,'crisp-edges':cairo.FILTER_NEAREST})
B=Path('/Users/taki/Desktop/Panvascular_Content_Alignment_ST_Audit_20260919');O=Path(__file__).resolve().parents[1];A=Path('/Users/taki/Desktop/Codex任务1/spatial_visium_geo_local/extracted')
NS='http://www.w3.org/2000/svg';N={'s':NS};E.register_namespace('',NS);E.register_namespace('xlink','http://www.w3.org/1999/xlink');MM=72/25.4
mpl.rcParams.update({'font.family':'Arial','font.size':6,'text.color':'#222222','svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.4})
params=json.load(open(B/'audit/Figure5_display_and_alignment.json'))['spatial_display']
# User-requested purple-blue-cyan-green-yellow-red sequence inspired by supplied screenshot.
hexes=['#482878','#3B528B','#2196C4','#35B779','#FDE725','#F68B32','#B2182B'];cm=LinearSegmentedColormap.from_list('user_reference_spatial',hexes);zero='#8B8F95'
f=plt.figure(figsize=(183/25.4,170/25.4));records=[]
def ax(x,y,w,h):return f.add_axes([x/183,1-(y+h)/170,w/183,h/170])
def tx(x,y,t,**kw):return f.text(x/183,1-y/170,t,fontsize=6,va='top',**kw)
for i,p in enumerate(params):
 gene=p['gene'];gse='GSE248608' if gene=='LIPA' else 'GSE314851';stem='GSM7917895_CA2' if gene=='LIPA' else 'GSM9414344_LCA';z=pd.read_csv(B/'audit'/f'spatial_counts_{gene}_plotted_spots.tsv',sep='\t',float_precision='round_trip');v=z.plotted_normalised_value.to_numpy();x=z.plot_hires_pixel_x.to_numpy();y=z.plot_hires_pixel_y.to_numpy();flag=z.barcode.isin(p['flagged_barcodes']).to_numpy();assert flag.sum()==p['n_flagged_hollow']
 image=Image.open(A/gse/(stem+'_tissue_hires_image.png'));norm=Normalize(p['vmin'],p['vmax'],clip=True);left=6+47*i
 tx(left,4.7,gene+' · '+stem)
 hea=ax(left,7.7,42,20.5);hea.imshow(image,alpha=1,interpolation='nearest');hea.set_xlim(p['xlim']);hea.set_ylim(p['ylim']);hea.set_aspect('equal');hea.axis('off')
 sa=ax(left,30,42,20.5);sa.set_facecolor('#050708');sa.set_xlim(p['xlim']);sa.set_ylim(p['ylim']);sa.set_aspect('equal');sa.set_xticks([]);sa.set_yticks([]);sa.spines[:].set_visible(False)
 # Marker diameter is a display choice, independent of the estimated capture footprint in histology review.
 positive=np.isfinite(v)&(v>0);iszero=np.isfinite(v)&(v==0);missing=~np.isfinite(v)
 sa.scatter(x[iszero&~flag],y[iszero&~flag],s=.42,c=zero,lw=0,rasterized=False)
 sa.scatter(x[positive&~flag],y[positive&~flag],s=.55,c=v[positive&~flag],cmap=cm,norm=norm,lw=0,rasterized=False)
 ec=[cm(norm(val)) if val>0 else zero for val in v[flag]]
 sa.scatter(x[flag],y[flag],s=3.0,facecolors='none',edgecolors=ec,lw=.45,rasterized=False)
 if missing.any():sa.scatter(x[missing],y[missing],s=2,c='white',marker='x',lw=.4)
 # Both rows have exactly the same data window, direction, axes size, and aspect ratio.
 cba=ax(left+6,52.1,30,1.5);cb=f.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cm),cax=cba,orientation='horizontal');cb.set_ticks([p['vmin'],p['vmax']]);cb.set_ticklabels([f"≤{p['vmin']:.2g}",f"≥{p['vmax']:.2g}"]);cb.ax.tick_params(labelsize=6,pad=1,length=1.5)
 f.canvas.draw();data_he=(hea.get_position().bounds);data_sp=(sa.get_position().bounds)
 records.append({'gene':gene,'section':stem,'n_spots':len(z),'n_positive':int(positive.sum()),'n_zero':int(iszero.sum()),'n_missing':int(missing.sum()),'flagged_barcodes':p['flagged_barcodes'],'vmin':p['vmin'],'vmax':p['vmax'],'xlim':p['xlim'],'ylim':p['ylim'],'same_H_and_E_spot_scale':np.allclose(data_he[2:],data_sp[2:]),'positive_size_pt2':.55,'zero_size_pt2':.42,'uncertain_hollow_size_pt2':3.0,'uncertain_linewidth_pt':.45,'capture_radius_used_as_main_marker_radius':False,'spot_background':'#050708','positive_palette':hexes,'zero_colour':zero,'separate_HE_and_spots':True,'expression_smoothing':False,'expression_interpolation':False,'allocated_region_mm':[left,4,42,52], 'actual_HE_data_region_mm':[data_he[0]*183,(1-data_he[1]-data_he[3])*170,data_he[2]*183,data_he[3]*170],'actual_spot_data_region_mm':[data_sp[0]*183,(1-data_sp[1]-data_sp[3])*170,data_sp[2]*183,data_sp[3]*170]})
# Necessary row identifiers, not panel subtitles.
f.text(3.2/183,1-18/170,'H&E',fontsize=6,ha='center',va='center',rotation=90)
f.text(3.2/183,1-40/170,'Expression',fontsize=6,ha='center',va='center',rotation=90)
tx(6,57.2,'log1p(CP10k)')
f.legend(handles=[Line2D([],[],marker='o',ls='',mfc=zero,mec='none',label='Zero',ms=3.2),Line2D([],[],marker='o',ls='',mfc='none',mec='#444',label='Histology uncertain',ms=4)],loc='upper left',bbox_to_anchor=(33/183,1-56.6/170),ncol=2,frameon=False,borderaxespad=0,handletextpad=.3,columnspacing=.7,fontsize=6)
f.savefig(O/'review/Figure5a_new_components.svg',dpi=500);plt.close(f)
new=E.parse(O/'review/Figure5a_new_components.svg').getroot();idmap={el.get('id'):'separate_st_'+el.get('id') for el in new.iter() if el.get('id')}
for el in new.iter():
 for key,val in list(el.attrib.items()):
  if key=='id':el.set(key,idmap[val]);continue
  for old,n in idmap.items():val=val.replace('url(#'+old+')','url(#'+n+')')
  if val.startswith('#') and val[1:] in idmap:val='#'+idmap[val[1:]]
  el.set(key,val)
tree=E.parse(B/'editable_svg/Figure5.svg');root=tree.getroot();fg=root.find("s:g[@id='figure_1']",N)
remove=['axes_1','axes_2','axes_3','axes_4','text_133','text_134','text_135','legend_1']
for gid in remove:
 old=fg.find(f"s:g[@id='{gid}']",N);assert old is not None;fg.remove(old)
wrapper=E.SubElement(fg,'{'+NS+'}g',{'id':'Figure5a_separate_HE_spots'})
for el in new.find("s:g[@id='separate_st_figure_1']",N):
 if el.get('id')!='separate_st_patch_1':wrapper.append(copy.deepcopy(el))
# Add top-level clip definitions, avoiding duplicates already inside the imported group.
existing={el.get('id') for el in root.iter() if el.get('id')}
for defs in new.findall('s:defs',N):
 d=E.SubElement(root,'{'+NS+'}defs')
 for el in defs:
  if not el.get('id') or el.get('id') not in existing:d.append(copy.deepcopy(el))
# Matplotlib clip definitions are normally top-level; all referenced imports verified below.
svg=O/'editable_svg/Figure5.svg';tree.write(svg,encoding='utf-8',xml_declaration=True);cairosvg.svg2pdf(url=str(svg),write_to=str(O/'main/Figure5.pdf'))
(O/'audit/Figure5a_display.json').write_text(json.dumps({'status':'HISTOLOGY_UNCONFIRMED_AUTHOR_REVIEW','baseline':str(B),'panel_allocated_region_mm':[0,0,100,60],'data_regions':records,'rule':'existing in_tissue=1 and total UMI >=500; no new exclusions','palette_reason':'Explicit user request to match supplied example; approximate custom sequence, not claimed exact source palette','uncertainty_legend':'Histology uncertain','ROI_unchanged':True},indent=2))
for side,p in [('before',B/'main/Figure5.pdf'),('after',O/'main/Figure5.pdf')]:
 d=fitz.open(p);d[0].get_pixmap(matrix=fitz.Matrix(2,2)).save(O/'review'/f'Figure5_{side}.png')
print('Figure5 exported')
