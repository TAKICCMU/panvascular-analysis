from pathlib import Path
import os,sys,copy,re,json,hashlib,io
os.environ['DYLD_FALLBACK_LIBRARY_PATH']='/Users/taki/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/lib'
sys.path.insert(0,'/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Local_Touchup_20260918/tools')
import numpy as np,pandas as pd,xml.etree.ElementTree as ET
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cairosvg,cairosvg.image,fitz
from cairosvg.surface import cairo
cairosvg.image.IMAGE_RENDERING.update({'pixelated':cairo.FILTER_NEAREST,'crisp-edges':cairo.FILTER_NEAREST})
plt.rcParams.update({'font.family':'Arial','font.size':6.5,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.45,'xtick.major.width':.45,'ytick.major.width':.45,'text.color':'#222222','axes.labelcolor':'#222222'})
O=Path(__file__).resolve().parents[1]; B=Path('/Users/taki/Desktop/Panvascular_Main_Figure_Refinement_20260920/reproducibility_dependencies/Panvascular_Marker_Spatial_Review_20260919'); A=Path('/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Full_Layout_Correction_20260918'); S=Path('/Users/taki/Desktop/投稿包_最新版_20260919/Source_Data_Figures');N={'s':'http://www.w3.org/2000/svg'};NS=N['s'];q=lambda s:'{'+NS+'}'+s
ET.register_namespace('',NS);ET.register_namespace('xlink','http://www.w3.org/1999/xlink');mm=72/25.4
P=pd.read_csv(S/'Figure1/L1_palette_bound.tsv',sep='\t').sort_values('order');pal=P.set_index('L1')['hex'].to_dict()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def text(g,x,y,t,size=6.5,anchor='start',**kw):
 e=ET.SubElement(g,q('text'),{'x':str(x),'y':str(y),'text-anchor':anchor,'style':f"font:{size}px Arial;fill:#222222",**kw});e.text=t;return e

def export(root,name):
 svg=O/'editable_svg'/f'{name}.svg';ET.ElementTree(root).write(svg,encoding='utf-8',xml_declaration=True);pdf=O/'figures'/f'{name}.pdf';cairosvg.svg2pdf(url=str(svg),write_to=str(pdf));d=fitz.open(pdf);d[0].get_pixmap(matrix=fitz.Matrix(2,2)).save(O/'review'/f'{name}.png');return pdf

def embed_plot(root,fig,prefix):
 buf=io.BytesIO();fig.savefig(buf,format='svg',transparent=True);plt.close(fig);r=ET.fromstring(buf.getvalue());ids={e.get('id'):prefix+e.get('id') for e in r.iter() if e.get('id')}
 for e in r.iter():
  for k,v in list(e.attrib.items()):
   if k=='id':e.set(k,ids[v])
   else:
    for old,new in ids.items():
     v=v.replace('url(#'+old+')','url(#'+new+')')
     if v=='#'+old:v='#'+new
    e.set(k,v)
 for ch in r:
  if ch.tag in [q('defs'),q('g')]:root.append(ch)

def remove_ids(root,ids):
 for p in root.iter():
  for ch in list(p):
   if ch.get('id') in ids:p.remove(ch)
