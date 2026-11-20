"""Standalone editable panels from the adopted SVG object tree. No PDF clipping.
Native data geometry is retained; page coordinates are tight-cropped, scaled
isotropically for legibility, and panel letters given a fixed outside anchor.
Statistical plots replaced by data-backed redraws are marked in panel_manifest.
"""
from pathlib import Path
import os,copy,io,json,re,hashlib
os.environ['DYLD_FALLBACK_LIBRARY_PATH']='/Users/taki/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/lib'
from lxml import etree as E
import numpy as np
from PIL import Image
import cairosvg
import pymupdf as fitz
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919'); O=P/'独立小图_最新版'
NS='http://www.w3.org/2000/svg';Q=lambda t:'{'+NS+'}'+t
ns={'s':NS};mm=72/25.4
MAN={
'Figure1':{'a':['axes_1'],'b':['axes_2','axes_3','axes_4','text_148'],'c':['axes_5'],'d':['axes_7'],'e':['preserved_e'],'f':['axes_6'],'g':['adopted_B_g_figure_1']},
'Figure2':{'a':['axes_1'],'b':['axes_2','line2d_109','text_178'],'c':['axes_5','axes_6'],'d':['axes_7'],'e':['axes_8'],'f':['axes_11']},
'Figure3':{'a':['identity_a_figure_1'],'b':['figure3_b_content','legend_2'],'c':['c_magma','c_pseudobulk'],'d':['figure3d_compact_forest'],'e':['axes_6','figure3e_three_symbol_key'],'f':['figure3f_three_evidence_columns']},
'Figure4':{'a':['axes_1','axes_2','axes_3','axes_4','patch_89','patch_90','patch_91']+['text_'+str(n) for n in range(99,107)],'b':['axes_5','legend_1'],'c':['axes_6','text_110'],'d':['axes_7','text_111'],'e':['axes_8','axes_9','text_113'],'f':['axes_10','axes_11','axes_12','Figure4f_module_key_0921']},
'Figure5':{'a':['r2_axes_'+str(i) for i in range(1,11)]+['r2_text_'+str(i) for i in [24,25,26,28,29,30,31,32]]+['r2_legend_3'],'b':['r2_axes_11']+['r2_text_'+str(i) for i in range(33,39)]+['r2_legend_1','r2_legend_2','Figure5b_pooled_key_0921'],'c':['FigureYa6_native_forest','legend_2'],'d':['axes_8']+['text_'+str(n) for n in range(138,148)]+['legend_3'],'e':['axes_9','text_156'],'f':['axes_10']+['text_'+str(n) for n in range(150,156)]+['legend_4']},
'ExtendedData_Figure1':{'a':['axes_2','axes_3','axes_4','axes_5'],'b':['axes_7'],'c':['axes_9'],'d':['axes_11']},
'ExtendedData_Figure2':{'a':['axes_1','axes_2'],'b':['newED2_figure_1']},
'ExtendedData_Figure3':{'a':['axes_1','axes_2','axes_3','axes_4']+['text_'+str(i) for i in range(146,150)],'b':['newED3_figure_1']},
'ExtendedData_Figure4':{'a':['axes_2','newED4_text_71'],'b':['newED4_axes_1','newED4_axes_2','newED4_axes_3'],'c':['axes_8','axes_10','text_127','text_128'],'d':['newED4_axes_4','newED4_axes_5','newED4_axes_6','newED4_text_72','newED4_text_73']},
'ExtendedData_Figure5':{'a':['axes_1','axes_2','axes_3','axes_5']+['text_'+str(i) for i in range(137,141)],'b':['newED5_figure_1']},
'ExtendedData_Figure6':{'a':['axes_1'],'b':['axes_2'],'c':['axes_3','axes_4','axes_5','text_83','legend_2'],'d':['axes_6']},
'ExtendedData_Figure7':{'a':['axes_1'],'b':['axes_3'],'c':['axes_6'],'d':['axes_8'],'e':['axes_10'],'f':['axes_12','text_143']},
'ExtendedData_Figure8':{'a':['axes_1'],'b':['axes_2'],'c':['axes_4'],'d':['axes_6'],'e':['axes_8'],'f':['axes_10'],'g':['axes_12']},
'ExtendedData_Figure9':{'a':['axes_1','text_94','text_95','text_96'],'b':['axes_2'],'c':['axes_3'],'d':['axes_4','axes_5']},
'ExtendedData_Figure10':{'a':['axes_1'],'b':['axes_2'],'c':['axes_5'],'d':['axes_3','legend_1'],'e':['axes_6','axes_7','text_89','text_90','text_91'],'f':['axes_4']}}

def extract(root,ids):
 existing={x.get('id') for x in root.iter() if x.get('id')};missing=set(ids)-existing
 if missing: raise ValueError(('missing SVG objects',missing))
 out=E.Element(root.tag,attrib=dict(root.attrib),nsmap=root.nsmap)
 # Collect every defs block: references often sit in axes/collections rather than root.
 defs=E.SubElement(out,Q('defs'));seen=set()
 for block in root.iter(Q('defs')):
  for el in block:
   identity=el.get('id') or E.tostring(el)
   if identity not in seen:defs.append(copy.deepcopy(el));seen.add(identity)
 def prune(el,chosen=False):
  if el.tag in [Q('defs'),Q('metadata')]:return None
  chosen=chosen or el.get('id') in ids
  if chosen:
   new=copy.deepcopy(el)
   for de in new.xpath('.//s:defs',namespaces=ns):de.getparent().remove(de)
   return new
  new=E.Element(el.tag,attrib=dict(el.attrib))
  for ch in el:
   if not isinstance(ch.tag,str):continue
   got=prune(ch)
   if got is not None:new.append(got)
  return new if len(new) else None
 for ch in root:
  if not isinstance(ch.tag,str):continue
  got=prune(ch)
  if got is not None:out.append(got)
 return out

def content_bbox(root):
 b=cairosvg.svg2png(bytestring=E.tostring(root),scale=2,background_color='white')
 im=np.array(Image.open(io.BytesIO(b)).convert('RGB'));mask=im.min(axis=2)<247
 ys,xs=np.where(mask);v=list(map(float,root.get('viewBox').split()))
 return [v[0]+xs.min()/im.shape[1]*v[2],v[1]+ys.min()/im.shape[0]*v[3],v[0]+(xs.max()+1)/im.shape[1]*v[2],v[1]+(ys.max()+1)/im.shape[0]*v[3]]

def native_save(root,figure,letter,ids):
 outdir=O/figure;outdir.mkdir(exist_ok=True)
 # Bounding by rendered visible content, including labels/keys, not original blank canvas.
 x0,y0,x1,y1=content_bbox(root)
 pdf=cairosvg.svg2pdf(bytestring=E.tostring(root))
 doc=fitz.open(stream=pdf,filetype='pdf')
 sizes=[s['size'] for b in doc[0].get_text('dict')['blocks'] if 'lines'in b for l in b['lines'] for s in l['spans'] if s['text'].strip()]
 base=float(np.percentile(sizes,10)) if sizes else 7
 scale=max(1.,7./base)
 # No nonuniform scaling: frozen images and circles retain their aspect.
 pad=10.;ww=(x1-x0)*scale+2*pad;hh=(y1-y0)*scale+2*pad+12
 final=E.Element(Q('svg'),nsmap={None:NS,'xlink':'http://www.w3.org/1999/xlink'},width=f'{ww}pt',height=f'{hh}pt',viewBox=f'0 0 {ww} {hh}')
 for el in list(root):
  if el.tag==Q('defs'):final.append(el)
 g=E.SubElement(final,Q('g'),transform=f'translate({pad},{pad+12}) scale({scale}) translate({-x0},{-y0})')
 for el in list(root):
  if el.tag!=Q('defs'):g.append(el)
 txt=E.SubElement(final,Q('text'),x='3',y='11',style="font: bold 10px Arial;fill:#20252A");txt.text=letter
 target=outdir/(figure+letter)
 E.ElementTree(final).write(str(target)+'.svg',encoding='utf-8',xml_declaration=True)
 cairosvg.svg2pdf(bytestring=E.tostring(final),write_to=str(target)+'.pdf')
 doc=fitz.open(str(target)+'.pdf');doc[0].get_pixmap(matrix=fitz.Matrix(2,2)).save(str(target)+'.png')
 return dict(figure=figure,panel=letter,mode='native_vector_reframe',source_objects=ids,original_font_p10=base,font_scale=scale,width_mm=ww/mm,height_mm=hh/mm,source_svg=str(P/'图件编辑资料/editable_svg'/f'{figure}.svg'),scientific_geometry='preserved; isotropic transform only',change='standalone visible-content frame; outside letter; readable text size; no page clipping')

def main():
 records=[]
 for figure,panels in MAN.items():
  root=E.parse(str(P/'图件编辑资料/editable_svg'/f'{figure}.svg')).getroot()
  for letter,ids in panels.items():
   try:
    r=extract(root,ids);rec=native_save(r,figure,letter,ids);records.append(rec);print(figure+letter,round(rec['width_mm'],1),round(rec['height_mm'],1),flush=True)
   except Exception as e:print('ERROR',figure+letter,repr(e),flush=True);raise
 (O/'review/native_panel_manifest.json').write_text(json.dumps(records,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
