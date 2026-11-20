"""Exact local mirror for BioRender assembly; no statistical recalculation.
Letters are independent SVG text objects in the review only, not panel assets.
"""
from pathlib import Path
import json,copy,os
os.environ['DYLD_FALLBACK_LIBRARY_PATH']='/opt/homebrew/lib'
from lxml import etree as E
import cairosvg,pymupdf as fitz
R=Path(__file__).resolve().parent
NS='http://www.w3.org/2000/svg';Q=lambda x:'{'+NS+'}'+x
MM=72/25.4
manifest={r['panel']:r for r in json.loads((R/'panel_manifest.json').read_text())}
W,H=183,210
placements=[
 ('a',1.5,3.0,180),
 ('b',1.5,42.0,78),
 ('c',83,45.0,98.5),
 ('d',1.5,109.0,180),
 ('e',1.5,160.0,64),
 ('f',68.0,160.0,39),
 ('g',109.5,160.0,72),
]
root=E.Element(Q('svg'),nsmap={None:NS,'xlink':'http://www.w3.org/1999/xlink'},width=f'{W}mm',height=f'{H}mm',viewBox=f'0 0 {W*MM} {H*MM}')
E.SubElement(root,Q('rect'),x='0',y='0',width=str(W*MM),height=str(H*MM),fill='white')
record=[]
for l,x,y,w in placements:
 name='Figure2'+l;src=R/'panels'/f'{name}.svg';r=E.parse(str(src)).getroot();meta=manifest[name];h=w/meta['width_mm']*meta['height_mm']
 v=list(map(float,r.get('viewBox').split()))
 group=E.SubElement(root,Q('svg'),id=f'panel_{l}',x=str(x*MM),y=str(y*MM),width=str(w*MM),height=str(h*MM),viewBox=' '.join(map(str,v)),overflow='hidden',preserveAspectRatio='xMidYMid meet')
 # Namespace every local SVG id before embedding so clip paths never collide.
 ids={el.get('id'):'fig2'+l+'_'+el.get('id') for el in r.iter() if el.get('id')}
 for el in r.iter():
  for k,val in list(el.attrib.items()):
   if k=='id':el.set(k,ids[val])
   else:
    for old,new in ids.items():val=val.replace('url(#'+old+')','url(#'+new+')').replace('#'+old if val=='#'+old else '\0','#'+new)
    el.set(k,val)
 for child in r:group.append(copy.deepcopy(child))
 ly={'a':2.5,'b':41,'c':41,'d':108,'e':159,'f':159,'g':159}[l]
 lx=x
 label=E.SubElement(root,Q('text'),x=str(lx*MM),y=str(ly*MM),id='native_panel_letter_'+l,**{'font-family':'Arial','font-size':'8','font-weight':'700','fill':'#20252A'});label.text=l
 record.append(dict(panel=name,x_mm=x,y_mm=y,width_mm=w,height_mm=h,letter_x_mm=lx,letter_baseline_y_mm=ly,native_BioRender_letter_required=True))
assert max(r['y_mm']+r['height_mm'] for r in record)<H
out=R/'review/Figure2_current_preview'
E.ElementTree(root).write(str(out)+'.svg',encoding='UTF-8',xml_declaration=True)
cairosvg.svg2pdf(bytestring=E.tostring(root),write_to=str(out)+'.pdf')
d=fitz.open(str(out)+'.pdf');d[0].get_pixmap(matrix=fitz.Matrix(2.5,2.5)).save(str(out)+'.png')
(R/'panel_placement.json').write_text(json.dumps({'page_width_mm':W,'page_height_mm':H,'status':'REVIEW_LAYOUT_NOT_170MM_PRODUCTION_SIZE','panels':record},ensure_ascii=False,indent=2))
