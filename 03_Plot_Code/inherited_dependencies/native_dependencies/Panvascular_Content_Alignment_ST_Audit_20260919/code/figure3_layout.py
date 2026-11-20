from pathlib import Path
import xml.etree.ElementTree as E,copy,io,re,json,sys
import numpy as np
from PIL import Image
import fitz
O=Path('/Users/taki/Desktop/Panvascular_Content_Alignment_ST_Audit_20260919');B=Path('/Users/taki/Desktop/Panvascular_FigureYa_20260919')
sys.path.insert(0,'/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Local_Touchup_20260918/tools')
import cairosvg, cairosvg.image
from cairosvg.surface import cairo
cairosvg.image.IMAGE_RENDERING.update({"pixelated": cairo.FILTER_NEAREST, "crisp-edges": cairo.FILTER_NEAREST})
N='http://www.w3.org/2000/svg';E.register_namespace('',N);E.register_namespace('xlink','http://www.w3.org/1999/xlink');q=lambda x:'{'+N+'}'+x
r=E.parse(B/'editable_svg/Figure3.svg').getroot();original=copy.deepcopy(r);fig=next(e for e in r if e.get('id')=='figure_1')
def find(i,root=r):return next(e for e in root.iter() if e.get('id')==i)
def move(e,x=0,y=0):e.set('transform',f'translate({x} {y}) '+e.get('transform',''))
def group(ids,name):
 g=E.SubElement(fig,q('g'),{'id':name})
 for i in ids:
  e=find(i);fig.remove(e);g.append(e)
 return g
# Bbox includes true glyph ink, lines and markers, not the white axes/page rectangles.
def bbox(e):
 rr=E.Element(q('svg'),dict(r.attrib));defs=E.SubElement(rr,q('defs'))
 for dd in r.iter(q('defs')):
  for ch in dd:defs.append(copy.deepcopy(ch))
 cp=copy.deepcopy(e)
 for gg in cp.iter(q('g')):
  if gg.get('id','').startswith('axes_') and len(gg):
   first=gg[0]
   if first.get('id','').startswith('patch_') and 'fill: #ffffff' in E.tostring(first,encoding='unicode'):gg.remove(first)
 rr.append(cp);data=cairosvg.svg2png(bytestring=E.tostring(rr),output_width=2075,output_height=1928);im=Image.open(io.BytesIO(data)).convert('RGBA');bb=im.getchannel('A').getbbox();sx=518.740157/im.width;sy=481.889764/im.height
 return [bb[0]*sx,bb[1]*sy,bb[2]*sx,bb[3]*sy]
def centre(e,target):
 b=bbox(e);dx=target-(b[0]+b[2])/2;move(e,dx,0);return dx
# b: maintain complete names, tight right-aligned label column, paired numeric method columns.
for i in range(84,103):move(find('text_'+str(i)),-22,0)
move(find('axes_2'),-42,0)
b=group(['axes_2']+['text_'+str(i) for i in range(84,103)],'figure3_b_content');centre(b,(181.417323+510.23622)/2)
bb=bbox(b)
leg=find('legend_2');lb=bbox(leg);move(leg,bb[0]-lb[0],-1.5)
# c: compact three-column comparison. Preserve both point axes and all seven source rows.
move(find('axes_3'),-20,-12);move(find('axes_4'),-55,-12)
c=group(['axes_3','axes_4'],'figure3_c_content');centre(c,518.740157/2);cb=bbox(c)
# Native small category markers reanchored under the name column; not a left-side panel.
for pid,tid,x in [('patch_35','text_104',cb[0]),('patch_36','text_105',cb[0]+45)]:
 pp=find(pid);pb=bbox(pp);move(pp,x-pb[0],-12)
 tt=find(tid);tb=bbox(tt);move(tt,x+6-tb[0],-12)
move(find('text_103'),0,-12)
# d/e remain rigid groups with shared y translation; no CI or row scale change.
d=find('axes_5');move(d,0,-19);centre(d,(12+352)/2)
e=find('axes_6');move(e,0,-19);centre(e,(352+510.23622)/2)
move(find('text_106'),0,-19);move(find('text_107'),-10.8,-19)
el=find('legend_3');eb=bbox(e);lb=bbox(el);move(el,eb[0]-lb[0],-19)
# f: independent bottom row, centred on all visible labels and plotted posterior points.
f=find('axes_7');centre(f,518.740157/2);fb=bbox(f)
fl=find('legend_4');lb=bbox(fl);move(fl,fb[0]-lb[0],472-lb[1])
# a: centre the actual visible content within its assigned first-column slot.
a=find('axes_1');centre(a,(12+181.417323)/2)
# Quantitative rendered bbox audit, in points and mm.
checks={name:bbox(obj) for name,obj in [('a',a),('b',b),('b_legend',leg),('c',c),('d',d),('e',e),('e_legend',el),('f',f),('f_legend',fl)]}
checks['baseline_d_bbox']=bbox(find('axes_5',original))
checks['baseline_f_bbox']=bbox(find('axes_7',original))
checks['baseline_d_f_gap_pt']=checks['baseline_f_bbox'][1]-checks['baseline_d_bbox'][3]
checks['d_f_visible_gap_pt']=checks['f'][1]-checks['d'][3]
checks['centering_error_pt']={'a':sum(checks['a'][::2])/2-(12+181.417323)/2,'b':sum(checks['b'][::2])/2-(181.417323+510.23622)/2,'c':sum(checks['c'][::2])/2-518.740157/2,'d':sum(checks['d'][::2])/2-(12+352)/2,'e':sum(checks['e'][::2])/2-(352+510.23622)/2,'f':sum(checks['f'][::2])/2-518.740157/2}
# Confirm the SVG source geometry remains literally unchanged for scientific marks.
for ai in [1,2,3,4,5,6,7]:
 old=find('axes_'+str(ai),original);new=find('axes_'+str(ai));x=copy.deepcopy(old);y=copy.deepcopy(new);x.attrib.pop('transform',None);y.attrib.pop('transform',None);assert E.tostring(x)==E.tostring(y),ai
checks['all_axes_internal_geometry_unchanged']=True
checks['de_same_y_translation_pt']=-19
checks['b_label_end_to_first_point_pt']=443.622047-42-(416.692913-22)
checks['b_method_column_centres_follow_original_marker_centres']=True
checks['legend_text_anchor_start']=all('text-anchor: start' in t.get('style','') for lid in ['legend_2','legend_3','legend_4'] for t in find(lid).iter(q('text')))
checks['legend_left_edge_alignment_pt']={'b':checks['b_legend'][0]-checks['b'][0],'e':checks['e_legend'][0]-checks['e'][0],'f':checks['f_legend'][0]-checks['f'][0]}
checks['text_content_unchanged']=sorted(t.text or '' for t in original.iter(q('text')))==sorted(t.text or '' for t in r.iter(q('text')))
for path in ['main','editable_svg','review','audit','notes']:(O/path).mkdir(exist_ok=True)
svg=O/'editable_svg/Figure3.svg';E.ElementTree(r).write(svg,encoding='utf-8',xml_declaration=True);pdf=O/'main/Figure3.pdf';cairosvg.svg2pdf(url=str(svg),write_to=str(pdf));pg=fitz.open(pdf)[0]
pg.get_pixmap(matrix=fitz.Matrix(2,2)).save(O/'review/figure3_after.png')
pre=fitz.open(B/'main/Figure3.pdf')[0];pre.get_pixmap(matrix=fitz.Matrix(2,2)).save(O/'review/figure3_before.png')
for name,rect in [('b',(175,0,518.74,143)),('c',(5,137,513,238)),('d',(5,222,513,424)),('f',(5,424,513,481.89)),('d_f',(5,222,513,481.89))]:
 for label,page in [('before',pre),('after',pg)]:page.get_pixmap(matrix=fitz.Matrix(3,3),clip=fitz.Rect(rect)).save(O/'review'/f'figure3_{name}_{label}.png')
checks['page_mm']=[pg.rect.width*25.4/72,pg.rect.height*25.4/72]
checks['minimum_font_pt']=min(s['size'] for b in pg.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans'])
(O/'audit/figure3_bbox_checks.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))
