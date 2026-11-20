from pathlib import Path
import xml.etree.ElementTree as E,copy,re,io,json,sys,hashlib
from PIL import Image
import fitz,numpy as np
B=Path('/Users/taki/Desktop/Panvascular_Content_Alignment_ST_Audit_20260919');O=Path('/Users/taki/Desktop/Panvascular_Marker_Spatial_Review_20260919')
sys.path.insert(0,'/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Local_Touchup_20260918/tools')
import cairosvg,cairosvg.image
from cairosvg.surface import cairo
cairosvg.image.IMAGE_RENDERING.update({'pixelated':cairo.FILTER_NEAREST,'crisp-edges':cairo.FILTER_NEAREST})
N='http://www.w3.org/2000/svg';E.register_namespace('',N);E.register_namespace('xlink','http://www.w3.org/1999/xlink');q=lambda s:'{'+N+'}'+s
r=E.parse(B/'editable_svg/Figure3.svg').getroot();old=copy.deepcopy(r);fig=next(e for e in r if e.get('id')=='figure_1')
def find(i,root=r):return next(e for e in root.iter() if e.get('id')==i)
def xyshift(e):
 return np.array([sum(float(v[k]) for v in re.findall(r'translate\(([-\d.]+)[ ,]+([-\d.]+)\)',e.get('transform',''))) for k in [0,1]])
def data_box(e):
 p=next(t for t in e[0] if t.tag==q('path'));v=np.array([float(x) for x in re.findall(r'-?\d+(?:\.\d+)?',p.get('d'))]).reshape(-1,2)+xyshift(e);return [v[:,0].min(),v[:,1].min(),v[:,0].max(),v[:,1].max()]
def bbox(e):
 rr=E.Element(q('svg'),dict(r.attrib));defs=E.SubElement(rr,q('defs'))
 for dd in r.iter(q('defs')):
  for ch in dd:defs.append(copy.deepcopy(ch))
 cp=copy.deepcopy(e)
 for gg in cp.iter(q('g')):
  if gg.get('id','').startswith('axes_') and len(gg) and 'fill: #ffffff' in E.tostring(gg[0],encoding='unicode'):gg.remove(gg[0])
 rr.append(cp);im=Image.open(io.BytesIO(cairosvg.svg2png(bytestring=E.tostring(rr),output_width=2075,output_height=1928))).convert('RGBA');b=im.getchannel('A').getbbox();return [b[0]*518.740157/2075,b[1]*481.889764/1928,b[2]*518.740157/2075,b[3]*481.889764/1928]
d=find('axes_5');e=find('axes_6');f=find('axes_7');db0=data_box(d);eb=data_box(e);target_gap=4*72/25.4;dx=eb[0]-target_gap-db0[2];d.set('transform',f'translate({dx} 0) '+d.get('transform',''))
# Additional key explains the ORIGINAL categorical classification fill rule.
# Source: build_root_figures.py:80; classification == strong_colocalisation.
leg=E.SubElement(fig,q('g'),{'id':'figure3_f_original_fill_key'})
for y,filled,label in [(442,True,'Strong support'),(453,False,'Other classification')]:
 E.SubElement(leg,q('circle'),{'cx':'414','cy':str(y),'r':'1.5','fill':'#356b8c' if filled else 'white','stroke':'#356b8c','stroke-width':'.6'})
 E.SubElement(leg,q('rect'),{'x':'420','y':str(y-1.5),'width':'3','height':'3','fill':'#555555' if filled else 'white','stroke':'#555555','stroke-width':'.6'})
 t=E.SubElement(leg,q('text'),{'x':'428','y':str(y+2.1),'style':"fill:#222222;font:6px 'Arial';text-anchor:start"});t.text=label
# b/c unchanged including parent group transforms and every original panel-letter anchor.
for i in ['figure3_b_content','figure3_c_content','legend_2','axes_7','legend_4']+[f'text_{n}' for n in [82,83,103,106,107,108]]:assert E.tostring(find(i))==E.tostring(find(i,old)),i
assert xyshift(d)[1]==xyshift(find('axes_5',old))[1] and xyshift(e)[1]==xyshift(find('axes_6',old))[1]
for name in ['axes_5','axes_6']:
 aa=copy.deepcopy(find(name));bb=copy.deepcopy(find(name,old));aa.attrib.pop('transform',None);bb.attrib.pop('transform',None);assert E.tostring(aa)==E.tostring(bb)
for sub in ['main','editable_svg','review','audit','notes']:(O/sub).mkdir(exist_ok=True)
svg=O/'editable_svg/Figure3.svg';E.ElementTree(r).write(svg,encoding='utf-8',xml_declaration=True);pdf=O/'main/Figure3.pdf';cairosvg.svg2pdf(url=str(svg),write_to=str(pdf));p=fitz.open(pdf)[0];p.get_pixmap(matrix=fitz.Matrix(2,2)).save(O/'review/figure3_after.png')
checks={'baseline':str(B),'d_horizontal_translation_pt':dx,'actual_d_data_bbox_before_pt':db0,'actual_d_data_bbox_pt':data_box(d),'actual_e_data_bbox_pt':data_box(e),'actual_f_data_bbox_pt':data_box(f),'de_data_gap_mm':(data_box(e)[0]-data_box(d)[2])*25.4/72,'de_all29_absolute_y_unchanged':True,'panel_letters_unchanged':True,'b_c_unchanged':True,'f_original_objects_unchanged':True,'original_fill_rule':"classification == 'strong_colocalisation'",'rule_source':'Panvascular_Nature_Reconstruction_20260918/code/build_root_figures.py:80','new_threshold_used':False,'allocated_panel_bbox_pt':{'de_as_pair':[40,220,510,424],'d':[40,220,359,424],'e':[359,220,510,424],'f':[12,424,510,481.889764]},'complete_content_bbox_pt':{'d':bbox(d),'e':bbox(e),'e_legend':bbox(find('legend_3')),'f':bbox(f),'f_scenario_legend':bbox(find('legend_4')),'f_fill_legend':bbox(leg)}}
checks['d_f_visible_gap_mm']=(bbox(f)[1]-bbox(d)[3])*25.4/72
checks['sha256']={str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in [B/'editable_svg/Figure3.svg',svg,pdf]}
checks['minimum_font_pt']=min(s['size'] for b in p.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans'])
(O/'audit/figure3_columns_fill_checks.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))
