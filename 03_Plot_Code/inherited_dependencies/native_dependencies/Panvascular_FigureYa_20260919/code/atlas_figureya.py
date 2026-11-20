"""Display-only port of FigureYa278 geom_tile + geom_point(shape=1), and
243 fixed-order grouped heatmap white borders. Never executes R analysis.
Uses approved SVG's exact rendered scalar colours, circle paths and coordinates.
"""
from pathlib import Path
import sys,re,json,hashlib,os,xml.etree.ElementTree as E
os.environ.setdefault("DYLD_FALLBACK_LIBRARY_PATH","/Users/taki/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/lib")
sys.path.insert(0,'/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Local_Touchup_20260918/tools')
import cairosvg,fitz
from PIL import Image,ImageOps,ImageDraw
B=Path('/Users/taki/Desktop/投稿包_最新版_20260919'); O=Path('/Users/taki/Desktop/Panvascular_FigureYa_20260919'); NS='http://www.w3.org/2000/svg'; N={'s':NS}; E.register_namespace('',NS);E.register_namespace('xlink','http://www.w3.org/1999/xlink')
sha=lambda b:hashlib.sha256(b).hexdigest()
num=lambda s:list(map(float,re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?',s)))
audit=[]
for name in ['Figure1']+[f'ExtendedData_Figure{i}' for i in range(2,6)]:
 src=B/'图件编辑资料/editable_svg'/f'{name}.svg';tree=E.parse(src);r=tree.getroot();images=[sha(E.tostring(x)) for x in r.findall('.//s:image',N)]; modifications=[];count=0
 for col in list(r.findall('.//s:g',N)):
  if not col.get('id','').startswith('PathCollection_'):continue
  paths=list(col.findall('s:path',N)); vals=[]
  for p in paths:
   xy=num(p.get('d',''));xs=xy[::2];ys=xy[1::2]
   if not xs:continue
   vals.append((p,(min(xs)+max(xs))/2,(min(ys)+max(ys))/2))
  if not vals:continue
  xs=sorted(set(round(x,4) for p,x,y in vals));ys=sorted(set(round(y,4) for p,x,y in vals))
  pitchx=min([b-a for a,b in zip(xs,xs[1:]) if b-a>.01],default=12);pitchy=min([b-a for a,b in zip(ys,ys[1:]) if b-a>.01],default=6.2)
  under=E.Element('{'+NS+'}g',{'id':'FigureYa278_tiles_'+col.get('id')})
  for p,x,y in vals:
   color=re.search(r'fill:\s*([^;]+)',p.get('style','')); color=color.group(1) if color else '#ffffff'
   tile=E.SubElement(under,'{'+NS+'}rect',{'x':str(x-pitchx*.48),'y':str(y-pitchy*.46),'width':str(pitchx*.96),'height':str(pitchy*.92),'fill':color,'stroke':'white','stroke-width':'.25'})
   if p.get('clip-path'):tile.set('clip-path',p.get('clip-path'))
   rgb=tuple(int(color.lstrip('#')[i:i+2],16)/255 for i in (0,2,4))
   linear=[v/12.92 if v<=0.04045 else ((v+0.055)/1.055)**2.4 for v in rgb]
   luminance=sum(v*w for v,w in zip(linear,(.2126,.7152,.0722)))
   edge='#ffffff' if luminance<.28 else '#374151'
   p.set('style',f'fill: none; stroke: {edge}; stroke-width: 0.28')
  col.insert(0,under);count+=len(vals)
  modifications.append({'collection':col.get('id'),'n_cells':len(vals),'tile_width':pitchx*.96,'tile_height':pitchy*.92,'encoding':'existing mean colour moved to tile; exact detection circle path retained as outline','circle_edge_rule':'white if relative tile luminance < 0.28; dark gray otherwise; contrast only, no data semantics'})
 # Area legends must match hollow circles. No category keys or UMAP objects edited.
 if count:
  for g in r.findall('.//s:g',N):
   if not g.get('id','').startswith('legend_'):continue
   for definition in g.findall('.//s:defs/s:path',N):
    definition.set('style','fill:none;stroke:#374151;stroke-width:0.4;stroke-opacity:1')
   for x in g.findall('.//s:use',N):
    x.set('style','fill: none; stroke: #374151; stroke-width: 0.4')
  # Explain changed encoding within existing legend labels, without a panel subtitle.
  for t in r.findall('.//s:text',N):
   if t.text=='Detected cells (%)':t.text='Circle: detected cells (%)'
 # Figure1 is frozen approved composition; make common count rows easier to trace with fine rules.
 if name=='Figure1':
  groups={g.get('id'):g for g in r.findall('.//s:g',N)}
  labels=groups['axes_3'].findall('.//s:text',N)
  ys=[]
  for t in labels:
   if t.get('y'):ys.append(float(t.get('y'))-1.8)
  rules=E.Element('{'+NS+'}g',{'id':'FigureYa243_shared_lineage_rules'})
  for y in ys:
   E.SubElement(rules,'{'+NS+'}path',{'d':f'M 7 {y+3.1} L 308 {y+3.1}','style':'fill:none;stroke:#e8edf0;stroke-width:0.3'})
  # Use common-row matrix separators only when literal coordinates exist.
  groups['axes_3'].insert(0,rules);modifications.append({'panel':'c/f','shared_row_rules':len(ys),'reason':'fixed-order row tracking, approved cells/studies bars remain'})
 out=O/'editable_svg'/f'{name}.svg';tree.write(out,encoding='utf-8',xml_declaration=True)
 assert images==[sha(E.tostring(x)) for x in r.findall('.//s:image',N)]
 pdf=O/('main' if name=='Figure1' else 'extended_data')/f'{name}.pdf';cairosvg.svg2pdf(url=str(out),write_to=str(pdf))
 doc=fitz.open(pdf);page=doc[0];pix=page.get_pixmap(matrix=fitz.Matrix(1.5,1.5));pix.save(str(O/'review'/f'atlas_{name}.png'))
 base=fitz.open(B/('Figures' if name=='Figure1' else 'Extended_Data')/f'{name}.pdf');bp=base[0].get_pixmap(matrix=fitz.Matrix(1.5,1.5));bi=Image.frombytes('RGB',[bp.width,bp.height],bp.samples);ai=Image.open(O/'review'/f'atlas_{name}.png').convert('RGB');canvas=Image.new('RGB',(bi.width+ai.width+30,max(bi.height,ai.height)+35),'white');canvas.paste(bi,(0,35));canvas.paste(ai,(bi.width+30,35));d=ImageDraw.Draw(canvas);d.text((10,10),'Approved baseline',fill='black');d.text((bi.width+40,10),'FigureYa display adaptation',fill='black');canvas.save(O/'review'/f'atlas_{name}_before_after.png')
 audit.append({'figure':name,'svg_source':str(src),'source_sha256':sha(src.read_bytes()),'output_sha256':sha(out.read_bytes()),'unchanged_embedded_images':len(images),'matrix_glyphs':count,'changes':modifications,'page_mm':[round(page.rect.width*25.4/72,3),round(page.rect.height*25.4/72,3)]})
 print(name,count,flush=True)
(O/'review/atlas_integrity.json').write_text(json.dumps(audit,indent=2))
