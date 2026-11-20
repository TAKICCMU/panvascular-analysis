from pathlib import Path
import os,sys,io,copy,json,re
os.environ.setdefault('DYLD_FALLBACK_LIBRARY_PATH','/Users/taki/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/lib')
sys.path.insert(0,'/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Local_Touchup_20260918/tools')
import cairosvg,fitz,numpy as np
import cairosvg.image
from cairosvg.surface import cairo
cairosvg.image.IMAGE_RENDERING.update({"pixelated":cairo.FILTER_NEAREST,"crisp-edges":cairo.FILTER_NEAREST})
from PIL import Image,ImageDraw
import xml.etree.ElementTree as E
B=Path('/Users/taki/Desktop/Panvascular_FigureYa_20260919');O=Path('/Users/taki/Desktop/Panvascular_Content_Alignment_ST_Audit_20260919');N={'s':'http://www.w3.org/2000/svg'};S=N['s'];E.register_namespace('',S);E.register_namespace('xlink','http://www.w3.org/1999/xlink');MM=72/25.4
names=['Figure1','Figure2','Figure4']+[f'ExtendedData_Figure{i}' for i in range(1,11)]
P={
'Figure1': [('a',0,0,183,36),('b',0,36,183,85),('c+f',0,85,112,138),('d',112,85,183,138),('e',0,138,183,170)],
'Figure2':[('a',0,0,183,22),('b',0,22,77,78),('c',77,22,183,78),('d',0,78,183,130),('e',0,130,95,170),('f',95,130,183,170)],
'Figure4':[('a',0,0,183,27),('b',0,27,183,56),('c',0,56,91,94),('d',91,56,183,94),('e',0,94,183,125),('f',0,125,183,170)],
'ExtendedData_Figure1':[('a',0,0,180,50),('b',0,50,180,96),('c',0,96,119,170),('d',119,96,180,170)],
'ExtendedData_Figure2':[('a',0,0,180,74),('b',0,74,180,170)],
'ExtendedData_Figure3':[('a',0,0,180,44),('b',0,44,180,170)],
'ExtendedData_Figure4':[('a',0,0,62,86),('b',62,0,180,86),('c',0,86,62,170),('d',62,86,180,170)],
'ExtendedData_Figure5':[('a',0,0,180,55),('b',0,55,180,170)],
'ExtendedData_Figure6':[('a',0,0,90,68),('b',90,0,180,68),('c',0,68,180,118),('d',0,118,180,170)],
'ExtendedData_Figure7':[('a',0,0,90,52),('b',90,0,180,52),('c',0,52,90,105),('d',90,52,180,105),('e',0,105,77,170),('f',77,105,180,170)],
'ExtendedData_Figure8':[('a',0,0,180,18),('b',0,18,90,75),('c',90,18,180,75),('d',0,75,90,134),('e',90,75,180,134),('f',0,134,78,170),('g',78,134,180,170)],
'ExtendedData_Figure9':[('a',0,0,180,76),('b',0,76,90,139),('c',90,76,180,139),('d',0,139,180,170)],
'ExtendedData_Figure10':[('a',0,0,36,85),('b',36,0,68,85),('c+e',68,0,180,170),('d',0,85,68,137),('f',0,137,68,170)]}

def bbox(root,element):
 rr=E.Element(root.tag,root.attrib)
 # Definitions may be local to other groups (SVG use references); collect all definitions.
 defs=E.SubElement(rr,'{'+S+'}defs')
 for d in root.findall('.//s:defs',N):
  for ch in d:defs.append(copy.deepcopy(ch))
 rr.append(copy.deepcopy(element))
 png=cairosvg.svg2png(bytestring=E.tostring(rr),scale=1,background_color='white')
 ar=np.asarray(Image.open(io.BytesIO(png)).convert('RGB'));mask=ar.min(axis=2)<245
 yy,xx=np.where(mask)
 return [float(xx.min())*.75,float(yy.min())*.75,float(xx.max()+1)*.75,float(yy.max()+1)*.75] if len(xx) else None

def union(boxes):return [min(x[0] for x in boxes),min(x[1] for x in boxes),max(x[2] for x in boxes),max(x[3] for x in boxes)]
def shift(g,x,y):g.set('transform',f'translate({x:.5f} {y:.5f}) '+g.get('transform',''))
def assigned(b,panels):
 x=(b[0]+b[2])/2/MM;y=(b[1]+b[3])/2/MM
 for p in panels:
  if p[1]<=x<=p[3] and p[2]<=y<=p[4]:return p[0]
 return min(panels,key=lambda p:abs(x-(p[1]+p[3])/2)+abs(y-(p[2]+p[4])/2))[0]
audit=[];legendaudit=[]
for name in names:
 tree=E.parse(B/'editable_svg'/f'{name}.svg');root=tree.getroot();fig=root.find("s:g[@id='figure_1']",N);panels=P[name];boxes={};members={p[0]:[] for p in panels};letters=[]
 for el in list(fig):
  if el.get('id')=='patch_1':continue
  texts=el.findall('.//s:text',N)
  if len(texts)==1 and texts[0].text in list('abcdefghi') and '700 8' in texts[0].get('style',''):letters.append(el);continue
  b=bbox(root,el)
  if b:boxes[el]=b;members[assigned(b,panels)].append(el)
 # Legend texts already mostly left anchored; normalize all and record exact original attributes.
 for legend in root.findall('.//s:g',N):
  if not legend.get('id','').startswith('legend_'):continue
  bb=bbox(root,legend)
  if not bb:continue
  panel=assigned(bb,panels);target=union([boxes[x] for x in members[panel]])[0]
  beforetexts=[]
  for t in legend.findall('.//s:text',N):
   st=t.get('style','');m=re.search('text-anchor: ([^;]+)',st);old=m.group(1) if m else 'start';beforetexts.append([t.text,old])
   if t.text and t.text.startswith('Circle:'):
    tb=bbox(root,t)
    if tb: t.set('transform',f'translate({bb[0]-tb[0]:.5f} 0) '+t.get('transform',''))
   if old!='start':t.set('style',re.sub(r'text-anchor: [^;]+','text-anchor: start',st))
  dx=target-bb[0]
  if name in ['Figure2','ExtendedData_Figure2','ExtendedData_Figure8']:
   dx=0 # Keep compound legend clear of existing detached title / in-panel text.
  # Horizontal legends can be left aligned safely without changing their vertical placement.
  shift(legend,dx,0)
  legendaudit.append(dict(figure=name,legend=legend.get('id'),panel=panel,before_bbox_pt=bb,translation_pt=[dx,0],internal_anchors_before=beforetexts,internal_anchor_after='start',reason='legend block left edge aligned to visible labelled data block; labels left anchored'))
 # Re-measure after legend alignment, then center visible content as a single linked group.
 for p in panels:
  key,x0,y0,x1,y1=p;els=members[key]
  if not els:continue
  bs=[bbox(root,el) for el in els];bb=union([b for b in bs if b]);target=[(x0+3)*MM,(y0+6)*MM,(x1-3)*MM,(y1-2)*MM]
  dx=(target[0]+target[2]-bb[0]-bb[2])/2;dy=(target[1]+target[3]-bb[1]-bb[3])/2
  pagew,pageh=map(float,root.get('viewBox').split()[2:]);dx=max(2*MM-bb[0],min(dx,pagew-2*MM-bb[2]));dy=max(2*MM-bb[1],min(dy,pageh-2*MM-bb[3]))
  reason='visible ink bbox centered in allocated panel; labels and legend travel together'
  # Keep reviewed workflows and full-width tiny summary strip completely intact.
  if key=='a' and (name.startswith('Figure') or name=='ExtendedData_Figure8'):dx=dy=0;reason='approved workflow/summary structure retained; no need to disrupt existing balanced content'
  # Preserve shared geometry and dense panels; translation does not change scale.
  for el in els:shift(el,dx,dy)
  after=[bb[0]+dx,bb[1]+dy,bb[2]+dx,bb[3]+dy]
  audit.append(dict(figure=name,panel=key,original_bbox_pt=union([boxes[x] for x in els]),bbox_after_legend_pt=bb,new_bbox_pt=after,allocated_mm=[x0,y0,x1,y1],translation_pt=[dx,dy],reason=reason,members=[el.get('id') for el in els]))
 for letter in letters:
  fig.remove(letter);fig.append(letter)
 svg=O/'editable_svg'/f'{name}.svg';tree.write(svg,encoding='utf-8',xml_declaration=True)
 pdf=O/('main' if name.startswith('Figure') else 'extended_data')/f'{name}.pdf';cairosvg.svg2pdf(url=str(svg),write_to=str(pdf))
 d=fitz.open(pdf);pm=d[0].get_pixmap(matrix=fitz.Matrix(1.3,1.3));im=Image.frombytes('RGB',[pm.width,pm.height],pm.samples);im.save(O/'review'/f'alignment_{name}.png')
 d0=fitz.open(B/('main' if name.startswith('Figure') else 'extended_data')/f'{name}.pdf');pm=d0[0].get_pixmap(matrix=fitz.Matrix(1.3,1.3));im0=Image.frombytes('RGB',[pm.width,pm.height],pm.samples);out=Image.new('RGB',(im.width+im0.width+20,max(im.height,im0.height)+25),'white');out.paste(im0,(0,25));out.paste(im,(im0.width+20,25));dr=ImageDraw.Draw(out);dr.text((5,5),'BEFORE',fill='black');dr.text((im0.width+25,5),'AFTER',fill='black');out.save(O/'review'/f'alignment_{name}_comparison.png');print(name,flush=True)
 (O/'audit/alignment_panels.json').write_text(json.dumps(audit,indent=2));(O/'audit/alignment_legends.json').write_text(json.dumps(legendaudit,indent=2))
