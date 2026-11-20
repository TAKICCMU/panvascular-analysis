"""Figure 4a display-only schematic. All drawn motifs are schematic, not data.
Reuses panels b–d byte-for-byte; no scientific analysis or data manipulation.
"""
from pathlib import Path
from lxml import etree as E
import copy,json,hashlib,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Circle, FancyBboxPatch
import pymupdf as fitz

O=Path(__file__).resolve().parents[1]
BASE=Path('/Users/taki/Desktop/泛血管稿件_批注修订_20261116')
OLD=BASE/'图件编辑资料/Figure4/panels'
MM=72/25.4
INK='#26313A';GRAY='#64717B';GUIDE='#D7E0E4';TEAL='#248C85';VIOLET='#8D7298';SLATE='#526D8C'
plt.rcParams.update({'font.family':'Arial','font.size':6.5,'svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42,'text.color':INK,'svg.hashsalt':'Figure4a_20261117','savefig.facecolor':'white'})
f=plt.figure(figsize=(183/25.4,25/25.4));a=f.add_axes([0,0,1,1]);a.set(xlim=(0,183),ylim=(25,0));a.axis('off')
def tx(x,y,s,size=6.5,**kw):
 kw.setdefault('ha','center');kw.setdefault('va','center');kw.setdefault('zorder',5);return a.text(x,y,s,fontsize=size,**kw)
def line(x0,y0,x1,y1,c=GRAY,lw=.7,z=2):a.plot([x0,x1],[y0,y1],c=c,lw=lw,solid_capstyle='round',zorder=z)
def arrow(x0,x1,y=13):a.add_patch(FancyArrowPatch((x0,y),(x1,y),arrowstyle='-|>',mutation_scale=7,lw=.8,color=GRAY,zorder=4))
def rect(x,y,w,h,fc,ec='none',lw=.5):a.add_patch(Rectangle((x,y),w,h,fc=fc,ec=ec,lw=lw))
def roundbox(x,y,w,h,fc,ec,lw=.5,r=1.4):
 a.add_patch(FancyBboxPatch((x,y),w,h,boxstyle=f'round,pad=0,rounding_size={r}',fc=fc,ec=ec,lw=lw,zorder=0))
def dot(x,y,r,col,ec='white',lw=.25,z=4):a.add_patch(Circle((x,y),r,fc=col,ec=ec,lw=lw,zorder=z))
# Two independent branches. Motifs encode workflow, not measured values.
# Each branch contains a larger input illustration, explicit modelling operation
# and a graphically differentiated output. Panel letter remains outside this SVG.
roundbox(3,1.1,86.3,22.8,'#F4F8F7','#D2E1DE',.55)
roundbox(92.3,1.1,89.7,22.8,'#F4F6F9','#D5DDE7',.55)
tx(6,4.1,'CoVarNet',7.8,fontweight='bold',ha='left')
tx(95.3,4.1,'Gene expression in tissue',7.8,fontweight='bold',ha='left')
# Cell proportions: rows = samples, columns = states. Conceptual matrix only.
mx,my=8,8.3
pals=[['#C6DED9','#A8C9C5','#E0E8EE','#CCD9E4','#E5DBEB','#CCBBD6'],
      ['#BCD9D3','#94BEB8','#D5E0EA','#B6CADA','#D4C4DF','#BBA3C9'],
      ['#D6E8E3','#BBD6D0','#AABFD0','#829EB8','#DDCEE6','#CDB9DA'],
      ['#C5DED9','#A8CEC7','#BCD0DE','#9DB6CA','#C7B1D3','#A488B6']]
for r,row in enumerate(pals):
 for c,col in enumerate(row):rect(mx+c*2.9,my+r*2.0,2.6,1.7,col)
# Samples and states are shown as non-numerical row/column identity tabs.
for r in range(4):rect(mx-1.5,my+r*2,0.65,1.7,GRAY)
for c in range(6):rect(mx+c*2.9,my-1.3,2.6,.55,[TEAL,TEAL,SLATE,SLATE,VIOLET,VIOLET][c])
tx(16.5,18.8,'Cell proportions',6.6,fontweight='bold')
tx(16.5,22,'223 profiles × 89 states',6.0)
arrow(28,33,12.2)
# NMF: X(sample,state) ~= W(sample,module) H(module,state).
tx(47,7.3,'NMF',6.8,fontweight='bold')
for r in range(4):
 for c,col in enumerate([TEAL,SLATE,VIOLET]):rect(36.5+c*1.8,9.5+r*1.7,1.5,1.4,col)
tx(44,12.7,'×',8.5)
for r,col in enumerate([TEAL,SLATE,VIOLET]):
 for c in range(6):rect(47+c*1.7,10+r*2,1.4,1.6,col)
tx(47,19,'Sample scores ×',6.3)
tx(47,22,'state weights',6.3)
arrow(60,65,12.2)
# Three schematic groups illustrate shared state membership, not actual edges.
# The numerical result is the independently recorded 14 modules.
for cx,cy,col,coords in [(70.3,10.6,TEAL,[(-2,-1),(1.7,-1.6),(0,1.6)]),
                        (81,10.8,VIOLET,[(-2,-1.4),(1.7,-.7),(.1,1.7)]),
                        (76,15.2,SLATE,[(-2.2,-.5),(1.7,-.8),(.2,1.5)])]:
 for k in range(3):
  p0=coords[k];p1=coords[(k+1)%3];line(cx+p0[0],cy+p0[1],cx+p1[0],cy+p1[1],col,.65)
 for dx,dy in coords:dot(cx+dx,cy+dy,.72,col)
tx(76,19,'14 modules',6.8,fontweight='bold')
tx(76,22,'States varying together',6.0)
# Spatial section motifs: physical sections and array spots are schematic.
# No microscopy image, coordinates, counts or inferred abundance are fabricated.
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch
for dx,dy in [(1.7,-1),(0.8,-.45),(0,0)]:
 roundbox(97+dx,7.8+dy,14.5,9.1,'white','#B5C4D3',.5,.7)
verts=[(98.5,11.4),(99,8.2),(103.5,9.2),(106,8.4),(110,9.1),(109,12.5),(111,14.4),(108,16),(104,14.8),(101,16),(98.5,13.6),(98.5,11.4)]
codes=[MPath.MOVETO]+[MPath.CURVE3]*10+[MPath.CLOSEPOLY]
a.add_patch(PathPatch(MPath(verts,codes),fc='#E0E8EF',ec='#9BAFBE',lw=.5))
for x,y in [(100.2,10.8),(102.3,10.6),(104.5,10.3),(106.5,10.4),(108.4,10.4),(100.6,12.6),(102.6,12.5),(104.7,12.2),(106.7,12.3),(108.6,12.4),(102.6,14.3),(104.6,14.1),(106.7,14.3),(108.5,14.4)]:dot(x,y,.41,SLATE,ec='white',lw=.2)
tx(105,19,'49 sections',6.6,fontweight='bold')
tx(105,22,'29 candidate genes',6.0)
arrow(114.5,119,12.2)
# Display the scientific question, leaving the model equation in Methods.
# A nondirectional link denotes association, not cell-intrinsic expression
# or a causal relationship. Abundance is inferred from spatial mapping.
dot(128.5,10.2,3.2,'#E4EBF3',ec='none',lw=0,z=1)
dot(148.8,10.2,3.2,'#E2F0ED',ec='none',lw=0,z=1)
# RNA and cells are symbolic inputs, not invented expression measurements.
for dx in [0,1.5,3.0]:
 a.plot([126.5+dx,127.2+dx,126.5+dx,127.2+dx],
        [8.6,9.7,10.8,11.9],color=SLATE,lw=.85,solid_capstyle='round')
for x,y in [(147.4,9.1),(149.7,9.1),(146.7,11.1),(149,11.2),(151,10.9)]:
 dot(x,y,.79,'#75B5A9',ec=TEAL,lw=.4)
 dot(x+.1,y,.25,TEAL,ec='none',lw=0)
tx(138.6,10.2,'vs',7.0,color=GRAY)
tx(128.5,16.0,'Gene\nexpression',6.5,linespacing=1.03)
tx(148.8,16.0,'Estimated\ncell abundance',6.5,linespacing=1.03)
tx(138.5,22,'Compare within sections',6.6,fontweight='bold')
arrow(158.8,162.6,12.2)
# Unit effects converge into a pooled result. Conceptual aggregation glyph.
for yy in [8.5,11.4,14.3]:
 roundbox(164.8,yy-.85,3.2,1.7,'#CAD8E6','none',0,.4)
 line(168.5,yy,172,11.4,SLATE,.6)
a.add_patch(FancyArrowPatch((172,11.4),(174,11.4),arrowstyle='-|>',mutation_scale=6,lw=.7,color=SLATE))
a.add_patch(plt.Polygon([(174.2,11.4),(176.5,9.9),(178.8,11.4),(176.5,12.9)],fc=SLATE,ec='white',lw=.3))
tx(171.4,19,'Combine results',6.6,fontweight='bold')
tx(171.4,22,'35 analysis units',6.0)

f.canvas.draw(); renderer=f.canvas.get_renderer(); ex=[]
for t in a.texts:
 bb=t.get_window_extent(renderer);ex.append({'text':t.get_text(),'font_pt':t.get_fontsize(),'bbox_px':list(bb.extents)})
for ext in ['pdf','svg','png']:f.savefig(O/'panels'/f'Figure4a.{ext}',dpi=450)
plt.close(f)
# Copy protected numerical panels unchanged.
hashes={}
for letter in 'bcd':
 for ext in ['pdf','svg','png']:
  src=OLD/f'Figure4{letter}.{ext}';dst=O/'panels'/src.name;shutil.copy2(src,dst)
  hashes[str(src)]=hashlib.sha256(src.read_bytes()).hexdigest();assert dst.read_bytes()==src.read_bytes()
# Composition retains the exact b–d positions and dimensions.
placements={'a':['Figure4a',0,0,183,25],'b':['Figure4b',0,27,83,59],'c':['Figure4c',0,90,83,80],'d':['Figure4d',87,27,96,143]}
ns='http://www.w3.org/2000/svg';root=E.Element('{'+ns+'}svg',nsmap={None:ns},width='183mm',height='170mm',viewBox=f'0 0 {183*MM} {170*MM}')
E.SubElement(root,'{'+ns+'}rect',width='100%',height='100%',fill='white')
doc=fitz.open();pg=doc.new_page(width=183*MM,height=170*MM)
for l,(asset,x,y,w,h) in placements.items():
 src=fitz.open(O/'panels'/f'{asset}.pdf');pg.show_pdf_page(fitz.Rect(x*MM,y*MM,(x+w)*MM,(y+h)*MM),src,0)
 sr=E.parse(str(O/'panels'/f'{asset}.svg')).getroot();vb=list(map(float,sr.attrib['viewBox'].split()));sx=w*MM/vb[2];sy=h*MM/vb[3];assert abs(sx-sy)<1e-5
 group=E.SubElement(root,'{'+ns+'}g',id=f'panel_{l}',transform=f'translate({x*MM} {y*MM}) scale({sx})')
 ids={el.get('id'):l+'_'+el.get('id') for el in sr.iter() if el.get('id')}
 for el in sr.iter():
  if el.get('id'):el.set('id',ids[el.get('id')])
  for k,v in list(el.attrib.items()):
   for old,new in ids.items():
    v=v.replace(f'url(#{old})',f'url(#{new})')
    if v==f'#{old}':v='#'+new
   el.set(k,v)
 for child in sr:group.append(copy.deepcopy(child))
 label=E.SubElement(root,'{'+ns+'}text',x=str((x+.4)*MM),y=str((y+2.4)*MM),style='font-family:Arial;font-weight:bold;font-size:8px;fill:#222831');label.text=l
 pg.insert_text(((x+.4)*MM,(y+2.4)*MM),l,fontsize=8,fontname='hebo',color=(.13,.16,.19))
E.ElementTree(root).write(str(O/'review/Figure4.svg'),encoding='utf-8',xml_declaration=True)
doc.save(O/'review/Figure4.pdf',garbage=4,deflate=True)
pg.get_pixmap(matrix=fitz.Matrix(180/72,180/72),alpha=False).save(O/'review/Figure4.png')
(O/'review/placements.json').write_text(json.dumps(placements,indent=2))
(O/'review/source_hashes.json').write_text(json.dumps(hashes,ensure_ascii=False,indent=2))
(O/'review/text_extents.json').write_text(json.dumps(ex,indent=2))
print('Generated panel a and assembly. All matrix, network and arrow motifs are schematic.')
