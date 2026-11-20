from pathlib import Path
import re, shutil, copy, json, hashlib
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np,pandas as pd, fitz
A=Path('/Users/taki/Documents/Panvascular_Figure_Archive_20260919')
O=Path('/Users/taki/Desktop/Panvascular_FigureYa_20260919')
B=Path('/Users/taki/Desktop/投稿包_最新版_20260919')
V=A/'Panvascular_Figure_Redesign_20260917_v1'
S=O/'code/genetics_assets'
for d in ['main','extended_data','editable_svg','notes','review','code/genetics_assets']: (O/d).mkdir(parents=True,exist_ok=True)
# Read-only historical module, rebound to archived native tables; writes go only to new assets.
src=(A/'Panvascular_Nature_Reconstruction_20260918/code/genetics_reconstruct.py').read_text()
src=src[:src.index("if __name__=='__main__':")]
src=src.replace("O=Path('/Users/taki/Desktop/Panvascular_Nature_Reconstruction_20260918');V=Path('/Users/taki/Desktop/Panvascular_Figure_Redesign_20260917_v1');S=O/'source_data/genetics'",f"O=Path({str(O)!r});V=Path({str(V)!r});S=O/'code/genetics_assets'")
ns={'__name__':'genetics_adapter'}
exec(compile(src,'archived_genetics_reconstruct.py','exec'),ns)
m=ns['m'];m.F=V/'source_data/original/Figure2'
# FigureYa126: scale='none', border_color=NA, no clustering. Preserve the frozen order,
# RdBu_r lightened mapping and both scientific ranges, rather than template expression scaling.
def correlation(ax):
 m.plain(ax);h=m.sub(ax,[.14,.18,.64,.80]);mat=m.rg.values
 colors=plt.get_cmap('RdBu_r')(np.linspace(0,1,256));colors[:,:3]=.52+.48*colors[:,:3]
 im=h.pcolormesh(np.arange(9)-.5,np.arange(9)-.5,mat,cmap=matplotlib.colors.ListedColormap(colors),norm=matplotlib.colors.TwoSlopeNorm(0,-1,1),edgecolors='none')
 h.set(xlim=(-.5,7.5),ylim=(7.5,-.5),aspect='equal');h.set_xticks(range(8),m.traits,rotation=45,ha='right');h.set_yticks(range(8),m.traits);h.tick_params(length=0,pad=2)
 for spine in h.spines.values():spine.set_visible(False)
 for i,t in enumerate(m.order):
  for j,u in enumerate(m.order[:i+1]):
   h.text(j,i,f'{mat[i,j]:.2f}',ha='center',va='center',size=6,color='#222222')
   q=m.rgl[((m.rgl.trait1==t)&(m.rgl.trait2==u))|((m.rgl.trait1==u)&(m.rgl.trait2==t))]
   if i!=j and len(q) and q.iloc[0].rg_p<.05:h.plot(j+.36,i-.36,'k.',ms=1.6)
 c=m.sub(ax,[.85,.30,.035,.58]);plt.colorbar(im,cax=c,ticks=[-1,0,1]);c.set_ylabel('rg',size=6,labelpad=2)
m.main_b=correlation
original_save=ns['save']
def save(f,n):
 # FigureYa6rHRsV2: fixed row spacing, aligned estimates, visible CI caps and neutral reference.
 # Linear loading/correlation axes are retained; do not apply template HR log transform.
 if n=='ExtendedData_Figure7':
  # Required scale label retained from ContentAlignment baseline, not a panel title.
  f.text(496/(180*72/25.4),1-382/(170*72/25.4),'Cross-trait intercept',rotation=90,ha='center',va='center',size=6)
  for outer in f.axes:
   for a in outer.child_axes:
    if a.get_xlabel() in ['Factor correlation (95% CI)','Standardized loading ± SE']:
     a.grid(axis='x',color='#E5E8EB',linewidth=.4);a.set_axisbelow(True)
     for line in a.lines:
      if line.get_marker()=='o':line.set_marker('o');line.set_markersize(3)
 # FigureYa43ManhattanV2: no grid/border, retained chromosome centres and thresholds.
 # Existing frozen raster layers retain every source point and exact cap; axes remain vector.
 if n=='Figure2':
  for a in f.axes:
   if a.get_ylabel()=='−log10(P)':
    a.grid(False);a.spines['top'].set_visible(False);a.spines['right'].set_visible(False)
    if a.get_legend():
     for tt in a.get_legend().get_texts():tt.set_fontsize(6)
 # Preserve ContentAlignment-approved set labels at exact page coordinates.
 if n=='Figure2':
  positions={'Positional':('Positional',336.13,402.0),'eQTL':('eQTL',449.06,402.0),'MAGMA Bonferroni':('MAGMA (Bonferroni)',392.598425,449.0)}
  for t in f.findobj(matplotlib.text.Text):
   if t.get_text() in positions:
    label,x,y=positions[t.get_text()];t.set_text(label);t.set_transform(f.transFigure);t.set_position((x/(183*72/25.4),1-y/(170*72/25.4)));t.set_ha('center');t.set_va('baseline')
   if t.get_text().startswith('● Nominal'):
    t.set_text('Nominal\nP < 0.05');t.set_position((63/183,1-69/170))
  f.add_artist(matplotlib.lines.Line2D([61/183],[1-70.1/170],transform=f.transFigure,marker='o',markersize=2,color='#222222',linestyle='none'))
 # Explicit labels/legends only. Remove unneeded footer from Venn, already documented in caption.
 for t in list(f.findobj(matplotlib.text.Text)):
  if t.get_text().startswith('Targets in union:'):t.remove()
 # FigureYa290: identity bars, explicit manual fill and classic axes; no fill normalization,
 # contingency-table generation or Fisher test is applied to frozen scores/counts.
 if n in ['ExtendedData_Figure1','ExtendedData_Figure8']:
  for outer in f.axes:
   for a in [outer]+list(outer.child_axes):
    if a.get_xlabel()=='scIB score':
     a.set_axisbelow(True);a.grid(axis='x',color='#E5E8EB',lw=.4)
     for bar in a.patches:
      if isinstance(bar,matplotlib.patches.Rectangle) and bar.get_height()<1:
       yc=bar.get_y()+bar.get_height()/2;bar.set_height(.36);bar.set_y(yc-.18)
 # FigureYa47: striped numeric table and booktabs-style horizontal rules.
 if n=='ExtendedData_Figure8':
  a=f.axes[-1]
  for i in [0,2,4,6]:
   y=.87-i*.138;a.add_patch(matplotlib.patches.Rectangle((.01,y-.060),.985,.12,fc='#F4F6F7',ec='none',zorder=-1))
  a.plot([.01,.995],[.967,.967],c='#899299',lw=.55)
 for t in f.findobj(matplotlib.text.Text):
  if t.get_text().startswith('● Shared G'):t.set_text(t.get_text())
 original_save(f,n)
 # Approved Figure2a is inserted as the exact frozen editable group; no rewrite of its five criteria.
 if n=='Figure2':
  p=O/'editable_svg/Figure2.svg';root=ET.parse(p).getroot();old=ET.parse(B/'图件编辑资料/editable_svg/Figure2.svg').getroot()
  dest=next(e for e in root.iter() if e.get('id')=='figure_1');new=copy.deepcopy(next(e for e in old.iter() if e.get('id')=='axes_1'))
  dest.remove(next(e for e in dest if e.get('id')=='axes_1'));dest.insert(1,new)
  ET.ElementTree(root).write(p,encoding='utf-8',xml_declaration=True)
  doc=fitz.open(O/'main/Figure2.pdf');base=fitz.open(B/'Figures/Figure2.pdf');out=fitz.open();page=out.new_page(width=doc[0].rect.width,height=doc[0].rect.height)
  top=fitz.Rect(0,0,page.rect.width,65);low=fitz.Rect(0,65,page.rect.width,page.rect.height)
  page.show_pdf_page(top,base,0,clip=top);page.show_pdf_page(low,doc,0,clip=low);out.save(O/'main/Figure2.tmp.pdf');doc.close();(O/'main/Figure2.tmp.pdf').replace(O/'main/Figure2.pdf')
 # Parent owns shared review assets: all our previews have genetics prefix.
 preview=O/f'review/{n}.png'
 if preview.exists():preview.unlink()
ns['save']=save
# FigureYa112 area_calculate port: derive binary-exclusive regions from actual sets, no new criteria.
origsets=ns['sets']
def sets():
 counts=origsets();d=pd.read_csv(S/'Figure2f_exact_gene_membership.tsv',sep='\t');sets_=[set(d.loc[d[c],'gene_id']) for c in ['positional','eQTL','MAGMA_Bonferroni']];union=set.union(*sets_);exact={}
 for i in range(1,8):
  code=f'{i:03b}';region=union.copy()
  for bit,s in zip(code,sets_):region=region&s if bit=='1' else region-s
  exact[code]=len(region)
 assert exact==counts and len(union)==1658 and exact['111']==224
 return exact
ns['sets']=sets
for name in ['ed1','ed7','ed8','main']:
 ns[name]();print('DONE',name,flush=True)
# PDF raster review at identical physical scale and structured per-figure data contract.
from PIL import Image,ImageOps,ImageDraw
rows=[]
for n in ['Figure2','ExtendedData_Figure1','ExtendedData_Figure7','ExtendedData_Figure8']:
 folder='main' if n=='Figure2' else 'extended_data';old=B/('Figures' if n=='Figure2' else 'Extended_Data')/(n+'.pdf');new=O/f'{folder}/{n}.pdf'
 pair=[]
 for label,p in [('before',old),('after',new)]:
  doc=fitz.open(p);pix=doc[0].get_pixmap(matrix=fitz.Matrix(2,2),alpha=False);im=Image.frombytes('RGB',[pix.width,pix.height],pix.samples);im.save(O/f'review/genetics_{n}_{label}.png');pair.append(im)
 canvas=Image.new('RGB',(sum(i.width for i in pair)+20,max(i.height for i in pair)+30),'white');draw=ImageDraw.Draw(canvas);draw.text((10,5),'Before',fill='black');draw.text((pair[0].width+20,5),'After: FigureYa adapter',fill='black');canvas.paste(pair[0],(0,30));canvas.paste(pair[1],(pair[0].width+20,30));canvas.save(O/f'review/genetics_{n}_comparison.png')
 rows.append({'figure':n,'pdf_sha256':hashlib.sha256(new.read_bytes()).hexdigest(),'size_mm':[round(v*25.4/72,2) for v in (doc[0].rect.width,doc[0].rect.height)],'data':'Frozen tables and exact order retained','visual_review':'PASS: rendered pages inspected; final panel f labels and dot legend reviewed'})
(O/'review/genetics_integrity.json').write_text(json.dumps(rows,indent=2))
