"""Author-requested layout only: Table10 values unchanged, no models refitted."""
from pathlib import Path
import hashlib,json,shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent
BK=R/'backup_before_efg_row'
protected={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (R/'panels').glob('Figure2[a-f].svg')}
for rel in ['panels/Figure2g.svg','panels/Figure2g.pdf','panels/Figure2g.png','upload/Figure2g.svg','panel_manifest.json','panel_placement.json','assemble_review.py','review/Figure2_current_preview.svg','review/Figure2_current_preview.pdf','review/Figure2_current_preview.png','review/current_biorender_status.json','README.md']:
 p=R/rel;q=BK/rel
 if p.exists() and not q.exists():
  q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
source=R/'alternatives/source_data/option_B_final_Table10_rows.tsv'
source_hash=hashlib.sha256(source.read_bytes()).hexdigest()
d=pd.read_csv(source,sep='\t')
assert d['Lead variant'].tolist()==['rs7525267','rs75135846','rs78051210','rs333929','rs7236320']
W,H=72,48
INK='#222B35'; MUTED='#697480'; BLUE='#355A80'; TEAL='#278B88'
plt.rcParams.update({'font.family':'Arial','font.size':6.2,'axes.labelsize':6.2,'xtick.labelsize':6,
 'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'axes.edgecolor':'#68737D',
 'axes.linewidth':.55,'svg.fonttype':'none','pdf.fonttype':42,'savefig.facecolor':'white',
 'svg.hashsalt':'Figure2g_single_row_author_reference'})
f=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
def ax(x,y,w,h):return f.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def text(x,y,s,**kw):return f.text(x/W,1-y/H,s,va=kw.pop('va','top'),**kw)
text(23,1,'Shared GWAS',fontsize=6.8,fontweight='bold',ha='center')
text(62,1,'Variant PIP',fontsize=6.8,fontweight='bold',ha='center')
ap=ax(31,8,15,30);ai=ax(54,8,16,30)
y=np.arange(5);p=d['P'].to_numpy(dtype=float);pip=d['Primary posterior inclusion probability (PIP)'].to_numpy(dtype=float)
for a in [ap,ai]:
 a.set_ylim(4.6,-.6);a.set_yticks([])
 for s in ['top','right','left']:a.spines[s].set_visible(False)
 a.tick_params(axis='x',length=2.5,width=.5,pad=2)
 for j in [0,2,4]:a.axhspan(j-.5,j+.5,color='#F6F8FA',lw=0,zorder=0)
ap.set_xlim(0,10);ap.set_xticks([0,5,10]);ap.set_xlabel('−log10(P)',labelpad=3)
ap.axvline(-np.log10(5e-8),ls=(0,(2,2)),lw=.6,color=MUTED,zorder=2)
ap.barh(y,-np.log10(p),color=BLUE,height=.30,zorder=3)
ai.set_xlim(-.02,1.04);ai.set_xticks([0,.5,1],['0.0','0.5','1.0']);ai.set_xlabel('PIP',labelpad=3)
ai.hlines(y,0,pip,color='#AED5D3',lw=1.5,zorder=2)
ai.scatter(pip,y,s=18,c=TEAL,edgecolors='none',zorder=3)
genes=['C1orf21 / COLGALT2','LINC01105 / SOX11','EPB41L2','MAP2K6 region','RAB27B']
for j,(rs,g) in enumerate(zip(d['Lead variant'],genes)):
 yy=8+(j+.6)/5.2*30
 text(1,yy-1.4,rs,fontsize=6.3,fontweight='bold',va='center')
 text(1,yy+1.3,g,fontsize=5.8,color=MUTED,va='center')
 text(50,yy,f'{pip[j]:.3f}',fontsize=6,ha='center',va='center',color=TEAL)
f.canvas.draw()
renderer=f.canvas.get_renderer()
over=[]
for t in f.findobj(matplotlib.text.Text):
 if not t.get_visible() or not t.get_text():continue
 b=t.get_window_extent(renderer)
 if b.x0<-.5 or b.y0<-.5 or b.x1>f.bbox.width+.5 or b.y1>f.bbox.height+.5:over.append(t.get_text())
assert not over,over
for ext in ['svg','pdf','png']:f.savefig(R/'panels'/f'Figure2g.{ext}',dpi=300)
plt.close(f)
shutil.copy2(R/'panels/Figure2g.svg',R/'upload/Figure2g.svg')
m=json.loads((R/'panel_manifest.json').read_text())
for item in m:
 if item['panel']=='Figure2g':item.update(width_mm=W,height_mm=H,change='Author screenshot style; compact same-row e/f/g layout; all five original P and PIP values unchanged; no panel letter baked',svg_sha256=hashlib.sha256((R/'panels/Figure2g.svg').read_bytes()).hexdigest())
(R/'panel_manifest.json').write_text(json.dumps(m,ensure_ascii=False,indent=2)+'\n')
assert source_hash==hashlib.sha256(source.read_bytes()).hexdigest()
assert protected=={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (R/'panels').glob('Figure2[a-f].svg')}
(R/'review/efg_row_checks.json').write_text(json.dumps({'scope':'layout only','source_sha256':source_hash,'all_five_original_records':True,'a_to_f_assets_byte_identical':True,'g_dimensions_mm':[W,H],'g_panel_letter_baked':False,'text_outside_page':over,'min_text_pt':5.8,'screenshot_threshold_line_preserved':True,'threshold_explanation':'existing figure legend, not an extra bottom key'},indent=2)+'\n')
print('G rendered 72 x 48 mm; original five records preserved; a-f assets unchanged.')
