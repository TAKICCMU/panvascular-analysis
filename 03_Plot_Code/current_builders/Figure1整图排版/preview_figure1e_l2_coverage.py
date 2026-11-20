"""A review-only alternative for Figure 1e; no scientific analysis is run.

One dot per named L2 state; x is the stored integer study count exactly.
Vertical offsets only separate duplicate x values within the parent category.
The adopted e panel and BioRender canvas are not touched.
"""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

O=Path(__file__).parent;P=O.parents[2];R=O/'review';R.mkdir(exist_ok=True)
sources=[P/'Source_Data_Figures/ExtendedData3_4_5/ED345_all_named_state_index.tsv',
 P/'Source_Data_Figures/Figure1/L2_cross_study_recurrence.tsv',O/'L1_display_palette.tsv']
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
hashes={str(p):sha(p) for p in sources}
adopted=[O/'panels/Figure1e.pdf',O/'panels/Figure1e.svg',O/'panels/Figure1e.png',O/'Figure1.pdf',O/'Figure1.svg']
adopted_hash={str(p):sha(p) for p in adopted}
named=pd.read_csv(sources[0],sep='\t');coverage=pd.read_csv(sources[1],sep='\t')
coverage[['level1_annotation','level2_annotation']]=coverage.state.str.split(' | ',regex=False,expand=True)
for frame in [named,coverage]:
 for c in ['level1_annotation','level2_annotation']:frame[c]=frame[c].str.strip()
z=named.merge(coverage,on=['level1_annotation','level2_annotation'],how='left',validate='one_to_one',suffixes=('_named',''))
assert len(z)==116 and z.studies.notna().all() and z.n_cells_named.eq(z.n_cells).all()
assert z.studies.between(0,26).all() and z.studies.eq(z.studies.astype(int)).all()
palette=pd.read_csv(sources[2],sep='\t').sort_values('order')
order=[s for s in palette.L1 if s in set(z.level1_annotation)]
assert len(order)==14
colors=palette.set_index('L1').hex.to_dict()

# Explicit millimetres retain print-size legibility. No panel letter is embedded.
mm=72/25.4;width=122.;dot_diameter=.85;dot_step=1.02
row_heights={}
for parent in order:
 max_ties=int(z[z.level1_annotation.eq(parent)].groupby('studies').size().max())
 row_heights[parent]=max(3.25,(max_ties-1)*dot_step+dot_diameter+1.5)
top=7.5;data_height=sum(row_heights.values());bottom=10.;height=top+data_height+bottom
mpl.rcParams.update({'font.family':'Arial','font.size':6.4,'text.color':'#222831',
 'axes.labelcolor':'#222831','xtick.color':'#48515A','svg.fonttype':'none','pdf.fonttype':42,
 'axes.linewidth':.5,'savefig.facecolor':'white'})
f=plt.figure(figsize=(width/25.4,height/25.4))
left=27.8;plot_width=79.0
ax=f.add_axes([left/width,bottom/height,plot_width/width,data_height/height])
ax.set_xlim(0,27);ax.set_ylim(data_height,0)
ax.set_xticks([0,5,10,15,20,26]);ax.set_yticks([])
ax.tick_params(axis='x',length=2.2,width=.45,pad=2.1,labelsize=6.2)
ax.set_xlabel('Number of studies',fontsize=6.5,labelpad=3.5)
ax.spines[['top','left','right']].set_visible(False)
ax.spines['bottom'].set_color('#87939D')
for x in [5,10,15,20,26]:ax.axvline(x,color='#E5E9EC',lw=.45,zorder=0)

overlay=f.add_axes([0,0,1,1],frameon=False)
overlay.set(xlim=(0,width),ylim=(height,0));overlay.set_axis_off()
overlay.text(25.8,3.2,'Cell type',ha='right',va='center',fontsize=6.5,weight='bold')
overlay.scatter([left+2.0],[3.2],s=(dot_diameter*mm)**2,color='#667783',edgecolors='none')
overlay.text(left+4.0,3.2,'One dot = one named L2 state',ha='left',va='center',fontsize=6.3)
overlay.text(115.6,3.2,'L2 states',ha='center',va='center',fontsize=6.4,weight='bold')

records=[];offset=0
for i,parent in enumerate(order):
 rh=row_heights[parent];cy=offset+rh/2
 if i%2==0:
  ax.axhspan(offset,offset+rh,fc='#F6F8F9',ec='none',zorder=-1)
 yy=top+cy
 overlay.text(25.8,yy,parent,ha='right',va='center',fontsize=6.4)
 overlay.add_patch(Rectangle((26.65,yy-.65),.42,1.3,facecolor=colors[parent],edgecolor='none'))
 group=z[z.level1_annotation.eq(parent)].sort_values(['studies','level2_annotation'])
 overlay.text(115.6,yy,str(len(group)),ha='center',va='center',fontsize=6.6,weight='bold')
 for count,tied in group.groupby('studies',sort=True):
  dy=(np.arange(len(tied))-(len(tied)-1)/2)*dot_step
  py=cy+dy
  ax.scatter(np.full(len(tied),count),py,s=(dot_diameter*mm)**2,c=colors[parent],
   edgecolors='none',linewidths=0,zorder=3)
  for (_,r),draw_y in zip(tied.iterrows(),py):
   row=r.to_dict();row['plot_y_mm']=float(draw_y);row['plot_x_studies']=int(count)
   row['display_hex']=colors[parent];records.append(row)
 offset+=rh

# Validate text bounds, dot count, exact x values and non-overlap before saving.
f.canvas.draw();renderer=f.canvas.get_renderer();bad=[]
for a in f.axes:
 visible_text=list(a.texts)
 if a.axison:visible_text += [a.xaxis.label]+list(a.get_xticklabels())
 for t in visible_text:
  if not t.get_text() or not t.get_visible():continue
  bb=t.get_window_extent(renderer).transformed(f.transFigure.inverted())
  if bb.x0<0 or bb.y0<0 or bb.x1>1 or bb.y1>1:bad.append(t.get_text())
assert not bad,bad
display=pd.DataFrame(records)
assert len(display)==116 and display.plot_x_studies.eq(display.studies).all()
xy=np.c_[display.plot_x_studies.to_numpy()*plot_width/27,display.plot_y_mm.to_numpy()]
dist=np.linalg.norm(xy[:,None,:]-xy[None,:,:],axis=2);np.fill_diagonal(dist,np.inf)
assert float(dist.min())>dot_diameter
prefix=R/'Figure1e_L2_coverage_prototype'
for ext in ['pdf','svg']:f.savefig(prefix.with_suffix('.'+ext))
f.savefig(prefix.with_suffix('.png'),dpi=300)
plt.close(f)
display.to_csv(R/'Figure1e_L2_coverage_plot_data.tsv',sep='\t',index=False)
check={'status':'PROTOTYPE_ONLY_NOT_ADOPTED','width_mm':width,'height_mm':round(height,3),
 'named_states':116,'parent_categories':14,'study_count_range':[int(z.studies.min()),int(z.studies.max())],
 'state_counts_match':True,'all_x_coordinates_equal_stored_counts':True,
 'vertical_stacking':'display-only separation of states with the same parent and study count',
 'minimum_dot_center_distance_mm':float(dist.min()),'dot_diameter_mm':dot_diameter,
 'font_pt':6.4,'panel_letter_embedded':False,'label_outside_page':bad,
 'interpretation':'Recorded study coverage; not an independent test of state reproducibility.',
 'source_sha256':hashes,'adopted_files_unchanged':all(sha(Path(p))==h for p,h in adopted_hash.items())}
assert check['adopted_files_unchanged'] and all(sha(Path(p))==h for p,h in hashes.items())
(R/'Figure1e_L2_coverage_prototype_check.json').write_text(json.dumps(check,indent=2))
print(json.dumps(check,indent=2))
