"""Review-only Figure 1e alternative, drawn from stored named-state summaries.

No expression, annotation, integration or inference is recomputed.  The left
plot retains one dot per L2 state, with the measured study count on the y axis.
Equal values spread only along the categorical x axis. The right plot shows
all stored cell counts, ranked for display. No adopted figure is replaced.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
import pymupdf

O=Path(__file__).parent
P=O.parents[2]
R=O/'review'
R.mkdir(exist_ok=True)
MM=72/25.4
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
sources=[P/'Source_Data_Figures/ExtendedData3_4_5/ED345_all_named_state_index.tsv',
         P/'Source_Data_Figures/Figure1/L2_cross_study_recurrence.tsv',
         O/'L1_display_palette.tsv']
before={str(p):sha(p) for p in sources}
protected=list((O/'panels').glob('Figure1*.*'))+[O/'Figure1.pdf',O/'Figure1.svg',O/'panel_placement.json']
protected_hash={str(p):sha(p) for p in protected if p.is_file()}
named=pd.read_csv(sources[0],sep='\t')
coverage=pd.read_csv(sources[1],sep='\t')
coverage[['level1_annotation','level2_annotation']]=coverage.state.str.split(' | ',regex=False,expand=True)
keys=['level1_annotation','level2_annotation']
for table in (named,coverage):
    for key in keys:table[key]=table[key].str.strip()
z=named.merge(coverage,on=keys,how='left',validate='one_to_one',suffixes=('_named',''))
assert len(z)==116 and z.studies.notna().all() and z.n_cells.eq(z.n_cells_named).all()
assert (z.studies==z.studies.astype(int)).all() and z.n_cells.gt(0).all()
palette=pd.read_csv(sources[2],sep='\t').sort_values('order')
order=[k for k in palette.L1 if k in set(z.level1_annotation)]
colors=palette.set_index('L1').hex.to_dict()
display={'Endothelial cell':'Endothelial','Fibroblast':'Fibroblast',
 'Smooth muscle cell':'Smooth muscle','Fibromyocyte':'Fibromyocyte',
 'Pericyte':'Pericyte','Macrophage':'Macrophage','Monocyte':'Monocyte',
 'Neutrophil':'Neutrophil','Dendritic cell':'Dendritic','T cell':'T cell',
 'Natural killer cell':'Natural killer','B cell':'B cell','Plasma cell':'Plasma cell','Mast cell':'Mast cell'}
assert len(order)==14 and len(set(display.values()))==14

mpl.rcParams.update({'font.family':'Arial','font.size':7,'text.color':'#222831',
 'axes.labelcolor':'#222831','axes.edgecolor':'#87939D','axes.linewidth':.5,
 'xtick.color':'#48515A','ytick.color':'#48515A','pdf.fonttype':42,
 'svg.fonttype':'none','savefig.facecolor':'white'})
W,H=183.,55.
f=plt.figure(figsize=(W/25.4,H/25.4))
def axes(box):
    x,y,w,h=box
    return f.add_axes([x/W,(H-y-h)/H,w/W,h/H])
left_box=[10.5,8.5,117.,25.5]
right_box=[144.,8.5,37.,25.5]
a=axes(left_box)
b=axes(right_box)
a.set_xlim(-.5,13.5)
a.set_ylim(5,27.5)
a.set_yticks([5,10,15,20,26])
a.set_ylabel('Studies per L2 state',fontsize=7,labelpad=3)
a.set_xticks(np.arange(14),[display[v] for v in order],rotation=55,ha='right',rotation_mode='anchor')
a.tick_params(axis='x',length=0,pad=3.5,labelsize=7)
a.tick_params(axis='y',length=2,pad=2,labelsize=6.7)
a.spines[['top','right']].set_visible(False)
a.spines['bottom'].set_visible(False)
a.grid(axis='y',color='#E5E9EC',lw=.45,zorder=0)
dot_diameter=1.05
dot_step=1.12
category_mm=left_box[2]/14
records=[]
for i,parent in enumerate(order):
    if i%2==0:a.axvspan(i-.5,i+.5,fc='#F7F8F9',ec='none',zorder=-1)
    g=z[z.level1_annotation.eq(parent)]
    for count,tied in g.groupby('studies',sort=True):
        tied=tied.sort_values('level2_annotation')
        xx=i+(np.arange(len(tied))-(len(tied)-1)/2)*dot_step/category_mm
        assert np.max(np.abs(xx-i))<.5
        a.scatter(xx,np.full(len(tied),count),s=(dot_diameter*MM)**2,
                  c=colors[parent],edgecolors='none',linewidths=0,zorder=3)
        for (_,row),x in zip(tied.iterrows(),xx):
            r=row.to_dict();r.update(display_label=display[parent],plot_x_category=float(x),
                plot_y_studies=int(count),display_hex=colors[parent]);records.append(r)
    a.plot([i-.34,i+.34],[5.15,5.15],color=colors[parent],lw=1.35,solid_capstyle='butt',clip_on=True)

ranked=z.sort_values(['n_cells','level1_annotation','level2_annotation'],ascending=[False,True,True]).copy()
ranked['rank_for_display']=np.arange(1,len(ranked)+1)
b.set_yscale('log')
b.set_xlim(-2,119)
b.set_ylim(100,100000)
b.set_xticks([1,60,116])
b.set_yticks([100,1000,10000,100000],['100','1k','10k','100k'])
b.yaxis.set_minor_locator(NullLocator())
b.set_ylabel('Cells per L2 state',fontsize=7,labelpad=2)
b.set_xlabel('L2 states ranked by size',fontsize=6.7,labelpad=3)
b.spines[['top','right']].set_visible(False)
b.tick_params(length=2,pad=2,labelsize=6.7)
b.grid(axis='y',color='#E5E9EC',lw=.45,zorder=0)
b.scatter(ranked.rank_for_display,ranked.n_cells,s=(.88*MM)**2,
    c=ranked.level1_annotation.map(colors),edgecolors='none',linewidths=0,zorder=3)

o=axes([0,0,W,H]);o.set(xlim=(0,W),ylim=(H,0));o.set_axis_off()
o.text(1.0,3.0,'L2 states',fontsize=6.7,va='center')
for i,parent in enumerate(order):
    x=left_box[0]+(i+.5)*category_mm
    o.text(x,3.0,str(int(z.level1_annotation.eq(parent).sum())),fontsize=7.4,
           fontweight='bold',ha='center',va='center')
o.scatter([74.],[51.6],s=(dot_diameter*MM)**2,c='#61717C',edgecolors='none')
o.text(76.,51.6,'One dot = one named L2 state',ha='left',va='center',fontsize=6.7)

f.canvas.draw();renderer=f.canvas.get_renderer();outside=[]
for axis in f.axes:
    texts=list(axis.texts)
    if axis.axison:texts += [axis.xaxis.label,axis.yaxis.label]+axis.get_xticklabels()+axis.get_yticklabels()
    for t in texts:
        if not t.get_visible() or not t.get_text():continue
        bb=t.get_window_extent(renderer).transformed(f.transFigure.inverted())
        if min(bb.x0,bb.y0)<0 or max(bb.x1,bb.y1)>1:outside.append(t.get_text())
assert not outside,outside
plotdata=pd.DataFrame(records).merge(ranked[keys+['rank_for_display']],on=keys,validate='one_to_one')
assert len(plotdata)==116 and plotdata.plot_y_studies.eq(plotdata.studies).all()
xy=np.c_[(plotdata.plot_x_category+.5).to_numpy()*category_mm,
         (plotdata.plot_y_studies-5).to_numpy()*left_box[3]/22.5]
dist=np.linalg.norm(xy[:,None,:]-xy[None,:,:],axis=2);np.fill_diagonal(dist,np.inf)
assert dist.min()>dot_diameter,dist.min()
prefix=R/'Figure1e_compact_coverage_and_size'
for ext in ['pdf','svg']:f.savefig(prefix.with_suffix('.'+ext))
f.savefig(prefix.with_suffix('.png'),dpi=300)
plt.close(f)
plotdata.to_csv(R/'Figure1e_compact_plot_data.tsv',sep='\t',index=False)
pd.DataFrame([{'full_parent_identity':k,'display_label':display[k],'hex':colors[k]} for k in order]).to_csv(
    R/'Figure1e_compact_display_labels.tsv',sep='\t',index=False)

# Same 183 x 70 mm review pages: exact old allocated area versus the new strip.
# This is a size comparison, not an official assembled Figure 1.
comp=pymupdf.open()
fontfile='/System/Library/Fonts/Supplemental/Arial.ttf'
for title,source,box in [
 ('Previous prototype in the current 122 x 35 mm allocation',R/'Figure1e_L2_coverage_prototype.pdf',[2,15,122,35]),
 ('Proposed full-width strip at native size: 183 x 55 mm',prefix.with_suffix('.pdf'),[0,11,183,55])]:
    page=comp.new_page(width=183*MM,height=70*MM)
    page.insert_font(fontname='ReviewArial',fontfile=fontfile)
    page.insert_text((2*MM,6*MM),title,fontsize=8,fontname='ReviewArial',color=(.13,.16,.19))
    x,y,w,h=box
    rect=pymupdf.Rect(x*MM,y*MM,(x+w)*MM,(y+h)*MM)
    with pymupdf.open(source) as src:page.show_pdf_page(rect,src,0,keep_proportion=True)
    if source.name.endswith('_prototype.pdf'):page.draw_rect(rect,color=(.7,.73,.76),width=.4)
comp.save(R/'Figure1e_actual_size_comparison.pdf',garbage=4,deflate=True)
comp.close()
assert all(sha(Path(k))==v for k,v in before.items())
assert all(sha(Path(k))==v for k,v in protected_hash.items())
audit={'status':'REVIEW_PROTOTYPE_ONLY_NOT_ADOPTED','native_size_mm':[W,H],
 'old_allocated_size_mm':[122,35],'previous_prototype_size_mm':[122,82.56],
 'previous_font_pt':6.4,'previous_font_when_fit_to_old_slot_pt':round(6.4*35/82.56,2),
 'native_label_font_pt':7,'native_coverage_dot_diameter_mm':dot_diameter,
 'previous_dot_diameter_mm':.85,'states_in_each_plot':116,'parent_categories':14,
 'minimum_coverage_dot_center_distance_mm':float(dist.min()),'text_outside_page':outside,
 'figure_panel_letters_embedded':False,
 'study_values_preserved':True,'cell_count_values_preserved':True,
 'layout_rule':'Exact study count on y; display-only horizontal separation within the categorical parent.',
 'size_plot_rule':'All 116 stored state cell counts; descending rank only; logarithmic y axis.',
 'interpretation':'Coverage and state size summaries, not an independent stability or annotation validation.',
 'integration_status':'Needs more space than the old e slot; no shrinking or changes to adopted panels performed.',
 'source_sha256':before,'adopted_files_unchanged':True}
(R/'Figure1e_compact_layout_check.json').write_text(json.dumps(audit,indent=2))
print(json.dumps({k:v for k,v in audit.items() if k!='source_sha256'},indent=2))
