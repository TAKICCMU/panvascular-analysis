"""Bounded visual updates. Read frozen 21 September backup; never change science."""
from f1_common import *
BACK=O/'merge_records/before_20260921_adjustments'
# Rebuild the two raster point layers from the exact frozen coordinate table.
script=(BACK/'code/figure1.py').read_text()
script=script.replace("'#CB7689'","'#A3966A'")
script=script.replace("ax.scatter(D.UMAP1,D.UMAP2,c=D[col].map(p)","shown=D.loc[D.level1_annotation.ne('Unknown')];ax.scatter(shown.UMAP1,shown.UMAP2,c=shown[col].map(p)")
script=script.replace("xy=mainax.transData.transform", "anchors=anchors.loc[anchors.L1.ne('Unknown')].reset_index(drop=True)\nxy=mainax.transData.transform")
script=script.replace('range(16)','range(len(anchors))')
# Preserve exactly the existing 15 labels' peripheral positions; removing Unknown is display-only.
script=script.replace("pos=np.array([s[:2] for s in slots]);ri,assign=linear_sum_assignment(cdist(xy,pos));", "pos=np.array([s[:2] for s in slots]);previous=pd.read_csv(O/'source_data/Figure1b_all16_direct_labels.tsv',sep='\\t').set_index('L1');assign=np.array([np.argmin(np.linalg.norm(pos-previous.loc[k,['new_label_x_mm','new_label_y_mm']].to_numpy(dtype=float),axis=1)) for k in anchors.L1]);")
script=script.replace('for step in range(100):','for step in range(0):')
# Fix loop over full 16-class summaries (unchanged).
script=script.replace('ax.set_yticks(range(len(anchors)),','ax.set_yticks(range(16),').replace('np.arange(len(anchors))','np.arange(16)')
# Avoid replacing full evidence/provenance tables while regenerating only b.
script=script.replace("anchors.to_csv(O/'source_data/Figure1b_all16_direct_labels.tsv'", "anchors.to_csv(O/'source_data/Figure1b_display15_direct_labels_20260921.tsv'")
script=script.replace("export(r,'Figure1')", "rebuilt=r")
script=script[:script.index("(O/'notes/Figure1_checks.json')")]
ns={};exec(compile(script,str(BACK/'code/figure1.py'),'exec'),ns)
root=ET.parse(BACK/'editable_svg/Figure1.svg').getroot()
# Swap only the two maps and their label overlay. Everything else is original SVG.
for ident in ['axes_2','axes_3','axes_4']:
 replacement=copy.deepcopy(next(e for e in ns['rebuilt'].iter() if e.get('id')==ident))
 # namespacing definitions avoids marker-id collision with immutable original panels.
 ids={e.get('id'):'b0921_'+e.get('id') for e in replacement.iter() if e.get('id') and e.get('id')!=ident}
 for el in replacement.iter():
  for k,v in list(el.attrib.items()):
   if k=='id' and v in ids:el.set(k,ids[v])
   elif k!='id':
    for a,b in ids.items():v=v.replace('url(#'+a+')','url(#'+b+')');v=('#'+b) if v=='#'+a else v
    el.set(k,v)
 for parent in root.iter():
  for j,child in enumerate(list(parent)):
   if child.get('id')==ident:parent.remove(child);parent.insert(j,replacement)
# align the three headings with actual frozen image rectangles, not legacy axes canvas.
for ai,ti in [('f1e_axes_8','f1e_text_186'),('f1e_axes_9','f1e_text_187'),('f1e_axes_10','f1e_text_188')]:
 a=next(e for e in root.iter() if e.get('id')==ai);im=a.find(q('image'));cx=float(im.get('x'))+float(im.get('width'))/2
 t=next(e for e in root.iter() if e.get('id')==ti).find(q('text'));t.set('x',str(cx));t.set('style',t.get('style').replace('text-anchor: start','text-anchor: middle'));t.attrib.pop('transform',None)
export(root,'Figure1')
D=ns['D'];report={'Figure1b_hidden_Unknown':int(D.level1_annotation.eq('Unknown').sum()),'Figure1b_shown':int(D.level1_annotation.ne('Unknown').sum()),'full_atlas_N':len(D),'hidden_in_both_panels':True,'annotation_and_all_other_panels_unchanged':True,'tissue_color_changed':{'Aortic valve':{'old':'#CB7689','new':'#A3966A'}},'heading_aligned_to_frozen_image_rectangles':True}
# Figure2: edit existing vector objects only.
r=ET.parse(BACK/'editable_svg/Figure2.svg').getroot();axes=next(e for e in r.iter() if e.get('id')=='axes_7');lg=next(e for e in r.iter() if e.get('id')=='legend_1')
# Horizontal legend in the dedicated gap directly above Manhattan plotting region.
for i,(lineid,textid,newx) in enumerate([('line2d_95','text_147',40),('line2d_96','text_148',180),('line2d_97','text_149',355)]):
 line=next(e for e in lg.iter() if e.get('id')==lineid);t=next(e for e in lg.iter() if e.get('id')==textid)
 oldy=[239.760632,248.145007,256.529382][i]
 for el in (line,t):el.set('transform',f'translate({newx-404.009033} {222.5-oldy})')
# the parent axes translation puts legend at y235.722 =83.16 mm, above the plot.
venn=next(e for e in r.iter() if e.get('id')=='axes_11');oldtexts={''.join(t.itertext()) for t in venn.findall('.//'+q('text'))}
# Three equal circles, schematic rather than area proportional. Seven exclusive counts frozen.
for child in list(venn):venn.remove(child)
venn.set('transform','')
centres=[(377,409),(423,409),(400,441)];rad=32
for (cx,cy),color in zip(centres,['#354f78','#278d87','#c68c35']):
 ET.SubElement(venn,q('circle'),{'cx':str(cx),'cy':str(cy),'r':str(rad),'fill':color,'fill-opacity':'.14','stroke':color,'stroke-width':'1.2'})
positions={'Positional':(371,391),'eQTL':(429,391),'MAGMA (Bonferroni)':(400,455),'65':(363,413),'1113':(437,413),'31':(400,446),'207':(400,403),'10':(381,433),'8':(419,433),'224':(400,424)}
assert set(positions)==oldtexts
for value,(x,y) in positions.items():
 if value=='MAGMA (Bonferroni)':
  text(venn,x,y,'MAGMA',size=6,anchor='middle');text(venn,x,y+8,'(Bonferroni)',size=6,anchor='middle')
 else:text(venn,x,y,value,size=6 if value in ['Positional','eQTL'] else 6.5,anchor='middle')
export(r,'Figure2')
report['Figure2f']={'circle_radius_pt':rad,'outline_pt':1.2,'unchanged_region_counts':[65,1113,31,207,10,8,224],'area_proportional':False};report['Figure2d_legend']='horizontal row above unchanged Manhattan data'
(O/'notes/Figure12_visual_20260921.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
