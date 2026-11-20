from common_main import *
from matplotlib.lines import Line2D
find=lambda r,i:next(e for e in r.iter() if e.get('id')==i)
base=O/'merge_records/before_20260921_adjustments/editable_svg/Figure3.svg'
r=ET.parse(base).getroot();old=copy.deepcopy(r)
# Rebuild only panel a at native font size, same points/statistics/axis limits.
aold=find(r,'axes_1');remove_ids(r,{'axes_1'})
f=plt.figure(figsize=(183/25.4,170/25.4));a=f.add_axes([15/183,1-36/170,63/183,27/170])
d=pd.read_csv(S/'Figure3/Fig3a_L1.tsv',sep='\t',float_precision='round_trip')
for _,v in d.iterrows():
 sig=v.FDR_BH<.05 and v['primary_r0.1__assoc_mcp_BH']<.05
 a.scatter(v.BETA_STD,v['primary_r0.1__assoc_mcz'],s=25+18*sig,facecolors=pal.get(v.group,'#777') if sig else 'white',edgecolors=pal.get(v.group,'#777'),lw=.8,zorder=3)
 if sig:
  lab=v.group.replace('Endothelial cell','EC').replace('Fibromyocyte','FMC')
  a.annotate(lab,(v.BETA_STD,v['primary_r0.1__assoc_mcz']),xytext=(-3,2) if lab=='Pericyte' else (3,-3),textcoords='offset points',ha='right' if lab=='Pericyte' else 'left',fontsize=6)
a.axhline(0,c='#aaa',lw=.55);a.axvline(0,c='#aaa',lw=.55)
a.set_xlim(-.0215,.040);a.set_ylim(-3.3,6.85)
a.set_xticks([-.02,0,.02,.04]);a.set_yticks([-2,0,2,4,6]);a.tick_params(labelsize=6)
a.set_xlabel('MAGMA standardised beta',labelpad=2,fontsize=6);a.set_ylabel('scDRS association Z',labelpad=2,fontsize=6)
a.legend(handles=[Line2D([],[],ls='',marker='o',ms=3.5,mfc='#888',mec='#888',label='Both FDR < 0.05')],loc='upper left',bbox_to_anchor=(-.03,-.33),frameon=False,fontsize=6,borderaxespad=0,handletextpad=.4)
embed_plot(r,f,'aligned_a_')
# Keep b objects/sizes; translate the complete table and its size legend together.
b=find(r,'figure3_b_content');b.set('transform','translate(65.8245708951751 0)') # old x +25 mm
lg=find(r,'legend_2');lg.set('transform','translate(129.1150364701148 -1.5)')
for t in r.iter(q('text')):
 if t.text=='b' and '700' in t.get('style',''):
  t.set('x',str(92*mm));t.set('transform',f'rotate(0, {92*mm}, {t.get("y")})')
# e: expand x positions only to the common right edge; preserve every y centre.
mat=find(r,'FigureYa_vector_matrix_Figure3e');ea=find(r,'axes_6')
rects=list(mat);left=min(float(x.get('x')) for x in rects);right=max(float(x.get('x'))+float(x.get('width')) for x in rects)
dx=-5.373607649156611;target=180*mm-dx;sx=(target-left)/(right-left)
for e in rects:
 e.set('x',str(left+(float(e.get('x'))-left)*sx));e.set('width',str(float(e.get('width'))*sx))
# Native rectangle gaps already provide the grid; remove old grid and label objects,
# then restore the five original headers at their exact new column centres.
remove_ids(ea,{'matplotlib.axis_11','matplotlib.axis_12'})
for parent in ea.iter():
 if parent.get('clip-path')=='url(#p2086020ded)':parent.attrib.pop('clip-path')
for col,lines in enumerate([['Primary','GWAS'],['Primary','PB'],['Leave-','CAD'],['Regional','ABF'],['Spatial']]):
 cell=next(e for e in mat if e.get('data-row')=='1' and e.get('data-col')==str(col+1))
 cx=float(cell.get('x'))+float(cell.get('width'))/2
 for j,label in enumerate(lines):text(ea,cx,255.760256-(len(lines)-1-j)*6.416438,label,size=6,anchor='middle')
# f uses the original saved 10 rows: two scenarios for each of five frozen contexts.
remove_ids(r,{'figure3f_three_evidence_columns'})
g=ET.SubElement(r,q('g'),{'id':'figure3f_three_evidence_columns'})
data=pd.read_csv(O/'source_data/Figure3f_same5_contexts_both_scenarios.tsv',sep='\t',float_precision='round_trip')
contexts=data[['gene','context']].drop_duplicates();assert len(data)==10 and len(contexts)==5
Y=np.linspace(156.5,166.1,5)*mm
cols=[(55,108,'PP.H4',0,1),(120,146,'PP.H3',0,.4),(158,180,'Shared SNPs',0,4500)]
for xa,xb,label,mn,mx in cols:
 text(g,(xa+xb)/2*mm,153.65*mm,label,size=6,anchor='middle')
 ticks=[0,.5,1] if label=='PP.H4' else [0,.2,.4] if label=='PP.H3' else [0,2000,4000]
 for v in ticks:
  x=(xa+(v-mn)/(mx-mn)*(xb-xa))*mm
  ET.SubElement(g,q('path'),{'d':f'M {x} {Y[0]-3} L {x} {Y[-1]+3}','stroke':'#e4e8eb','stroke-width':'.4','fill':'none'})
  text(g,x,169*mm,str(v) if label!='Shared SNPs' else ('0' if v==0 else str(v//1000)+'k'),size=5.7,anchor='middle')
for i,(_,v) in enumerate(contexts.iterrows()):
 cy=Y[i];lab=v.gene+' | '+v.context.replace('_',' ').replace('Cells Cultured fibroblasts','Cultured fibroblasts')
 text(g,52.5*mm,cy+2,lab,size=6,anchor='end')
 rows=data[(data.gene==v.gene)&(data.context==v.context)]
 for j,(xa,xb,label,mn,mx) in enumerate(cols):
  XX=lambda z:(xa+(float(z)-mn)/(mx-mn)*(xb-xa))*mm
  if j<2:
   pair=rows.set_index('scenario');ET.SubElement(g,q('path'),{'d':f'M {XX(pair.loc["primary",label])} {cy} L {XX(pair.loc["leave_one_CAD",label])} {cy}','stroke':'#bbc3c8','stroke-width':'.6','fill':'none'})
   for _,rr in rows.iterrows():
    xx=XX(rr[label]);col='#356b8c' if rr.scenario=='primary' else '#555555';fill=col if rr.classification=='strong_colocalisation' else 'white';attr={'fill':fill,'stroke':col,'stroke-width':'.7'}
    if rr.scenario=='primary':ET.SubElement(g,q('circle'),{'cx':str(xx),'cy':str(cy),'r':'1.8',**attr})
    else:ET.SubElement(g,q('rect'),{'x':str(xx-1.7),'y':str(cy-1.7),'width':'3.4','height':'3.4',**attr})
  else:
   assert rows.n_overlap.nunique()==1
   ET.SubElement(g,q('rect'),{'x':str(XX(0)),'y':str(cy-1.5),'width':str(XX(rows.n_overlap.iloc[0])-XX(0)),'height':'3','fill':'#9ba9b3'})
for y,typ,label,fill in [(156,'o','Primary','white'),(159,'s','Leave-CAD','white'),(162,'o','Strong','#555555'),(165,'o','Other','white')]:
 cx=4.8*mm;cy=y*mm;attr={'stroke':'#356b8c' if label=='Primary' else '#555555','fill':fill,'stroke-width':'.7'}
 if typ=='o':ET.SubElement(g,q('circle'),{'cx':str(cx),'cy':str(cy),'r':'1.7',**attr})
 else:ET.SubElement(g,q('rect'),{'x':str(cx-1.6),'y':str(cy-1.6),'width':'3.2','height':'3.2',**attr})
 text(g,cx+4,cy+1.9,label,size=5.7)
export(r,'Figure3')
# Evidence of unchanged values and geometry where shared.
assert ET.tostring(find(r,'figure3d_compact_forest'))==ET.tostring(find(old,'figure3d_compact_forest'))
for id in ['c_magma','c_pseudobulk']:assert ET.tostring(find(r,id))==ET.tostring(find(old,id))
assert [(e.get('y'),e.get('height'),e.get('data-source-value')) for e in mat]==[(e.get('y'),e.get('height'),e.get('data-source-value')) for e in find(old,'FigureYa_vector_matrix_Figure3e')]
(O/'notes/Figure3_alignment_0921.json').write_text(json.dumps({'scientific_analysis_rerun':False,'a_same15_source_rows':True,'b_native_objects_translation_only':True,'c_source_objects_identical':True,'d_source_objects_identical':True,'de_29_row_y_geometry_unchanged':True,'f_saved10_rows_same5_contexts':True,'f_fill_original_classification':True,'f_legend_left_vertical':True,'outer_e_f_right_edge_mm':180,'width_mm':183,'height_mm':170},indent=2))
print('Figure3 exported')
