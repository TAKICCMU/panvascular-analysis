from common_main import *
BASE=Path('/Users/taki/Desktop/Panvascular_Main_Figure_Refinement_20260920/reproducibility_dependencies/Panvascular_Figure3_Group_Spacing_20260920')
MC=Path('/Users/taki/Desktop/Panvascular_Main_Figure_Refinement_20260920/reproducibility_dependencies/Panvascular_Marker_Communication_Enrichment_20260920')
find=lambda r,i:next(e for e in r.iter() if e.get('id')==i)
# Figure4: replace exactly ten lollipop glyphs with zero-origin bars of identical length.
r=ET.parse(MC/'editable_svg/Figure4.svg').getroot();barcheck=[]
for ai in [6,7]:
 a=find(r,f'axes_{ai}')
 for g in a:
  circ=g.find(q('circle'));path=g.find(q('path'))
  if circ is None or path is None:continue
  vals=list(map(float,re.findall(r'-?\d+(?:\.\d+)?',path.get('d'))));assert len(vals)==4
  x0,y,x1,y1=vals; assert abs(y-float(circ.get('cy')))<1e-8 and y==y1
  for ch in list(g):g.remove(ch)
  ET.SubElement(g,q('rect'),{'x':str(x0),'y':str(y-3.1),'width':str(x1-x0),'height':'6.2','fill':circ.get('fill'),'clip-path':path.get('clip-path')})
  barcheck.append({'axis':ai,'x0':x0,'x1':x1,'y':y,'bar_length_unchanged':True})
assert len(barcheck)==10
export(r,'Figure4');pd.DataFrame(barcheck).to_csv(O/'notes/Figure4_bar_geometry.tsv',sep='\t',index=False)
# Figure3 accepted source: only named areas change.
r=ET.parse(BASE/'editable_svg/Figure3.svg').getroot();original=copy.deepcopy(r)
# Databases form one common, left-aligned column in each enrichment block.
groups=[]
for bid,name_x in [('c_magma',16),('c_pseudobulk',13)]:
 block=find(r,bid); rows=[t for t in block.iter(q('text')) if 'font: 6px' in t.get('style','') and 'text-anchor: end' in t.get('style','')]
 assert len(rows)==9
 # Native cropped axes include a 0.576pt origin offset, retained in all objects.
 origin=.576
 for t in rows:
  t.set('x',str(origin+name_x*mm));t.set('style',t.get('style').replace('text-anchor: end','text-anchor: start'));t.attrib.pop('transform',None)
 line_x=origin+(name_x-1.2)*mm
 paths=[p for p in block.iter(q('path')) if 'stroke: #aaaaaa' in p.get('style','')];assert len(paths)==3
 for p in paths:
  v=list(map(float,re.findall(r'-?\d+(?:\.\d+)?',p.get('d'))));p.set('d',f'M {line_x} {v[1]} L {line_x} {v[3]}')
 for t in block.iter(q('text')):
  if t.text in ['GO BP','KEGG','Reactome']:
   t.set('x',str(line_x-mm));t.attrib.pop('transform',None)
 groups.append({'block':bid,'group_line_x_pt':line_x,'all_three_lines_aligned':True,'name_start_gap_mm':1.2})
# Replace forest drawing at original exact row centres, expand leftwards to remove unused margin.
C=pd.read_csv(S/'Figure3/Fig3d_e_all29.tsv',sep='\t',float_precision='round_trip')
a=find(r,'axes_5'); labels=[t.text for t in find(a,'matplotlib.axis_10').iter(q('text'))];assert len(labels)==29
mat=find(r,'FigureYa_vector_matrix_Figure3e');rects=[e for e in mat if e.get('data-col')=='1']
ycent=np.array([float(e.get('y'))+float(e.get('height'))/2-19 for e in rects])
assert len(ycent)==29
oldedgepath=find(a,'patch_12').find(q('path'));v=list(map(float,re.findall(r'-?\d+(?:\.\d+)?',oldedgepath.get('d'))));
xhi=348.661417+37.83591751150942-31.87094183783131
xlo=29*mm;top=257.952756-19;bottom=416.125984-19
xmin=C.meta_CI_low.min()-.12;xmax=C.meta_CI_high.max()+.12
X=lambda val:xlo+(float(val)-xmin)/(xmax-xmin)*(xhi-xlo)
remove_ids(r,{'axes_5'})
g=ET.SubElement(r,q('g'),{'id':'figure3d_compact_forest'})
for i,(idx,row) in enumerate(C.iterrows()):
 cy=ycent[i];pitch=ycent[1]-ycent[0]
 if i%2==0:ET.SubElement(g,q('rect'),{'x':str(xlo),'y':str(cy-pitch/2),'width':str(xhi-xlo),'height':str(pitch),'fill':'#f5f6f7'})
for tick in [-2,-1,0,1,2]:
 x=X(tick);ET.SubElement(g,q('path'),{'d':f'M {x} {top} L {x} {bottom+3}','stroke':'#aeb5ba' if tick==0 else '#e3e7ea','stroke-width':'.7' if tick==0 else '.4','fill':'none'})
 text(g,x,bottom+11,str(tick).replace('-','−'),size=6,anchor='middle')
for i,(_,row) in enumerate(C.iterrows()):
 cy=ycent[i];col=pal.get(row.cell_state,'#777777')
 text(g,xlo-3.5,cy+2.0,labels[i],size=6,anchor='end')
 lo,hi,xx=X(row.meta_CI_low),X(row.meta_CI_high),X(row.meta_logFC)
 ET.SubElement(g,q('path'),{'d':f'M {lo} {cy} L {hi} {cy} M {lo} {cy-1.2} L {lo} {cy+1.2} M {hi} {cy-1.2} L {hi} {cy+1.2}','stroke':col,'stroke-width':'.85','fill':'none'})
 rr=2.35;ET.SubElement(g,q('path'),{'d':f'M {xx} {cy-rr} L {xx+rr} {cy} L {xx} {cy+rr} L {xx-rr} {cy} Z','fill':col})
ET.SubElement(g,q('path'),{'d':f'M {xlo} {top} L {xlo} {bottom} L {xhi} {bottom}','stroke':'#222222','stroke-width':'.65','fill':'none'})
text(g,(xlo+xhi)/2,414.9,'Meta log2 fold change (95% CI)',size=6,anchor='middle')
# Same evidence matrix, only one-row legend; no not-evaluated entries currently present.
vals=[int(e.get('data-source-value')) for e in mat];assert 0 not in vals
remove_ids(r,{'legend_3'})
lg=ET.SubElement(r,q('g'),{'id':'figure3e_three_symbol_key'})
for x,label,col in [(128,'Support','#278d87'),(143,'Inverse','#d18a32'),(157,'No support','#dce1e5')]:
 ET.SubElement(lg,q('rect'),{'x':str(x*mm),'y':str(143*mm),'width':'3.8','height':'3.8','fill':col});text(lg,x*mm+5.5,143*mm+3.9,label,size=6)
# Expand f into three compact quantitative columns, same five contexts, no new model.
cf=pd.read_csv(S/'Figure3/Fig3f_strong5_contexts.tsv',sep='\t',float_precision='round_trip')
cp=Path('/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Figure_Redesign_20260917_v1/source_data/extra/all29_complete_cis_coloc_results.tsv')
co=pd.read_csv(cp,sep='\t',float_precision='round_trip');frozen=cf[['gene','context']].merge(co,on=['gene','context'],validate='one_to_many',sort=False)
assert len(frozen)==10 and set(frozen.scenario)=={'primary','leave_one_CAD'}
for _,v in cf.iterrows():
 a0=frozen[(frozen.gene==v.gene)&(frozen.context==v.context)&(frozen.scenario=='primary')].iloc[0]
 for field in ['PP.H3','PP.H4','n_overlap']:
  assert np.isclose(a0[field],v[field],rtol=0,atol=1e-15)
  frozen.loc[(frozen.gene==v.gene)&(frozen.context==v.context)&(frozen.scenario=='primary'),field]=v[field]
frozen.to_csv(O/'source_data/Figure3f_same5_contexts_both_scenarios.tsv',sep='\t',index=False)
remove_ids(r,{'axes_7','legend_4','figure3_f_original_fill_key'})
g=ET.SubElement(r,q('g'),{'id':'figure3f_three_evidence_columns'})
Y=np.linspace(156.5,166.1,5)*mm
cols=[(65,104,'PP.H4',0,1),(115,143,'PP.H3',0,.4),(155,178,'Shared SNPs',0,4500)]
# Column names centred over their own data areas.
for xa,xb,label,mn,mx in cols:
 text(g,(xa+xb)/2*mm,153.65*mm,label,size=6,anchor='middle')
 ticks=[0,.5,1] if label=='PP.H4' else [0,.2,.4] if label=='PP.H3' else [0,2000,4000]
 for t in ticks:
  x=(xa+(t-mn)/(mx-mn)*(xb-xa))*mm
  ET.SubElement(g,q('path'),{'d':f'M {x} {Y[0]-3} L {x} {Y[-1]+3}','stroke':'#e4e8eb','stroke-width':'.4','fill':'none'})
  text(g,x,169*mm,('0' if t==0 else str(t)) if label!='Shared SNPs' else ('0' if t==0 else str(int(t/1000))+'k'),size=5.7,anchor='middle')
for i,(_,v) in enumerate(cf.iterrows()):
 cy=Y[i];lab=v.gene+' | '+v.context.replace('_',' ').replace('Cells Cultured fibroblasts','Cultured fibroblasts')
 text(g,62.5*mm,cy+2,lab,size=6,anchor='end')
 rows=frozen[(frozen.gene==v.gene)&(frozen.context==v.context)]
 for j,(xa,xb,label,mn,mx) in enumerate(cols):
  XX=lambda z:(xa+(float(z)-mn)/(mx-mn)*(xb-xa))*mm
  if j<2:
   field=label;pair=rows.set_index('scenario')
   ET.SubElement(g,q('path'),{'d':f'M {XX(pair.loc["primary",field])} {cy} L {XX(pair.loc["leave_one_CAD",field])} {cy}','stroke':'#bbc3c8','stroke-width':'.6','fill':'none'})
   for _,rr in rows.iterrows():
    xx=XX(rr[field]);col='#356b8c' if rr.scenario=='primary' else '#555555';fill=col if rr.classification=='strong_colocalisation' else 'white'
    attr={'fill':fill,'stroke':col,'stroke-width':'.7'}
    if rr.scenario=='primary':ET.SubElement(g,q('circle'),{'cx':str(xx),'cy':str(cy),'r':'1.8',**attr})
    else:ET.SubElement(g,q('rect'),{'x':str(xx-1.7),'y':str(cy-1.7),'width':'3.4','height':'3.4',**attr})
  else:
   assert rows.n_overlap.nunique()==1
   val=rows.n_overlap.iloc[0];ET.SubElement(g,q('rect'),{'x':str(XX(0)),'y':str(cy-1.5),'width':str(XX(val)-XX(0)),'height':'3','fill':'#9ba9b3'})
# Compact scenario and fill symbols in left part of the column-header line.
for x,typ,label,fill in [(5,'o','Primary','white'),(21,'s','Leave-CAD','white'),(43,'o','Strong','#555555'),(54,'o','Other','white')]:
 cx=x*mm;cy=153*mm
 if typ=='o':ET.SubElement(g,q('circle'),{'cx':str(cx),'cy':str(cy),'r':'1.7','stroke':'#356b8c' if x==5 else '#555555','fill':fill,'stroke-width':'.7'})
 else:ET.SubElement(g,q('rect'),{'x':str(cx-1.6),'y':str(cy-1.6),'width':'3.2','height':'3.2','stroke':'#555555','fill':fill,'stroke-width':'.7'})
 text(g,cx+4,cy+1.9,label,size=5.7)
# Full old figure f is still traceable; legend explicitly preserves the stored classification.
export(r,'Figure3')
C.to_csv(O/'source_data/Figure3d_all29_unchanged.tsv',sep='\t',index=False)
pd.DataFrame({'gene':C.gene_symbol,'forest_y_pt':ycent,'evidence_y_pt':ycent,'delta_pt':0}).to_csv(O/'notes/Figure3_shared_row_check.tsv',sep='\t',index=False)
(O/'notes/Figure3_checks.json').write_text(json.dumps({'database_lines':groups,'all_29_forest_rows_unchanged':True,'evidence_matrix_SVG_unchanged':ET.tostring(mat)==ET.tostring(find(original,'FigureYa_vector_matrix_Figure3e')),'no_not_evaluated_cells_in_matrix':True,'f_existing_rows':10,'f_shared_SNP_counts_equal_between_scenarios':True,'new_models_or_statistics':False},indent=2))
(O/'notes/Figure3_4_legend_suggestions.md').write_text('Figure3c：GO BP、KEGG、Reactome 采用统一左对齐分组列，通路名左对齐；原MAGMA9条及pseudobulk9×7统计值不变。\n\nFigure3d/e：29个候选行严格共用原页面y中心；d从左侧增加有效森林图区。e仅删Not evaluated图例（现有145单元无此类），其余3类排成一行。\n\nFigure3f待采用图注：原五个gene–tissue contexts的Primary与Leave-CAD结果。新增PP.H3代表两个性状由不同因果变异解释的假设后验，与PP.H4区分；Shared SNPs为两个场景共同的n_overlap，不是独立样本数。圆形/方形分别为Primary/Leave-CAD，实心沿用原strong_colocalisation分类，空心为其他分类；PP.H4较高不自动等于强支持。新增列均读取原结果，不引入新模型或阈值。\n\nFigure4c/d：将NMF权重棒棒糖改为零起点水平柱形；原权重、颜色、顺序和横轴不变。ORA继续PENDING。\n')
print('Figure3/4 exported')
