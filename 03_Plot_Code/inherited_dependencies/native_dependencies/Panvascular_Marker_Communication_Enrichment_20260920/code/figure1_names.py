from common import *
src=Path('/Users/taki/Desktop/Panvascular_Figure1b_Complete_Labels_20260920/figures/Figure1.svg');r=ET.parse(src).getroot();images=[ET.tostring(x) for x in r.iter(q('image'))]
remove_ids(r,{'Figure1b_all16_annotations'})
# Remove old number-to-name legend and its swatch objects, not the study legend.
parents={c:p for p in r.iter() for c in p}
for e in list(r.iter(q('text'))):
 if re.match(r'^\d+  ',e.text or '') or (e.text=='●' and 188<float(e.get('x','-1'))<198 and 117<float(e.get('y','-1'))<233):
  g=parents[e];parents[g].remove(g)
# Existing L1 swatches are circle/patches at this frozen legend position.
for g in list(r.iter(q('g'))):
 if not g.get('id','').startswith('patch_'):continue
 paths=list(g.iter(q('path')))
 if not paths:continue
 v=np.array(re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?',paths[0].get('d','')),float)
 if len(v)>4:
  xy=v.reshape(-1,2)
  if 188<xy[:,0].mean()<198 and 117<xy[:,1].mean()<232:
   if g in parents:parents[g].remove(g)
a=pd.read_csv('/Users/taki/Desktop/Panvascular_Figure1b_Complete_Labels_20260920/notes/Figure1b_16_label_mapping.tsv',sep='\t').sort_values('anchor_page_y_pt').reset_index(drop=True)
g=ET.SubElement(r,q('g'),{'id':'Figure1b_direct_cell_names'})
for i,row in a.iterrows():
 y=123+i*7.35; x=174
 ET.SubElement(g,q('path'),{'d':f'M {row.anchor_page_x_pt} {row.anchor_page_y_pt} L 163 {y-2} L {x-4} {y-2}','fill':'none','stroke':'#7b7b7b','stroke-width':'.3','stroke-opacity':'.75'})
 ET.SubElement(g,q('circle'),{'cx':str(row.anchor_page_x_pt),'cy':str(row.anchor_page_y_pt),'r':'.65','fill':row.hex})
 text(g,x,y,row.L1,6.5)
 a.loc[i,'label_x_pt']=x;a.loc[i,'label_y_pt']=y
assert images==[ET.tostring(x) for x in r.iter(q('image'))]
export(r,'Figure1');a.to_csv(O/'source_data/Figure1b_direct_names.tsv',sep='\t',index=False)
print('Figure1 names: 16; all image objects identical')
