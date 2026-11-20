from pathlib import Path
from openpyxl import load_workbook
from collections import Counter
import json
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');O=P/'图件编辑资料/补充表与方法同步_20260921'
w=load_workbook(P/'Supplementary_Tables.xlsx');out=[]
for s in w:
 nonempty={}
 formulas=[];links=[];comments=[];types=Counter()
 for row in s:
  vals=[(c.column,c.value) for c in row if c.value is not None]
  if vals:nonempty[row[0].row]=vals
  for c in row:
   if c.value is None:continue
   types[c.data_type]+=1
   if c.data_type=='f':formulas.append([c.coordinate,c.value])
   if c.hyperlink:links.append([c.coordinate,c.hyperlink.target,c.hyperlink.location])
   if c.comment:comments.append([c.coordinate,c.comment.text])
 # Need all structural separators, not only first rows.
 structure=[]
 for r,vs in nonempty.items():
  if r<=8 or len(vs)==1 or r-1 not in nonempty or (len(vs)>1 and len(nonempty.get(r-1,[]))==1):
   structure.append([r,vs])
 obj={'sheet':s.title,'rows':s.max_row,'cols':s.max_column,'structure':structure,'merged':list(map(str,s.merged_cells.ranges)),'types':dict(types),'formulas':formulas,'links':links,'comments':comments,'hidden_rows':[k for k,v in s.row_dimensions.items() if v.hidden],'hidden_cols':[k for k,v in s.column_dimensions.items() if v.hidden],'filter':s.auto_filter.ref,'freeze':s.freeze_panes,'validations':len(s.data_validations.dataValidation),'tables':list(s.tables)}
 out.append(obj)
 print(s.title,s.max_row,s.max_column,'MERGED',obj['merged'],'TYPES',dict(types))
 for r,v in structure:
  print(' ',r,repr(v)[:2400])
(O/'workbook_structure.json').write_text(json.dumps(out,ensure_ascii=False,indent=2,default=str))
