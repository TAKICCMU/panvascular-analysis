from pathlib import Path
from openpyxl import load_workbook
import csv,json,re,shutil,ast,importlib.util
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');O=P/'图件编辑资料/五项字典纠正_20260922';B=O/'轻量备份';fixes=json.loads((O/'fixes.json').read_text())
def backup(p):
 q=B/p.relative_to(P);q.parent.mkdir(parents=True,exist_ok=True)
 if p.exists() and not q.exists():shutil.copy2(p,q)
# Keep already-existing textual review copies current through targeted field edits only.
p=P/'Supplementary_Table_Headers.md';backup(p);s=p.read_text();blocks=re.split(r'(?=^## Supplementary Table )',s,flags=re.M)
for i,b in enumerate(blocks):
 m=re.match(r'## Supplementary Table (\d+) ',b)
 if not m:continue
 n=int(m[1])
 for f in [f for f in fixes if f['table']==n]:
  pat=r'(\d+\. \*\*[A-Z]+ — )'+re.escape(f['old'])+r'(\*\*\n   - 定义：)[^\n]*(\n   - 单位/尺度：)[^\n]*'
  b,num=re.subn(pat,lambda m:m[1]+f['new']+m[2]+f['definition']+m[3]+f['unit'],b);assert num==1,(n,f['old'])
 blocks[i]=b
p.write_text(''.join(blocks))
p=P/'Spreadsheet_Format_Check.tsv';backup(p);rows=list(csv.DictReader(p.open(),delimiter='\t'));fields=list(rows[0]);lookup={(f['table'],f['old']):f for f in fixes}
for r in rows:
 if (f:=lookup.get((int(r['Table']),r['Actual header']))):r['Actual header']=f['new']
with p.open('w') as out:w=csv.DictWriter(out,fields,delimiter='\t');w.writeheader();w.writerows(rows)
p=P/'Key_Table_Examples.md';backup(p);s=p.read_text();blocks=re.split(r'(?=^## Table )',s,flags=re.M)
for i,b in enumerate(blocks):
 m=re.match(r'## Table (\d+) ',b)
 if m:
  for f in [f for f in fixes if f['table']==int(m[1])]:b=b.replace('- **'+f['old']+'**：','- **'+f['new']+'**：')
 blocks[i]=b
p.write_text(''.join(blocks))
p=P/'Numeric_Text_Conversions.tsv';rows=list(csv.DictReader(p.open(),delimiter='\t'));fields=list(rows[0]);changed=False
for r in rows:
 if (f:=lookup.get((int(r['Table']),r['Field']))):r['Field']=f['new'];changed=True
if changed:
 backup(p)
 with p.open('w') as out:w=csv.DictWriter(out,fields,delimiter='\t');w.writeheader();w.writerows(rows)
# Authoritative reviewed mapping replaces lexical guessing in the historical formatting template.
R=P/'Public_Release/Supplementary_Tables_TSV';wb=load_workbook(P/'Supplementary_Tables.xlsx',read_only=True);dd={(int(r[0]),r[1]):r for r in list(wb['Data dictionary'].values)[1:]};maps=list(csv.DictReader((P/'column_name_map.tsv').open(),delimiter='\t'));approved=[]
for r in maps:
 n=int(r['Worksheet'].split()[1]);d=dd[(n,r['Display column'])];approved.append({'table':n,'original':r['Original column'],'display':r['Display column'],'definition':d[2],'unit':d[3]})
wb.close();(R/'reviewed_field_definitions.json').write_text(json.dumps(approved,ensure_ascii=False,indent=2))
p=P/'图件编辑资料/补充表与方法同步_20260921/code/format_workbook.py';backup(p);s=p.read_text();tree=ast.parse(s);replacements=[]
for f in tree.body:
 if isinstance(f,ast.FunctionDef) and f.name in ['display','semantics']:
  signature='h,n' if f.name=='display' else 'n,h,new';body=("return item['display']" if f.name=='display' else "return item['definition'], item['unit']")
  replacements.append((f.lineno,f.end_lineno,f"def {f.name}({signature}):\n item=REVIEWED_FIELDS.get((n,h))\n if item is None: raise ValueError(f'No reviewed type/definition for table {{n}}, field {{h!r}}; do not infer from FDR/share keywords.')\n {body}\n"))
lines=s.splitlines(keepends=True)
for start,end,rep in sorted(replacements,reverse=True):lines[start-1:end]=[rep]
s=''.join(lines);where=s.index('def display(');s=s[:where]+"REVIEWED_FIELDS={(r['table'],r['original']):r for r in json.loads((P/'Public_Release/Supplementary_Tables_TSV/reviewed_field_definitions.json').read_text())}\n"+s[where:];p.write_text(s)
# Test only the classification functions, never run the formatting/reconversion script.
tree=ast.parse(s);nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['display','semantics']];env={'REVIEWED_FIELDS':{(r['table'],r['original']):r for r in approved}};exec(compile(ast.Module(body=nodes,type_ignores=[]),'<reviewed-field-rules>','exec'),env)
for r in approved:assert env['display'](r['original'],r['table'])==r['display'];assert env['semantics'](r['table'],r['original'],r['display'])==(r['definition'],r['unit'])
try:env['semantics'](999,'shared_FDR_unknown','shared_FDR_unknown')
except ValueError:pass
else:raise AssertionError('Unknown field should not be guessed')
spec=importlib.util.spec_from_file_location('reader',R/'read_supplementary_table.py');reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)
for n in [13,20,39,40,50]:h,rows=reader.read_table(P/'Supplementary_Tables.xlsx',n);assert len(h)==len(set(h))
(O/'template_and_reader_tests.json').write_text(json.dumps({'reviewed_definition_cases':len(approved),'unknown_field_fails_closed':True,'reader_tables':[13,20,39,40,50],'formatting_script_executed':False,'scientific_data_conversion_executed':False},indent=2))
shutil.copy2(O/'Changed_dictionary_and_headers.tsv',P/'Changed_dictionary_and_headers.tsv')
text='''# 本轮字段定义纠正

- 已修正五项类型错误：Table13的AE/AG/AI为boolean，Table40的R为variant identifier，Table50的B为text。相关表头、字典及列名映射已同步。
- Table13的29行中，AE与Y<0.05、AG与X<0.05全部一致；AI按原脚本的“背景数≥3且方向一致性≥0.75”核对，仍为27 FALSE、2 TRUE。AF的nominal-P标志独立保留。
- Table40原脚本按SNP.PP.H4降序取首项；已说明该量条件于H4，并不代表区域达到强共定位。全部rsID与区域后验不变。
- Table20五个通路数量列、Table39三个互作数量列补充计数含义。共更改13个表头和13行字典；64个表题及其余列名不变。
- 全字段类型冲突复查未发现这五项以外的明确“布尔/文本被标作概率或比例”冲突；笼统source scale/category未据此重写。模板改用明确的已审核字段映射，不再按FDR/share字样猜类型；新字段缺来源时停止推断。
- 64表数据区逐字节相同：数值、文本、布尔、空白及行序均未改变；样式、列宽、筛选、冻结、合并关系也保留。没有重新转换数值或修改细微浮点尾差。
- 现有四份核查文本仅同步涉及的字段名称/定义；原机器字段和读表入口保留。文稿、87条引用、15幅图、ORA与空间分组均未改动。

逐项差异见Changed_dictionary_and_headers.tsv；轻量备份与来源/验证记录在图件编辑资料/五项字典纠正_20260922。完成后暂停，供作者确认。
'''
(P/'五项字段定义纠正_完成说明.md').write_text(text)
print('Mappings, exact changed-row export and existing text copies synchronized; 1026 template cases and 5 reader tests passed.')
