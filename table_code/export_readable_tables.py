"""Read-only, lossless value export of the adopted workbook and DOCX context."""
from pathlib import Path
import csv, json, datetime, zipfile, xml.etree.ElementTree as ET, gzip
import openpyxl

WORK=Path('/Users/taki/Desktop/泛血管稿件_批注修订_20261116/图件编辑资料/重作图数据包')
R=WORK/'package'
wb=openpyxl.load_workbook(R/'02_Plot_Data/Supplementary_Tables.xlsx',read_only=True,data_only=False)
OUT=R/'02_Plot_Data/Supplementary_Tables_Text'
OUT.mkdir(exist_ok=True)
catalog={int(row[0]):row for row in list(wb['Information'].values)[1:] if isinstance(row[0],(int,float))}
def safe(v):
    if isinstance(v,(datetime.datetime,datetime.date)):return v.isoformat()
    return v
index=[];formulas=[]
for ws in wb:
    is_table=ws.title.startswith('Table ')
    n=int(ws.title.split()[-1]) if is_table else None
    stem=('Table_%02d'%n) if n else ws.title.replace(' ','_')
    header_row=2 if is_table else 1
    row_count=0;nonempty=0
    with (OUT/(stem+'.tsv')).open('w',newline='') as f, gzip.open(OUT/(stem+'.rows.jsonl.gz'),'wt',encoding='utf-8') as j:
        w=csv.writer(f,delimiter='\t')
        for rowno,cells in enumerate(ws.iter_rows(),1):
            vals=[safe(c.value) for c in cells]
            for c in cells:
                if c.data_type=='f':formulas.append({'sheet':ws.title,'cell':c.coordinate,'formula':c.value})
            j.write(json.dumps({'excel_row':rowno,'values':vals},ensure_ascii=False,allow_nan=False)+'\n')
            if rowno>=header_row:w.writerow(vals)
            if rowno>header_row:
                row_count+=1;nonempty+=int(any(v is not None for v in vals))
    info=catalog.get(n,[])
    index.append({'table':n or '', 'sheet':ws.title,'title':info[1] if info else ws.title,
                  'row_unit':info[3] if info else '', 'excel_rows':ws.max_row,'excel_columns':ws.max_column,
                  'header_excel_row':header_row,'data_rows_including_blanks':row_count,'nonempty_data_rows':nonempty,
                  'TSV':str((OUT/(stem+'.tsv')).relative_to(R)),
                  'typed_rows_JSONL_GZ':str((OUT/(stem+'.rows.jsonl.gz')).relative_to(R))})
with (R/'06_Maps/TABLE_INDEX.tsv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(index[0]),delimiter='\t');w.writeheader();w.writerows(index)
(OUT/'formula_cells.json').write_text(json.dumps(formulas,ensure_ascii=False,indent=2))
(OUT/'README.md').write_text('''# 现行工作簿的只读文本导出

唯一来源是包内 `02_Plot_Data/Supplementary_Tables.xlsx`，60张表，另有Information、Abbreviations和Data dictionary。
没有重新计算、筛行、排序、改名、修改数值或数字格式。TSV第一行为实际列名，表格标题见TABLE_INDEX.tsv；原工作簿前两行中的标题行不混入数据。
`*.rows.jsonl.gz`逐行保存完整Excel值及原行号，包含标题行。JSON中的null对应空单元格，0、FALSE和字面字符串NA保持不同；TSV的空字段对应空单元格。
JSON保留Excel存储值的数值/文本/布尔类型，显示精度由原XLSX保留。仅将日期序列化为ISO文本。
跨记录类型的空字段不能自动补0；请同时读Data_dictionary.tsv中的字段定义、尺度和缺失语义。
这些文件便于阅读，不替代原始Source Data或原XLSX。未执行公式；公式单元格（若有）在formula_cells.json中逐项列出。
''')

NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
for src in sorted((R/'05_Manuscript_Context').glob('*.docx')):
    with zipfile.ZipFile(src) as z:
        doc=ET.fromstring(z.read('word/document.xml'))
    lines=['# '+src.name+'（只读正文文本导出）','',
           '引用域、表格版式及格式以原DOCX为准。本文件只方便检索，不用于替换文稿。','']
    for i,p in enumerate(doc.findall('.//w:p',NS),1):
        text=''.join(p.itertext()) if False else ''.join(t.text or '' for t in p.findall('.//w:t',NS))
        if text.strip():lines.append(f'[{i}] {text}')
    src.with_suffix('.md').write_text('\n\n'.join(lines))
(WORK/'text_export_check.json').write_text(json.dumps({'table_count':sum(x['table']!='' for x in index),'sheets':len(index),'formula_cell_count':len(formulas),'data_rows':sum(x['nonempty_data_rows'] for x in index)},indent=2))
print((WORK/'text_export_check.json').read_text())
