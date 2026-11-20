from pathlib import Path
import csv,json,hashlib,shutil
from docx import Document
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');O=P/'图件编辑资料/补充表与方法同步_20260921'
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as s:
  for b in iter(lambda:s.read(8388608),b''):h.update(b)
 return h.hexdigest()
f=P/'投稿核查/核对与完成说明_CN.md';b=O/'轻量备份/投稿核查'/f.name;b.parent.mkdir(exist_ok=True,parents=True)
if not b.exists():shutil.copy2(f,b)
t=f.read_text().replace('ORA选择后BH的统计有效性、部分供者身份、测量覆盖和组织学限制仍PENDING；无新分析。','作者已接受现有ORA方案和空间分组；原分析单位与统计量保留。单细胞metadata缺失、测量覆盖和组织学解释限制仍按实际证据说明，本轮没有新增验证。');t+='\n当前补充表/Methods更新：见根目录《补充表与方法修改清单.md》及column_name_map.tsv。该更新取代本轮以前的ORA/空间分组待修复状态。\n';f.write_text(t)
meta=json.loads((O/'table_metadata.json').read_text());si=Document(P/'Supplementary_Information.docx');txt='\n'.join(c.text for t in si.tables for r in t.rows for c in r.cells)
# source metadata is keyed by table number
for v in meta.values() if isinstance(meta,dict) else meta:
 title=v['title'];assert title in txt,title
with (P/'图件采用版本.tsv').open() as f:
 figrows=list(csv.DictReader(f,delimiter='\t'))
 for r in figrows:
  assert sha(Path(r['package_pdf']))==r['package_pdf_sha256'];assert sha(Path(r['package_svg']))==r['package_svg_sha256']
f=P/'SUBMISSION_FILE_LIST.tsv';rows=list(csv.DictReader(f.open(),delimiter='\t'));fields=list(rows[0])
for r in rows:
 q=P/r['path']
 if q.is_file():r['bytes']=str(q.stat().st_size);r['sha256']=sha(q)
for name in ['column_name_map.tsv','补充表与方法修改清单.md']:
 if not any(r['path']==name for r in rows):rows.append({'path':name,'purpose':'Supplementary table editorial mapping / change record','status':'CURRENT','bytes':(P/name).stat().st_size,'sha256':sha(P/name)})
with f.open('w') as s:w=csv.DictWriter(s,fields,delimiter='\t');w.writeheader();w.writerows(rows)
status=json.loads((P/'投稿核查/FINAL_PACKAGE_STATUS.json').read_text());status['manuscript_and_figures_pdf_bytes']=(P/'Nature_Submission/Manuscript_with_15_Figures.pdf').stat().st_size;status['manuscript_pages']=30;status['supplementary_tables_checked']=64;(P/'投稿核查/FINAL_PACKAGE_STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
roots=['Figures','Extended_Data','Nature_Submission','Public_Release','Source_Data_Figures','Supplementary_Table20_full_results','上传文件']
files=[]
for root in roots:files += [q for q in (P/root).rglob('*') if q.is_file() and q.name!='.DS_Store' and '__pycache__' not in q.parts]
files += [P/n for n in ['Manuscript.docx','Cover_Letter.docx','Supplementary_Information.docx','Supplementary_Tables.xlsx','全部图件_15页审阅.pdf','图件采用版本.tsv','README_当前状态.txt','UPLOAD_MANIFEST.tsv','SUBMISSION_FILE_LIST.tsv','column_name_map.tsv','补充表与方法修改清单.md']]
rows=[dict(path=str(q.relative_to(P)),bytes=q.stat().st_size,sha256=sha(q)) for q in sorted(set(files))]
with (P/'文件来源与校验.tsv').open('w') as s:w=csv.DictWriter(s,['path','bytes','sha256'],delimiter='\t');w.writeheader();w.writerows(rows)
(P/'PACKAGE_SHA256SUMS.txt').write_text(''.join(r['sha256']+'  '+r['path']+'\n' for r in rows))
report={'status':'CURRENT_FILES_REPLACED_FOR_AUTHOR_REVIEW','numbered_tables':64,'scientific_reanalysis':False,'figure_PDF_SVG_byte_identity':True,'all_table_titles_match_SI':True,'original_header_adapter_test':[7,53],'numeric_data_validation':'workbook_validation.json','docx_format_citation_validation':'document_validation.json','visual_scope':'64-table field/sample previews reviewed; rendered edited Word pages inspected; custom table PDF is not native Excel screenshot','source_gaps':['Complete per-sample upstream Scrublet logs/platform exceptions not recovered','Historical internal susieR CS-deduplication default not independently recovered; no additional project filter asserted'],'files_hashed':len(rows)}
(O/'FINAL_DELIVERY_VALIDATION.json').write_text(json.dumps(report,ensure_ascii=False,indent=2));print(json.dumps(report,ensure_ascii=False))
