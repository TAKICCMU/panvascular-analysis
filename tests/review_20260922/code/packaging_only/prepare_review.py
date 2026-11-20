from pathlib import Path
import shutil,hashlib,json
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919'); O=P/'图件编辑资料/全流程代码审核_20260922'; B=O/'Panvascular_code_review'; B.mkdir(exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
protected=[P/'Manuscript.docx',P/'Supplementary_Information.docx',P/'Supplementary_Tables.xlsx',P/'全部图件_15页审阅.pdf']+list((P/'Figures').glob('*.pdf'))+list((P/'Extended_Data').glob('*.pdf'))
(O/'protected_before.json').write_text(json.dumps({str(x):sha(x) for x in protected},indent=2))
for src,dst in [('Public_Release/Analysis_Code','analysis_code'),('Public_Release/Figure_Code','figure_code'),('Source_Data_Figures','evidence/source_data'),('Public_Release/Supplementary_Tables_TSV','evidence/supplementary_tables'),('Supplementary_Table20_full_results','evidence/table20_full_results')]:
 shutil.copytree(P/src,B/dst,dirs_exist_ok=True,ignore=shutil.ignore_patterns('.DS_Store','__pycache__','*.pyc'))
for d in ['audit','audit/logs','audit/code','recovered_sources','evidence/current_figures','evidence/release_contracts']:(B/d).mkdir(parents=True,exist_ok=True)
for src in ['Manuscript.docx','Supplementary_Information.docx','图件采用版本.tsv','Supplementary_Table_Headers.md','Changed_dictionary_and_headers.tsv']:
 shutil.copy2(P/src,B/'evidence'/Path(src).name)
for src in ['Figures','Extended_Data']:
 for x in (P/src).glob('*.pdf'):shutil.copy2(x,B/'evidence/current_figures'/x.name)
for x in (P/'Public_Release/Atlas_Metadata').glob('*'):
 if x.is_file() and x.stat().st_size<1000000:shutil.copy2(x,B/'evidence/release_contracts'/x.name)
shutil.copy2('/Users/taki/Desktop/Panvascular_Main_Figure_Refinement_20260920/merge_records/PATH_MAP.json',B/'recovered_sources/original_figure_PATH_MAP.json')
print(B)
print('FILES',len([x for x in B.rglob('*') if x.is_file()]))
