from pathlib import Path
import csv,json,hashlib
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
f=P/'SUBMISSION_FILE_LIST.tsv';rows=list(csv.DictReader(f.open(),delimiter='\t'));fields=list(rows[0]);names=['Supplementary_Table_Headers.md','Spreadsheet_Format_Check.tsv','Numeric_Text_Conversions.tsv','Key_Table_Examples.md','Changed_dictionary_and_headers.tsv','五项字段定义纠正_完成说明.md','Public_Release/Supplementary_Tables_TSV/reviewed_field_definitions.json']
for n in names:
 if not any(r['path']==n for r in rows):rows.append(dict(path=n,purpose='Read-only Excel review export',status='CURRENT',bytes='',sha256=''))
for r in rows:
 q=P/r['path']
 if q.is_file():r['bytes']=q.stat().st_size;r['sha256']=sha(q)
with f.open('w') as s:w=csv.DictWriter(s,fields,delimiter='\t');w.writeheader();w.writerows(rows)
f=P/'文件来源与校验.tsv';rows=list(csv.DictReader(f.open(),delimiter='\t'));fields=list(rows[0])
for n in names:
 if not any(r['path']==n for r in rows):rows.append(dict(path=n,bytes='',sha256=''))
for r in rows:
 q=P/r['path'];assert q.exists();r['bytes']=q.stat().st_size;r['sha256']=sha(q)
with f.open('w') as s:w=csv.DictWriter(s,fields,delimiter='\t');w.writeheader();w.writerows(rows)
(P/'PACKAGE_SHA256SUMS.txt').write_text(''.join(r['sha256']+'  '+r['path']+'\n' for r in rows))
print('Current manifest updated:',len(rows),'files')
