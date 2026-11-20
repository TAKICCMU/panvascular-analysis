from pathlib import Path
import json,csv,gzip,math
from openpyxl import load_workbook
B=Path(__file__).resolve().parents[2];T=B/'evidence/supplementary_tables'
structure=json.loads((T/'table_structure.json').read_text());mapping=list(csv.DictReader((T/'column_name_map.tsv').open(),delimiter='\t'))
w=load_workbook(B/'evidence/Supplementary_Tables.xlsx',read_only=True,data_only=False);rows=[]
def equivalent(x,y):
 if x is None:return y==''
 if isinstance(x,bool):return str(x).lower()==str(y).lower()
 if isinstance(x,(int,float)):
  try:return math.isclose(float(x),float(y),rel_tol=1e-12,abs_tol=1e-14)
  except:return False
 return str(x)==y
for s in structure:
 name=s['sheet'];ws=w[name];m=[v for v in mapping if v['Worksheet']==name]
 header=list(next(ws.iter_rows(min_row=3,max_row=3,values_only=True)))
 with gzip.open(T/(name+'.tsv.gz'),'rt') as f:export=list(csv.reader(f,delimiter='\t'))
 original=export[s['data_start']-1:s['data_end']];cells=0;mismatches=[]
 for rn,(old,new) in enumerate(zip(original,ws.iter_rows(min_row=s['data_start'],max_row=s['data_end'],values_only=True)),s['data_start']):
  for cn,(a,b) in enumerate(zip(old,new),1):
   cells+=1
   if not equivalent(b,a):mismatches.append([rn,cn])
 rows.append(dict(table=name,rows=len(original),cells=cells,header_matches=header==[r['Display column'] for r in m],machine_headers_match=export[2]==[r['Original column'] for r in m],data_mismatches=len(mismatches),mismatch_cells=str(mismatches[:30])))
with (B/'audit/WORKBOOK_MACHINE_DATA_CHECK.tsv').open('w') as f:
 writer=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t');writer.writeheader();writer.writerows(rows)
print(json.dumps({'tables':len(rows),'cells':sum(x['cells'] for x in rows),'data_mismatches':sum(x['data_mismatches'] for x in rows),'header_mismatches':sum(not x['header_matches'] for x in rows)}))
