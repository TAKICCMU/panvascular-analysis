from pathlib import Path
import subprocess,json,hashlib,csv
B=Path(__file__).resolve().parents[2]; cache={};rows=[]
for root in [B/'analysis_code',B/'figure_code',B/'recovered_sources']:
 for p in root.rglob('*'):
  if not p.is_file() or p.suffix.lower() not in ['.r','.sh']:continue
  h=hashlib.sha256(p.read_bytes()).hexdigest();key=(p.suffix.lower(),h)
  if key not in cache:
   if p.suffix.lower()=='.r':cmd=['/usr/local/bin/Rscript','--vanilla','-e','invisible(parse(commandArgs(TRUE)[1]))',str(p)]
   else:cmd=['bash','-n',str(p)]
   r=subprocess.run(cmd,capture_output=True,text=True,timeout=30);cache[key]=(r.returncode,r.stderr[-1500:])
  code,err=cache[key];rows.append({'path':str(p.relative_to(B)),'language':'R' if p.suffix.lower()=='.r' else 'Shell','status':'PASS_PARSE' if code==0 else 'PARSE_REVIEW','detail':err,'sha256':h})
with (B/'audit/R_SHELL_SYNTAX.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t');w.writeheader();w.writerows(rows)
print(json.dumps({'files':len(rows),'unique_source_contents':len(cache),'parse_failures':[x for x in rows if x['status']!='PASS_PARSE']}))
