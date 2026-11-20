import time,json,csv,subprocess,os,hashlib
from pathlib import Path
D=Path(__file__).parent.resolve();progress=D/'NEW_QTL_PROGRESS.json';deadline=time.time()+7200
while not (D/'NEW_QTL_COMPLETE.json').exists():
 if (D/'NEW_QTL_BLOCKED.json').exists():
  (D/'NEW_COLOC_STATUS.json').write_text(json.dumps({'status':'BLOCKED_RESOURCE_SERVICE','details':'NEW_QTL_BLOCKED.json'}));raise SystemExit(2)
 if time.time()>deadline:
  (D/'NEW_COLOC_STATUS.json').write_text(json.dumps({'status':'WAIT_TIMEOUT_NOT_COMPLETE'}));raise SystemExit(3)
 time.sleep(15)
R=D/'coloc_new12';(R/'00_metadata').mkdir(parents=True,exist_ok=True)
(R/'00_metadata/qtl_download_manifest.tsv').write_bytes((D/'new_qtl_manifest.tsv').read_bytes())
source=D.parent/'Tier1_full_colocalisation_v1/scripts/02_run_full_coloc.R';s=source.read_text().replace('ROOT <- "/Users/taki/Desktop/Codex任务1/Tier1_full_colocalisation_v1"','ROOT <- '+json.dumps(str(R),ensure_ascii=False));s=s.replace('library(data.table)','library(data.table)\n  setDTthreads(1)');s=s.split('# Plot gene-level maxima')[0]+'\nwriteLines(capture.output(sessionInfo()), file.path(ROOT, "00_metadata", "R_sessionInfo.txt"))\n';script=D/'run_new12_coloc.R';script.write_text(s)
(D/'NEW_COLOC_STATUS.json').write_text(json.dumps({'status':'RUNNING','method':'unchanged original ABF/prior grid; 1 thread'}))
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
with open(D/'coloc_new12.log','w') as log:p=subprocess.run(['Rscript',str(script)],env=env,stdout=log,stderr=subprocess.STDOUT)
status={'status':'COMPLETE_PENDING_VALIDATION' if p.returncode==0 else 'FAILED','exit_code':p.returncode}
(D/'NEW_COLOC_STATUS.json').write_text(json.dumps(status,indent=2))

if p.returncode==0:
 subprocess.run(["python3",str(D/"validate_and_integrate_all29.py")],check=True)
