import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ['CUDA_VISIBLE_DEVICES']=''
import json,signal,subprocess,hashlib,time,fcntl,csv,sys
from pathlib import Path
B=Path(__file__).resolve().parents[1];W=B.parent
sys.path.insert(0,str(W/'code'));from common import processes
PREP_PID=int(sys.argv[1]);proc=Path('/proc')/str(PREP_PID)
assert str(B/'code/prepare.py') in proc.joinpath('cmdline').read_bytes().replace(b'\0',b' ').decode()
start=proc.joinpath('stat').read_text().split()[21]
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for z in iter(lambda:f.read(8<<20),b''):h.update(z)
 return h.hexdigest()
assert json.loads((B/'tests/native_v1/PASS.json').read_text())['status']=='PASS'
(B/'results').mkdir(exist_ok=True);early=[]
# Stop only this branch's new export process while running early native fits.
os.kill(PREP_PID,signal.SIGSTOP)
(B/'EARLY_PREP_PAUSE.json').write_text(json.dumps({'pid':PREP_PID,'start_ticks':start,'reason':'Bounded early CKB fits using already-completed L1 input; resume in finally','time':time.time()}))
try:
 for state in ['Macrophage','Fibroblast','Smooth muscle cell']:
  for ctx in ['G1','G2','G3']:
   name='L1_'+ctx+'_'+hashlib.sha256(state.encode()).hexdigest()[:12];p=B/'inputs'/name
   if not (p/'units.tsv').exists():continue
   keys=[B/'code/fit.R',B/'tests/native_v1/PASS.json',B/'genes.tsv',p/'counts.mtx.gz',p/'units.tsv',W/'disease_background_3groups_review_v2/sample_background_sidecar_v2.tsv']
   bound={str(q):sha(q) for q in keys}
   receipt={'authorization':'User approved three-background analysis and CKB focus','model':name,'cell_state':state,'context':ctx,'input_source':'hash-verified frozen sample-level raw pseudobulk; completed L1 export','bound_files':bound,'prep_pid_suspended':PREP_PID,'scope':'one complete all-gene edgeR test, CKB extracted after full-gene BH'}
   rp=B/'logs'/(name+'_EARLY_RECEIPT.json');rp.write_text(json.dumps(receipt,indent=2))
   handle=None
   while not handle:
    ps=processes();n=sum(x['threads'] for x in ps)
    if n+2<=40:
     for i in range(4):
      h=(W/'locks'/('cpu_heavy_'+str(i)+'.lock')).open('a')
      try:fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB);handle=h;break
      except BlockingIOError:h.close()
    if not handle:time.sleep(5)
   cmd=['/usr/bin/Rscript','--vanilla',str(B/'code/fit.R'),str(p),str(B/'genes.tsv'),str(B/'results'/name)]
   with (B/'logs'/(name+'.log')).open('x') as log:
    child=subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    (B/'EARLY_PROGRESS.json').write_text(json.dumps({'status':'RUNNING_NATIVE_EDGER','controller_pid':os.getpid(),'model':name,'cell_state':state,'context':ctx,'pid':child.pid,'command':cmd,'completed':early},indent=2));print('START',name,child.pid,flush=True);rc=child.wait()
   handle.close();res=B/'results'/name
   if rc!=0:raise RuntimeError('Early model failed '+name)
   model=json.loads((res/'MODEL.json').read_text());entry={'model_id':name,'receipt':str(rp),'receipt_sha256':sha(rp),'result_status':model['status'],'MODEL_sha256':sha(res/'MODEL.json'),'DE_sha256':sha(res/'DE.tsv') if (res/'DE.tsv').exists() else None}
   early.append(entry);(B/'EARLY_VALIDATED_MODELS.json').write_text(json.dumps(early,indent=2))
   if (res/'CKB_result.tsv').exists():print((res/'CKB_result.tsv').read_text(),flush=True)
 finally_dummy=False
finally:
 if proc.exists() and proc.joinpath('stat').read_text().split()[21]==start:
  os.kill(PREP_PID,signal.SIGCONT);print('PREPARATION_RESUMED',PREP_PID,flush=True)
 (B/'EARLY_PREP_RESUMED.json').write_text(json.dumps({'pid':PREP_PID,'start_ticks':start,'time':time.time()}))
