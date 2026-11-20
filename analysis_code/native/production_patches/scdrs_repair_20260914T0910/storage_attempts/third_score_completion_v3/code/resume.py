#!/usr/bin/env python3
"""Resume only the completed third score validation and native L1/L2 tests."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ['CUDA_VISIBLE_DEVICES']=''
os.environ['MALLOC_CONF']='background_thread:false'
from pathlib import Path
import sys,json,time,subprocess,fcntl,traceback,signal
W=Path(__file__).resolve().parents[1]
C=json.loads((W/'config.json').read_text())
sys.path.insert(0,str(Path(C['previous_attempt'])/'code'))
import coordinator as legacy
from storage_utils import sha,write,proc,same,send,mount_budget,cpu_budget
D=Path(C['output'])
def admitted(other,stage):
 return other+1+stage<=40
def verify():
 legacy.verify()
 receipt=json.loads((W/'RECEIPT.json').read_text())
 for x in receipt['bound_files']:
  if sha(x['path'])!=x['sha256']:raise ValueError('Bound content changed '+x['path'])
 for x in receipt['score_stat_bindings']:
  s=Path(x['path']).stat()
  for k in ['st_size','st_mtime_ns','st_ino','st_dev']:
   if getattr(s,k)!=x[k]:raise ValueError('Reused score changed')
def resources(stage,target=None):
 own={os.getpid()}|({target['pid']} if target else set())
 cpu=cpu_budget(legacy.C,own)
 actual=cpu['actual_other_threads']+proc(os.getpid())['threads']
 if target and same(target):actual+=proc(target['pid'])['threads']
 b=mount_budget(legacy.C)
 out=dict(time=time.time(),actual_project_threads=actual,other_threads=cpu['actual_other_threads'],stage_reservation=stage,projected=cpu['actual_other_threads']+1+stage,pass_gate=admitted(cpu['actual_other_threads'],stage),budget=b,others=cpu['others'])
 write(W/'RESOURCE_STATUS.json',out)
 return out
def run(cmd,log,native=False):
 stage=13 if native else 1
 while True:
  verify();r=resources(stage)
  if r['pass_gate'] and not r['budget']['reasons']:break
  write(W/'COORD_STATUS.json',dict(status='WAITING_ACTUAL_RESOURCE',resource=r));time.sleep(15)
 env=os.environ.copy()
 if native:
  for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:env[k]='12'
 with open(log,'x') as f:
  p=subprocess.Popen(list(map(str,cmd)),stdout=f,stderr=subprocess.STDOUT,env=env,stdin=subprocess.DEVNULL,start_new_session=True)
  t=proc(p.pid)
  write(W/'CURRENT_CHILD.json',dict(target=t,cmd=list(map(str,cmd)),log=str(log),native=native))
  write(W/'COORD_STATUS.json',dict(status='RUNNING',pid=p.pid,cmd=list(map(str,cmd))))
  try:
   while p.poll() is None:
    r=resources(stage,t)
    if r['actual_project_threads']>40 or r['budget']['reasons']:raise RuntimeError('Actual project resource guard exceeded')
    time.sleep(15)
  except BaseException:
   if same(t):os.killpg(p.pid,signal.SIGSTOP)
   raise
 record=dict(cmd=list(map(str,cmd)),exit_code=p.returncode,log=str(log),completed=time.time())
 history=json.loads((W/'COMMANDS.json').read_text()) if (W/'COMMANDS.json').exists() else []
 history.append(record);write(W/'COMMANDS.json',history)
 if p.returncode:raise RuntimeError('Native/validator failed '+str(record))
def main():
 with open(legacy.C['project_lock'],'a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  verify()
  if (W/'STARTED.json').exists():raise FileExistsError('No automatic replay')
  write(W/'STARTED.json',dict(pid=os.getpid(),time=time.time()))
  validator=Path(C['previous_attempt'])/'code/validate_native.py'
  py='/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python'
  run([py,validator,'score','--run',D],D/'score_validation.log')
  schema=json.loads((D/'SCORE_VALIDATION.json').read_text())
  score=list((D/'scores').glob('*.full_score.gz'))
  assert len(score)==1
  for key in legacy.OLD['group_keys']:
   group=D/key;group.mkdir(exist_ok=False)
   run([legacy.OLD['scdrs_executable'],'perform-downstream','--h5ad-file',legacy.OLD['derived'],'--score-file',score[0],'--out-folder',group,'--group-analysis',key,'--flag-filter-data','True','--flag-raw-count','True'],D/(key+'.log'),True)
   run([py,validator,'group','--run',D,'--key',key],D/(key+'_validation.log'))
   schema[key]=json.loads((D/(key+'_VALIDATION.json')).read_text())
  schema['trait']='leave_one_CAD_r0.1'
  write(D/'SCHEMA_VALIDATION.json',schema)
  verify()
  manifests={}
  for trait in ['primary_r0.1','primary_r0.6']:
   source=Path(legacy.C['output'])/trait/'SCHEMA_VALIDATION.json'
   manifests[trait]=dict(path=str(source),sha256=sha(source),reused=True)
  manifests['leave_one_CAD_r0.1']=dict(path=str(D/'SCHEMA_VALIDATION.json'),sha256=sha(D/'SCHEMA_VALIDATION.json'),score_recomputed=False)
  write(W/'COMPLETE.json',dict(status='COMPLETE_NATIVE_3_SCENARIOS_VALIDATED',scenarios=manifests,score_recomputed=False,time=time.time()))
  write(W/'COORD_STATUS.json',dict(status='COMPLETE',time=time.time()))
if __name__=='__main__':
 try:main()
 except BaseException:
  write(W/'COORD_EXIT.json',dict(status='FAILED_REVIEW_REQUIRED',error=traceback.format_exc(),time=time.time()))
  raise
