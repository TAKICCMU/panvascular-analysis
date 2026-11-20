#!/usr/bin/env python3
"""Adopt existing score PID, then continue original CLI sequence on /data; no score replay."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ['CUDA_VISIBLE_DEVICES']='';os.environ['MALLOC_CONF']='background_thread:false'
from pathlib import Path
import json,time,signal,subprocess,fcntl,shutil,traceback,sys
from storage_utils import *
W=Path(__file__).resolve().parents[1];C=json.loads((W/'config.json').read_text());OLD=json.loads(Path(C['original_config']).read_text());D=Path(C['output']);COMMANDS=[]

def verify():
 for file in [W/'RECEIPT.json',W/'COORDINATOR_RECEIPT.json',Path(C['original_receipt'])]:
  r=json.loads(file.read_text())
  for x in r['bound_files']:
   if sha(x['path'])!=x['sha256']:raise ValueError('Bound input/code changed '+x['path'])
 # The new original attempt certificate binds both full hashes to this current boot/stat.
 for e in OLD['authenticated_files']:
  cert=json.loads(Path(e['certificate']).read_text());s=Path(e['path']).stat()
  if cert['observed_sha256']!=e['sha256'] or not cert['full_sha_verified'] or cert['boot_id']!=Path('/proc/sys/kernel/random/boot_id').read_text().strip():raise ValueError('Certified identity changed')
  for k,v in [('size',s.st_size),('inode',s.st_ino),('device',s.st_dev),('mtime_ns',s.st_mtime_ns),('nlink',s.st_nlink)]:
   if cert[k]!=v:raise ValueError('Certified file stat changed '+k)
 if sha(OLD['covariates'])!=OLD['covariates_sha256']:raise ValueError('Covariates changed')
 for s in OLD['scenarios']:
  if sha(s['geneset'])!=s['sha256']:raise ValueError('Gene-set changed')
 if OLD['native_threads']!=12 or OLD['n_ctrl']!=1000 or OLD['graph_policy']!='frozen_global_n30':raise ValueError('Original parameters differ')

def refresh(target=None,paused=False):
 b=mount_budget(C);cpu=cpu_budget(C,{os.getpid()}|({target['pid']} if target else set()))
 own_threads=proc(os.getpid())['threads']+(proc(target['pid'])['threads'] if target and same(target) and proc(target['pid'])['state'] not in ['T','t','Z'] else 0)
 danger=any(x['free_GiB']<x['floor_GiB']+2 for x in b['filesystems'].values()) or cpu['actual_other_threads']+own_threads>40
 if target and same(target):
  state=proc(target['pid'])['state']
  if danger and not paused:
   if state in ['T','t']:raise ValueError('External pause not owned by coordinator')
   send(target,signal.SIGSTOP);paused=True;write(W/'RUNTIME_PAUSE.json',dict(time=time.time(),target=target,budget=b,cpu=cpu))
  elif paused and not b['reasons'] and cpu['pass_gate']:
   send(target,signal.SIGCONT);paused=False;write(W/'RUNTIME_RESUME.json',dict(time=time.time(),target=target,budget=b,cpu=cpu))
 write(W/'RESOURCE_STATUS.json',dict(time=time.time(),budget=b,cpu=cpu,paused_owned=paused))
 return b,cpu,paused

def gate():
 while True:
  verify();b,cpu,_=refresh()
  if not b['reasons'] and cpu['pass_gate']:return
  write(W/'COORD_STATUS.json',dict(status='WAITING_MOUNT_OR_CPU_RESOURCE',time=time.time(),budget=b,cpu=cpu));time.sleep(20)

def call(cmd,log,native=False):
 gate();env=os.environ.copy()
 if native:
  for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:env[k]='12'
 if Path(log).exists():raise FileExistsError('Do not overwrite command log')
 with open(log,'w') as f:
  child=subprocess.Popen([str(x) for x in cmd],stdout=f,stderr=subprocess.STDOUT,env=env,stdin=subprocess.DEVNULL)
  target=proc(child.pid);write(W/'CURRENT_CHILD.json',dict(target=target,cmd=[str(x) for x in cmd],log=str(log),time=time.time()));paused=False
  while child.poll() is None:
   try:_,_,paused=refresh(target,paused)
   except Exception:
    if same(target) and proc(target['pid'])['state'] not in ['T','t','Z']:send(target,signal.SIGSTOP)
    raise
   time.sleep(15)
  rc=child.returncode
 COMMANDS.append(dict(cmd=[str(x) for x in cmd],exit_code=rc,log=str(log),completed=time.time()));write(D/'COMMANDS.json',COMMANDS)
 if rc:raise RuntimeError('Command failed '+str(rc)+' '+str(log))

def validate_score(run):call(['/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python',W/'code/validate_native.py','score','--run',run],run/'score_validation.log')

def adopt_first():
 native=C['native'];parent=C['old_controller'];paused=False
 while same(native,zombie_ok=True):
  p=proc(native['pid'])
  if p['state']=='Z':
   code=decode_exit(p['exit_wait_status']);break
  if not same(parent) or proc(parent['pid'])['state'] not in ['T','t']:raise ValueError('Old controller pause/identity lost')
  _,_,paused=refresh(native,paused);write(W/'COORD_STATUS.json',dict(status='ADOPTED_NATIVE_COMPUTING_NO_RESTART',native=proc(native['pid']),controller_paused=parent['pid'],time=time.time()));time.sleep(15)
 else:raise ValueError('Native disappeared without capture of real exit status')
 event=dict(native=p,exit_code=code,time=time.time(),returncode_source='Linux zombie /proc/stat exit_code; original parent held stopped')
 write(W/'ADOPTED_NATIVE_EXIT.json',event)
 if code:raise RuntimeError('Adopted native failed; old controller remains paused, no automatic replay')
 run=D/'primary_r0.1';source=Path(C['old_attempt'])/'outputs/primary_r0.1/compute-score.log';dest=run/'compute-score.log'
 if dest.exists():raise FileExistsError('Final compute log copy already exists')
 before=sha(source);shutil.copy2(source,dest)
 if sha(source)!=before or sha(dest)!=before:raise ValueError('Compute log was not stable after native exit')
 validate_score(run);validation=json.loads((run/'SCORE_VALIDATION.json').read_text())
 write(W/'FIRST_SCORE_ACCEPTED.json',dict(status='VALIDATED_EXISTING_NATIVE_SCORE_NO_RECOMPUTATION',event=event,validation=validation,original_log_sha256=before))
 # Scientific child has exited and full score/log/schema have passed before supervisor retirement.
 if not same(parent) or proc(parent['pid'])['state'] not in ['T','t']:raise ValueError('Old controller identity/state changed')
 send(parent,signal.SIGKILL)
 until=time.monotonic()+5
 while proc(parent['pid']) and proc(parent['pid'])['state']!='Z' and time.monotonic()<until:time.sleep(.02)
 if proc(parent['pid']) and proc(parent['pid'])['state']!='Z':raise RuntimeError('Old controller did not retire')
 write(W/'OLD_CONTROLLER_RETIRED.json',dict(pid=parent['pid'],start_ticks=parent['start_ticks'],signal='SIGKILL stopped supervisor only AFTER native exited and full score validated',time=time.time(),native_exit_code=code,native_not_killed=True))
 COMMANDS.append(dict(adopted_native_command=native['cmd'],pid=native['pid'],start_ticks=native['start_ticks'],exit_code=code,score_recomputed=False,completed=time.time()));write(D/'COMMANDS.json',COMMANDS)
 return run

def remainder(first):
 with open(C['project_lock'],'a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  for i,trait in enumerate(OLD['scenarios']):
   run=D/trait['name']
   if i:
    run.mkdir(exist_ok=False);scores=run/'scores';scores.mkdir()
    write(D/'CURRENT_STAGE.json',dict(trait=trait['name'],status='COMPUTE_SCORE_START',time=time.time()))
    cmd=[OLD['scdrs_executable'],'compute-score','--h5ad-file',OLD['derived'],'--h5ad-species','human','--gs-file',trait['geneset'],'--gs-species','human','--cov-file',OLD['covariates'],'--ctrl-match-opt','mean_var','--weight-opt','vs','--flag-filter-data','True','--flag-raw-count','True','--n-ctrl','1000','--flag-return-ctrl-raw-score','False','--flag-return-ctrl-norm-score','True','--out-folder',scores]
    call(cmd,run/'compute-score.log',True);validate_score(run)
   schema=json.loads((run/'SCORE_VALIDATION.json').read_text());files=list((run/'scores').glob('*.full_score.gz'))
   if len(files)!=1:raise ValueError('Expected one full score')
   for key in OLD['group_keys']:
    group=run/key;group.mkdir(exist_ok=False);write(D/'CURRENT_STAGE.json',dict(trait=trait['name'],status='PERFORM_DOWNSTREAM',group_key=key,time=time.time()))
    cmd=[OLD['scdrs_executable'],'perform-downstream','--h5ad-file',OLD['derived'],'--score-file',files[0],'--out-folder',group,'--group-analysis',key,'--flag-filter-data','True','--flag-raw-count','True']
    call(cmd,run/(key+'.log'),True)
    call(['/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python',W/'code/validate_native.py','group','--run',run,'--key',key],run/(key+'_validation.log'))
    schema[key]=json.loads((run/(key+'_VALIDATION.json')).read_text())
   schema['trait']=trait['name'];write(run/'SCHEMA_VALIDATION.json',schema);write(D/'CURRENT_STAGE.json',dict(trait=trait['name'],status='TRAIT_COMPLETE',time=time.time()))
  verify();files={str(p.relative_to(D)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in D.rglob('*') if p.is_file()};write(W/'OUTPUT_HASHES.json',files)
  write(W/'COMPLETE.json',dict(status='COMPLETE_NATIVE_3_SCENARIOS_ON_DATA',data_root=str(D),first_score_recomputed=False,derived_rebuilt=False,output_hash_manifest_sha256=sha(W/'OUTPUT_HASHES.json'),time=time.time()))

def main():
 with open(W/'.coordinator.lock','a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  if (W/'COORD_STARTED.json').exists() or (W/'COMPLETE.json').exists():raise FileExistsError('Coordinator already started')
  verify();migration=json.loads((W/'MIGRATION_COMPLETE.json').read_text())
  if Path(C['old_scores']).resolve()!=Path(C['data_scores']).resolve():raise ValueError('Data alias changed')
  guard=json.loads((W/'COORDINATOR_RECEIPT.json').read_text())['guard_to_replace']
  s=json.loads((W/'MOUNT_GUARD_STATUS.json').read_text())
  if s['status']!='WATCHING_DATA_SCORE' or not same(guard):raise ValueError('Standalone guard is not safely replaceable')
  b=mount_budget(C)
  write(W/'COORD_STARTED.json',dict(pid=os.getpid(),time=time.time(),guard_to_retire=guard,first_score_recomputed=False))
  send(guard,signal.SIGTERM)
  write(D/'CURRENT_STAGE.json',dict(trait='primary_r0.1',status='ADOPTED_EXISTING_COMPUTE_SCORE_PID',pid=C['native']['pid'],time=time.time()))
  first=adopt_first();remainder(first)
if __name__=='__main__':
 try:main()
 except Exception:
  traceback.print_exc()
  if (W/'COORD_STARTED.json').exists() and json.loads((W/'COORD_STARTED.json').read_text())['pid']==os.getpid() and same(C['native']) and proc(C['native']['pid'])['state'] not in ['T','t','Z']:
   with open(W/'emergency_guard.log','a') as f:
    g=subprocess.Popen(['/usr/bin/python3','-u',str(W/'code/mount_guard.py')],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
   write(W/'EMERGENCY_GUARD_PID.json',dict(pid=g.pid,time=time.time()))
  if (W/'COORD_STARTED.json').exists() and json.loads((W/'COORD_STARTED.json').read_text())['pid']==os.getpid():write(W/'COORD_EXIT.json',dict(status='FAILED_REVIEW_REQUIRED_NO_REPLAY',error=traceback.format_exc(),time=time.time()))
  raise
