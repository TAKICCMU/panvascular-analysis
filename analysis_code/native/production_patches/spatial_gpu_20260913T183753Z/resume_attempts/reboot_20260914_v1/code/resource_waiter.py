"""One bounded-thread GPU queue. Waiting is not scientific completion."""
import os,sys,json,time,subprocess,shutil,fcntl,hashlib
from pathlib import Path
from resource_coordination import coordination,project_processes
import signal
JOB=Path(__file__).resolve().parents[1];cfg=json.loads((JOB/'configs/config.json').read_text())
lock=open(JOB/'GPU_JOB.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def status(state,**kw):
 d={'status':state,'time_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'waiter_pid':os.getpid(),'science_started':False,**kw};(JOB/'RUN_STATUS.json').write_text(json.dumps(d,indent=2));print(json.dumps(d),flush=True)
def resources():
 raw=subprocess.check_output(['nvidia-smi','--query-gpu=uuid,memory.free,utilization.gpu','--format=csv,noheader,nounits'],text=True).strip().splitlines()
 gpu=next(x.split(', ') for x in raw if x.startswith(cfg['gpu_uuid']))
 proc=subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,gpu_uuid','--format=csv,noheader,nounits'],text=True).strip().splitlines()
 foreign=[x for x in proc if cfg['gpu_uuid'] in x]
 free=shutil.disk_usage(JOB).free/(1024**3)
 return {'gpu_uuid':gpu[0],'gpu_free_MiB':int(gpu[1]),'gpu_util_pct':int(gpu[2]),'foreign_compute_processes':foreign,'disk_free_GiB':free,'coordination':coordination()}
# Configured 16GiB GPU reservation + 10GiB first CPU phase + 20GiB absolute floor.
stable=0
while True:
 try:r=resources();ready=not r['foreign_compute_processes'] and r['gpu_util_pct']<10 and r['gpu_free_MiB']>16384 and r['disk_free_GiB']>46 and r['coordination']['pass_gate']
 except Exception as e:status('WAITING_RESOURCE',reason=repr(e));time.sleep(60);continue
 if not ready:stable=0;status('WAITING_RESOURCE',reason='Need no existing GPU compute process, util<10%, freeGPU>16GiB and disk>46GiB',resource=r);time.sleep(60);continue
 stable+=1
 if stable<3:status('WAITING_RESOURCE',reason='Three consecutive idle observations required',resource=r,stable=stable);time.sleep(60);continue
 aux=open(cfg['shared_auxiliary_lock'],'a');fcntl.flock(aux,fcntl.LOCK_EX)
 # Recheck after waiting for the shared auxiliary lease.
 r=resources()
 if r['foreign_compute_processes'] or r['gpu_util_pct']>=10 or r['disk_free_GiB']<=46 or not r['coordination']['pass_gate']:
  fcntl.flock(aux,fcntl.LOCK_UN);aux.close();stable=0
  status('WAITING_RESOURCE',reason='Resource changed during auxiliary lease wait; returning to queue',resource=r);time.sleep(60);continue
 break
py=cfg['gpu_python'];env=dict(os.environ,CUDA_VISIBLE_DEVICES=cfg['gpu_uuid'],OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MKL_NUM_THREADS='4',NUMBA_NUM_THREADS='4',NUMEXPR_NUM_THREADS='4',MPLBACKEND='Agg',PYTHONUNBUFFERED='1')
for stage in ['smoke','production']:
 if stage=='production' and not (JOB/'smoke/PASS.json').exists():status('FAILED',reason='Native smoke receipt absent');sys.exit(1)
 status('RUNNING',stage=stage,science_started=(stage=='production'))
 cmd=[py,'-u','-B',str(JOB/'code/spatial_gpu.py'),stage]
 with open(JOB/f'logs/{stage}.log','a') as log:
  p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
  (JOB/f'{stage}_PID.json').write_text(json.dumps({'pid':p.pid,'command':cmd,'GPU_UUID':cfg['gpu_uuid'],'cpu_threads':4,'science_started':stage=='production'},indent=2))
  while p.poll() is None:
   live=project_processes();actual=sum(x['threads'] for x in live)
   (JOB/'NATIVE_CPU_STATUS.json').write_text(json.dumps({'actual_project_threads':actual,'limit':40,'stage_pid':p.pid,'processes':live},indent=2))
   if actual>40:
    os.killpg(p.pid,signal.SIGTERM)
    try:p.wait(timeout=10)
    except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
    status('FAILED_RESOURCE_LIMIT',stage=stage,actual_project_threads=actual,science_started=stage=='production');sys.exit(2)
   time.sleep(2)
  rc=p.returncode
 if rc:status('FAILED',stage=stage,exit_code=rc,science_started=stage=='production');sys.exit(rc)
final=json.loads((JOB/'RUN_STATUS.json').read_text())
if final.get('status')!='COMPLETE':
 status('PARTIAL_FAILED',science_started=True,detail=final);sys.exit(1)
status('COMPLETE',science_started=True,scientific_review_required=True)
