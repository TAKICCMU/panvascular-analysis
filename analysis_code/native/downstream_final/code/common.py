import os
THREAD_KEYS=['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS','POLARS_MAX_THREADS','RAYON_NUM_THREADS','ARROW_NUM_THREADS','BLOSC_NTHREADS','VECLIB_MAXIMUM_THREADS']
for k in THREAD_KEYS:os.environ[k]='1'
os.environ['MALLOC_CONF']='background_thread:false'
from pathlib import Path
import json,hashlib,time,subprocess,fcntl,signal
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908')
D=Path('/data/anzhen_vhd/scRNAseq/codex/panvascular_downstream_storage_v1')
W=Path(__file__).resolve().parents[1]
C=R/'downstream_CPU_L1L2_v2_1_release_compat_20260914'
AUDIT=D/'covarnet_pseudobulk_launch_review_20260915/updated_scope_audit_v2'
PY='/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python'
RS='/home/anzhen_vhd/miniconda3/envs/communication_r/bin/Rscript'
ATLAS=R/'panvascular_CPU_theta2_L1L2_final_v1.h5ad'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def dump(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);t=p.with_name(p.name+'.tmp.'+str(os.getpid()));t.write_text(json.dumps(x,indent=2,ensure_ascii=False,default=str));t.replace(p)
def processes():
 out=[]
 for p in Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   cmd=(p/'cmdline').read_bytes().replace(b'\0',b' ').decode()
   if not any(x in cmd for x in [str(R),str(D)]):continue
   if any(x in cmd for x in ['python3 -','sshd:','/bin/sh -c','/bin/bash -c','scp -']):continue
   z={a.split(':',1)[0]:a.split(':',1)[1].strip() for a in (p/'status').read_text().splitlines() if ':' in a}
   out.append(dict(pid=int(p.name),ppid=int(z['PPid']),threads=int(z['Threads']),state=z['State'],cmd=cmd))
  except (OSError,UnicodeError,KeyError):pass
 return out
def verify():
 rec=json.loads((W/'RECEIPT.json').read_text())
 for p,h in rec['bound_files'].items():
  if sha(p)!=h:raise ValueError('Bound content changed '+p)
 st=ATLAS.stat()
 assert [st.st_dev,st.st_ino,st.st_size,st.st_mtime_ns]==rec['atlas_stat']
 assert rec['scope_authorized_by_user']=='2026-09-16 expanded downstream execution request'
 return rec
def run_job(name,cmd,log,reserve=4,env=None):
 verify();handles=[]
 (W/'locks').mkdir(exist_ok=True)
 while True:
  active=processes();n=sum(x['threads'] for x in active)
  for i in range(4):
   h=open(W/'locks'/('cpu_heavy_'+str(i)+'.lock'),'a')
   try:
    fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if n+reserve>40:h.close();break
    handles=[h];break
   except BlockingIOError:h.close()
  if handles:break
  dump(W/'status'/(name+'.json'),dict(status='WAITING_SHARED_40_THREAD_OR_4_JOB_BUDGET',actual_threads=n,reservation=reserve,time=time.time()))
  time.sleep(20)
 e=os.environ.copy();e.update(env or {});e['CUDA_VISIBLE_DEVICES']=''
 with open(log,'x') as f:
  p=subprocess.Popen([str(x) for x in cmd],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,env=e,start_new_session=True)
  dump(W/'status'/(name+'.json'),dict(status='RUNNING',pid=p.pid,command=cmd,log=str(log),started=time.time()))
  code=p.wait()
 dump(W/'status'/(name+'.json'),dict(status='COMPLETE' if code==0 else 'FAILED',pid=p.pid,command=cmd,exit_code=code,log=str(log),ended=time.time()))
 for h in handles:h.close()
 if code:raise RuntimeError(name+' failed: '+str(log))
