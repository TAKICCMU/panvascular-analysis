from pathlib import Path
import os,json,time,subprocess,fcntl,traceback
R=Path(__file__).resolve().parents[1];PY='/home/anzhen_vhd/miniconda3/envs/scib_GPU/bin/python'
lock=open(R/'CONTROLLER.lock','a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def status(stage,**kw):
 p=R/'STATUS.json';tmp=R/'STATUS.tmp';tmp.write_text(json.dumps({'stage':stage,'pid':os.getpid(),'updated_unix':time.time(),**kw},indent=2));tmp.replace(p)
env=os.environ.copy();env.pop('LD_LIBRARY_PATH',None);env.pop('CUDA_VISIBLE_DEVICES',None)
env.update(XLA_PYTHON_CLIENT_PREALLOCATE='false',XLA_PYTHON_CLIENT_MEM_FRACTION='.30',OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8',MKL_NUM_THREADS='8',NUMBA_NUM_THREADS='8')
cpuenv=env.copy();cpuenv.update(CUDA_VISIBLE_DEVICES='',JAX_PLATFORMS='cpu')
def run(script,args,env,name):
 with (R/f'logs/{name}.log').open('a') as log:
  subprocess.run([PY,'-u',str(R/'code'/script),*args],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
try:
 while True:
  s=json.loads((R/'STATUS.json').read_text())
  if s['stage']=='CPU_HARMONY_COMPLETE':break
  if s['stage']!='CPU_HARMONY':raise RuntimeError('Unexpected integration status '+str(s))
  if not Path('/proc/'+str(s['pid'])).exists():raise RuntimeError('Integration process stopped before completion')
  time.sleep(5)
 status('CPU_NEIGHBOR_PREPARATION',jobs=8,backend='PyNNDescent CPU',Leiden='not required; original KMeans metrics')
 run('evaluate_cpu.py',['--method','harmony_cpu','--seed','0','--jobs','8','--prepare-only'],cpuenv,'neighbors_cpu')
 # Independent CPU visualization tasks use exactly the trusted evaluation set.
 plots=[]
 for m in ['harmony_cpu','harmony_gpu']:
  log=(R/f'logs/visualize_{m}.log').open('a');p=subprocess.Popen([PY,'-u',str(R/'code/visualize.py'),m],env=cpuenv,stdout=log,stderr=subprocess.STDOUT);plots.append((p,log,m))
 status('SCIB_METRICS_WAITING_GPU_LOCK',neighbors_backend='CPU',UMAP_backend='CPU',Leiden='not run')
 gpu=open(R.parents[1]/'GPU_HEAVY_STAGE.lock','a');fcntl.flock(gpu,fcntl.LOCK_EX)
 status('SCIB_METRICS',integration_backend='CPU',neighbors_backend='CPU',metric_backend='historical JAX GPU',UMAP_backend='CPU')
 run('evaluate_cpu.py',['--method','harmony_cpu','--seed','0','--jobs','8'],env,'metrics_cpu_harmony')
 fcntl.flock(gpu,fcntl.LOCK_UN);gpu.close()
 run('evaluate_cpu.py',['--aggregate','--seed','0','--jobs','8'],cpuenv,'aggregate_cpu')
 status('WAITING_CPU_UMAP_REVIEW_FIGURES')
 for p,log,m in plots:
  rc=p.wait();log.close()
  if rc:raise RuntimeError(m+' UMAP failed '+str(rc))
 status('FINALIZING_REVIEW')
 run('finalize.py',[],cpuenv,'finalize')
except Exception:
 status('FAILED_REQUIRES_INSPECTION',error=traceback.format_exc());raise
