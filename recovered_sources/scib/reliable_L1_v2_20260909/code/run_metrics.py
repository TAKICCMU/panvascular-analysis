from pathlib import Path
import os, subprocess, json, time, fcntl
O=Path(__file__).resolve().parents[1]
env=os.environ.copy();env.pop('LD_LIBRARY_PATH',None)
env.update(XLA_PYTHON_CLIENT_PREALLOCATE='false',XLA_PYTHON_CLIENT_MEM_FRACTION='.30',OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8',NUMBA_NUM_THREADS='8')
py='/home/anzhen_vhd/miniconda3/envs/scib_GPU/bin/python'
lock=open(O/'metrics/CONTROLLER.lock','a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def status(stage,**extra):(O/'metrics/STATUS.json').write_text(json.dumps(dict(stage=stage,pid=os.getpid(),updated_unix=time.time(),**extra),indent=2))
status('waiting_supervised_embeddings')
while not all((O/f'embeddings/{m}.json').exists() for m in ['scanvi','scpoli']):time.sleep(15)
gpu=open(O.parents[1]/'GPU_HEAVY_STAGE.lock','a')
# O = takeover/scib/version; O.parents[1] is takeover.
fcntl.flock(gpu,fcntl.LOCK_EX)
for method in ['pca','harmony','scvi','scanorama','liger','scalex','scanvi','scpoli']:
    status('official_metrics',method=method)
    with (O/f'logs/metrics_{method}.log').open('a') as f:
        rc=subprocess.run([py,str(O/'code/evaluate_official.py'),'--method',method,'--seed','0','--jobs','8'],env=env,stdout=f,stderr=subprocess.STDOUT).returncode
    if rc:status('failed_requires_review',method=method,exit_code=rc);raise SystemExit(rc)
status('aggregating')
rc=subprocess.run([py,str(O/'code/evaluate_official.py'),'--aggregate','--seed','0','--jobs','8'],env=env).returncode
status('complete_requires_review' if rc==0 else 'failed_aggregation',exit_code=rc)
raise SystemExit(rc)
