import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS','POLARS_MAX_THREADS','RAYON_NUM_THREADS','ARROW_NUM_THREADS','BLOSC_NTHREADS']:
    os.environ[k]='1'
os.environ['CUDA_VISIBLE_DEVICES']=''
from pathlib import Path
import sys,json,hashlib,subprocess,time,signal
B=Path(__file__).resolve().parent
N=B.parent;P=N/'communication';L=P/'L2'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):
    t=p.with_suffix('.tmp');t.write_text(json.dumps(x,indent=2));t.replace(p)
def check(row):
    out=L/'outputs'/row['run_id'];z=json.loads((out/'COMPLETE.json').read_text())
    assert z['status']=='COMPLETE_NATIVE_LIANA_CELLCHAT' and z['sample']==row['sample_uid']
    for k,v in dict(min_cells=50,max_cells_per_type=200,n_perms=1000,nboot=100,liana_seed=1337,threads_CellChat=1,liana_n_jobs=1).items():assert z[k]==v,(row['run_id'],k)
    assert z['sample_seed']==int(row['seed']) and not z['raw_source_modified']
    return dict(run_id=row['run_id'],files={str(p):sha(p) for p in out.iterdir() if p.is_file()})
def verify():
    r=json.loads((B/'RECEIPT.json').read_text())
    for p,h in r['bound_files'].items():assert sha(p)==h,p
    sys.path.insert(0,str(P/'code'))
    from verify import verify as original_verify
    original_verify()
    return json.loads((B/'PLAN.json').read_text())
def identity(pid,start):
    try:return Path('/proc/'+str(pid)+'/stat').read_text().split()[21]==str(start)
    except FileNotFoundError:return False
if __name__=='__main__':
    plan=verify();mode=sys.argv[1]
    if mode=='lane':
        lane=int(sys.argv[2]);rows=plan['lanes'][lane];done=[]
        for row in rows:
            out=L/'outputs'/row['run_id'];scratch=L/'scratch'/row['run_id']
            if not (out/'COMPLETE.json').exists():
                assert not out.exists() and not scratch.exists(),'Unexpected existing partial run '+row['run_id']
                affinity=plan['cpu_pairs'][lane];assert set(affinity)<=os.sched_getaffinity(0),'CPU affinity changed; review required'
                cmd=plan['native_command_prefix']+['--run-id',row['run_id'],'--scratch',str(scratch),'--out',str(out),'--threads','1','--receipt',str(P/'RECEIPT.json')]
                with (L/'logs'/(row['run_id']+'.log')).open('x') as log:
                    child=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,preexec_fn=lambda:os.sched_setaffinity(0,affinity))
                    dump(B/('lane_'+str(lane)+'_ACTIVE.json'),dict(pid=child.pid,run_id=row['run_id'],command=cmd,cpu_pair=affinity))
                    assert child.wait()==0,'Native worker failed; inspect sample log'
            done.append(check(row));dump(B/('lane_'+str(lane)+'_PROGRESS.json'),dict(completed=len(done),expected=len(rows),last=row['run_id'],updated=time.time()))
        dump(B/('lane_'+str(lane)+'_COMPLETE.json'),dict(status='COMPLETE_VALIDATED',results=done))
    elif mode=='resume':
        for i in range(3):
            c=json.loads((B/('lane_'+str(i)+'_COMPLETE.json')).read_text());assert c['status']=='COMPLETE_VALIDATED'
            for r in c['results']:
                for p,h in r['files'].items():assert sha(p)==h,p
        current=plan['current_sample'];cp=L/'outputs'/current['run_id']/'COMPLETE.json'
        while not cp.exists():
            assert identity(plan['current_worker_pid'],plan['current_worker_start']),'Original worker disappeared without COMPLETE'
            time.sleep(10)
        for row in plan['all_samples']:check(row)
        pid=plan['original_controller_pid'];assert identity(pid,plan['original_controller_start'])
        assert Path('/proc/'+str(pid)+'/status').read_text().split('State:')[1].lstrip().startswith('T')
        os.kill(pid,signal.SIGCONT)
        dump(B/'RESUMED_ORIGINAL_CONTROLLER.json',dict(status='ALL_SAMPLES_VALIDATED_ORIGINAL_CONTROLLER_RESUMED',samples=len(plan['all_samples']),pid=pid,time=time.time()))
    else:raise ValueError(mode)
