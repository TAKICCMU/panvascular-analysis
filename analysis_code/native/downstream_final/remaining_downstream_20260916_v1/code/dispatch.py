import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS','POLARS_MAX_THREADS','RAYON_NUM_THREADS','ARROW_NUM_THREADS','BLOSC_NTHREADS']:
    os.environ[k]='1'
os.environ['CUDA_VISIBLE_DEVICES']=''
import json,subprocess,time,hashlib,fcntl
from pathlib import Path
N=Path(__file__).resolve().parents[1]
W=N.parent
def dump(p,x):
    t=p.with_suffix('.tmp');t.write_text(json.dumps(x,indent=2,ensure_ascii=False));t.replace(p)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def actual_threads():
    total=0
    for p in Path('/proc').iterdir():
        if not p.name.isdigit():continue
        try:
            c=(p/'cmdline').read_bytes().replace(b'\0',b' ').decode()
            if str(W.parent) not in c:continue
            if any(x in c for x in ['python3 -','bash -c','sh -c','sshd:','scp ']):continue
            s=dict(x.split(':',1) for x in (p/'status').read_text().splitlines() if ':' in x)
            if s['State'].strip().startswith(('T','Z')):continue
            total+=int(s['Threads'])
        except (OSError,UnicodeError,KeyError):pass
    return total
(N/'logs').mkdir(exist_ok=True)
lock=open(N/'dispatcher.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
state_path=N/'QUEUE_STATUS.json'
states=json.loads(state_path.read_text())['jobs'] if state_path.exists() else {}
if any(x.get('status')=='RUNNING' for x in states.values()):raise RuntimeError('Previous RUNNING jobs require process verification before restart')
active={}
while True:
    config=json.loads((N/'QUEUE.json').read_text())
    for jid,(p,log,h) in list(active.items()):
        rc=p.poll()
        if rc is not None:
            states[jid].update(status='PROCESS_EXITED_OK_PENDING_RESULT_REVIEW' if rc==0 else 'FAILED',exit_code=rc,ended=time.time())
            log.close();h.close();del active[jid]
    for j in config['jobs']:
        jid=j['id']
        if jid in states and states[jid]['status'] not in ['QUEUED','WAITING_DEPENDENCY','WAITING_RESOURCES']:continue
        if any(states.get(d,{}).get('status')!='PROCESS_EXITED_OK_PENDING_RESULT_REVIEW' for d in j.get('depends_on',[])):
            states[jid]={'status':'WAITING_DEPENDENCY','dependencies':j['depends_on']};continue
        n=actual_threads()
        if len(active)>=4 or n+j.get('threads_reserved',6)>40:
            states[jid]={'status':'WAITING_RESOURCES','actual_project_threads':n};continue
        handle=None
        for i in range(4):
            h=open(W/'locks'/('cpu_heavy_'+str(i)+'.lock'),'a')
            try:fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB);handle=h;break
            except BlockingIOError:h.close()
        if handle is None:continue
        for f,hsh in j.get('bound_files',{}).items():
            if sha(f)!=hsh:raise RuntimeError('Changed launch input '+f)
        logpath=N/'logs'/(jid+'.log')
        log=logpath.open('x')
        p=subprocess.Popen(j['command'],stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True,env=os.environ.copy())
        active[jid]=(p,log,handle)
        states[jid]={'status':'RUNNING','pid':p.pid,'command':j['command'],'log':str(logpath),'started':time.time(),'threads_reserved':j.get('threads_reserved',6)}
        print('START',jid,p.pid,flush=True)
    dump(state_path,{'updated':time.time(),'jobs':states,'active_heavy_jobs':len(active),'actual_project_threads':actual_threads(),'max_heavy_jobs':4,'max_project_threads':40})
    if config.get('sealed') and not active and all(states.get(j['id'],{}).get('status') in ['PROCESS_EXITED_OK_PENDING_RESULT_REVIEW','FAILED'] for j in config['jobs']):break
    time.sleep(10)
