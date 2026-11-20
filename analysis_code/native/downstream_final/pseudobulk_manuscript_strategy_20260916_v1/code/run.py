import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import subprocess,json,time,csv,sys,hashlib,fcntl
B=Path(__file__).resolve().parents[1];W=B.parent;sys.path.insert(0,str(W/'code'))
from common import processes
R='/usr/bin/Rscript';PY=sys.executable

def dump(name,x):
 p=B/name;t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2));t.replace(p)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def run(cmd,log):
 with (B/'logs'/log).open('x') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
def budget():
 ps=processes();n=sum(p['threads'] for p in ps)
 # Explicitly count this new branch if legacy resource detector does not yet match it.
 seen={p['pid'] for p in ps}
 for d in Path('/proc').glob('[0-9]*'):
  try:
   cmd=(d/'cmdline').read_bytes().replace(b'\0',b' ').decode()
   if str(B) in cmd and int(d.name) not in seen:
    n+=int(next(l.split(':')[1] for l in (d/'status').read_text().splitlines() if l.startswith('Threads:')))
  except (OSError,StopIteration,UnicodeDecodeError):pass
 return n

def queue(tasks,stage):
 todo=list(tasks);active=[];done=[]
 while todo or active:
  for q in active[:]:
   rc=q['p'].poll()
   if rc is None:continue
   q['log'].close();q['lock'].close();active.remove(q)
   if rc:raise RuntimeError(f"{stage} {q['id']} failed code {rc}")
   done.append(q['id'])
  n=budget()
  if todo and len(active)<3 and n+3<=40:
   lock=None
   for i in range(4):
    f=(W/f'locks/cpu_heavy_{i}.lock').open('a')
    try:fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB);lock=f;break
    except BlockingIOError:f.close()
   if lock:
    ident,cmd=todo.pop(0);f=(B/'logs'/f'{stage}_{ident}.log').open('x');p=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT)
    active.append(dict(id=ident,p=p,log=f,lock=lock));print(stage,'START',ident,p.pid,flush=True)
  dump('PROGRESS.json',dict(stage=stage,completed=len(done),expected=len(tasks),pending=len(todo),active=[dict(id=q['id'],pid=q['p'].pid) for q in active],project_threads=n,time=time.time()))
  if todo or active:time.sleep(2)
try:
 (B/'logs').mkdir(exist_ok=True);(B/'tests').mkdir(exist_ok=True)
 run([R,str(B/'code/native.R'),str(B),'test'],'native_test.log')
 assert json.loads((B/'tests/PASS.json').read_text())['status']=='PASS'
 run([PY,str(B/'code/prepare.py')],'prepare.log')
 receipt=json.loads((B/'RECEIPT.json').read_text())
 for p,h in receipt['bound_files'].items():assert sha(p)==h,p
 rows=list(csv.DictReader((B/'tables/model_roster.tsv').open(),delimiter='\t'));rows=[r for r in rows if r['status']=='READY']
 # Prioritize the two manuscript CKB populations, then the rest, without changing any test.
 rows.sort(key=lambda r:(r['cell_state'] not in ['Fibroblast','Macrophage'],r['level'],r['cell_state'],r['context']))
 queue([(r['model_id'],[R,str(B/'code/native.R'),str(B),'fit',r['model_id']]) for r in rows],'DE_and_LOS')
 (B/'meta').mkdir();unique=sorted({(r['level'],r['cell_state']) for r in rows},key=lambda x:(x[1] not in ['Fibroblast','Macrophage'],x))
 tasks=[]
 for level,state in unique:
  ident=level+'_'+hashlib.sha256(state.encode()).hexdigest()[:12];tasks.append((ident,[R,str(B/'code/native.R'),str(B),'meta',state,level,str(B/'meta'/f'{ident}.tsv')]))
 queue(tasks,'REML_target_set')
 run([PY,str(B/'code/report.py')],'report.log')
 print('COMPLETE',flush=True)
except Exception as e:
 dump('FAILED.json',dict(error=repr(e),time=time.time()));raise
