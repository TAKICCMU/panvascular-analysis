"""Import-safe stdlib identity, mount and resource checks for lossless scDRS storage handoff."""
from pathlib import Path
import os,json,hashlib,shutil,time,signal,ctypes

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def write(p,d):
 p=Path(p);q=p.with_name(p.name+'.partial');q.write_text(json.dumps(d,indent=2,default=str));q.replace(p)
def proc(pid):
 p=Path('/proc')/str(pid)
 try:
  s=(p/'stat').read_text().rsplit(')',1)[1].split();lines=dict(x.split(':',1) for x in (p/'status').read_text().splitlines() if ':' in x)
  return dict(pid=pid,start_ticks=s[19],state=s[0],ppid=int(s[1]),threads=int(lines['Threads']),cmd=(p/'cmdline').read_bytes().replace(b'\0',b' ').decode(),exit_wait_status=int(s[49]) if len(s)>49 else None)
 except (OSError,ValueError,KeyError):return None
def same(saved,zombie_ok=False):
 p=proc(saved['pid'])
 return bool(p and p['start_ticks']==saved['start_ticks'] and ((zombie_ok and p['state']=='Z') or p['cmd']==saved['cmd']))
def send(saved,sig):
 if not same(saved):raise ValueError('Process identity changed; refuse signal')
 os.kill(saved['pid'],sig)
def wait_state(saved,states,seconds=5):
 t=time.monotonic()
 while time.monotonic()-t<seconds:
  p=proc(saved['pid'])
  if p and p['start_ticks']==saved['start_ticks'] and p['state'] in states:return p
  time.sleep(.01)
 raise RuntimeError('Expected process state not reached')
def fds(saved):
 if not same(saved):raise ValueError('PID changed')
 rows=[];p=Path('/proc')/str(saved['pid'])
 for f in (p/'fd').iterdir():
  try:
   info=dict(x.split(':',1) for x in (p/'fdinfo'/f.name).read_text().splitlines() if ':' in x)
   rows.append(dict(fd=int(f.name),target=os.readlink(f),flags=int(info['flags'].strip(),8)))
  except OSError:continue
 return rows
def assert_no_score_writer(records,scores):
 prefix=str(scores)+'/'
 if any((x['target']==str(scores) or x['target'].startswith(prefix)) and (x['flags']&os.O_ACCMODE)!=os.O_RDONLY for x in records):raise ValueError('Open score writer; migration refused')
def exchange(a,b):
 libc=ctypes.CDLL(None,use_errno=True);f=libc.renameat2;f.argtypes=[ctypes.c_int,ctypes.c_char_p,ctypes.c_int,ctypes.c_char_p,ctypes.c_uint];f.restype=ctypes.c_int
 if f(-100,os.fsencode(a),-100,os.fsencode(b),2):raise OSError(ctypes.get_errno(),os.strerror(ctypes.get_errno()))
def mount_budget(c):
 if sha(c['ledger'])!=c['ledger_sha256']:raise ValueError('Mount resource ledger changed')
 ledger=json.loads(Path(c['ledger']).read_text());out={};reasons=[]
 for name,m in ledger['filesystems'].items():
  path=Path(m['mount_path'])
  if not path.exists() or path.stat().st_dev!=m['device']:raise ValueError('Absent/changed mount '+name)
  used={};remaining={};over=[]
  for key,a in m['reservations'].items():
   n=reservation_bytes(a['paths'],m['device'])
   growth=max(0,n-a['baseline_bytes'])/1024**3
   used[key]=n;remaining[key]=max(0,a['allocation_gib']-growth)
   if growth>a['allocation_gib']:over.append(key)
  free=shutil.disk_usage(path).free/1024**3;needed=m['floor_gib']+sum(remaining.values())
  out[name]=dict(device=m['device'],free_GiB=free,required_free_GiB=needed,floor_GiB=m['floor_gib'],remaining_GiB=remaining,used_bytes=used,over_allocation=over,pass_gate=free>=needed and not over)
  if not out[name]['pass_gate']:reasons.append(name+'_future_disk_budget')
 return dict(filesystems=out,reasons=reasons)
def cpu_budget(c,own):
 ledger=json.loads(Path(c['ledger']).read_text());ps={}
 for p in Path('/proc').iterdir():
  if p.name.isdigit():
   v=proc(int(p.name))
   if v:ps[v['pid']]=v
 targets={p for p,v in ps.items() if any(root in v['cmd'] for root in ledger['project_roots'])};owners=set(own)
 for _ in range(30):
  prev=(len(targets),len(owners));targets|={p for p,v in ps.items() if v['ppid'] in targets};owners|={p for p,v in ps.items() if v['ppid'] in owners}
  if prev==(len(targets),len(owners)):break
 other=[ps[p] for p in targets-owners if ps[p]['state'] not in ['T','t','Z']];actual=sum(p['threads'] for p in other)
 comm=sum(p['threads'] for p in other if c['communication_root'] in p['cmd']);gpu=sum(p['threads'] for p in other if c['gpu_root'] in p['cmd'])
 projected=actual+14+max(0,10-comm)+max(0,4-gpu)
 return dict(actual_other_threads=actual,projected=projected,pass_gate=projected<=ledger['cpu_total_max'],others=other)

def decode_exit(status):
 return (status >> 8) & 255 if (status & 127)==0 else -(status & 127)

def reservation_bytes(roots,device):
 # Snapshot each file stat once. Atomic status-file replacement may remove
 # a listed entry before stat; only ENOENT is harmless. Other errors propagate.
 import stat
 seen={}
 def onerror(error):
  if not isinstance(error,FileNotFoundError):raise error
 for root in roots:
  for directory,dirs,files in os.walk(root,onerror=onerror):
   for filename in files:
    p=Path(directory)/filename
    try:
     v=p.stat()
     identity=str(p.resolve())
    except FileNotFoundError:
     continue
    if not stat.S_ISREG(v.st_mode):continue
    if v.st_dev!=device:raise ValueError('Reservation crosses mount '+str(p))
    seen[identity]=v.st_size
 return sum(seen.values())
