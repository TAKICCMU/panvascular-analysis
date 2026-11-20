#!/usr/bin/env python3
"""Scoring-only reboot attempt. Reuses certified derived input; old outputs remain immutable."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ['CUDA_VISIBLE_DEVICES']=''
import json,hashlib,fcntl,time,sys,shutil,traceback,argparse,importlib.util
from pathlib import Path
W=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def write(p,d):
 p=Path(p);q=p.with_name(p.name+'.partial');q.write_text(json.dumps(d,indent=2,default=str));q.replace(p)
def check_stat(record,path,expected):
 if str(record['path'])!=str(path) or record.get('sha256',record.get('expected_sha256'))!=expected:raise ValueError('Wrong certified file/hash')
 if record.get('full_sha_verified',record.get('full_sha256_verified')) is not True:raise ValueError('Certificate lacks full hash verification')
 if 'observed_sha256' in record and record['observed_sha256']!=expected:raise ValueError('Observed complete SHA differs')
 if 'boot_id' in record and record['boot_id']!=Path('/proc/sys/kernel/random/boot_id').read_text().strip():raise ValueError('Certificate belongs to a different boot')
 saved=record.get('stat',record.get('stat_after',record));s=Path(path).stat()
 for k,v in [('size',s.st_size),('mtime_ns',s.st_mtime_ns),('inode',s.st_ino),('device',s.st_dev),('nlink',s.st_nlink)]:
  if saved.get(k)!=v:raise ValueError('Reboot certificate stat changed '+k)
 return True

def validate_parameters(c):
 old=json.loads(Path(c['original_config']).read_text())
 if c['scenarios']!=old['scenarios'] or c['native_threads']!=12 or c['n_ctrl']!=1000 or c['graph_policy']!='frozen_global_n30':raise ValueError('Scientific parameters changed')
 if c['scdrs_executable']!=old['scdrs_executable'] or c['rebuild_derived'] is not False:raise ValueError('Forbidden native executable or derived rebuild')
 if c['group_keys']!=['level1_annotation_scdrs_test','level2_annotation_scdrs_test']:raise ValueError('Group keys changed')
 if Path(c['output']).exists() or (W/'STARTED.json').exists() or (W/'EXIT_STATUS.json').exists():raise FileExistsError('Existing attempt output/start/exit: no overwrite or duplicate restart')

def verify(c,receipt=None,check_new=True):
 for x in c['bound_files']:
  if sha(x['path'])!=x['sha256']:raise ValueError('Bound file changed '+x['path'])
 for entry in c['authenticated_files']:
  cert=json.loads(Path(entry['certificate']).read_text());record=cert
  if cert.get('full_sha_verified') is not True or cert.get('observed_sha256')!=entry['sha256'] or cert.get('boot_id')!=Path('/proc/sys/kernel/random/boot_id').read_text().strip():raise ValueError('Root complete SHA/boot certificate not verified')
  for key in entry['selector']:record=record[key]
  check_stat(record,entry['path'],entry['sha256'])
 old=json.loads(Path(c['original_receipt']).read_text())
 for x in old['bound_files']:
  if sha(x['path'])!=x['sha256']:raise ValueError('Historical method/input changed '+x['path'])
 if sha(c['covariates'])!=c['covariates_sha256']:raise ValueError('Derived covariate hash changed')
 derived=json.loads(Path(c['derived_report']).read_text())
 if derived['derived_h5ad_sha256']!=c['derived_sha256'] or derived['covariates_sha256']!=c['covariates_sha256'] or derived['source_symbol_only_features']!=735 or derived['group_keys']!=c['group_keys']:raise ValueError('Derived semantics changed')
 for x in c['scenarios']:
  if sha(x['geneset'])!=x['sha256']:raise ValueError('Gene-set changed')
 if receipt:
  r=json.loads(Path(receipt).read_text())
  if r['mode']!='approved_execution' or r['allowed_modules']!=['scdrs_scoring_and_native_downstream'] or r['config_sha256']!=sha(W/'config.json'):raise ValueError('Receipt mismatch')
  for x in r['bound_files']:
   if sha(x['path'])!=x['sha256']:raise ValueError('Receipt file changed '+x['path'])
 if check_new:validate_parameters(c)

def process_budget(c,procs,own):
 others=[p for p in procs if p['pid'] not in own and p['state'] not in ['T','t']]
 actual=sum(p['threads'] for p in others)
 comm=sum(p['threads'] for p in others if c['communication_root'] in p['cmd'])
 gpu=sum(p['threads'] for p in others if c['gpu_root'] in p['cmd'])
 reserved_comm=max(0,10-comm);reserved_gpu=max(0,4-gpu)
 projected=actual+14+reserved_comm+reserved_gpu
 return dict(actual_other_threads=actual,scdrs_reserved=14,remaining_communication_threads=reserved_comm,remaining_gpu_threads=reserved_gpu,projected=projected,pass_gate=projected<=40,others=others)

def resource_snapshot(c):
 ps={}
 for path in Path('/proc').iterdir():
  if not path.name.isdigit():continue
  try:
   fields=dict(x.split(':',1) for x in (path/'status').read_text().splitlines() if ':' in x)
   ps[int(path.name)]=dict(pid=int(path.name),ppid=int(fields['PPid']),threads=int(fields['Threads']),state=fields['State'].strip()[0],cmd=(path/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace'))
  except (OSError,ValueError,KeyError):pass
 targets={p for p,x in ps.items() if c['root'] in x['cmd']};own={os.getpid()}
 for _ in range(30):
  old=(len(targets),len(own));targets|={p for p,x in ps.items() if x['ppid'] in targets};own|={p for p,x in ps.items() if x['ppid'] in own}
  if old==(len(targets),len(own)):break
 cpu=process_budget(c,[ps[p] for p in targets],own)
 if sha(c['resource_ledger'])!=c['resource_ledger_sha256']:raise ValueError('Shared resume resource ledger changed')
 ledger=json.loads(Path(c['resource_ledger']).read_text())
 remaining={};usage={};over=[]
 for name,a in ledger['disk_reservations'].items():
  paths=set()
  for root in a['paths']:paths|={p.resolve() for p in Path(root).rglob('*') if p.is_file()}
  current=sum(p.stat().st_size for p in paths);growth=max(0,current-a['baseline_bytes'])/1024**3
  usage[name]=current;remaining[name]=max(0,a['allocation_gib']-growth)
  if growth>a['allocation_gib']:over.append(name)
 free=shutil.disk_usage(c['root']).free/1024**3;needed=ledger['unreserved_floor_gib']+sum(remaining.values())
 mem={x.split(':')[0]:int(x.split()[1]) for x in Path('/proc/meminfo').read_text().splitlines() if ':' in x};avail=mem['MemAvailable']/1024**2
 reasons=[]
 if not cpu['pass_gate']:reasons.append('project_CPU40_future_communication10_GPU4')
 if free<needed or over:reasons.append('future_disk_reservations_plus_20GiB_floor')
 if avail<220:reasons.append('MemAvailable_below220GiB')
 return dict(status='WAITING_RESOURCE' if reasons else 'RESOURCE_GATE_PASSED',time=time.time(),reasons=reasons,cpu=cpu,free_GiB=free,required_free_GiB=needed,remaining_reservations_GiB=remaining,current_branch_bytes=usage,over_allocation=over,mem_available_GiB=avail)

def wait_resources(c):
 while True:
  s=resource_snapshot(c);write(W/'RESOURCE_STATUS.json',s)
  if not s['reasons']:return
  time.sleep(45)

def main():
 p=argparse.ArgumentParser();p.add_argument('--receipt',required=True);p.add_argument('--execute',action='store_true');a=p.parse_args()
 c=json.loads((W/'config.json').read_text());verify(c,a.receipt)
 if not a.execute:print('PASS scoring-only reboot contract');return
 with open(c['lock'],'a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  wait_resources(c);verify(c,a.receipt)
  write(W/'STARTED.json',dict(pid=os.getpid(),time=time.time(),reused_derived=c['derived'],derived_sha256=c['derived_sha256'],rebuild_derived=False,old_native_outputs_preserved=c['old_native_outputs']))
  # Import same immutable scientific adapter; call only native_sequence, never build_derived.
  spec=importlib.util.spec_from_file_location('unchanged_scdrs_adapter',c['adapter']);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
  write(W/'STAGE_STATUS.json',dict(status='SCORING_RESTART_FROM_CERTIFIED_DERIVED',pid=os.getpid(),time=time.time()))
  module.native_sequence(c['derived'],c['covariates'],c['scenarios'],c['output'],c['scdrs_executable'],c['group_keys'],12,n_ctrl=1000,resource_check=lambda:wait_resources(c))
  verify(c,a.receipt,check_new=False)
  outputs={str(p.relative_to(W)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in Path(c['output']).rglob('*') if p.is_file()}
  write(W/'OUTPUT_HASHES.json',outputs);write(W/'EXIT_STATUS.json',dict(exit_code=0,status='COMPLETE_NATIVE_3_SCENARIOS',time=time.time(),rebuild_derived=False,output_hash_manifest_sha256=sha(W/'OUTPUT_HASHES.json')))
if __name__=='__main__':
 try:main()
 except Exception:
  traceback.print_exc()
  if '--execute' in sys.argv and (W/'STARTED.json').exists() and json.loads((W/'STARTED.json').read_text()).get('pid')==os.getpid():write(W/'EXIT_STATUS.json',dict(exit_code=1,status='FAILED',time=time.time(),error=traceback.format_exc()))
  raise
