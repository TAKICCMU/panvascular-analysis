"""Single-thread control process; an explicit sealed receipt is required for full execution."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[key]='1'
os.environ['CUDA_VISIBLE_DEVICES']='';os.environ['PYTHONDONTWRITEBYTECODE']='1'
import json,hashlib,argparse,time,shutil,fcntl,traceback,sys,datetime
from pathlib import Path
W=Path(__file__).resolve().parents[1];R=W.parents[1]
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def write(p,d):
 tmp=Path(str(p)+'.tmp');tmp.write_text(json.dumps(d,indent=2));tmp.replace(p)
def verify(receipt):
 d=json.loads(Path(receipt).read_text())
 assert d['mode']=='approved_execution' and d['allowed_modules']==['scdrs'] and d['approval_source']=='user_authorized_traits_repair_20260914'
 for entry in d['bound_files']:
  if sha(entry['path'])!=entry['sha256']:raise ValueError('Bound file changed: '+entry['path'])
 cfg=json.loads(Path(d['config']).read_text());assert [x['name'] for x in cfg['scenarios']]==['primary_r0.1','primary_r0.6','leave_one_CAD_r0.1'];assert cfg['graph_policy']=='frozen_global_n30' and cfg['n_ctrl']==1000 and cfg['native_threads']==12
 cert=json.loads(Path(cfg['atlas_stat']).read_text());s=Path(cert['path']).stat()
 for key,value in [('size',s.st_size),('mtime_ns',s.st_mtime_ns),('inode',s.st_ino),('device',s.st_dev),('nlink',s.st_nlink)]:
  if cert[key]!=value:raise ValueError('Certified atlas stat changed: '+key)
 assert cert['expected_sha256']==cfg['atlas_sha256']=='00fc74c1ddaa5e4b32d8664f134e0ec5a0c762f5e9a9fb45b96acbd40804a310'
 manifest=json.loads(Path(cfg['manifest']).read_text());assert manifest['release_status']=='FINAL_APPROVED' and manifest['ready_for_downstream'] and manifest['level1_frozen'] and manifest['level2_frozen']
 assert manifest['h5ad_sha256']==cfg['atlas_sha256'] and manifest['shape_verified']==[1175115,37336]
 for name,digest in manifest['release_sidecar_sha256'].items():
  if sha(Path(manifest['release_support_dir'])/name)!=digest:raise ValueError('Release sidecar changed: '+name)
 smoke=json.loads(Path(cfg['smoke_validation']).read_text());assert smoke['status']=='PASS' and smoke['shared_adapter_sha256']==sha(W/'code/scdrs_adapter.py')
 return cfg

def directory_bytes(paths):
 return sum(p.stat().st_size for root in paths for p in Path(root).rglob('*') if p.is_file())
def remaining_disk(free,allocations):
 remaining={};over=[]
 for name,a in allocations.items():
  growth=max(0,a['current_bytes']-a['baseline_bytes'])/1024**3
  remaining[name]=max(0,a['allocation_gib']-growth)
  if growth>a['allocation_gib']:over.append(name)
 required=20+sum(remaining.values())
 return {'required_free_gib':required,'remaining_reservation_gib':remaining,'over_allocation':over,'disk_pass':free>=required and not over}
def resource_snapshot(cfg):
 procs={}
 for path in Path('/proc').iterdir():
  if not path.name.isdigit():continue
  try:
   cmd=(path/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace');lines=(path/'status').read_text().splitlines();fields={x.split(':',1)[0]:x.split(':',1)[1].strip() for x in lines if ':' in x}
   procs[int(path.name)]={'pid':int(path.name),'ppid':int(fields['PPid']),'threads':int(fields['Threads']),'cmd':cmd}
  except (OSError,ValueError,KeyError):pass
 targets={pid for pid,p in procs.items() if str(R) in p['cmd']};own={os.getpid()}
 changed=True
 while changed:
  extra={pid for pid,p in procs.items() if p['ppid'] in targets};extraown={pid for pid,p in procs.items() if p['ppid'] in own};changed=not extra<=targets or not extraown<=own;targets|=extra;own|=extraown
 others=[procs[p] for p in sorted(targets-own)];threads=sum(p['threads'] for p in others)
 free=shutil.disk_usage(R).free/1024**3;mem={x.split(':')[0]:int(x.split()[1]) for x in Path('/proc/meminfo').read_text().splitlines() if ':' in x};avail=mem['MemAvailable']/1024**2
 used=sum(p.stat().st_size for p in (W/'outputs').rglob('*') if p.is_file())/1024**3 if (W/'outputs').exists() else 0
 reasons=[]
 # Conservative native library may have 12 workers plus main/control; GPU4 future reserved.
 projected=threads+14+4
 if projected>40:reasons.append('live_project_threads_plus_native14_plus_GPU4_exceed40')
 allocations={name:dict(entry,current_bytes=directory_bytes(entry['paths'])) for name,entry in cfg['disk_reservations'].items()}
 disk=remaining_disk(free,allocations)
 if not disk['disk_pass']:reasons.append('remaining_disk_reservation_unmet_or_branch_over_budget')
 if avail<220:reasons.append('MemAvailable_below220GiB')
 return {'status':'WAITING_RESOURCE' if reasons else 'RESOURCE_GATE_PASSED','timestamp':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reasons':reasons,'existing_project_threads':threads,'projected_conservative_threads':projected,'other_project_processes':others,'free_GiB':free,'own_output_GiB':used,'mem_available_GiB':avail,'disk_accounting':disk,'disk_allocations':allocations}

def wait_resources(cfg):
 while True:
  d=resource_snapshot(cfg);write(W/'RESOURCE_STATUS.json',d)
  if not d['reasons']:return
  time.sleep(60)

def main():
 p=argparse.ArgumentParser();p.add_argument('--receipt',required=True);p.add_argument('--execute',action='store_true');args=p.parse_args();cfg=verify(args.receipt)
 if not args.execute:print('PASS bound approval/config/code/inputs/smoke/atlas-stat/release-sidecars');return
 lockpath=R/'production_launches/execution_control_20260914/.scdrs_full_project.lock'
 with open(lockpath,'a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  if (W/'EXIT_STATUS.json').exists() or (W/'outputs').exists():raise FileExistsError('Existing production output; inspect before a new versioned attempt')
  wait_resources(cfg);cfg=verify(args.receipt);write(W/'STAGE_STATUS.json',{'status':'BUILDING_DERIVED','pid':os.getpid(),'time':time.time()});print('RESOURCE_GATE_PASSED; BUILDING_DERIVED; pid='+str(os.getpid()),flush=True)
  import anndata as ad,pandas as pd,gc
  from scdrs_adapter import build_derived,native_sequence
  a=ad.read_h5ad(cfg['atlas'],backed='r');assert list(a.shape)==[1175115,37336]
  assert not a.obs.level2_annotation.astype(str).str.contains('Pending',case=False).any()
  gm=pd.read_csv(cfg['gene_map'],sep='\t');cov=pd.read_csv(cfg['covariates'],sep='\t',index_col=0)
  report=build_derived(a,gm,cov,W/'outputs/derived/scdrs_raw_symbols.h5ad');a.file.close();del a,gm,cov;gc.collect()
  write(W/'STAGE_STATUS.json',{'status':'DERIVED_COMPLETE','derived':report,'time':time.time()});print('DERIVED_COMPLETE '+json.dumps(report),flush=True)
  # A new child loads one derived object at a time; release parent expression/graphs first.
  native_sequence(W/'outputs/derived/scdrs_raw_symbols.h5ad',W/'outputs/derived/native_covariates.cov.gz',cfg['scenarios'],W/'outputs/native',cfg['scdrs_executable'],report['group_keys'],12,resource_check=lambda:wait_resources(cfg))
  write(W/'EXIT_STATUS.json',{'exit_code':0,'status':'PAUSED_FOR_SCDRS_TRAITS_REVIEW','completed':time.time()})
if __name__=='__main__':
 try:main()
 except Exception:
  traceback.print_exc()
  if '--execute' in sys.argv:write(W/'EXIT_STATUS.json',{'exit_code':1,'error':traceback.format_exc(),'time':time.time()})
  raise
