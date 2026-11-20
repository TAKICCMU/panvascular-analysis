from pathlib import Path
import json,subprocess,datetime,os,shutil,sys,time,csv
ROOT=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908')
BASE=ROOT/'production_launches';ADMIN=BASE/'execution_control_20260914'
CPU=BASE/'phase1_/20260913T183720_920345Z'
def command(args):
 p=subprocess.run(args,capture_output=True,text=True);return p.stdout.strip() or p.stderr.strip()
def snapshot():
 runs=list((CPU/'outputs/runs').glob('*/RUN_STATUS.json'))
 statuses={str(p.parent):json.loads(p.read_text()) for p in runs}
 pid=json.loads((CPU/'PID.json').read_text())['pid']
 procs={}
 for p in Path('/proc').glob('[0-9]*'):
  try:
   d=dict(x.split(':',1) for x in (p/'status').read_text().splitlines() if ':' in x)
   procs[int(p.name)]={'pid':int(p.name),'ppid':int(d['PPid']),'threads':int(d['Threads']),'rss':d.get('VmRSS','').strip(),'cmd':(p/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace').strip()}
  except (OSError,ValueError,KeyError):pass
 own={pid};changed=True
 while changed:
  nxt=own|{p for p,x in procs.items() if x['ppid'] in own};changed=nxt!=own;own=nxt
 mem=dict(x.split(':',1) for x in Path('/proc/meminfo').read_text().splitlines())
 result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'cpu_modules':statuses,'cpu_exit':json.loads((CPU/'CPU_EXIT_STATUS.json').read_text()) if (CPU/'CPU_EXIT_STATUS.json').exists() else None,'cpu_processes':[procs[x] for x in sorted(own) if x in procs],'all_project_processes':[v for v in procs.values() if str(BASE) in v['cmd']],'mem_available_kib':int(mem['MemAvailable'].split()[0]),'disk_free_gib':shutil.disk_usage(ROOT).free/1024**3,'load_average':os.getloadavg(),'gpu':command(['nvidia-smi','--query-gpu=uuid,utilization.gpu,memory.used,memory.total','--format=csv,noheader']),'gpu_processes':command(['nvidia-smi','--query-compute-apps=pid,process_name,used_memory','--format=csv,noheader']),'other_branch_statuses':{str(p):p.read_text()[:5000] for p in BASE.glob('*/RUN_STATUS.json')}}
 target=ADMIN/'CURRENT_PROGRESS.json';tmp=target.with_suffix('.'+str(os.getpid())+'.tmp');tmp.write_text(json.dumps(result,indent=2));tmp.replace(target)
 with (ADMIN/'RESOURCE_PROGRESS.jsonl').open('a') as f:f.write(json.dumps(result)+'\n')
 return result
if __name__=='__main__':
 while True:
  r=snapshot()
  if not any(x in sys.argv for x in ['--watch','--watch-all']):print(json.dumps(r,indent=2));break
  if '--watch-all' in sys.argv:
   active=[x for x in r['all_project_processes'] if 'check_progress.py' not in x['cmd'] and 'execution_control_20260914' not in x['cmd']]
   if not active:break
  elif r['cpu_exit'] is not None:break
  time.sleep(60)
