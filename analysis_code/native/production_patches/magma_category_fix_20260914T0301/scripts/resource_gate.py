from pathlib import Path
import time,os,json,datetime,shutil
W=Path(__file__).resolve().parents[1];R=W.parents[1];prior=[R/'production_launches/phase1_/20260913T183720_920345Z',R/'production_launches/traits_20260914T023900Z']
while True:
 reasons=[];active=[]
 for job in prior:
  marker=job/('CPU_EXIT_STATUS.json' if 'phase1_' in str(job) else 'MAGMA_EXIT_STATUS.json')
  if not marker.exists():reasons.append('Prior job exit marker absent: '+str(marker))
  else:
   try:assert isinstance(json.loads(marker.read_text())['exit_code'],int)
   except Exception:reasons.append('Invalid exit marker: '+str(marker))
 for p in Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   cmd=(p/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
   if any(str(job) in cmd for job in prior):active.append({'pid':int(p.name),'command':cmd})
  except (PermissionError,FileNotFoundError,ProcessLookupError):pass
 if active:reasons.append('Prior scientific job process still active')
 free=shutil.disk_usage(R).free/1024**3
 if free<40:reasons.append('Require40GiB: floor20+GPU16+MAGMA4')
 mem={x.split(':')[0]:int(x.split()[1]) for x in Path('/proc/meminfo').read_text().splitlines() if ':' in x};available=mem['MemAvailable']/1024**2
 if available<96:reasons.append('Require96GiB available RAM')
 report={'status':'WAITING_RESOURCE' if reasons else 'RESOURCE_GATE_PASSED','time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reasons':reasons,'active_prior_jobs':active,'free_gib':free,'mem_available_gib':available,'reservations_gib':{'floor':20,'GPU':16,'MAGMA':4}}
 (W/'RESOURCE_GATE.json').write_text(json.dumps(report,indent=2))
 if not reasons:print(json.dumps(report),flush=True);break
 print(json.dumps(report),flush=True);time.sleep(60)
