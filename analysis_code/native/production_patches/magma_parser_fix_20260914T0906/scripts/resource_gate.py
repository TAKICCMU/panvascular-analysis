from pathlib import Path
import os,json,time,datetime,shutil
W=Path(__file__).resolve().parents[1];R=W.parents[1];previous=R/'production_launches/magma_category_fix_20260914T0301'
while True:
 reasons=[];active=[];exitfile=previous/'MAGMA_EXIT_STATUS.json'
 if not exitfile.exists():reasons.append('Prior MAGMA has not exited')
 else:
  try:assert isinstance(json.loads(exitfile.read_text())['exit_code'],int)
  except Exception:reasons.append('Invalid prior exit marker')
 for d in Path('/proc').iterdir():
  if not d.name.isdigit():continue
  try:
   cmd=(d/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
   if str(previous) in cmd:active.append({'pid':int(d.name),'command':cmd})
  except (PermissionError,ProcessLookupError,FileNotFoundError):pass
 if active:reasons.append('Prior MAGMA process still active')
 free=shutil.disk_usage(R).free/1024**3
 if free<40:reasons.append('Disk needs 20 floor+16 GPU+4 MAGMA GiB')
 mem={x.split(':')[0]:int(x.split()[1]) for x in Path('/proc/meminfo').read_text().splitlines() if ':' in x};available=mem['MemAvailable']/1024**2
 if available<96:reasons.append('RAM available <96GiB')
 status={'status':'WAITING_RESOURCE' if reasons else 'RESOURCE_GATE_PASSED','time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reasons':reasons,'active_previous':active,'disk_free_GiB':free,'mem_available_GiB':available,'CPU_scientific_threads':12};(W/'RESOURCE_GATE.json').write_text(json.dumps(status,indent=2));print(json.dumps(status),flush=True)
 if not reasons:break
 time.sleep(60)
