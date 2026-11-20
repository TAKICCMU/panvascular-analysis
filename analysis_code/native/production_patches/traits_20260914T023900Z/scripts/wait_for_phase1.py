import os,time,json,shutil,datetime
from pathlib import Path
W=Path(__file__).resolve().parents[1];R=W.parents[1];P=R/'production_launches/phase1_/20260913T183720_920345Z'
status=W/'MAGMA_QUEUE_STATUS.json';once=False
while True:
 reasons=[];done=P/'CPU_EXIT_STATUS.json';phase1=None
 if not done.is_file():reasons.append('Exact phase1 completion marker absent')
 else:
  try:phase1=json.loads(done.read_text());assert isinstance(phase1['exit_code'],int)
  except Exception:reasons.append('Invalid phase1 completion marker')
 active=[]
 for d in Path('/proc').iterdir():
  if not d.name.isdigit():continue
  try:
   cmd=(d/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
   if str(P) in cmd and int(d.name)!=os.getpid():active.append({'pid':int(d.name),'command':cmd})
  except (PermissionError,ProcessLookupError,FileNotFoundError):pass
 if active:reasons.append('Phase1 job process still active')
 free=shutil.disk_usage(R).free/1024**3
 if free<34:reasons.append('Free disk below 34 GiB including 20 GiB floor, 10 GiB phase1 reserve, 4 GiB MAGMA reserve')
 mem={x.split(':')[0]:int(x.split()[1]) for x in Path('/proc/meminfo').read_text().splitlines() if ':' in x};available=mem['MemAvailable']/1024**2
 if available<96:reasons.append('MemAvailable below 96 GiB')
 # This is a queue gate only. Native authorization and preflight run afterward without bypass.
 payload={'status':'WAITING_RESOURCE' if reasons else 'RESOURCE_GATE_PASSED','time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reasons':reasons,'active_phase1':active,'phase1_exit':phase1,'free_gib':free,'mem_available_gib':available,'scientific_work_started':False}
 status.write_text(json.dumps(payload,indent=2))
 if not reasons:print(json.dumps(payload),flush=True);break
 if not once:print(json.dumps(payload),flush=True);once=True
 time.sleep(60)
