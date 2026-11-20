from pathlib import Path
import json,os,shutil
ROOT=Path("/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908")
LEDGER=ROOT/"maintenance_resume_20260914/RESOURCE_RESUME_LEDGER.json"
def project_processes():
 result=[]
 for p in Path("/proc").iterdir():
  if not p.name.isdigit():continue
  try:
   args=(p/"cmdline").read_bytes().split(bytes([0]));words=[x.decode(errors="replace") for x in args if x]
   # Python controller/worker or R program with explicit project file arguments.
   if not any(x.startswith(str(ROOT)+"/") for x in words):continue
   if any("python3 - <<" in x for x in words):continue
   st=dict(x.split(":",1) for x in (p/"status").read_text().splitlines() if ":" in x)
   if st["State"].strip()[0] in "TZt":continue
   result.append(dict(pid=int(p.name),ppid=int(st["PPid"]),threads=int(st["Threads"]),cmd=words))
  except (FileNotFoundError,PermissionError,ProcessLookupError):pass
 return result
def disk_bytes(paths):
 total=0
 for x in paths:
  p=Path(x)
  if p.exists():total+=sum(f.stat().st_size for f in p.rglob("*") if f.is_file())
 return total
def decision(free,remaining,threads,extra=5):
 required=20+sum(remaining.values())
 return dict(pass_gate=free>=required and threads+extra<=40,required_free_gib=required,projected_threads=threads+extra)
def coordination(extra=5):
 d=json.loads(LEDGER.read_text());remaining={};exceeded=[]
 for k,v in d["disk_reservations"].items():
  used=max(0,disk_bytes(v["paths"])-v["baseline_bytes"])/1024**3
  remaining[k]=max(0,v["allocation_gib"]-used)
  if used>v["allocation_gib"]:exceeded.append(k)
 proc=project_processes();t=sum(x["threads"] for x in proc)
 free=shutil.disk_usage(ROOT).free/1024**3;r=decision(free,remaining,t,extra)
 r.update(free_gib=free,remaining_gib=remaining,processes=proc,threads=t,over_budget=exceeded)
 r["pass_gate"]=r["pass_gate"] and not exceeded
 return r
