"""Replacement guard for currently running adopted native. Does not resume old controller."""
from pathlib import Path
import os,json,time,signal
from storage_utils import *
W=Path(__file__).resolve().parents[1]
c=json.loads((W/'config.json').read_text());target=c['native'];paused=False
while same(target,zombie_ok=True) and proc(target['pid'])['state']!='Z':
 try:
  b=mount_budget(c);danger=any(x['free_GiB']<x['floor_GiB']+2 for x in b['filesystems'].values());state=proc(target['pid'])['state']
  if danger and not paused:
   if state in ['T','t']:write(W/'MOUNT_GUARD_STATUS.json',dict(status='EXTERNAL_PAUSE_NOT_OWNED',time=time.time()));break
   send(target,signal.SIGSTOP);paused=True;write(W/'MOUNT_GUARD_PAUSE.json',dict(time=time.time(),target=target,budget=b))
  elif paused and not b['reasons']:
   send(target,signal.SIGCONT);paused=False;write(W/'MOUNT_GUARD_RESUME.json',dict(time=time.time(),target=target,budget=b))
  write(W/'MOUNT_GUARD_STATUS.json',dict(status='PAUSED_FOR_MOUNT_SAFETY' if paused else 'WATCHING_DATA_SCORE',pid=os.getpid(),time=time.time(),target=target,budget=b))
 except Exception as e:
  if same(target) and proc(target['pid'])['state'] not in ['T','t']:send(target,signal.SIGSTOP);paused=True
  write(W/'MOUNT_GUARD_STATUS.json',dict(status='PAUSED_MOUNT_IDENTITY_OR_GUARD_ERROR',error=repr(e),pid=os.getpid(),time=time.time()));break
 time.sleep(15)
else:write(W/'MOUNT_GUARD_STATUS.json',dict(status='NATIVE_EXITED_NO_FURTHER_ACTION',pid=os.getpid(),time=time.time(),target=proc(target['pid'])))
