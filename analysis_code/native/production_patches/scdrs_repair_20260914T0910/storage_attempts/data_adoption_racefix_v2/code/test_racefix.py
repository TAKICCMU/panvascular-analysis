from pathlib import Path
import sys,json,tempfile,os,threading,time,py_compile
from unittest.mock import patch
W=Path(__file__).resolve().parents[1];sys.path.insert(0,str(W/'code'))
import storage_utils as u
passed=[]
with tempfile.TemporaryDirectory(prefix='scdrs_race_fixture_') as t:
 root=Path(t);p=root/'status.tmp';p.write_bytes(b'1234');dev=root.stat().st_dev
 assert u.reservation_bytes([root],dev)==4;passed.append('normal_bytes')
 real=Path.stat
 def vanish(self,*args,**kwargs):
  if self==p:
   self.unlink(missing_ok=True);raise FileNotFoundError(str(self))
  return real(self,*args,**kwargs)
 # Python3.8 Path.unlink missing_ok is supported.
 with patch.object(Path,'stat',vanish):assert u.reservation_bytes([root],dev)==0
 passed.append('deterministic_disappearing_file_between_listing_and_stat')
 p.write_bytes(b'1234')
 with patch.object(Path,'stat',side_effect=PermissionError('fixture')):
  try:u.reservation_bytes([root],dev)
  except PermissionError:passed.append('permission_failure_not_swallowed')
  else:raise AssertionError('PermissionError ignored')
 try:u.reservation_bytes([root],dev+1)
 except ValueError:passed.append('wrong_device_rejected')
 else:raise AssertionError('Wrong device ignored')
 p.unlink();stop=threading.Event()
 def replace():
  for i in range(3000):
   tmp=root/'status.json.tmp';tmp.write_bytes(b'x'*32);tmp.replace(root/'status.json')
 thread=threading.Thread(target=replace);thread.start()
 for i in range(3000):assert u.reservation_bytes([root],dev)>=0
 thread.join();passed.append('3000_concurrent_atomic_replacements')
C=json.loads((W/'config.json').read_text())
assert not u.mount_budget(C)['reasons'];passed.append('real_three_mount_budget_pass')
for p in (W/'code').glob('*.py'):py_compile.compile(str(p),doraise=True)
passed.append('code_compiles')
(W/'TEST_RACEFIX.json').write_text(json.dumps({'status':'PASS','tests':passed,'scientific_calculation_repeated':False},indent=2))
print(passed)
