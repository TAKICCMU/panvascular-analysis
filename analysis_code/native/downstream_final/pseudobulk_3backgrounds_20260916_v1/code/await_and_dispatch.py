import os,json,time,subprocess,sys
from pathlib import Path
B=Path(__file__).resolve().parents[1]
while not (B/'PREPARATION.json').exists():
 # A disappeared exporter without its pass receipt is a failure, never READY.
 p=Path('/proc/1218024/cmdline')
 if not p.exists():raise RuntimeError('Exporter ended without PREPARATION.json; inspect preparation.log')
 time.sleep(15)
assert json.loads((B/'PREPARATION.json').read_text())['status']=='PASS_READY_FOR_MODEL_EXECUTION'
cmd=['/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python','-u',str(B/'code/dispatch.py')]
print('EXECUTING',cmd,flush=True);os.execv(cmd[0],cmd)
