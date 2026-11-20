import os,json,hashlib,subprocess,sys,time
from pathlib import Path
P=Path(__file__).resolve().parents[1];W=P.parents[1]
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
side=W/'disease_background_3groups_review_v2/sample_background_sidecar_v2.tsv';root=W/'covarnet'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(side)=='852c08353c53a1b4d32b0b5928c9d1ac6a0278e365748912adfaf2ff30cf2ce4'
os.environ['COVARNET_SIDECAR']=str(side)
mode=sys.argv[1];assert mode in ['refits','associations']
inputs=list((P/'code').glob('*'))+[side]+list(root.glob('*.tsv'))+list((root/'stable_NMF').glob('*K14*'))+[root/'stable_NMF/NMF_input_minmax.tsv']
binding={str(x):sha(x) for x in inputs if x.is_file()}
out=P/mode;out.mkdir(exist_ok=True)
receipt=out/'INPUT_BINDING.json'
if receipt.exists():assert json.loads(receipt.read_text())==binding,'Input/code changed: fresh run output required'
else:receipt.write_text(json.dumps(binding,indent=2))
rs='/usr/bin/Rscript' if mode=='associations' else '/home/anzhen_vhd/miniconda3/envs/communication_r/bin/Rscript'
cmd=[rs,'--vanilla',str(P/'code'/('remaining_refits.R' if mode=='refits' else 'associations.R')),str(root)]
if mode=='associations':cmd.append(str(side))
cmd.append(str(out))
status=out/'STATUS.json';status.write_text(json.dumps({'status':'RUNNING','started':time.time(),'command':cmd}))
r=subprocess.run(cmd)
status.write_text(json.dumps({'status':'COMPLETED' if r.returncode==0 else 'FAILED','returncode':r.returncode,'finished':time.time(),'inputs':binding,'outputs':{str(x.relative_to(out)):sha(x) for x in out.rglob('*') if x.is_file() and x!=status}},indent=2))
sys.exit(r.returncode)
