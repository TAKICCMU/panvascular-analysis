import os,json,hashlib,subprocess,sys,time,shutil
from pathlib import Path
P=Path(__file__).resolve().parent;C=P.parent;W=C.parents[1];old=C/'associations';root=W/'covarnet'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
for k in ['R_HOME','R_LIBS','R_LIBS_SITE','R_LIBS_USER']:os.environ.pop(k,None)
old_binding=json.loads((old/'INPUT_BINDING.json').read_text())
for p,h in old_binding.items():assert sha(p)==h,'Original content changed: '+p
reused=['aligned_latest_metadata.tsv','cross_study_donor_overlap.tsv','disease_same_tissue_cross_study_meta.tsv','disease_within_study_tissue_background.tsv','raw_H_descriptive_by_group.tsv','scope_audit.tsv']
inputs={str(old/x):sha(old/x) for x in reused};inputs.update({str(p):sha(p) for p in P.glob('*') if p.is_file()});inputs['original_binding_sha']=sha(old/'INPUT_BINDING.json')
out=P/'results';out.mkdir(exist_ok=True);r=out/'INPUT_BINDING.json'
if r.exists():assert json.loads(r.read_text())==inputs,'Patch binding changed'
else:r.write_text(json.dumps(inputs,indent=2))
for x in reused:shutil.copy2(old/x,out/x)
(out/'REUSED_DISEASE_RECEIPT.json').write_text(json.dumps({x:sha(old/x) for x in reused},indent=2))
status=out/'STATUS.json';status.write_text(json.dumps({'status':'RUNNING','started':time.time()}))
cmd=['/usr/bin/Rscript','--vanilla',str(P/'association_patch_v2.R'),str(root),str(old),str(out)]
a=subprocess.run(cmd)
status.write_text(json.dumps({'status':'COMPLETED' if a.returncode==0 else 'FAILED','returncode':a.returncode,'finished':time.time(),'inputs':inputs,'outputs':{str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file() and p!=status}},indent=2))
sys.exit(a.returncode)
