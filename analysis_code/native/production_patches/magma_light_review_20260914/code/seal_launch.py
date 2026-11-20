import os,json,hashlib,subprocess,time
from pathlib import Path
W=Path(__file__).resolve().parents[1];R=W.parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
source=R/'production_launches/magma_recovery_roundtrip_fix_20260914T0930/outputs/runs/20260914T013548_608658_d62bd9f120de/magma';out=Path('/data/anzhen_vhd/scRNAseq/codex/panvascular_downstream_storage_v1/scdrs_storage_v1/magma_completed_review_v1')
assert not out.exists() and not (W/'LAUNCH.json').exists()
c={'source':str(source),'output':str(out),'atlas_sha256':'00fc74c1ddaa5e4b32d8664f134e0ec5a0c762f5e9a9fb45b96acbd40804a310','scenarios':['primary_r0.1','primary_r0.6','leave_one_CAD_r0.1'],'levels':['L1','L2'],'weights':['cell_weighted','study_balanced'],'project_roots':[str(R),'/data/anzhen_vhd/scRNAseq/codex/panvascular_downstream_storage_v1']};(W/'config.json').write_text(json.dumps(c,indent=2))
assert json.loads((W/'TEST_RESULTS.json').read_text())['status']=='PASS'
inputs=[source.parent/n for n in ['RUN_IDENTITY.json','RUN_STATUS.json','magma.receipt.json']]+[source/n for n in ['COMPLETE.json','RECOVERY_PROGRESS.json','L1_variables.tsv','L2_variables.tsv']]
for s in c['scenarios']:
 for l in c['levels']:
  for w in c['weights']:
   family=f'{s}_{l}_{w}'
   inputs.extend(source/(family+x) for x in ['_with_BH.tsv','_model_mapping.tsv','.gsa.out'])
identity=json.loads((source.parent/'RUN_IDENTITY.json').read_text());assert identity['atlas_sha256']==c['atlas_sha256']
for p,h in identity['files'].items():assert sha(p)==h;inputs.append(Path(p))
inputs+=list((W/'code').glob('*.py'))+[W/'config.json',W/'TEST_RESULTS.json',Path('/data/anzhen_vhd/scRNAseq/codex/panvascular_downstream_storage_v1/communication_data_v1/code/mount_resources.py')]
receipt={'scope':'magma_completed_output_review_only','authorization':'Main thread relays explicit user authorization for one-thread summary/plots only; no new MAGMA test','bound_files':[{'path':str(p),'sha256':sha(p)} for p in sorted(set(inputs))]};(W/'RECEIPT.json').write_text(json.dumps(receipt,indent=2))
cmd=['/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python','-B',str(W/'code/run_review.py')];env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:env[k]='1'
with open(W/'run.log','w') as f:p=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT,start_new_session=True,env=env)
(W/'LAUNCH.json').write_text(json.dumps({'pid':p.pid,'command':cmd,'receipt_sha256':sha(W/'RECEIPT.json'),'time':time.time()},indent=2));print(json.dumps({'pid':p.pid,'receipt_sha256':sha(W/'RECEIPT.json')}))
