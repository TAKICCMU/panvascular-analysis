"""Prepare a content-bound scDRS-only receipt after actual native smoke passes. Does not launch."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
import json,hashlib,datetime,subprocess,sys
from pathlib import Path
W=Path(__file__).resolve().parents[1];R=W.parents[1]
from production_runner import sha,directory_bytes
for dirname in ['configs','approvals']: (W/dirname).mkdir(exist_ok=True)
validation=W/'tests/FINAL_SMOKE_VALIDATION.json';l2=W/'tests/actual_l2_obs_smoke_v3/SMOKE_RESULT.json';s=json.loads(l2.read_text());assert s['status']=='PASS_NATIVE_SHARED_ADAPTER_ACTUAL_OBS'
subprocess.run([sys.executable,str(W/'tests/validate_native_outputs.py')],check=True,stdout=subprocess.DEVNULL)
for folder in ['actual_obs_smoke_v1','actual_l2_obs_smoke_v3']:
 for file in (W/'tests'/folder/'native').glob('*/SCHEMA_VALIDATION.json'):
  q=json.loads(file.read_text());assert q['n_control_columns']==1000 and q['n_scored_cells']==1400
  for key in ['level1_annotation_scdrs_test','level2_annotation_scdrs_test']:
   assert q[key]['nonfinite_assoc']==q[key]['nonfinite_hetero']==0
  if folder.endswith('v3'):assert q['level2_annotation_scdrs_test']['n_eligible_groups']>=2
subprocess.run([sys.executable,str(W/'tests/test_resource_accounting.py')],check=True)
subprocess.run([sys.executable,str(W/'tests/test_exact_counts.py')],check=True)
# Actual numeric control matrix checks are deliberately small-input only.
import pandas as pd,numpy as np
for folder in ['actual_obs_smoke_v1','actual_l2_obs_smoke_v3']:
 for file in (W/'tests'/folder/'native').glob('*/scores/*.full_score.gz'):
  z=pd.read_csv(file,sep='\t',index_col=0);assert np.isfinite(z.to_numpy()).all();assert z.index.is_unique
scope=W/'approvals/USER_SCOPE.md';scope.write_text('''用户已明确授权修正并重跑 MAGMA/scDRS。此回执仅覆盖 scDRS 原生 compute-score 与 perform-downstream 三个既定遗传场景，n_ctrl=1000；不修改标签、图谱、旧代码/结果、GPU分支或其他模块。\n当前负责人已授权 READY map 经真实小测后，在严格资源门控与完整内容绑定后排队或启动，不需要重复请求泛泛确认。\n图策略已明确批准保留 frozen_global_n30；不声称等价于原生无图默认PCA20/n15，hetero_mcp依赖该图而assoc_mcp不依赖邻接。\n原始详细范围见 parallel_execution_handoff_20260914/CODEX_START_PARALLEL_EXECUTION_PROMPT.md。\n''')
manifest=R/'panvascular_CPU_theta2_L1L2_final_v1.acceptance_r1.manifest.json';stat=R/'production_launches/execution_control_20260914/AUTHENTICATED_ATLAS_STAT.json';ledger=R/'production_launches/execution_control_20260914/RESOURCE_LEDGER_traits_repair_v2.json';m=json.loads((R/'enrichment/scdrs_label_independent/genesets/GENESETS_FROZEN.json').read_text())
assert sha(ledger)=='dd9d55c3d83a4bc867eb3d879b3afadb6a44d6748da75c6f9f4f27f049286c28'
external={'GPU':{'paths':[str(R/'production_launches/spatial_gpu_20260913T183753Z/outputs'),str(R/'production_launches/spatial_gpu_20260913T183753Z/smoke')],'allocation_gib':16},'MAGMA':{'paths':[str(R/'production_launches/magma_recovery_roundtrip_fix_20260914T0930/outputs')],'allocation_gib':4},'scDRS':{'paths':[str(W/'outputs')],'allocation_gib':32}}
for record in external.values():record['baseline_bytes']=directory_bytes(record['paths'])
assert external['scDRS']['baseline_bytes']==0
cfg={'atlas':str(R/'panvascular_CPU_theta2_L1L2_final_v1.h5ad'),'atlas_sha256':'00fc74c1ddaa5e4b32d8664f134e0ec5a0c762f5e9a9fb45b96acbd40804a310','atlas_stat':str(stat),'manifest':str(manifest),'gene_map':str(W/'inputs/gene_map_scDRS_READY.tsv'),'covariates':str(R/'production_launches/traits_20260914T023900Z/inputs/scDRS_raw_verified.cov.gz'),'smoke_validation':str(validation),'graph_policy':'frozen_global_n30','n_ctrl':1000,'native_threads':12,'controller_threads':1,'scdrs_executable':'/home/anzhen_vhd/miniconda3/envs/scdrs_env/bin/scdrs','scenarios':[{'name':v['scenario'],'geneset':v['gs_file'],'sha256':v['gs_sha256']} for v in m['scenarios']],'disk_reservations':external,'reservation_baseline_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'resource_ledger':str(ledger)}
assert sha(cfg['gene_map'])=='55378df49e46c988dd641a633c041951281af5189447ca952fc805bfad85e24c'
assert sha(cfg['covariates'])=='826c7d919686301f548decf35685417d0d41e4b25ad67ca715dbdb4d87305037'
config=W/'configs/scdrs_APPROVED.json';receipt=W/'approvals/scdrs_APPROVAL_RECEIPT.json'
if config.exists() or receipt.exists():raise FileExistsError('Sealed config/receipt already exists')
config.write_text(json.dumps(cfg,indent=2))
files=[config,scope,manifest,stat,ledger,validation,l2,W/'evidence/NATIVE_VERSIONS_AND_SOURCE.json',Path(cfg['gene_map']),Path(cfg['covariates']),R/'parallel_execution_handoff_20260914/CODEX_START_PARALLEL_EXECUTION_PROMPT.md']
files+=list((W/'audit_scripts').glob('*.py'))+list((W/'code').glob('*.py'))+list((W/'tests').glob('*.py'))+list((W/'inputs').glob('*.json'))
files+=[Path(v['geneset']) for v in cfg['scenarios']]
cert=json.loads(stat.read_text());files.append(Path(cert['preflight']))
versions=json.loads((W/'evidence/NATIVE_VERSIONS_AND_SOURCE.json').read_text());files += [Path(x['path']) for x in versions['source']]
record={'mode':'approved_execution','allowed_modules':['scdrs'],'approval_source':'user_authorized_traits_repair_20260914','config':str(config),'bound_files':[{'path':str(p),'sha256':sha(p)} for p in sorted(set(files))],'generated_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'no_production_started_by_prepare':True}
receipt.write_text(json.dumps(record,indent=2));print(str(receipt))
