import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='2'
os.environ['PYTHONDONTWRITEBYTECODE']='1';os.environ['CUDA_VISIBLE_DEVICES']=''
import sys,json,hashlib,datetime,shutil,shlex
from pathlib import Path
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908');C=R/'downstream_CPU_L1L2_v2_1_release_compat_20260914';W=R/'production_launches/traits_20260914T023900Z';P=R/'production_launches/phase1_/20260913T183720_920345Z'
sys.path.insert(0,str(C));import yaml
from src.dispatcher import code_digest,authorization
from src.review_dependencies import reviewed_dependencies
from src.core import sha
for d in ['configs','approvals','outputs']:(W/d).mkdir(exist_ok=True)
assert code_digest(C)=='7ecd30a62e3528a4c3c6ec8329124c8f2754ff90f65c8b28e75fcfc0f7e695df'
audit=json.loads((W/'evidence/MAGMA_NATIVE_SMOKE.json').read_text());assert audit['returncode']==0 and audit['schema_pass']
resources=json.loads((W/'evidence/TRAIT_RESOURCES_VERIFIED.json').read_text());scenario=[{'name':r['trait'],'gene_results':r['path'],'sha256':r['sha256']} for r in resources if r.get('path') and r['module']=='magma']
for s in scenario:assert sha(s['gene_results'])==s['sha256']
c=yaml.safe_load((C/'configs/phase1_config_PROPOSED.yaml').read_text());c['execution'].update(mode='approved_execution',run_full_analysis=True,cpu_only=True,total_threads=40,max_workers=1,per_worker_threads=12,approval_receipt=str(W/'approvals/magma_approval.json'),enabled_modules=['cohort','magma'],overwrite=False,resume=False);c['outputs']['root']=str(W/'outputs');c['covarnet']={'run_nmf':False,'approve_named_state_frequency_rule':False,'workers':1,'threads_per_worker':12};c['external']={'gene_map':str(W/'inputs/gene_map_MAGMA_READY.tsv'),'magma':{'executable':'/home/anzhen_vhd/scRNAseq/magma','scenarios':scenario}};c.setdefault('pseudobulk',{})['approve_legacy_ql_se_approximation']=False
cf=W/'configs/magma_config_APPROVED.yaml';cf.write_text(json.dumps(c,indent=2,ensure_ascii=False));dep,manifest=reviewed_dependencies(c,cf)
scope=W/'approvals/USER_EXECUTION_SCOPE.md';scope.write_text('MAGMA execution scope: the three registered scenarios after native input and mapping verification; cohort and MAGMA only, one 12-thread worker within the 40-thread project limit; do not overwrite accepted outputs.\n')
r={'approved':True,'approval_basis':'Explicit current user request limited to three registered MAGMA scenarios after native smoke and mapping verification','approved_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'atlas_sha256':c['atlas']['input_sha256'],'config_sha256':sha(cf),'code_sha256':code_digest(C),'reviewed_dependency_digest':dep,'reviewed_dependency_manifest':manifest,'user_approval_evidence':str(scope),'user_approval_evidence_sha256':sha(scope),'allowed_modules':['cohort','magma'],'scdrs_approved_in_this_receipt':False,'clinical_model_approval':False,'cpu_worker_threads':12,'max_workers':1,'scheduling_constraint':'Only after exact phase1 job CPU_EXIT_STATUS and no live process references its job path'}
(W/'approvals/magma_approval.json').write_text(json.dumps(r,indent=2,ensure_ascii=False));authorization(c,cf)
cmd=['/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python','-B',str(C/'scripts/run_downstream.py'),'--config',str(cf),'--execute','--modules','cohort','magma']
sh='''#!/usr/bin/env bash
set -euo pipefail
JOB=__WORK__
exec 9>"$JOB/.magma_queue.lock"
flock -n 9 || exit 73
if [[ -e "$JOB/MAGMA_EXIT_STATUS.json" ]]; then exit 74; fi
export CUDA_VISIBLE_DEVICES=""
export PYTHONDONTWRITEBYTECODE=1
export MPLBACKEND=Agg
export OMP_NUM_THREADS=12 OPENBLAS_NUM_THREADS=12 MKL_NUM_THREADS=12 NUMBA_NUM_THREADS=12 NUMEXPR_NUM_THREADS=12
/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python -B "$JOB/scripts/wait_for_phase1.py"
date -Is > "$JOB/MAGMA_STARTED.txt"
set +e
__CMD__
rc=$?
set -e
printf '{"exit_code":%s,"completed_at":"%s"}\\n' "$rc" "$(date -Is)" > "$JOB/MAGMA_EXIT_STATUS.json"
exit "$rc"
'''.replace('__WORK__',shlex.quote(str(W))).replace('__CMD__',shlex.join(cmd))
(W/'start_magma_after_phase1.sh').write_text(sh);(W/'start_magma_after_phase1.sh').chmod(0o750)
for p in [cf,W/'approvals/magma_approval.json',scope]:p.chmod(0o440)
(W/'MAGMA_LAUNCH_MANIFEST.json').write_text(json.dumps({'status':'WAITING_RESOURCE','actual_science_started':False,'command':cmd,'CPU_worker_threads':12,'max_workers':1,'dispatcher_threads':12,'GPU_CPU_reserved':4,'estimated_magma_output_GiB':4,'minimum_unreserved_free_GiB':20,'phase1_output_reserve_GiB':10,'required_free_GiB':34,'config_sha256':sha(cf),'code_sha256':code_digest(C),'approval_receipt':str(W/'approvals/magma_approval.json'),'phase1_job':str(P)},indent=2));print(str(W))
