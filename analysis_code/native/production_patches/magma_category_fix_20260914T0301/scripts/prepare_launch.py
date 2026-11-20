import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='2'
os.environ['CUDA_VISIBLE_DEVICES']='';os.environ['PYTHONDONTWRITEBYTECODE']='1'
from pathlib import Path
import sys,json,datetime,shutil,shlex,hashlib,fcntl
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908');O=R/'production_launches/traits_20260914T023900Z';W=R/'production_launches/magma_category_fix_20260914T0301';C=W/'magma_adapter_v1';tests=json.loads((W/'evidence/CATEGORY_TEST_RESULTS.json').read_text());assert tests['tests_passed']==6
lock=open(R/'production_launches/execution_control_20260914/.auxiliary_compute.lock','a');fcntl.flock(lock,fcntl.LOCK_EX)
sys.path.insert(0,str(C));from src.dispatcher import code_digest,authorization
from src.review_dependencies import reviewed_dependencies
from src.core import sha
assert code_digest(C)==tests['new_code_digest'];c=json.loads((O/'configs/magma_config_APPROVED_r2.yaml').read_text());c['execution']['approval_receipt']=str(W/'approvals/magma_approval.json');c['outputs']['root']=str(W/'outputs');cf=W/'configs/magma_config_APPROVED.yaml';cf.write_text(json.dumps(c,indent=2));dep,dm=reviewed_dependencies(c,cf);scope=W/'approvals/USER_EXECUTION_SCOPE.md';shutil.copyfile(O/'approvals/USER_EXECUTION_SCOPE.md',scope)
command=['/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python','-B',str(C/'scripts/run_downstream.py'),'--config',str(cf),'--execute','--modules','cohort','magma']
sh='''#!/usr/bin/env bash
set -euo pipefail
JOB=__W__
exec 9>"$JOB/.magma_launch.lock"
flock -n 9 || exit 73
[[ ! -e "$JOB/MAGMA_EXIT_STATUS.json" ]] || exit 74
export CUDA_VISIBLE_DEVICES="" PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
export OMP_NUM_THREADS=12 OPENBLAS_NUM_THREADS=12 MKL_NUM_THREADS=12 NUMBA_NUM_THREADS=12 NUMEXPR_NUM_THREADS=12
/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python -B "$JOB/scripts/resource_gate.py"
date -Is > "$JOB/MAGMA_STARTED.txt"
set +e
__CMD__
rc=$?
set -e
printf '{"exit_code":%s,"completed_at":"%s"}\\n' "$rc" "$(date -Is)" > "$JOB/MAGMA_EXIT_STATUS.json"
exit "$rc"
'''.replace('__W__',shlex.quote(str(W))).replace('__CMD__',shlex.join(command));launch=W/'start_magma.sh';launch.write_text(sh);launch.chmod(0o750)
r={'approved':True,'approval_basis':'Current explicit bounded MAGMA execution authorization, necessary independent identity-code compatibility repair only','approved_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'atlas_sha256':c['atlas']['input_sha256'],'config_sha256':sha(cf),'code_sha256':code_digest(C),'reviewed_dependency_digest':dep,'reviewed_dependency_manifest':dm,'allowed_modules':['cohort','magma'],'clinical_model_approval':False,'scdrs_approval':False,'prior_failed_run':str(O),'user_approval_evidence':str(scope),'user_approval_evidence_sha256':sha(scope),'validation':str(W/'evidence/CATEGORY_TEST_RESULTS.json'),'validation_sha256':sha(W/'evidence/CATEGORY_TEST_RESULTS.json'),'exact_diff_sha256':sha(W/'evidence/EXACT_SOURCE_CHANGE.diff'),'launch_control_sha256':{str(p):sha(p) for p in [launch,W/'scripts/resource_gate.py']},'disk_reservations_GiB':{'floor':20,'GPU':16,'MAGMA':4},'workers':1,'per_worker_threads':12}
rp=W/'approvals/magma_approval.json';rp.write_text(json.dumps(r,indent=2));authorization(c,cf)
for p in [cf,rp,scope]:p.chmod(0o440)
(W/'LAUNCH_MANIFEST.json').write_text(json.dumps({'command':command,'config':str(cf),'receipt':str(rp),'code_dir':str(C),'code_sha256':code_digest(C),'tests_passed':6,'science_started':False,'prior_failure_preserved':str(O),'output_root':str(W/'outputs')},indent=2));fcntl.flock(lock,fcntl.LOCK_UN);print(str(launch))
