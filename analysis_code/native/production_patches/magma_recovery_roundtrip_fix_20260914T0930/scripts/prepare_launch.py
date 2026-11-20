import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='2'
os.environ['PYTHONDONTWRITEBYTECODE']='1';os.environ['CUDA_VISIBLE_DEVICES']=''
from pathlib import Path
import sys,json,shutil,datetime,shlex,hashlib
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908');OLD=R/'production_launches/magma_category_fix_20260914T0301';W=R/'production_launches/magma_recovery_roundtrip_fix_20260914T0930';C=W/'magma_adapter_v3';RUN=OLD/'outputs/runs/20260913T190814_808850_7ccf71427e40';S=RUN/'magma';sys.path.insert(0,str(C));from src.dispatcher import code_digest,authorization
from src.review_dependencies import reviewed_dependencies
from src.core import sha
assert json.loads((W/'evidence/PARSER_TEST_RESULTS.json').read_text())['new_code_digest']==code_digest(C)
assert json.loads((W/'evidence/REAL_ROUNDTRIP_TEST_RESULTS.json').read_text())['code_digest']==code_digest(C)
oldcf=OLD/'configs/magma_config_APPROVED.yaml';oldc=json.loads(oldcf.read_text());auth=json.loads((S/'EXTERNAL_AUTHENTICATED.json').read_text());assert auth['config_sha256']==sha(oldcf);assert auth['code_sha256']==code_digest(OLD/'magma_adapter_v1');assert auth['atlas_sha256']==oldc['atlas']['input_sha256'];assert auth['preflight_sha256']==sha(S/'preflight.json')
commands=[json.loads(x) for x in (S/'commands.jsonl').read_text().splitlines()];assert len(commands)==1;cmd=commands[0];scenarios=oldc['external']['magma']['scenarios'];first=next(s for s in scenarios if s['name']=='primary_r0.1');assert cmd==[oldc['external']['magma']['executable'],'--gene-results',first['gene_results'],'--gene-covar',str(S/'L1_cell_weighted_specificity.tsv'),'--model','direction=pos','condition=Average','--out',str(S/'primary_r0.1_L1_cell_weighted')]
assert sha(first['gene_results'])==first['sha256'];assert 'End time is' in (S/'primary_r0.1_L1_cell_weighted.log').read_text()
files={}
for p in [S/x for x in ['L1_cell_weighted_specificity.tsv','L1_variables.tsv','L1_eligibility.tsv','L1_gene_mapping_audit.tsv','primary_r0.1_L1_cell_weighted.gsa.out','primary_r0.1_L1_cell_weighted.log','commands.jsonl','EXTERNAL_AUTHENTICATED.json','preflight.json']]+[oldcf,R/'production_launches/traits_20260914T023900Z/inputs/gene_map_MAGMA_READY.tsv',OLD/'magma_adapter_v1/src/external_modules.py',Path('/home/anzhen_vhd/scRNAseq/magma')]:files[p.name]={'path':str(p),'sha256':sha(p)}
cert={'atlas_sha256':oldc['atlas']['input_sha256'],'source_config_sha256':sha(oldcf),'source_code_digest':auth['code_sha256'],'specificity_source_sha256':sha(C/'src/external_modules.py'),'annotation_and_study_keys':{k:oldc['atlas'].get(k,{'level1_key':'level1_annotation','level2_key':'level2_annotation','study_key':'study'}[k]) for k in ['level1_key','level2_key','study_key']},'specificity_parameters':{'normalization':'full raw library CP1e6 natural log1p','min_cells':500,'min_studies':3,'min_study_cells':20,'cell_weighted_and_study_balanced_separate_specifications':True},'files':files,'reuse_policy':'Copy certified L1 cell-weighted table and first complete raw native gsa/log; read raw for unsaved L1 study-balanced and L2; validate recomputed L1 equivalence at serialization tolerance; retain all outputs and all conditional rows'}
mp=W/'evidence/RECOVERY_INPUT_MANIFEST.json';mp.write_text(json.dumps(cert,indent=2));c=oldc;c['outputs']['root']=str(W/'outputs');c['execution']['approval_receipt']=str(W/'approvals/magma_approval.json');c['execution']['max_workers']=1;c['execution']['per_worker_threads']=12;c['external']['magma']['recovery_manifest']=str(mp);c['external']['magma']['recovery_manifest_sha256']=sha(mp);cf=W/'configs/magma_config_APPROVED.yaml';cf.write_text(json.dumps(c,indent=2));dep,dm=reviewed_dependencies(c,cf)
scope=W/'approvals/USER_EXECUTION_SCOPE.md';shutil.copyfile(OLD/'approvals/USER_EXECUTION_SCOPE.md',scope)
command=['/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python','-B',str(W/'scripts/run_magma_light.py'),'--config',str(cf),'--execute','--modules','cohort','magma']
sh='''#!/usr/bin/env bash
set -euo pipefail
JOB=__W__
exec 9>"$JOB/.magma_execution.lock"
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
r={'approved':True,'approval_basis':'User explicitly requested repair and rerun MAGMA/scDRS; this receipt only authorizes MAGMA existing three scenarios with validated native conditional-output parser and certified reuse','approved_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'atlas_sha256':c['atlas']['input_sha256'],'config_sha256':sha(cf),'code_sha256':code_digest(C),'reviewed_dependency_digest':dep,'reviewed_dependency_manifest':dm,'allowed_modules':['cohort','magma'],'clinical_model_approval':False,'scdrs_approval_in_this_receipt':False,'user_approval_evidence':str(scope),'user_approval_evidence_sha256':sha(scope),'parser_tests_sha256':sha(W/'evidence/PARSER_TEST_RESULTS.json'),'actual_roundtrip_tests_sha256':sha(W/'evidence/REAL_ROUNDTRIP_TEST_RESULTS.json'),'recovery_manifest_sha256':sha(mp),'launch_control_sha256':{str(p):sha(p) for p in [launch,W/'scripts/resource_gate.py',W/'scripts/run_magma_light.py']},'disk_reservations_GiB':{'floor':20,'GPU':16,'MAGMA':4},'workers':1,'per_worker_threads':12};rp=W/'approvals/magma_approval.json';rp.write_text(json.dumps(r,indent=2));authorization(c,cf)
for p in [cf,rp,scope,mp]:p.chmod(0o440)
(W/'LAUNCH_MANIFEST.json').write_text(json.dumps({'command':command,'config':str(cf),'receipt':str(rp),'code_dir':str(C),'code_sha256':code_digest(C),'output_root':str(W/'outputs'),'prior_failed_run':str(RUN),'native_outputs_to_reuse':1,'expected_total_outputs':12,'expected_new_native_calls':11,'tests':20,'science_started':False},indent=2));print(str(launch))
