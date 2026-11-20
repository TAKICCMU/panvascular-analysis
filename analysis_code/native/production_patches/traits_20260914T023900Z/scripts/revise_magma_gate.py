import os,signal,json,hashlib,datetime,shlex,sys
from pathlib import Path
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908');W=R/'production_launches/traits_20260914T023900Z';C=R/'downstream_CPU_L1L2_v2_1_release_compat_20260914'
assert not (W/'MAGMA_STARTED.txt').exists(),'Science already started; do not revise'
terminated=[]
for pid,fragment in [(34069,str(W/'start_magma_after_phase1.sh')),(34071,str(W/'scripts/wait_for_phase1.py'))]:
 p=Path(f'/proc/{pid}/cmdline')
 if p.exists():
  cmd=p.read_bytes().replace(b'\0',b' ').decode();assert fragment in cmd;os.kill(pid,signal.SIGTERM);terminated.append(pid)
(W/'MAGMA_QUEUE_R1_SUPERSEDED.json').write_text(json.dumps({'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'previous_pids':terminated,'scientific_work_started':False,'reason':'Prestart GPU reservation updated; old files preserved, new gate requires 40 GiB'},indent=2))
old=(W/'scripts/wait_for_phase1.py').read_text();new=old.replace('if free<34:', 'if free<40:').replace("'Free disk below 34 GiB including 20 GiB floor, 10 GiB phase1 reserve, 4 GiB MAGMA reserve'", "'Free disk below 40 GiB including 20 GiB floor, 16 GiB GPU reservation, 4 GiB MAGMA reservation'").replace("'scientific_work_started':False}", "'scientific_work_started':False,'fixed_reservations_GiB':{'floor':20,'GPU':16,'MAGMA':4},'current_resource_ledger':json.loads((R/'production_launches/execution_control_20260914/RESOURCE_LEDGER.json').read_text())}")
newgate=W/'scripts/wait_for_phase1_r2.py';newgate.write_text(new);compile(new,str(newgate),'exec')
c=json.loads((W/'configs/magma_config_APPROVED.yaml').read_text());c['execution']['approval_receipt']=str(W/'approvals/magma_approval_r2.json');cf=W/'configs/magma_config_APPROVED_r2.yaml';cf.write_text(json.dumps(c,indent=2,ensure_ascii=False))
sh=(W/'start_magma_after_phase1.sh').read_text().replace('scripts/wait_for_phase1.py','scripts/wait_for_phase1_r2.py').replace('magma_config_APPROVED.yaml','magma_config_APPROVED_r2.yaml');launch=W/'start_magma_after_phase1_r2.sh';launch.write_text(sh);launch.chmod(0o750)
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='2'
sys.path.insert(0,str(C));from src.core import sha
from src.dispatcher import code_digest,authorization
from src.review_dependencies import reviewed_dependencies
dep,manifest=reviewed_dependencies(c,cf);r=json.loads((W/'approvals/magma_approval.json').read_text());r.update(approved_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),config_sha256=sha(cf),reviewed_dependency_digest=dep,reviewed_dependency_manifest=manifest,launch_control_sha256={str(p):sha(p) for p in [launch,newgate]},resource_policy={'disk_floor_GiB':20,'GPU_reservation_GiB':16,'MAGMA_reservation_GiB':4,'minimum_free_GiB':40,'live_ledger':'Resource states only; not immutable scientific input'});rp=W/'approvals/magma_approval_r2.json';rp.write_text(json.dumps(r,indent=2));cf.chmod(0o440);rp.chmod(0o440);authorization(c,cf)
(W/'MAGMA_LAUNCH_MANIFEST_r2.json').write_text(json.dumps({'status':'WAITING_RESOURCE','science_started':False,'replaces_prestart_queue':str(W/'MAGMA_QUEUE_R1_SUPERSEDED.json'),'config':str(cf),'receipt':str(rp),'launch_script':str(launch),'code_sha256':code_digest(C),'command':['/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python','-B',str(C/'scripts/run_downstream.py'),'--config',str(cf),'--execute','--modules','cohort','magma'],'resource_policy':r['resource_policy']},indent=2));print(str(launch))
