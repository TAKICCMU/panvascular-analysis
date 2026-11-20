"""Real native preflight and metadata-only cohort entrance, no scientific workers."""
from pathlib import Path
import sys,os,json,time,traceback,datetime,shutil,importlib.metadata as im
sys.dont_write_bytecode=True
B=Path(__file__).resolve().parents[1];R=B.parent;sys.path.insert(0,str(B))
from src.core import preflight,sha,atomic_json
from src.cohort import cohort_from_obs,apply_unit_review,read_unit_review
from src.anchors import final_anchor_panel
from src.state_identity import build_state_identity_map
from src.core import read_frame
from src.dispatcher import authorization
c=json.loads((B/'configs/phase1_config_PROPOSED.yaml').read_text());p=Path(c['atlas']['path']);before=p.stat();t=time.time()
r={'status':'BLOCKED','h5ad_path':str(p),'h5ad_expected_sha256':c['atlas']['input_sha256'],'full_h5ad_hash_stream_bytes':before.st_size,'read_only':True,'scientific_workers_started':False,'production_receipt_issued':False}
try:
 obs,ont,pf=preflight(c,B/'evidence/NATIVE_PREFLIGHT.json');r['native_preflight']=pf;r['preflight_seconds']=time.time()-t
 print('Native preflight passed',flush=True)
 assert str(obs.donor_uid.dtype)=='category' and str(obs.biological_sample_uid.dtype)=='category'
 cm=cohort_from_obs(obs);review=read_unit_review(c['statistical_units']['cohort_design_file']);cohort=apply_unit_review(cm,review)
 r['actual_cohort']={'status':'PASS','native_input_dtypes':{x:str(obs[x].dtype) for x in ['donor_uid','biological_sample_uid','sample_uid']},'samples':len(cohort),'aggregate_approved':int(cohort.aggregate_approved.sum()),'inference_eligible':int(cohort.inference_eligible.sum()),'identity_values_unchanged':True}
 # Evidence only, not counts/marker analysis and no metadata writeback.
 genes=read_frame(p,'var');panel=final_anchor_panel(ont,genes);state=build_state_identity_map(ont)
 panel.to_csv(B/'evidence/ACTUAL_ANCHORS.tsv',sep='\t',index=False)
 r['actual_ontology']={'named_states':len(state),'anchors':len(panel),'special_counts':obs[c['atlas']['level2_key']].value_counts().reindex(['Unknown','Not_assessed']).to_dict(),'sample_uid':obs.sample_uid.nunique()}
 try:authorization(c,B/'configs/phase1_config_PROPOSED.yaml');raise RuntimeError('Unexpectedly authorized')
 except Exception as exc:
  if str(exc)!='CODE_READY_FOR_REVIEW_NOT_EXECUTED':raise
 r['execution_gate']='CODE_READY_FOR_REVIEW_NOT_EXECUTED';r['status']='PASS_FINAL_RELEASE_ACCEPTANCE';r['stop_status']='READY_FOR_PHASE1_EXECUTION_APPROVAL'
except Exception as exc:
 r['error']=str(exc);r['traceback']=traceback.format_exc();print(r['traceback'],flush=True)
after=p.stat();r['file_stat_unchanged']=(before.st_size,before.st_mtime_ns,before.st_ino)==(after.st_size,after.st_mtime_ns,after.st_ino)
r['input_stat']={'size':after.st_size,'mtime_ns':after.st_mtime_ns,'inode':after.st_ino,'nlink':after.st_nlink}
protected=json.loads((B/'evidence/PROTECTED_SMALL_FILES_BEFORE.json').read_text());changes=[p for p,h in protected.items() if not Path(p).is_file() or sha(p)!=h];r['old_source_manifest_support_unchanged']=not changes;r['changed_protected_files']=changes
if changes or not r['file_stat_unchanged']:r['status']='BLOCKED';r['stop_status']='BLOCKED'
r['elapsed_seconds']=time.time()-t;r['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();r['config_sha256']=sha(B/'configs/phase1_config_PROPOSED.yaml')
atomic_json(B/'FINAL_RELEASE_ACCEPTANCE_RECHECK.json',r)
env={'timestamp':r['finished_utc'],'cpu_count':os.cpu_count(),'load_average':os.getloadavg(),'memory_kib':{l.split(':')[0]:int(l.split()[1]) for l in Path('/proc/meminfo').read_text().splitlines() if l.startswith(('MemTotal:','MemAvailable:'))},'disk':shutil.disk_usage(R)._asdict(),'versions':{x:im.version(x) for x in ['anndata','pandas','numpy','scipy','h5py','omicverse']},'threads_for_this_process':{x:os.getenv(x) for x in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']}}
atomic_json(B/'evidence/ACTUAL_ENVIRONMENT.json',env);print(json.dumps({k:v for k,v in r.items() if k not in ['native_preflight','traceback']},default=str),flush=True)
if r['status']!='PASS_FINAL_RELEASE_ACCEPTANCE':sys.exit(2)
