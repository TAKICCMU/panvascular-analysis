from common import *
import importlib.util,sys,traceback
OLD=R/'production_launches/scdrs_repair_20260914T0910/storage_attempts/third_score_completion_v3'
spec=importlib.util.spec_from_file_location('prior_verified_resume',OLD/'code/resume.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
OUT=Path(old.C['output'])
def main():
 old.verify()
 target=json.loads((OLD/'CURRENT_CHILD.json').read_text())['target'];pid=target['pid']
 p=Path('/proc')/str(pid)
 stat=(p/'stat').read_text().rsplit(')',1)[1].split();assert stat[19]==str(target['start_ticks'])
 cmd=(p/'cmdline').read_bytes().replace(b'\0',b' ').decode()
 assert 'perform-downstream' in cmd and str(OUT) in cmd and 'level1_annotation_scdrs_test' in cmd
 with open(old.legacy.C['project_lock'],'a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  assert sum(x['threads'] for x in processes())<40
  os.kill(pid,signal.SIGCONT)
  dump(W/'status/scdrs_adopted.json',dict(status='RESUMED_EXISTING_NATIVE_L1',pid=pid,original_start_ticks=target['start_ticks'],original_score_recomputed=False))
  while p.exists():
   try:
    s=(p/'stat').read_text().rsplit(')',1)[1].split()
    if s[0]=='Z':break
    if s[19]!=str(target['start_ticks']):break
   except OSError:break
   time.sleep(20)
  validator=Path(old.C['previous_attempt'])/'code/validate_native.py'
  key='level1_annotation_scdrs_test'
  run_job('scdrs_third_L1_validation',[PY,validator,'group','--run',OUT,'--key',key],W/'logs/scdrs_third_L1_validation.log',reserve=1)
  key='level2_annotation_scdrs_test';group=OUT/key;group.mkdir(exist_ok=False)
  score=list((OUT/'scores').glob('*.full_score.gz'));assert len(score)==1
  cmd=[old.legacy.OLD['scdrs_executable'],'perform-downstream','--h5ad-file',old.legacy.OLD['derived'],'--score-file',score[0],'--out-folder',group,'--group-analysis',key,'--flag-filter-data','True','--flag-raw-count','True']
  run_job('scdrs_third_L2',cmd,W/'logs/scdrs_third_L2.log',reserve=3)
  run_job('scdrs_third_L2_validation',[PY,validator,'group','--run',OUT,'--key',key],W/'logs/scdrs_third_L2_validation.log',reserve=1)
  schema=json.loads((OUT/'SCORE_VALIDATION.json').read_text())
  for key in old.legacy.OLD['group_keys']:schema[key]=json.loads((OUT/(key+'_VALIDATION.json')).read_text())
  schema['trait']='leave_one_CAD_r0.1';dump(OUT/'SCHEMA_VALIDATION.json',schema)
  dump(W/'status/scdrs_adopted.json',dict(status='THIRD_SCENARIO_COMPLETE_VALIDATED',original_L1_exit_code='UNAVAILABLE_ORPHAN_ADOPTED; outputs independently validated',score_recomputed=False,schema_sha256=sha(OUT/'SCHEMA_VALIDATION.json')))
if __name__=='__main__':
 try:main()
 except BaseException:
  dump(W/'status/scdrs_adoption_error.json',dict(status='FAILED',error=traceback.format_exc()));raise
