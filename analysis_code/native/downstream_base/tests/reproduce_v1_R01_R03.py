"""Independent reproduction from reviewer prose; original reviewer scripts absent."""
import sys,os,json,tempfile
from pathlib import Path
from unittest.mock import patch
from types import SimpleNamespace
sys.path.insert(0,sys.argv[1])
from src import dispatcher
from src.core import sha,atomic_json
from tests.test_core import CoreTests
from src.cohort import cohort_from_obs
from src.covarnet import prepare
import pandas as pd
results={'source_package':sys.argv[1],'evidence_type':'INDEPENDENT_REPRODUCTION_NOT_ORIGINAL_REVIEWER_SCRIPT'}
with tempfile.TemporaryDirectory() as td:
 r=Path(td);(r/'atlas.bin').write_bytes(b'synthetic');atomic_json(r/'manifest.json',{'release_status':'FINAL_APPROVED','level2_frozen':True,'ready_for_downstream':True});(r/'cohort.tsv').write_text('sample_uid\taggregate_approved\ns0\tFalse\n')
 c={'atlas':{'path':str(r/'atlas.bin'),'manifest':str(r/'manifest.json'),'input_sha256':sha(r/'atlas.bin')},'execution':{'mode':'approved_execution','run_full_analysis':True,'approval_receipt':str(r/'receipt.json'),'cpu_only':True,'max_workers':1,'per_worker_threads':20},'statistical_units':{'cohort_design_file':str(r/'cohort.tsv')}}
 cp=r/'config.json';atomic_json(cp,c);atomic_json(r/'receipt.json',{'approved':True,'atlas_sha256':c['atlas']['input_sha256'],'config_sha256':sha(cp),'code_sha256':dispatcher.code_digest(Path(sys.argv[1]))})
 dispatcher.authorization(c,cp);(r/'cohort.tsv').write_text('sample_uid\taggregate_approved\ns0\tTrue\n');dispatcher.authorization(c,cp)
 results['R01']={'reproduced':True,'old_receipt_accepted_after_legal_cohort_change':True}
 def fakeprocess(cmd,**kw):
  run=Path(cmd[cmd.index('--run-dir')+1]);module=cmd[cmd.index('--worker')+1];atomic_json(run/(module+'.receipt.json'),{'status':'COMPLETE_LEGACY_APPROXIMATION_EXPLICIT' if module=='meta' else 'COMPLETE'});return SimpleNamespace(returncode=0)
 c['outputs']={'root':str(r)}
 with patch.object(dispatcher,'authorization'),patch.object(dispatcher,'cache_identity',return_value=('abc',{})),patch('src.core.preflight',return_value=(None,None,{})),patch.object(dispatcher,'DAG',{'meta':[],'targets':['meta']}),patch.object(dispatcher.subprocess,'run',side_effect=fakeprocess):
  _,states=dispatcher.run(c,cp,selected=['targets'])
 results['R02']={'reproduced':states['meta']['status']=='PREPARED_NOT_EXECUTED','statuses':states}
 fixture=CoreTests();fixture.setUp()
 try:
  cfg=fixture.make_release();run=fixture.root/'run';run.mkdir();(run/'cohort').mkdir();cohort=cohort_from_obs(fixture.obs);cohort.to_csv(run/'cohort/cohort_design.tsv',sep='\t',index=False)
  ont=fixture.ont.loc[~fixture.ont.level2_annotation.isin(['Unknown','Not_assessed'])].copy();ont['covarnet_eligible']=True;ont.to_csv(fixture.root/'support/level2_ontology.tsv',sep='\t',index=False)
  atomic_json(run/'preflight.json',{});atomic_json(run/'RUN_IDENTITY.json',{})
  with patch.object(dispatcher,'cache_identity',return_value=('x',{})):dispatcher.worker(cfg,'counts',run)
  counts=pd.read_csv(run/'counts/sample_state_counts.tsv',sep='\t');outont=pd.read_csv(run/'counts/ontology.tsv',sep='\t');cohort['condition']=cohort.comparison_group
  try:prepare(counts,cohort,outont,{'min_profile_cells':1,'min_node_cells':1,'min_node_profiles':1,'min_node_studies':1});results['R03']={'reproduced':False}
  except ValueError as e:results['R03']={'reproduced':'Explicit boolean' in str(e),'error':str(e),'missing_eligibility':int(outont.covarnet_eligible.isna().sum())}
 finally:fixture.tearDown()
Path(sys.argv[2]).write_text(json.dumps(results,indent=2));print(json.dumps(results,indent=2))
assert all(results[k]['reproduced'] for k in ['R01','R02','R03'])
