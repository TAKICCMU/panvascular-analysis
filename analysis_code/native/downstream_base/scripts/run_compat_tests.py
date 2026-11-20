from pathlib import Path
import sys,unittest,json,time,hashlib
sys.dont_write_bytecode=True
B=Path(__file__).resolve().parents[1];sys.path.insert(0,str(B))
names=['tests.test_release_compat','tests.test_core','tests.test_revision_interfaces.InterfaceRevisionTests.test_approval_binds_legal_cohort_contrast_model_and_transitive_gene_content','tests.test_revision_interfaces.InterfaceRevisionTests.test_relative_named_panels_and_resource_tables_bound','tests.test_revision_interfaces.InterfaceRevisionTests.test_missing_dependency_rejects']
suite=unittest.defaultTestLoader.loadTestsFromNames(names)
def flatten(s):
 for x in s:
  if isinstance(x,unittest.TestSuite):yield from flatten(x)
  else:yield x.id()
ids=list(flatten(suite));assert len(ids)==len(set(ids));start=time.time()
with (B/'logs/COMPAT_TESTS_FINAL.log').open('w') as log:r=unittest.TextTestRunner(stream=log,verbosity=2).run(suite)
failed={t.id():e for t,e in r.errors+r.failures};skips={t.id():e for t,e in r.skipped}
records=[{'test_id':x,'category':'NEW_SYNTHETIC_NATIVE' if x.startswith('tests.test_release_compat.') else 'HISTORICAL_RELEVANT_RERUN','status':'FAIL' if x in failed else 'SKIP' if x in skips else 'PASS','details':failed.get(x,skips.get(x,''))} for x in ids]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
out={'tests':records,'total':r.testsRun,'passed':sum(x['status']=='PASS' for x in records),'new_unique':sum(x['category']=='NEW_SYNTHETIC_NATIVE' for x in records),'historical_relevant_rerun':sum(x['category']=='HISTORICAL_RELEVANT_RERUN' for x in records),'historical_full_suite_count_not_rerun':123,'seconds':time.time()-start,'actual_models_on_real_data':False,'test_receipts':'Temporary synthetic authorization fixtures only, no production receipt','code_hashes':{str(p.relative_to(B)):sha(p) for d in ['src','scripts','tests'] for p in (B/d).glob('*.py')}}
(B/'evidence/COMPAT_TEST_RESULTS.json').write_text(json.dumps(out,ensure_ascii=False,indent=2));print({k:v for k,v in out.items() if k not in ['tests','code_hashes']},flush=True)
if not r.wasSuccessful():sys.exit(1)
