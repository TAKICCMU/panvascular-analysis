"""Auditable suite runner; only synthetic fixtures and static checks."""
from pathlib import Path
import unittest,sys,json,datetime,subprocess,os,importlib.metadata,ast
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root))
class AuditResult(unittest.TextTestResult):
 def __init__(self,*a,**k):super().__init__(*a,**k);self.records=[]
 def addSuccess(self,test):super().addSuccess(test);self.records.append({'test':test.id(),'status':'PASS','scope':'synthetic'})
 def addError(self,test,err):super().addError(test,err);self.records.append({'test':test.id(),'status':'ERROR','detail':self._exc_info_to_string(err,test)})
 def addFailure(self,test,err):super().addFailure(test,err);self.records.append({'test':test.id(),'status':'FAIL','detail':self._exc_info_to_string(err,test)})
 def addSkip(self,test,reason):super().addSkip(test,reason);self.records.append({'test':test.id(),'status':'NOT_RUN','reason':reason})
suite=unittest.TestSuite()
for module in ['tests.test_core','tests.test_models','tests.test_covarnet','tests.test_external_modules','tests.test_external_statistics','tests.test_revision_interfaces','tests.test_revision_covarnet','tests.test_revision_spatial','tests.test_revision_statistics','tests.test_external_gate','tests.test_revision_composition_scope','tests.test_v2_1_projection','tests.test_v2_1_sensitivity']:
 suite.addTests(unittest.defaultTestLoader.loadTestsFromName(module))
with open(root/'tests/CORE_TEST_LOG.txt','w') as log:r=unittest.TextTestRunner(stream=log,verbosity=2,resultclass=AuditResult).run(suite)
records=r.records
for folder in ['src','scripts','tests']:
 for p in (root/folder).rglob('*.py'):
  try:ast.parse(p.read_text());records.append({'test':'syntax:'+str(p.relative_to(root)),'status':'PASS','scope':'static_only'})
  except SyntaxError as e:records.append({'test':'syntax:'+str(p.relative_to(root)),'status':'FAIL','detail':str(e)})
for p in list((root/'scripts').glob('*.R'))+list((root/'src').glob('*.R')):
 code='parse(file='+json.dumps(str(p))+');cat("PARSE_OK\\n")'
 rr=subprocess.run(['/usr/bin/Rscript','-e',code],capture_output=True,text=True)
 records.append({'test':'R_parse:'+str(p.relative_to(root)),'status':'PASS' if rr.returncode==0 else 'FAIL','scope':'static_only','detail':rr.stderr[-2000:]})
rr=subprocess.run([sys.executable,str(root/'scripts/run_downstream.py'),'--config',str(root/'configs/analysis_config.yaml')],capture_output=True,text=True)
records.append({'test':'default_entrypoint_does_not_execute','status':'PASS' if rr.returncode==0 and any(s in rr.stdout for s in ['WAITING_FOR_FINAL_L2_ATLAS','CODE_READY_FOR_REVIEW_NOT_EXECUTED']) and json.loads(rr.stdout).get('executed') is False else 'FAIL','scope':'read_only_status','stdout':rr.stdout,'stderr':rr.stderr[-2000:]})
for module in ['final_release_full_preflight','full_scale_pseudobulk','real_cohort_design_fit','full_scale_CoVarNet','MAGMA_scDRS_native_full','CellChat_LIANA_native_full','spatial_cell2location_full','full_GWAS_QTL_revalidation','full_figures2to5']:
 records.append({'test':module,'status':'NOT_RUN','reason':'Final atlas/execution approval or external resource review missing; this turn code and synthetic tests only'})
versions={}
for name in ['omicverse','scanpy','anndata','numpy','pandas','scipy','h5py','statsmodels','PyYAML']:
 try:versions[name]=importlib.metadata.version(name)
 except importlib.metadata.PackageNotFoundError:versions[name]='NOT_INSTALLED'
report={'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'python':sys.executable,'versions':versions,'records':records,'synthetic_test_count':r.testsRun,'synthetic_success':r.wasSuccessful(),'no_real_atlas_analysis':True,'static_is_not_native_validation':True,'additional_reports':['v2_1_sensitivity_results/TEST_RESULTS.json','v2_1_baseline_reproduction/targeted_findings.json','revision_logs/covarnet_server_native_tests_final.log','revision_logs/covarnet_native_nmf_final.log','revision_evidence/BASELINE_R04_R05.json','BASELINE_R01_R03.json','revision_statistics_results/STATISTICS_TEST_FINAL.log']}
(root/'tests/TEST_RESULTS.json').write_text(json.dumps(report,indent=2,ensure_ascii=False))
print(json.dumps({'synthetic_test_count':r.testsRun,'synthetic_success':r.wasSuccessful(),'total_failures':sum(x['status'] in ['FAIL','ERROR'] for x in records),'n_records':len(records)}))
sys.exit(0 if all(x['status'] not in ['FAIL','ERROR'] for x in records) else 1)
