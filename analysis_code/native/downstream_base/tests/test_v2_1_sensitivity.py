"""V2_S02 real control-flow regressions and native synthetic partial GLMM."""
import unittest,tempfile,json
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd
from src.core import atomic_json
from src.sensitivity import run_same_models,execute as sensitivity_execute
from src.composition import execute as composition_execute

class SensitivityFailurePaths(unittest.TestCase):
 def fixture(self,t,primary=True):
  p=t/'primary';s=t/'scope';(p/'counts').mkdir(parents=True);(s/'counts').mkdir(parents=True);(p/'composition_models').mkdir()
  pd.DataFrame({'biological_sample_uid':['b'],'n_cells':[2]}).to_csv(p/'counts/biological_state_counts.tsv',sep='\t',index=False)
  if primary:pd.DataFrame({'model_id':['A','B'],'term':['age','age'],'estimate':[1.,2.]}).to_csv(p/'composition_models/all_effects.tsv',sep='\t',index=False)
  atomic_json(p/'composition_models/MODEL_AUDIT.json',[{'model_id':'A','status':'COMPLETE'},{'model_id':'B','status':'COMPLETE'}])
  atomic_json(t/'models.json',[{'approved':True,'model_id':k,'tested_terms':['age']} for k in ['A','B']])
  cfg={'composition':{'models_file':str(t/'models.json')},'sensitivity':{'model_modules':['composition']}}
  return p,s,cfg
 def scope_output(self,s,mode):
  d=s/'composition_models';d.mkdir();ids=['A','B'] if mode=='success' else ['A'] if mode=='partial' else []
  pd.DataFrame({'model_id':ids,'term':['age']*len(ids),'estimate':[.5]*len(ids)}).to_csv(d/'all_effects.tsv',sep='\t',index=False)
  audit=[{'model_id':k,'status':'COMPLETE' if k in ids else 'BLOCKED','reason':'' if k in ids else 'Insufficient independent units after complete-case: donor support below threshold'} for k in ['A','B']]
  atomic_json(d/'MODEL_AUDIT.json',audit)
  return {'status':{'success':'COMPLETE','partial':'PARTIAL','blocked':'BLOCKED'}[mode],'models':audit}
 def check_mode(self,mode,expected):
  with tempfile.TemporaryDirectory() as z:
   p,s,cfg=self.fixture(Path(z))
   with patch('src.composition.execute',side_effect=lambda *args:self.scope_output(s,mode)):
    report=run_same_models(cfg,p,s,['b'])
   self.assertEqual(report[0]['status'],expected);d=pd.read_csv(s/'composition_primary_vs_scope.tsv',sep='\t');self.assertEqual(len(d),2)
   if mode!='success':
    b=d.set_index('model_id').loc['B'];self.assertEqual(b.comparison_status,'NOT_ESTIMABLE_IN_SCOPE');self.assertTrue(pd.isna(b.effect_difference));self.assertIn('Insufficient independent units',b.scope_reason)
   if mode!='blocked':self.assertEqual(d.set_index('model_id').loc['A','effect_difference'],-.5)
 def test_all_success(self):self.check_mode('success','COMPLETE')
 def test_partial_preserves_missing_reason(self):self.check_mode('partial','PARTIAL')
 def test_all_blocked_preserves_primary_rows(self):self.check_mode('blocked','BLOCKED')
 def test_missing_primary_keeps_scope_only(self):
  with tempfile.TemporaryDirectory() as z:
   p,s,cfg=self.fixture(Path(z),primary=False)
   with patch('src.composition.execute',side_effect=lambda *args:self.scope_output(s,'success')):report=run_same_models(cfg,p,s,['b'])
   self.assertEqual(report[0]['status'],'BLOCKED');d=pd.read_csv(s/'composition_primary_vs_scope.tsv',sep='\t');self.assertTrue(d.comparison_status.eq('NO_PRIMARY_EFFECT').all());self.assertTrue(d.primary_reason.eq('MISSING_PRIMARY_RESULTS').all())
 def test_program_exception_is_failed_not_unestimable(self):
  with tempfile.TemporaryDirectory() as z:
   p,s,cfg=self.fixture(Path(z))
   with patch('src.composition.execute',side_effect=RuntimeError('SYNTHETIC_PROGRAM_BUG')):report=run_same_models(cfg,p,s,['b'])
   self.assertEqual(report[0]['status'],'FAILED');self.assertEqual(report[0]['executor_receipt']['exception_type'],'RuntimeError')
   d=pd.read_csv(s/'composition_primary_vs_scope.tsv',sep='\t');self.assertEqual(len(d),2);self.assertTrue(d.scope_reason.eq('SYNTHETIC_PROGRAM_BUG').all());self.assertTrue(d.scope_failure_kind.eq('PROGRAM_EXCEPTION').all())
 def test_de_equivalent_partial_gap_preserved(self):
  with tempfile.TemporaryDirectory() as z:
   t=Path(z);p=t/'primary';s=t/'scope';p.mkdir();s.mkdir();(p/'de/de_000001').mkdir(parents=True)
   d=pd.DataFrame({'level':['L2','L2'],'cell_state':['A','B'],'gene_id':['g','g'],'contrast_id':['c','c'],'logFC':[1.,2.]});d.to_csv(p/'de/de_000001/DE.tsv',sep='\t',index=False)
   def fake(*args):
    out=s/'de/de_000001';out.mkdir(parents=True);d.iloc[:1].assign(logFC=.5).to_csv(out/'DE.tsv',sep='\t',index=False)
    audit=[{'level':'L2','cell_state':'A','contrast_id':'c','status':'COMPLETE'},{'level':'L2','cell_state':'B','contrast_id':'c','status':'BLOCKED','reason':'Insufficient independent donors: reference'}];atomic_json(s/'de/MODEL_AUDIT.json',audit);return {'status':'PARTIAL','models':audit}
   with patch('src.sensitivity.filter_pseudobulk'),patch('src.inference.execute_de',side_effect=fake):report=run_same_models({'sensitivity':{'model_modules':['de']}},p,s,[])
   self.assertEqual(report[0]['status'],'PARTIAL');out=pd.read_csv(s/'de_primary_vs_scope.tsv',sep='\t').set_index('cell_state');self.assertEqual(out.loc['A','effect_difference'],-.5);self.assertIn('Insufficient independent donors',out.loc['B','scope_reason']);self.assertTrue(pd.isna(out.loc['B','effect_difference']))

class DEExecutionAudit(unittest.TestCase):
 def execute_fixture(self,t,eligibility,exception=None):
  from src.inference import execute_de
  from types import SimpleNamespace
  def prepared(config,source,out):
   out=Path(out);out.mkdir();pd.DataFrame(eligibility).to_csv(out/'MODEL_ELIGIBILITY.tsv',sep='\t',index=False)
   return [] if all(x['status']=='NOT_SELECTED' for x in eligibility) else [{'id':'de_000001','argv':['synthetic'],'level':'L2','cell_state':'A','contrast_id':'c'}]
  with patch('src.inference.prepare_de_tasks',side_effect=prepared),patch('src.inference.subprocess.run',side_effect=exception,return_value=SimpleNamespace(returncode=0)):
   return execute_de({},t)
 def test_de_eligibility_gap_propagates_partial(self):
  with tempfile.TemporaryDirectory() as z:
   rows=[dict(level='L2',cell_state='A',contrast_id='c',status='READY_FOR_APPROVED_EXECUTION',reason=''),dict(level='L2',cell_state='B',contrast_id='c',status='BLOCKED_MODEL',reason='Insufficient complete independent donor pairs')]
   result=self.execute_fixture(Path(z),rows);self.assertEqual(result['status'],'PARTIAL');self.assertIn('Insufficient complete independent donor pairs',[x.get('reason','') for x in result['models']])
 def test_de_program_exception_not_model_block(self):
  with tempfile.TemporaryDirectory() as z:
   result=self.execute_fixture(Path(z),[dict(level='L2',cell_state='A',contrast_id='c',status='READY_FOR_APPROVED_EXECUTION',reason='')],OSError('native executable failure'))
   self.assertEqual(result['status'],'FAILED');self.assertEqual(result['tasks'][0]['exception_type'],'OSError');self.assertEqual(result['tasks'][0]['failure_kind'],'PROGRAM_EXCEPTION')
 def test_de_all_unapproved_remains_blocked(self):
  with tempfile.TemporaryDirectory() as z:
   result=self.execute_fixture(Path(z),[dict(level='L2',cell_state='A',contrast_id='c',status='NOT_SELECTED',reason='Contrast not human approved')]);self.assertEqual(result['status'],'BLOCKED');self.assertEqual(result['tasks'],[])

class NativePartialComposition(unittest.TestCase):
 def test_native_glmm_partial_keeps_comparison_and_status(self):
  with tempfile.TemporaryDirectory() as z:
   t=Path(z);(t/'cohort').mkdir();(t/'counts').mkdir();rng=np.random.default_rng(14);n=120
   age=rng.uniform(25,75,n);az=(age-age.mean())/age.std();st=np.repeat(np.arange(6),20);re=rng.normal(0,.6,6);p=1/(1+np.exp(-(-1+az*.7+re[st])));prob=rng.beta(p*25,(1-p)*25);y=rng.binomial(200,prob)
   m=pd.DataFrame({'sample_uid':[f's{i}' for i in range(n)],'biological_sample_uid':[f'b{i}' for i in range(n)],'donor_uid':[f'd{i}' for i in range(n)],'study':[f'S{i}' for i in st],'vascular_bed':['v']*n,'assay_type':np.where(st==5,'snRNA-seq','scRNA-seq'),'age':age,'comparison_group':['reference','case']*60,'donor_resolved':[True]*n,'inference_eligible':[True]*n});m.to_csv(t/'cohort/cohort_design.tsv',sep='\t',index=False)
   d=pd.DataFrame({'biological_sample_uid':m.biological_sample_uid,'level1_annotation':['synthetic']*n,'level2_annotation':['synthetic_state']*n,'n_cells':y,'n_all_qc':[200]*n,'n_parent':[200]*n});d.to_csv(t/'counts/biological_state_counts.tsv',sep='\t',index=False);d.assign(sample_uid=m.sample_uid).to_csv(t/'counts/sample_state_counts.tsv',sep='\t',index=False)
   base={'approved':True,'age_center':float(age.mean()),'age_scale':float(age.std()),'formula':'cbind(n_state,n_parent-n_state) ~ age_z + comparison_group + (1|study)','tested_terms':['age_z'],'FDR_family':'synthetic','state':'synthetic','level':'L1','min_denominator':100,'min_units':8,'reference_levels':{'comparison_group':'reference'}}
   atomic_json(t/'models.json',[dict(base,model_id='A',min_donors=8),dict(base,model_id='B',min_donors=110)])
   cfg={'composition':{'rscript':'/home/anzhen_vhd/miniconda3/envs/proportion_r/bin/Rscript','models_file':str(t/'models.json')},'sensitivity':{'run':True,'model_modules':['composition'],'assay_column':'assay_type','scRNA_only_values':['scRNA-seq'],'exclude_studies':['S5']}}
   primary=composition_execute(cfg,t);self.assertEqual(primary['status'],'COMPLETE');receipt=sensitivity_execute(cfg,t)
   self.assertEqual(receipt['status'],'PARTIAL');self.assertEqual(receipt['scopes'][0]['status'],'PARTIAL');self.assertEqual(receipt['scopes'][1]['status'],'PARTIAL');self.assertEqual(receipt['scopes'][1]['phase'],'REUSED_IDENTICAL_SCOPE')
   out=pd.read_csv(t/'sensitivity/scRNA_only/composition_primary_vs_scope.tsv',sep='\t').set_index('model_id');self.assertEqual(out.loc['A','comparison_status'],'COMPARABLE');self.assertTrue(np.isfinite(out.loc['A','effect_difference']));self.assertTrue(pd.isna(out.loc['B','effect_difference']));self.assertIn('Insufficient independent units',out.loc['B','scope_reason'])
   self.assertFalse((t/'sensitivity/exclude_target_study/composition_models').exists())
   print(json.dumps({'status':'PASS_NATIVE_SYNTHETIC_PARTIAL','n_primary':120,'n_scope':100,'primary_status':primary['status'],'scope_status':receipt['status'],'retained_missing_reason':out.loc['B','scope_reason'],'same_scope_not_refit':True,'real_atlas_used':False}),flush=True)
if __name__=='__main__':unittest.main(verbosity=2)
