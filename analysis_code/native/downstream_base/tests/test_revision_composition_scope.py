"""Native beta-binomial scope refit, synthetic only, no atlas input.

The approved synthetic baseline age center/scale is fixed across scopes.
"""
import unittest,tempfile,json,subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from src.core import atomic_json
from src.composition import execute as composition_execute
from src.sensitivity import execute as sensitivity_execute

class RevisionCompositionScope(unittest.TestCase):
 def test_native_missing_age_transform_is_blocked(self):
  with tempfile.TemporaryDirectory() as z:
   t=Path(z);m=pd.DataFrame({'biological_sample_uid':[f'b{i}' for i in range(12)],'donor_uid':[f'd{i}' for i in range(12)],'age':np.arange(30,42),'study':['A']*6+['B']*6,'comparison_group':['reference','case']*6,'n_state':[20]*12,'n_parent':[200]*12});m.to_csv(t/'m.tsv',sep='\t',index=False)
   atomic_json(t/'s.json',{'approved':True,'formula':'cbind(n_state,n_parent-n_state) ~ age_z + comparison_group + (1|study)'})
   root=Path(__file__).resolve().parents[1]
   p=subprocess.run(['/home/anzhen_vhd/miniconda3/envs/proportion_r/bin/Rscript',str(root/'scripts/composition_glmm.R'),str(t/'m.tsv'),str(t/'s.json'),str(t/'out')],capture_output=True,text=True)
   self.assertNotEqual(p.returncode,0);self.assertIn('age_z requires explicitly approved',p.stderr)
 def test_native_same_formula_scope_refit_and_dedup(self):
  with tempfile.TemporaryDirectory() as z:
   t=Path(z);(t/'cohort').mkdir();(t/'counts').mkdir()
   rng=np.random.default_rng(14);n=120;age=rng.uniform(25,75,n);az=(age-age.mean())/age.std();study=np.repeat(np.arange(6),20);re=rng.normal(0,.6,6)
   p=1/(1+np.exp(-(-1+az*.7+re[study])));prob=rng.beta(p*25,(1-p)*25);y=rng.binomial(200,prob)
   meta=pd.DataFrame({'sample_uid':[f's{i}' for i in range(n)],'biological_sample_uid':[f'b{i}' for i in range(n)],'donor_uid':[f'd{i}' for i in range(n)],'study':[f'S{v}' for v in study],'vascular_bed':['bed']*n,'assay_type':np.where(study==5,'snRNA-seq','scRNA-seq'),'age':age,'comparison_group':['reference','case']*60,'donor_resolved':[True]*n,'inference_eligible':[True]*n})
   meta.to_csv(t/'cohort/cohort_design.tsv',sep='\t',index=False)
   counts=pd.DataFrame({'biological_sample_uid':meta.biological_sample_uid,'level1_annotation':['synthetic']*n,'level2_annotation':['synthetic_state']*n,'n_cells':y,'n_all_qc':[200]*n,'n_parent':[200]*n})
   counts.to_csv(t/'counts/biological_state_counts.tsv',sep='\t',index=False)
   sample=counts.copy();sample.insert(0,'sample_uid',meta.sample_uid);sample.to_csv(t/'counts/sample_state_counts.tsv',sep='\t',index=False)
   spec={'approved':True,'age_center':float(age.mean()),'age_scale':float(age.std()),'formula':'cbind(n_state,n_parent-n_state) ~ age_z + comparison_group + (1|study)','tested_terms':['age_z'],'model_id':'synthetic_age','FDR_family':'synthetic','state':'synthetic','level':'L1','min_denominator':100,'min_units':8,'min_donors':8,'reference_levels':{'comparison_group':'reference'}}
   atomic_json(t/'models.json',[spec]);cfg={'composition':{'rscript':'/home/anzhen_vhd/miniconda3/envs/proportion_r/bin/Rscript','models_file':str(t/'models.json')},'sensitivity':{'run':True,'model_modules':['composition'],'assay_column':'assay_type','scRNA_only_values':['scRNA-seq'],'exclude_studies':['S5']}}
   primary=composition_execute(cfg,t);self.assertEqual(primary['status'],'COMPLETE')
   report=sensitivity_execute(cfg,t);self.assertEqual(report['status'],'COMPLETE')
   source=t/'composition_models';scope=t/'sensitivity/scRNA_only/composition_models'
   pmodel=json.loads((source/'synthetic_age/MODEL.json').read_text());smodel=json.loads((scope/'synthetic_age/MODEL.json').read_text())
   self.assertEqual(pmodel['formula'],spec['formula']);self.assertEqual(smodel['formula'],spec['formula']);self.assertEqual(pmodel['n_after'],120);self.assertEqual(smodel['n_after'],100)
   pdsgn=json.loads((source/'synthetic_age.json').read_text());sdsgn=json.loads((scope/'synthetic_age.json').read_text());self.assertEqual(pdsgn,sdsgn)
   result=pd.read_csv(t/'sensitivity/scRNA_only/composition_primary_vs_scope.tsv',sep='\t');self.assertEqual(len(result),1);self.assertEqual(result.term.iloc[0],'age_z');self.assertEqual(result.comparison_status.iloc[0],'COMPARABLE');self.assertGreater(result.estimate_primary.iloc[0],.01);self.assertGreater(result.estimate_sensitivity.iloc[0],.01)
   pin=pd.read_csv(source/'synthetic_age/included_units.tsv',sep='\t').set_index('biological_sample_uid');sin=pd.read_csv(scope/'synthetic_age/included_units.tsv',sep='\t').set_index('biological_sample_uid');np.testing.assert_allclose(pin.loc[sin.index,'age_z'],sin.age_z,rtol=0,atol=1e-12);self.assertEqual(pmodel['age_transform'],smodel['age_transform'])
   self.assertTrue(np.isfinite(result.effect_difference.iloc[0]));self.assertEqual(report['scopes'][1]['status'],'COMPLETE');self.assertEqual(report['scopes'][1]['phase'],'REUSED_IDENTICAL_SCOPE');self.assertFalse((t/'sensitivity/exclude_target_study/composition_models').exists())
   print(json.dumps({'test':'NATIVE_SYNTHETIC_COMPOSITION_SCOPE','n_primary':120,'n_scope':100,'studies_primary':6,'studies_scope':5,'same_model_spec':pdsgn==sdsgn,'formula':spec['formula'],'effect_scale':'log odds per fixed baseline age SD','primary_age_effect':float(result.estimate_primary.iloc[0]),'scope_age_effect':float(result.estimate_sensitivity.iloc[0]),'equivalent_scope_not_rerun':True,'production_data_used':False}),flush=True)
if __name__=='__main__':unittest.main(verbosity=2)
