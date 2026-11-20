import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np,pandas as pd
from external_statistics import communication_effects,spatial_estimators,spatial_module_projection
from external_modules import render_spatial_adapter
class StatisticalTests(unittest.TestCase):
 def test_inherited_communication_effects(self):
  d=[];long=[]
  for study in ['s1','s2']:
   for i in range(8):
    sample=f'{study}_{i}';d.append(dict(sample_unit=sample,donor_uid=sample,donor_verified=True,study_canonical=study,tissue_v2='artery',disease_context='CAD',condition_class='control' if i<4 else 'disease'))
    long.append(dict(sample_unit=sample,interaction_key='A|B|L|R',source='A',target='B',ligand='L',receptor='R',score=float(i%4+(2 if i>=4 else 0)),threshold_scheme='main50',method='CellChat'))
  e,m=communication_effects(pd.DataFrame(long),pd.DataFrame(d));self.assertEqual(len(e),2);np.testing.assert_allclose(e.effect_disease_minus_control,2);self.assertAlmostEqual(m.iloc[0].pooled_effect_disease_minus_control,2)
  d[1]['donor_uid']=d[0]['donor_uid']
  with self.assertRaises(ValueError):communication_effects(pd.DataFrame(long),pd.DataFrame(d))
 def test_spatial_permutation_reproducible(self):
  ns=spatial_estimators(n_permutations=9);coords=np.column_stack([np.arange(20),np.zeros(20)]);x=np.arange(20,dtype=float)
  a=ns['moran_permutation'](x,coords,np.random.default_rng(0));b=ns['moran_permutation'](x,coords,np.random.default_rng(0));self.assertEqual(a,b);self.assertGreaterEqual(a[1],.1)
 def test_spatial_patch_contract(self):
  p=Path(__file__).resolve().parents[1]/'external_legacy/spatial/02_spatial_v4_cell2location_GPU.py'
  s=render_spatial_adapter(p.read_text(),'/tmp/no-atlas','/tmp/no-run','/tmp/no-manifest','level1_annotation','level2_annotation','study',[],'blocked')
  self.assertIn('USE_GPU = False',s);self.assertNotIn('atlas.raw.',s);self.assertNotIn('LOCKED71',s)
  self.assertIn("'level2_annotation'",s);compile(s,'new_spatial','exec')
 def test_dynamic_module_projection(self):
  a=pd.DataFrame([[1.,2.,4.],[2.,3.,5.]],columns=['A','B','C'])
  w=pd.DataFrame({'one':[1.,0.,0.],'two':[0.,1.,2.]},index=['A','B','C'])
  scores,audit=spatial_module_projection(a,w)
  log=np.log(a.to_numpy()+np.median(a.to_numpy())*.01);clr=log-log.mean(axis=1,keepdims=True)
  np.testing.assert_allclose(scores['one'],clr[:,0]);np.testing.assert_allclose(scores['two'],(clr[:,1]+2*clr[:,2])/3)
  np.testing.assert_allclose(audit.mapped_weight_fraction,1)
if __name__=='__main__':unittest.main()
