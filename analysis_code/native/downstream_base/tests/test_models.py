"""Real R backend tests with artificial counts; not atlas inference."""
import unittest,tempfile,json,subprocess,gzip,os
from pathlib import Path
import numpy as np,pandas as pd
from scipy import sparse,io
from src.core import atomic_json
ROOT=Path(__file__).resolve().parents[1]
class RModelTests(unittest.TestCase):
 def test_edger_case_reference_direction(self):
  with tempfile.TemporaryDirectory() as t:
   t=Path(t);rng=np.random.default_rng(7);n=20;g=250
   count=rng.negative_binomial(20,.5,(g,n));count[:30,10:]*=5
   with gzip.open(t/'x.mtx.gz','wb') as f:io.mmwrite(f,sparse.csr_matrix(count),field='integer')
   meta=pd.DataFrame({'pseudobulk_id':[f'p{i}' for i in range(n)],'donor_uid':[f'd{i}' for i in range(n)],'comparison_group':['reference']*10+['case']*10})
   meta.to_csv(t/'units.tsv',sep='\t',index=False);pd.DataFrame({'gene_id':[f'g{i}' for i in range(g)],'gene_symbol':[f'G{i}' for i in range(g)]}).to_csv(t/'genes.tsv',sep='\t',index=False)
   spec={'reference':'reference','case':'case','formula':'~ group','contrast_id':'synthetic','cell_state':'test','level':'L2','meta_context':'fake','min_tested_genes':100,'rows':json.loads(meta.to_json(orient='records'))};atomic_json(t/'spec.json',spec)
   cmd=['/usr/bin/Rscript',str(ROOT/'scripts/pseudobulk_edgeR.R'),str(t/'x.mtx.gz'),str(t/'units.tsv'),str(t/'genes.tsv'),str(t/'spec.json'),str(t/'result')]
   r=subprocess.run(cmd,capture_output=True,text=True);self.assertEqual(r.returncode,0,r.stdout+r.stderr)
   d=pd.read_csv(t/'result/DE.tsv',sep='\t');self.assertGreater(d.iloc[:30].logFC.median(),1.0);self.assertTrue(d.SE_kind.str.contains('APPROXIMATE').all())
 def test_edger_singular_design_stops(self):
  with tempfile.TemporaryDirectory() as t:
   t=Path(t);n=12;g=120;count=np.random.default_rng(2).poisson(20,(g,n));
   with gzip.open(t/'x.mtx.gz','wb') as f:io.mmwrite(f,sparse.csr_matrix(count),field='integer')
   m=pd.DataFrame({'pseudobulk_id':[f'p{i}' for i in range(n)],'donor_uid':[f'd{i}' for i in range(n)],'comparison_group':['reference']*6+['case']*6,'study':['a']*6+['b']*6});m.to_csv(t/'u.tsv',sep='\t',index=False);pd.DataFrame({'gene_id':[f'g{i}' for i in range(g)]}).to_csv(t/'g.tsv',sep='\t',index=False)
   atomic_json(t/'spec.json',{'reference':'reference','case':'case','formula':'~ study + group','rows':json.loads(m.to_json(orient='records'))})
   r=subprocess.run(['/usr/bin/Rscript',str(ROOT/'scripts/pseudobulk_edgeR.R'),str(t/'x.mtx.gz'),str(t/'u.tsv'),str(t/'g.tsv'),str(t/'spec.json'),str(t/'out')],capture_output=True,text=True)
   self.assertNotEqual(r.returncode,0);self.assertIn('Non-identifiable',r.stderr)
 def test_meta_approximation_runs_when_explicit(self):
  with tempfile.TemporaryDirectory() as t:
   t=Path(t);pd.DataFrame({'level':['L2']*3,'cell_state':['fake']*3,'gene_id':['G']*3,'meta_context':['A','B','C'],'logFC':[1.,1.2,.8],'legacy_SE_approx':[.2,.25,.22]}).to_csv(t/'m.tsv',sep='\t',index=False)
   atomic_json(t/'parameters.json',{'approve_legacy_ql_se_approximation':True,'meta_min_contexts':2,'strict_shared_min_contexts':3,'direction_consistency':.75,'fdr_threshold':.05,'context_synthesis_strategy':'one_approved_model_per_context'})
   r=subprocess.run(['/usr/bin/Rscript',str(ROOT/'scripts/cross_context_meta.R'),str(t/'m.tsv'),str(t/'out.tsv'),str(t/'parameters.json')],capture_output=True,text=True);self.assertEqual(r.returncode,0,r.stderr)
   d=pd.read_csv(t/'out.tsv',sep='\t');self.assertAlmostEqual(d.effect.iloc[0],1.,delta=.15)
class CompositionRTests(unittest.TestCase):
 def test_beta_binomial_positive_control(self):
  with tempfile.TemporaryDirectory() as t:
   t=Path(t);rng=np.random.default_rng(14);n=120;age=rng.uniform(25,75,n);az=(age-age.mean())/age.std();studies=np.repeat(np.arange(6),20);re=rng.normal(0,.6,6);p=1/(1+np.exp(-(-1+az*.7+re[studies])));prob=rng.beta(p*25,(1-p)*25);y=rng.binomial(200,prob)
   m=pd.DataFrame({'biological_sample_uid':[f'b{i}' for i in range(n)],'donor_uid':[f'd{i}' for i in range(n)],'study':studies.astype(str),'age':age,'comparison_group':['reference','case']*60,'n_state':y,'n_parent':200});m.to_csv(t/'m.tsv',sep='\t',index=False)
   spec={'approved':True,'age_center':float(age.mean()),'age_scale':float(age.std()),'formula':'cbind(n_state,n_parent-n_state) ~ age_z + comparison_group + (1|study)','tested_terms':['age_z'],'model_id':'fake','FDR_family':'test','state':'fake','level':'L1','min_denominator':100,'min_units':8,'min_donors':8,'reference_levels':{'comparison_group':'reference'}};atomic_json(t/'s.json',spec)
   r=subprocess.run(['/home/anzhen_vhd/miniconda3/envs/proportion_r/bin/Rscript',str(ROOT/'scripts/composition_glmm.R'),str(t/'m.tsv'),str(t/'s.json'),str(t/'out')],capture_output=True,text=True)
   self.assertEqual(r.returncode,0,r.stdout+r.stderr);d=pd.read_csv(t/'out/effects.tsv',sep='\t');self.assertGreater(d.estimate.iloc[0],.2)
if __name__=='__main__':unittest.main(verbosity=2)
