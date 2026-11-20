"""Native synthetic revision tests. No production input or approvals changed."""
import unittest,tempfile,json,subprocess,os
from pathlib import Path
import numpy as np,pandas as pd
from scipy import sparse
from src.core import Blocked,atomic_json
from src.inference import meta_parameters,contrast_feasibility,execute_de
from src.sensitivity import complete_biosample_scope,compare_effects,execute
ROOT=Path(__file__).resolve().parents[1]
class RevisionStatistics(unittest.TestCase):
 def test_se_gate_and_threshold_validation(self):
  p={'meta_min_contexts':2,'strict_shared_min_contexts':3,'direction_consistency':.75,'approve_legacy_ql_se_approximation':False}
  with self.assertRaises(Blocked):meta_parameters({'pseudobulk':p})
  p['approve_legacy_ql_se_approximation']=True
  self.assertEqual(meta_parameters({'pseudobulk':p})['meta_min_contexts'],2)
  p['strict_shared_min_contexts']=1
  with self.assertRaises(Blocked):meta_parameters({'pseudobulk':p})
 def test_scope_partial_technical_unit_refused(self):
  m=pd.DataFrame({'sample_uid':['a','b'],'biological_sample_uid':['same','same']})
  with self.assertRaises(Blocked):complete_biosample_scope(m,m.iloc[:1])
 def test_comparison_keeps_missing(self):
  a=pd.DataFrame({'key':['a','b'],'logFC':[1.,2.]});b=pd.DataFrame({'key':['a'],'logFC':[.5]})
  d=compare_effects(a,b,['key'],'logFC')
  self.assertEqual(d.loc[d.key=='b','comparison_status'].iloc[0],'NOT_ESTIMABLE_IN_SCOPE')
 def test_actual_current_contrast_upper_bound(self):
  m=pd.read_csv(ROOT/'configs/cohort_design.tsv',sep='\t');c=pd.read_csv(ROOT/'configs/contrasts.tsv',sep='\t',keep_default_na=False)
  f=contrast_feasibility(m,c)
  self.assertEqual(len(f),38);self.assertEqual(int(f.both_arms_nonempty.sum()),11);self.assertEqual(list(f.loc[f.passes_donor_upper_bound,'contrast_id']),['contrast_013'])
  self.assertFalse(f.approved.any())
 def test_native_meta_configuration_applies(self):
  with tempfile.TemporaryDirectory() as z:
   t=Path(z);rows=[]
   for gene,effects in [('two',[1.,1.1]),('three',[1.,1.1,.9]),('discordant',[1.,-1.,1.]),('direction_boundary',[1.,1.1,.9,-.01])]:
    for i,e in enumerate(effects):rows.append(dict(level='L2',cell_state='S',gene_id=gene,meta_context=str(i),logFC=e,legacy_SE_approx=10. if gene=='direction_boundary' and i==3 else .1))
   pd.DataFrame(rows).to_csv(t/'in.tsv',sep='\t',index=False)
   p={'pseudobulk':dict(meta_min_contexts=2,strict_shared_min_contexts=3,direction_consistency=.75,approve_legacy_ql_se_approximation=True)}
   atomic_json(t/'p.json',meta_parameters(p))
   r=subprocess.run(['/usr/bin/Rscript',str(ROOT/'scripts/cross_context_meta.R'),str(t/'in.tsv'),str(t/'out.tsv'),str(t/'p.json')],capture_output=True,text=True)
   self.assertEqual(r.returncode,0,r.stderr)
   d=pd.read_csv(t/'out.tsv',sep='\t').set_index('gene_id');self.assertEqual(d.loc['two','status'],'COMPLETE');self.assertTrue(d.loc['two','passes_main_definition']);self.assertFalse(d.loc['two','passes_strict_definition']);self.assertTrue(d.loc['three','passes_strict_definition']);self.assertFalse(d.loc['discordant','passes_main_definition'])
   p['pseudobulk']['meta_min_contexts']=3;atomic_json(t/'p3.json',meta_parameters(p))
   r=subprocess.run(['/usr/bin/Rscript',str(ROOT/'scripts/cross_context_meta.R'),str(t/'in.tsv'),str(t/'out3.tsv'),str(t/'p3.json')],capture_output=True,text=True)
   self.assertEqual(r.returncode,0,r.stderr);d=pd.read_csv(t/'out3.tsv',sep='\t').set_index('gene_id');self.assertEqual(d.loc['two','status'],'NOT_TESTED_INSUFFICIENT_CONTEXTS')
   p['pseudobulk']['direction_consistency']=.8;atomic_json(t/'p8.json',meta_parameters(p))
   r=subprocess.run(['/usr/bin/Rscript',str(ROOT/'scripts/cross_context_meta.R'),str(t/'in.tsv'),str(t/'out8.tsv'),str(t/'p8.json')],capture_output=True,text=True);self.assertEqual(r.returncode,0,r.stderr)
   self.assertTrue(d.loc['direction_boundary','passes_main_definition']);d8=pd.read_csv(t/'out8.tsv',sep='\t').set_index('gene_id');self.assertFalse(d8.loc['direction_boundary','passes_main_definition'])
 def test_nmf_sensitivity_matching_is_explicit(self):
  from src.sensitivity import compare_nmf_outputs
  with tempfile.TemporaryDirectory() as z:
   t=Path(z);a=t/'primary/covarnet/stable_NMF';b=t/'scope/covarnet/stable_NMF';a.mkdir(parents=True);b.mkdir(parents=True)
   pd.DataFrame({'A':[.8,.1,.1],'B':[.7,.2,.1]},index=['s1','s2','s3']).to_csv(a/'module_node_weights_K02.tsv',sep='\t')
   pd.DataFrame({'X':[.9,.1],'Y':[.1,.9]},index=['s1','s2']).to_csv(b/'module_node_weights_K02.tsv',sep='\t')
   d=compare_nmf_outputs(t/'primary',t/'scope');self.assertTrue(d.many_to_one_duplicate.all());self.assertEqual(d.one_to_one_module.nunique(),2);self.assertTrue(d.evidence_status.eq('LIMITED_MISSING_WEIGHT').all())
 def test_native_same_model_sensitivity_and_scope_dedup(self):
  with tempfile.TemporaryDirectory() as z:
   t=Path(z);(t/'cohort').mkdir();(t/'counts').mkdir();(t/'pseudobulk').mkdir()
   n=24;g=140;rng=np.random.default_rng(76)
   m=pd.DataFrame({'sample_uid':[f's{i}' for i in range(n)],'biological_sample_uid':[f'b{i}' for i in range(n)],'donor_uid':[f'd{i}' for i in range(n)],'study':['A']*16+['B']*8,'vascular_bed':['v']*n,'assay_type':['scRNA-seq']*16+['snRNA-seq']*8,'cohort_disease_context':['disease']*n,'comparison_group':['reference','case']*12,'inference_eligible':[True]*n,'donor_resolved':[True]*n,'n_cells':[40]*n})
   m.to_csv(t/'cohort/cohort_design.tsv',sep='\t',index=False);m[['sample_uid']].assign(n_cells=40).to_csv(t/'counts/sample_state_counts.tsv',sep='\t',index=False)
   for level in ['L1','L2']:
    d=t/'pseudobulk'/level;d.mkdir();u=m.drop(columns='sample_uid').assign(pseudobulk_id=[f'{level}{i}' for i in range(n)],cell_state='state');u.to_csv(d/'units.tsv',sep='\t',index=False)
    count=rng.negative_binomial(20,.5,(n,g));count[1::2,:20]*=4;sparse.save_npz(d/'counts_units_by_genes.npz',sparse.csr_matrix(count));pd.DataFrame({'gene_id':[f'g{i}' for i in range(g)],'gene_symbol':[f'G{i}' for i in range(g)]}).to_csv(d/'genes.tsv',sep='\t',index=False)
   spec={'contrast_id':'syn','study':'','vascular_bed':'v','assay_type':'','cohort_disease_context':'disease','case':'case','reference':'reference','min_independent_per_arm':5,'paired':False,'approved':True,'formula':'~ group','covariates':'','meta_context':'disease'}
   pd.DataFrame([spec]).to_csv(t/'contrasts.tsv',sep='\t',index=False)
   cfg={'statistical_units':{'contrasts_file':str(t/'contrasts.tsv')},'pseudobulk':{'min_cells':20,'min_tested_genes':100,'rscript':'/usr/bin/Rscript'},'sensitivity':{'run':True,'model_modules':['de'],'assay_column':'assay_type','scRNA_only_values':['scRNA-seq'],'exclude_studies':['B']}}
   execute_de(cfg,t);receipt=execute(cfg,t);self.assertEqual(receipt['status'],'COMPLETE')
   self.assertEqual(receipt['scopes'][1]['status'],'COMPLETE');self.assertEqual(receipt['scopes'][1]['phase'],'REUSED_IDENTICAL_SCOPE');self.assertFalse((t/'sensitivity/exclude_target_study/de').exists())
   cmp=pd.read_csv(t/'sensitivity/scRNA_only/de_primary_vs_scope.tsv',sep='\t');self.assertTrue(cmp.comparison_status.eq('COMPARABLE').all());self.assertEqual(len(cmp),2*g)
   first=json.loads((t/'de/de_000000.json').read_text());second=json.loads((t/'sensitivity/scRNA_only/de/de_000000.json').read_text());self.assertEqual(first['formula'],second['formula']);self.assertEqual(first['min_tested_genes'],second['min_tested_genes']);self.assertEqual(len(first['rows']),24);self.assertEqual(len(second['rows']),16)
if __name__=='__main__':unittest.main(verbosity=2)
