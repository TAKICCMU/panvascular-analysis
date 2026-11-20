"""Real dispatcher/worker tests with synthetic release; no real atlas operations."""
import unittest,tempfile,json,os,copy
from pathlib import Path
from unittest.mock import patch
from types import SimpleNamespace
import pandas as pd,numpy as np
from src import dispatcher
from src.core import sha,atomic_json,Blocked,preflight
from src.review_dependencies import reviewed_dependencies
from src.status import receipt_state
from tests import test_core

class InterfaceRevisionTests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
 def tearDown(self):self.tmp.cleanup()
 def approved_config(self):
  r=self.root;(r/'atlas.bin').write_bytes(b'synthetic only')
  atomic_json(r/'manifest.json',{'release_status':'FINAL_APPROVED','level2_frozen':True,'ready_for_downstream':True})
  (r/'cohort.tsv').write_text('sample_uid\taggregate_approved\ns0\tFalse\n')
  (r/'contrasts.tsv').write_text('contrast_id\tformula\nx\t~ group\n')
  atomic_json(r/'models.json',[{'approved':False,'formula':'~ group'}]);(r/'genes.tsv').write_text('gene_id\tsymbol\ng1\tG\n')
  atomic_json(r/'resource.json',{'gene_map':str(r/'genes.tsv')})
  c={'atlas':{'path':str(r/'atlas.bin'),'manifest':str(r/'manifest.json'),'input_sha256':sha(r/'atlas.bin')},'execution':{'mode':'approved_execution','run_full_analysis':True,'approval_receipt':str(r/'receipt.json'),'cpu_only':True,'max_workers':1,'per_worker_threads':4},'statistical_units':{'cohort_design_file':str(r/'cohort.tsv'),'contrasts_file':str(r/'contrasts.tsv')},'composition':{'models_file':str(r/'models.json')},'covarnet':{'threads_per_worker':4},'external':{'manifest':str(r/'resource.json')}}
  cp=r/'config.json';atomic_json(cp,c)
  dg,payload=reviewed_dependencies(c,cp)
  atomic_json(r/'receipt.json',{'approved':True,'atlas_sha256':c['atlas']['input_sha256'],'config_sha256':sha(cp),'code_sha256':dispatcher.code_digest(Path(dispatcher.__file__).resolve().parents[1]),'reviewed_dependency_digest':dg})
  return c,cp
 def test_approval_binds_legal_cohort_contrast_model_and_transitive_gene_content(self):
  c,cp=self.approved_config();dispatcher.authorization(c,cp)
  for filename,a,b in [('cohort.tsv','False','True'),('contrasts.tsv','~ group','~ age + group'),('models.json','false','true'),('genes.tsv','\tG','\tH')]:
   p=self.root/filename;old=p.read_text();p.write_text(old.replace(a,b))
   with self.assertRaisesRegex(Blocked,'dependency digest'):dispatcher.authorization(c,cp)
   p.write_text(old);dispatcher.authorization(c,cp)
 def test_relative_named_panels_and_resource_tables_bound(self):
  c,cp=self.approved_config();(self.root/'anchors.tsv').write_text('gene_id\nG1\n')
  c['figures']={'l2_anchor_panel':'anchors.tsv','marker_panel':'genes.tsv'}
  c['external']['evidence_tables']={'arbitrary_role':'anchors.tsv'}
  first,_=reviewed_dependencies(c,cp);(self.root/'anchors.tsv').write_text('gene_id\nG2\n');second,_=reviewed_dependencies(c,cp)
  self.assertNotEqual(first,second)
 def test_missing_dependency_rejects(self):
  c,cp=self.approved_config();(self.root/'genes.tsv').unlink()
  with self.assertRaisesRegex(Blocked,'missing'):dispatcher.authorization(c,cp)
 def dispatch(self,reported,returncode=0):
  c={'outputs':{'root':str(self.root)},'execution':{'per_worker_threads':1,'max_workers':1}}
  atomic_json(self.root/'config',c)
  def fakeprocess(cmd,**kw):
   run=Path(cmd[cmd.index('--run-dir')+1]);module=cmd[cmd.index('--worker')+1]
   atomic_json(run/(module+'.receipt.json'),{'status':reported if module=='meta' else 'COMPLETE','method':'REML','SE_kind':'APPROXIMATE_FROM_QL_F_NOT_COEFFICIENT_SE'})
   return SimpleNamespace(returncode=returncode if module=='meta' else 0)
  with patch.object(dispatcher,'authorization'),patch.object(dispatcher,'cache_identity',return_value=('abc',{})),patch('src.core.preflight',return_value=(None,None,{})),patch.object(dispatcher,'DAG',{'meta':[],'targets':['meta']}),patch.object(dispatcher.subprocess,'run',side_effect=fakeprocess):
   return dispatcher.run(c,self.root/'config',selected=['targets'])[1]
 def test_real_dispatcher_meta_complete_unblocks_targets(self):
  r=self.dispatch('COMPLETE');self.assertEqual(r['meta']['status'],'COMPLETE');self.assertEqual(r['targets']['status'],'COMPLETE')
 def test_real_dispatcher_prepared_and_failed_block_targets(self):
  for state,rc in [('PREPARED_NOT_EXECUTED',0),('FAILED',0),('COMPLETE',1),('BLOCKED',3)]:
   r=self.dispatch(state,rc);self.assertNotEqual(r['meta']['status'],'COMPLETE');self.assertEqual(r['targets']['status'],'BLOCKED_DEPENDENCY')
 def test_unknown_or_absent_receipt_never_success(self):
  for r in [None,{}, {'status':'COMPLETE_LEGACY_APPROXIMATION_EXPLICIT'}]:self.assertEqual(receipt_state(0,r),'FAILED')
 def test_native_release_counts_worker_prepare_spatial_map(self):
  from src.cohort import cohort_from_obs
  from src.covarnet import prepare
  from src.state_identity import map_state_axis
  from src.external_statistics import spatial_module_projection
  fixture=test_core.CoreTests();fixture.setUp()
  try:
   cfg=fixture.make_release();root=fixture.root;cp=root/'config.json';run=root/'run';run.mkdir();(run/'cohort').mkdir()
   # Production named ontology has boolean eligibility; reserved rows are added by real counts worker.
   ont=fixture.ont.loc[~fixture.ont.level2_annotation.isin(['Unknown','Not_assessed'])].copy();ont['covarnet_eligible']=True;ont['level2_state_id']=['STATE_'+str(i) for i in range(len(ont))];ont.to_csv(root/'support/level2_ontology.tsv',sep='\t',index=False)
   manifest=json.load(open(root/'manifest.json'));manifest['release_sidecar_sha256']['level2_ontology.tsv']=sha(root/'support/level2_ontology.tsv');atomic_json(root/'manifest.json',manifest)
   cohort=cohort_from_obs(fixture.obs);cohort.to_csv(run/'cohort/cohort_design.tsv',sep='\t',index=False)
   cfg['execution']={'block_rows':3};cfg['statistical_units']={};atomic_json(cp,cfg)
   preflight(cfg,run/'preflight.json');_,payload=dispatcher.cache_identity(cfg,cp);atomic_json(run/'RUN_IDENTITY.json',payload)
   with patch.dict(os.environ,{'DOWNSTREAM_CONFIG_PATH':str(cp)}):result=dispatcher.worker(cfg,'counts',run)
   self.assertEqual(result['status'],'COMPLETE')
   outont=pd.read_csv(run/'counts/ontology.tsv',sep='\t');special=outont.level2_annotation.isin(['Unknown','Not_assessed']);self.assertFalse(outont.loc[special,'covarnet_eligible'].any());self.assertFalse(outont.covarnet_eligible.isna().any())
   counts=pd.read_csv(run/'counts/sample_state_counts.tsv',sep='\t');self.assertEqual(int(counts['count'].sum()),24)
   cohort['condition']=cohort.comparison_group
   prepared=prepare(counts,cohort,outont,{'min_profile_cells':1,'min_node_cells':1,'min_node_profiles':1,'min_node_studies':1})
   self.assertEqual(prepared['manifest']['n_cells'],24);self.assertEqual(set(prepared['counts'].index),{'STATE_0','STATE_1'})
   state_map=prepared['state_identity_map'];ids=prepared['counts'].index
   safe=map_state_axis(ids,state_map,'state_id','safe_label');self.assertEqual(list(map_state_axis(safe,state_map,'safe_label')),list(ids))
   weights=pd.DataFrame({'module':[.6,.4]},index=ids)
   abundance=pd.DataFrame([[2.,3.]],columns=safe)
   projection,audit=spatial_module_projection(abundance,weights,state_map=state_map,abundance_key='safe_label')
   self.assertEqual(projection.shape,(1,1))
  finally:fixture.tearDown()
 def test_l2_within_parent_removes_shared_lineage_marker(self):
  import anndata as ad
  from scipy import sparse
  from src.expression import marker_summary
  obs=pd.DataFrame({'level1_annotation':['Fib']*4+['Mac']*2,'level2_annotation':['A','A','B','B','C','C']},index=[f'c{i}' for i in range(6)])
  x=sparse.csr_matrix([[10,9],[10,9],[10,1],[10,1],[0,0],[0,0]])
  a=ad.AnnData(x,obs=obs,var=pd.DataFrame(index=['shared','state_specific']));a.layers['log1p']=x.astype(float);path=self.root/'marker.h5ad';a.write_h5ad(path)
  cfg={'atlas':{'path':str(path),'level1_key':'level1_annotation','level2_key':'level2_annotation'},'execution':{'block_rows':2}}
  marker_summary(cfg,obs,self.root/'markers')
  g=pd.read_csv(self.root/'markers/L2_all_gene_marker_evidence.tsv.gz',sep='\t');w=pd.read_csv(self.root/'markers/L2_within_parent_marker_evidence.tsv.gz',sep='\t')
  self.assertGreater(g.query('level2_annotation=="A" and gene_id=="shared"').mean_difference.iloc[0],0)
  self.assertEqual(w.query('level2_annotation=="A" and gene_id=="shared"').mean_difference.iloc[0],0)
  self.assertGreater(w.query('level2_annotation=="A" and gene_id=="state_specific"').mean_difference.iloc[0],0)
  self.assertTrue(w.query('level2_annotation=="C"').mean_difference.isna().all())
 def test_final_ontology_anchor_exact_resolution(self):
  from src.anchors import final_anchor_panel
  ont=pd.DataFrame({'level1_annotation':['Fib'],'level2_annotation':['StateA'],'marker_anchor':['G1;G2'],'approved':[True]})
  genes=pd.DataFrame({'gene_symbol':['G1','G2']},index=['id1','id2'])
  a=final_anchor_panel(ont,genes);self.assertEqual(list(a.gene_id),['id1','id2'])
  bad=genes.copy();bad.loc['id3']=['G1']
  with self.assertRaisesRegex(Blocked,'exactly once'):final_anchor_panel(ont,bad)
  ont.marker_anchor='not a gene'
  with self.assertRaises(Blocked):final_anchor_panel(ont,genes)
 def test_actual_l2_dotplot_uses_final_ontology_anchors(self):
  from src.expression import marker_summary
  from src.figures import dotplot
  fixture=test_core.CoreTests();fixture.setUp()
  try:
   cfg=fixture.make_release();cfg['execution']={'block_rows':3};cfg['figures']={'l2_anchor_panel':None}
   op=fixture.root/'support/level2_ontology.tsv';ont=pd.read_csv(op,sep='\t');ont['marker_anchor']='G0;G1';ont.to_csv(op,sep='\t',index=False)
   run=fixture.root/'review';run.mkdir();marker_summary(cfg,fixture.obs,run/'markers')
   result=dotplot(cfg,run,'L2');self.assertEqual(result['status'],'COMPLETE')
   source=pd.read_csv(run/'figures/L2_DotPlot_source.tsv',sep='\t');self.assertEqual(set(source.gene_id),{'G0','G1'})
   self.assertTrue((run/'figures/L2_DotPlot.png').is_file());self.assertTrue((run/'figures/resolved_final_L2_anchors.tsv').is_file())
   panel=pd.read_csv(run/'figures/resolved_final_L2_anchors.tsv',sep='\t');panel.loc[0,'gene_id']='MISSING_GENE';panel.to_csv(fixture.root/'bad_panel.tsv',sep='\t',index=False);cfg['figures']['l2_anchor_panel']=str(fixture.root/'bad_panel.tsv')
   with self.assertRaisesRegex(Blocked,'missing exact gene'):dotplot(cfg,run,'L2')
  finally:fixture.tearDown()
if __name__=='__main__':unittest.main(verbosity=2)
