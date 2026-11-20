"""Native AnnData release→authenticated real CLI→artifact verification→projection."""
import unittest,json,copy,sys,importlib.util
from pathlib import Path
from unittest.mock import patch
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from tests import test_core
from src.core import preflight,atomic_json,sha,Blocked
from src.dispatcher import code_digest
from src.review_dependencies import reviewed_dependencies
from src.spatial_projection_contract import weights_record,abundance_record,validate_projection_inputs,ontology_digest
from src.spatial_contract import digest
from src.spatial_renderer import effective_parameters
from src.state_identity import build_state_identity_map
ROOT=Path(__file__).resolve().parents[1]
class ProjectionV21Tests(unittest.TestCase):
 def setUp(self):
  self.fixture=test_core.CoreTests();self.fixture.setUp();self.r=self.fixture.root;self.c=self.fixture.make_release();self.ont=pd.read_csv(self.r/'support/level2_ontology.tsv',sep='\t')
  self.run=self.r/'source_run';self.run.mkdir();(self.run/'counts').mkdir();(self.run/'covarnet/stable_NMF').mkdir(parents=True)
  self.cp=self.r/'config.json';self.c.update(execution={'mode':'approved_execution','run_full_analysis':True,'cpu_only':True,'max_workers':1,'per_worker_threads':2,'approval_receipt':str(self.r/'approval.json')},statistical_units={},covarnet={'workers':1,'threads_per_worker':2},external={'spatial':{'reference_mode':'descriptive','excluded_targets':[],'forced_genes':[]}})
  preflight(self.c,self.run/'preflight.json');atomic_json(self.run/'RUN_IDENTITY.json',{'atlas_sha256':self.c['atlas']['input_sha256'],'synthetic_fixture':'not a real fitted NMF'})
  atomic_json(self.run/'RUN_STATUS.json',{'covarnet':{'status':'COMPLETE'}})
  self.ont.to_csv(self.run/'counts/ontology.tsv',sep='\t',index=False)
  atomic_json(self.run/'covarnet/INPUT_MANIFEST.json',{'source_identity':{'files':{str(self.run/'counts/ontology.tsv'):sha(self.run/'counts/ontology.tsv')}}})
  m=build_state_identity_map(self.ont);self.m=m
  self.w=self.run/'covarnet/stable_NMF/module_node_weights_K02.tsv';pd.DataFrame({'CM01':[.8,.2],'CM02':[.2,.8]},index=m.state_id).to_csv(self.w,sep='\t');self.complete_weights()
  self.post=self.r/'posterior.h5ad';spots=['p0','p1','p2'];ob=pd.DataFrame({'sample_id':['s']*3,'section_id':['section']*3,'biological_unit_id':['donor']*3},index=spots)
  a=ad.AnnData(sparse.csr_matrix(np.ones((3,1))),obs=ob,var=pd.DataFrame(index=['REF']));a.obsm['q05_cell_abundance_w_sf']=pd.DataFrame([[1.,4.],[4.,1.],[2.,3.]],index=spots,columns=m.safe_label);a.write_h5ad(self.post)
  self.ref=self.r/'reference.tsv';pd.DataFrame({s:[1.] for s in m.safe_label},index=['REF']).to_csv(self.ref,sep='\t');self.prov=self.r/'posterior.provenance.json'
  params=effective_parameters(self.c['external']['spatial']);params.update(reference_mode='descriptive',excluded_targets=[],forced_genes=[],only_gses=[])
  atomic_json(self.prov,{'schema':'spatial_provenance_v2','atlas_sha256':self.c['atlas']['input_sha256'],'ontology_sha256':ontology_digest(self.ont),'parameters_sha256':digest(params),'reference_artifact_sha256':sha(self.ref),'posterior_sha256':sha(self.post),'reference_genes_sha256':digest(['REF']),'spot_ids_sha256':digest(spots),'sample_ids_sha256':digest(['s']*3),'targets_sha256':digest([]),'spatial_units_sha256':digest([['section','donor']]*3),'reference_mode':'descriptive','excluded_targets':[]})
  self.ab=self.r/'abundance.tsv';pd.DataFrame([[1.,4.],[4.,1.],[2.,3.]],index=spots,columns=m.state_id).to_csv(self.ab,sep='\t')
  self.refresh_records();self.approve()
 def tearDown(self):self.fixture.tearDown()
 def complete_weights(self):
  selected=self.run/'covarnet/stable_NMF/selected_K.txt';selected.write_text('K_FINAL=2\n')
  nodes=self.run/'covarnet/CoVarNet_node_counts.tsv'
  if not nodes.exists():pd.DataFrame({'profile1':[30,20],'profile2':[10,25]},index=self.m.state_id).to_csv(nodes,sep='\t')
  inp=self.run/'covarnet/INPUT_MANIFEST.json';rec=json.loads(inp.read_text());rec['files']={nodes.name:sha(nodes)};atomic_json(inp,rec)
  atomic_json(self.run/'covarnet/COVARNET_COMPLETE.json',{'status':'COMPLETE','input_manifest_sha256':sha(self.run/'covarnet/INPUT_MANIFEST.json'),'files':{'stable_NMF/'+self.w.name:sha(self.w),'stable_NMF/selected_K.txt':sha(selected),'CoVarNet_node_counts.tsv':sha(nodes)}})
 def refresh_records(self):
  wm=weights_record(self.w,self.run,self.ont,self.c['atlas']['input_sha256']);am=abundance_record(self.ab,self.post,self.prov,self.ref,self.ont,self.c,'state_id')
  atomic_json(self.r/'weights.json',wm);atomic_json(self.r/'abundance.json',am)
  self.c['external']['spatial']['module_projection']={'approved':True,'weights_manifest':str(self.r/'weights.json'),'weights_manifest_sha256':sha(self.r/'weights.json'),'abundance_manifest':str(self.r/'abundance.json'),'abundance_manifest_sha256':sha(self.r/'abundance.json'),'covarnet_run_dir':str(self.run),'covarnet_run_identity_sha256':sha(self.run/'RUN_IDENTITY.json'),'posterior_path':str(self.post),'posterior_sha256':sha(self.post)}
 def approve(self):
  atomic_json(self.cp,self.c);d,_=reviewed_dependencies(self.c,self.cp);atomic_json(self.r/'approval.json',{'approved':True,'atlas_sha256':self.c['atlas']['input_sha256'],'config_sha256':sha(self.cp),'code_sha256':code_digest(ROOT),'reviewed_dependency_digest':d})
 def validate(self,c=None):return validate_projection_inputs(c or self.c,self.ont,self.ab,self.w,'state_id')
 def main(self,name):
  sp=importlib.util.spec_from_file_location('projection_actual_cli',ROOT/'scripts/external_spatial_modules.py');e=importlib.util.module_from_spec(sp);sp.loader.exec_module(e)
  args=['cli','--config',str(self.cp),'--weights',str(self.w),'--abundance',str(self.ab),'--abundance-key','state_id','--out',str(self.r/name)]
  with patch.object(sys,'argv',args):e.main()
 def test_native_release_authenticated_cli_projection(self):
  self.main('out');p=self.r/'out';v=json.loads((p/'PROJECTION_INPUT_VALIDATION.json').read_text());self.assertFalse(v['independent_spatial_validation']);self.assertEqual(pd.read_csv(p/'spatial_module_scores.tsv.gz',sep='\t').shape[0],3)
 def test_missing_manifests_or_selection_block(self):
  c=copy.deepcopy(self.c);c['external']['spatial']['module_projection']['weights_manifest']=None
  with self.assertRaisesRegex(ValueError,'Missing'):self.validate(c)
  c=copy.deepcopy(self.c);c['external']['spatial']['module_projection']['approved']=False
  with self.assertRaisesRegex(ValueError,'not approved'):self.validate(c)
 def test_wrong_atlas_ontology_or_run_block(self):
  c=copy.deepcopy(self.c);c['atlas']['input_sha256']='wrong'
  with self.assertRaisesRegex(ValueError,'atlas'):self.validate(c)
  old=self.ont.copy();self.ont.loc[0,'level2_annotation']='another_state'
  with self.assertRaisesRegex(ValueError,'ontology'):self.validate()
  self.ont=old;c=copy.deepcopy(self.c);c['external']['spatial']['module_projection']['covarnet_run_dir']=str(self.r/'other_run')
  with self.assertRaisesRegex(ValueError,'run'):self.validate(c)
 def test_changed_weights_rejected_without_new_approval(self):
  w=pd.read_csv(self.w,sep='\t',index_col=0);w['CM01']=[.1,.9];w.to_csv(self.w,sep='\t')
  with self.assertRaisesRegex(ValueError,'digest'):self.validate()
  self.complete_weights();self.refresh_records()
  # Valid regenerated artifacts still require a corresponding new approval receipt.
  atomic_json(self.cp,self.c)
  with self.assertRaises(Blocked):self.main('old_approval')
  self.approve();self.main('new_approval');self.assertTrue((self.r/'new_approval/spatial_module_scores.tsv.gz').exists())
 def test_abundance_values_spot_order_and_content_binding(self):
  original=self.ab.read_text();ab=pd.read_csv(self.ab,sep='\t',index_col=0);ab.iloc[0,0]=123.;ab.to_csv(self.ab,sep='\t')
  with self.assertRaisesRegex(ValueError,'values'):self.validate()
  self.ab.write_text(original);ab=pd.read_csv(self.ab,sep='\t',index_col=0).iloc[::-1];ab.to_csv(self.ab,sep='\t')
  with self.assertRaisesRegex(ValueError,'spot'):self.validate()
  self.ab.write_text(original);self.ab.write_text(original+'\n')
  with self.assertRaisesRegex(ValueError,'manifest'):self.validate()
 def test_actual_cli_other_table_not_config_bound_copy(self):
  other=self.r/'renamed_old_weights.tsv';other.write_bytes(self.w.read_bytes())
  with self.assertRaisesRegex(ValueError,'actual CLI'):validate_projection_inputs(self.c,self.ont,self.ab,other,'state_id')
 def test_posterior_or_provenance_changed_rejects(self):
  record=json.loads(self.prov.read_text());record['atlas_sha256']='old_locked71';atomic_json(self.prov,record)
  with self.assertRaises(ValueError):self.validate()
 def test_counts_table_cannot_masquerade_as_weights(self):
  with self.assertRaisesRegex(ValueError,'selected primary W'):weights_record(self.run/'covarnet/CoVarNet_node_counts.tsv',self.run,self.ont,self.c['atlas']['input_sha256'])
 def test_module_ids_and_normalized_weight_semantics(self):
  w=pd.read_csv(self.w,sep='\t',index_col=0);w.columns=['profile1','profile2'];w.to_csv(self.w,sep='\t');self.complete_weights()
  with self.assertRaisesRegex(ValueError,'module IDs'):self.refresh_records()
  w.columns=['CM01','CM02'];w*=2;w.to_csv(self.w,sep='\t');self.complete_weights()
  with self.assertRaisesRegex(ValueError,'normalized'):self.refresh_records()
 def test_embedding_is_not_posterior_abundance(self):
  with self.assertRaisesRegex(ValueError,'Unsupported'):abundance_record(self.ab,self.post,self.prov,self.ref,self.ont,self.c,'state_id','X_umap')
 def test_posterior_factor_names_conflict_rejected(self):
  a=ad.read_h5ad(self.post);a.uns['mod']={'factor_names':list(reversed(self.m.safe_label))};a.write_h5ad(self.post)
  prov=json.loads(self.prov.read_text());prov['posterior_sha256']=sha(self.post);atomic_json(self.prov,prov)
  with self.assertRaisesRegex(ValueError,'factor_names'):abundance_record(self.ab,self.post,self.prov,self.ref,self.ont,self.c,'state_id')
 def test_spot_ids_with_leading_zero_preserved_by_real_cli(self):
  import gzip
  a=ad.read_h5ad(self.post);spots=['001','002','003'];a.obs_names=spots;a.write_h5ad(self.post)
  ab=pd.read_csv(self.ab,sep='\t',index_col=0);ab.index=spots;ab.to_csv(self.ab,sep='\t')
  prov=json.loads(self.prov.read_text());prov['posterior_sha256']=sha(self.post);prov['spot_ids_sha256']=digest(spots);atomic_json(self.prov,prov)
  self.refresh_records();self.approve();self.main('leading_zero')
  with gzip.open(self.r/'leading_zero/spatial_module_scores.tsv.gz','rt') as f:self.assertTrue(f.read().splitlines()[1].startswith('001\t'))
 def test_partial_receipt_blocks_actual_dispatcher_targets(self):
  from tests import test_revision_interfaces
  fixture=test_revision_interfaces.InterfaceRevisionTests();fixture.setUp()
  try:
   states=fixture.dispatch('PARTIAL');self.assertEqual(states['meta']['status'],'PARTIAL');self.assertEqual(states['targets']['status'],'BLOCKED_DEPENDENCY')
  finally:fixture.tearDown()
if __name__=='__main__':unittest.main(verbosity=2)
