"""Synthetic spatial interface/provenance regression; no atlas/model training."""
import ast,copy,json,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np,pandas as pd
from state_identity import build_state_identity_map,map_state_axis
from external_statistics import spatial_module_projection
from external_modules import render_spatial_adapter
from spatial_renderer import effective_parameters,GLOBALS
from spatial_contract import digest,file_digest,validate_spatial_provenance,strict_fit_count_models,spatial_design_gate,matched_person_empirical
ROOT=Path(__file__).resolve().parents[1]
class RevisionSpatialTests(unittest.TestCase):
 def ontology(self):return pd.DataFrame({'level1_annotation':['Macrophage','T cell'],'level2_annotation':['C1QC/FOLR2 state','CD4+ memory']})
 def test_state_mapping_chain(self):
  m=build_state_identity_map(self.ontology());w=pd.DataFrame({'M1':[.4,.6]},index=m.state_id)
  a=pd.DataFrame([[1.,2.],[3.,4.]],columns=m.safe_label)
  s,audit=spatial_module_projection(a,w,m,'safe_label')
  self.assertEqual(s.shape,(2,1));self.assertTrue(audit.full_weight_coverage.all())
  for namespace in ['parent_qualified','safe_label','spatial_reference_state_original']:
   mapped=map_state_axis(m.state_id,m,'state_id',namespace)
   self.assertEqual(map_state_axis(mapped,m,namespace).tolist(),m.state_id.tolist())
 def test_explicit_ids_distinguish_same_display_name(self):
  o=self.ontology();o.level2_annotation='shared';o['state_id']=['MA_01','TC_01']
  m=build_state_identity_map(o);self.assertEqual(m.state_id.tolist(),['MA_01','TC_01'])
  with self.assertRaises(ValueError):map_state_axis(['shared'],m,'level2_annotation')
  with self.assertRaises(ValueError):build_state_identity_map(o.drop(columns='state_id'))
 def test_release_level2_state_id_priority(self):
  o=self.ontology();o['level2_state_id']=['mac_01','T_01'];m=build_state_identity_map(o)
  self.assertEqual(m.state_id.tolist(),['mac_01','T_01'])
  o['state_id']=['mac_wrong','T_01']
  with self.assertRaisesRegex(ValueError,'Conflicting'):build_state_identity_map(o)
 def test_missing_important_weights_blocked_and_audited(self):
  m=build_state_identity_map(self.ontology());w=pd.DataFrame({'M1':[.1,.9]},index=m.state_id)
  with self.assertRaises(ValueError) as cm:spatial_module_projection(pd.DataFrame([[1.]],columns=m.safe_label[:1]),w,m,'safe_label')
  self.assertAlmostEqual(cm.exception.audit.iloc[0].missing_weight_fraction,.9)
 def test_all_rendered_parameters_locked_and_perturbable(self):
  source=(ROOT/'external_legacy/spatial/02_spatial_v4_cell2location_GPU.py').read_text()
  for key,glob in GLOBALS.items():
   value=33.5 if key=='detection_alpha' else (effective_parameters()[key]+17)
   s=render_spatial_adapter(source,'/tmp/a','/tmp/new','/tmp/manifest','level1_annotation','level2_annotation','study',[],'blocked',training_config={key:value})
   tree=ast.parse(s);literal={}
   for n in tree.body:
    if isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name):
     try:literal[n.targets[0].id]=ast.literal_eval(n.value)
     except (ValueError,TypeError):pass
   self.assertEqual(literal[glob],value);self.assertEqual(literal['EFFECTIVE_PARAMETERS'][key],value)
   self.assertNotIn('os.environ.get',s);self.assertNotIn('atlas.raw.',s)
   self.assertIn('selected_genes.difference_update',s)
   self.assertIn('Existing posterior forbidden',s)
 def test_target_excluded_renderer_rejects_conflicting_forced_genes(self):
  source=(ROOT/'external_legacy/spatial/02_spatial_v4_cell2location_GPU.py').read_text()
  with self.assertRaises(ValueError):render_spatial_adapter(source,'a','b','c','l1','l2','study',['T'],'a',training_config={'reference_mode':'target_excluded','excluded_targets':['T']})
 def record(self,p):
  expected={k:k*2 for k in ['atlas_sha256','ontology_sha256','parameters_sha256','reference_artifact_sha256']}
  r=dict(expected,schema='spatial_provenance_v2',posterior_sha256=file_digest(p),reference_mode='target_excluded',excluded_targets=['T'],reference_genes_sha256=digest(['G']),targets_sha256=digest(['T']),spot_ids_sha256=digest(['s1','s2']),sample_ids_sha256=digest(['A','A']),spatial_units_sha256=digest([['S','D'],['S','D']]))
  return r,expected
 def test_provenance_positive_and_rejections(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'posterior';p.write_text('synthetic')
   r,e=self.record(p)
   def validate(r=r,e=e,genes=['G'],mode='target_validation'):
    return validate_spatial_provenance(r,e,genes,['T'],['s1','s2'],['A','A'],mode,p,[['S','D'],['S','D']])
   self.assertTrue(validate()['inferential_target_validation_allowed'])
   for key in ['atlas_sha256','ontology_sha256','parameters_sha256','reference_artifact_sha256','spot_ids_sha256','spatial_units_sha256','posterior_sha256']:
    bad=dict(r);bad[key]='wrong'
    with self.assertRaises(ValueError):validate(r=bad)
   bad=dict(r,source='locked71')
   with self.assertRaises(ValueError):validate(r=bad)
   with self.assertRaises(ValueError):validate(r={})
   bad=dict(r,reference_genes_sha256=digest(['G','T']))
   with self.assertRaises(ValueError):validate(r=bad,genes=['G','T'])
   self.assertFalse(validate(r=bad,genes=['G','T'],mode='descriptive')['inferential_target_validation_allowed'])
 def test_covariance_failure_never_hc1(self):
  import statsmodels.api as sm
  n=100;y=np.arange(n)%7+1;c=np.arange(n,dtype=float);xy=np.column_stack([np.arange(n)%10,np.arange(n)//10])
  with patch.object(sm.GLM,'fit',side_effect=RuntimeError('synthetic cluster failure')) as fit:
   with self.assertRaisesRegex(ValueError,'spatial_block_cluster failed'):strict_fit_count_models(y,c,np.ones(n)*100,xy)
   self.assertEqual(fit.call_count,1);self.assertEqual(fit.call_args.kwargs['cov_type'],'cluster')
 def test_support_boundaries_and_unapproved_science(self):
  for n in [39,40,99,100]:
   for min_n in [40,100]:
    for positive in [10,20]:
     y=np.zeros(n);y[:positive]=1;c=np.arange(n,dtype=float);xy=np.column_stack([np.arange(n)%10,np.arange(n)//10])
     for minimum in [10,20]:
      if n<min_n or positive<minimum:
       with self.assertRaisesRegex(ValueError,'Insufficient'):strict_fit_count_models(y,c,np.ones(n)*100,xy,minimum,min_n)
  with self.assertRaisesRegex(ValueError,'BLOCKED'):spatial_design_gate({})
  base={'approved':True,'min_spots':100,'min_positive_spots':20,'interaction_model':'none_additive_section','section_to_person':'mean_section_z_same_person_complete_support','bh_family':'all_target_state_pairs_person_empirical','variance_method':'spatial_block_cluster'}
  self.assertEqual(spatial_design_gate(base)['min_spots'],100)
  for key in ['interaction_model','section_to_person','bh_family','variance_method']:
   with self.assertRaisesRegex(ValueError,'BLOCKED'):spatial_design_gate(dict(base,**{key:'unapproved_other'}))
 def test_matched_null_identical_complete_sections(self):
  design={'approved':True,'min_spots':100,'min_positive_spots':20,'interaction_model':'none_additive_section','section_to_person':'mean_section_z_same_person_complete_support','bh_family':'all_target_state_pairs_person_empirical','variance_method':'spatial_block_cluster'}
  rows=[]
  for i in range(3):
   for section in range(2):
    for rank in [0,1,2]:
     rows.append(dict(pair_id='T|state',dataset=f'D{i}',section_id=f'S{i}_{section}',biological_unit_id=f'P{i}',match_rank=rank,is_target=rank==0,z=3-rank,tested_gene=f'gene_{i}_{rank}'))
  c=pd.DataFrame(rows);c.loc[(c.biological_unit_id=='P0')&(c.section_id=='S0_0')&(c.match_rank==1),'z']=np.nan
  out=matched_person_empirical(c,design)
  self.assertEqual(out.iloc[0].n_sections_complete,5);self.assertEqual(out.iloc[0].n_biological_units,3);self.assertEqual(out.iloc[0].n_control_ranks,2)
  self.assertTrue(np.isfinite(out.iloc[0].BH_q))
 def test_actual_spatial_statistics_entry_with_synthetic_posterior(self):
  import importlib.util,anndata as ad
  from scipy import sparse
  spec=importlib.util.spec_from_file_location('spatial_entry',ROOT/'scripts/external_spatial_statistics.py');entry=importlib.util.module_from_spec(spec);spec.loader.exec_module(entry)
  with tempfile.TemporaryDirectory() as td:
   d=Path(td);rng=np.random.default_rng(12);n=200;o=self.ontology();m=build_state_identity_map(o)
   obs=pd.DataFrame({'sample_id':['sample1']*100+['sample2']*100,'section_id':['section1']*100+['section2']*100,'biological_unit_id':['person1']*100+['person2']*100},index=[f'spot{i}' for i in range(n)])
   raw=ad.AnnData(sparse.csr_matrix(rng.poisson(4,(n,9))),obs=obs,var=pd.DataFrame(index=['T']+[f'G{i}' for i in range(8)]));raw.write_h5ad(d/'raw.h5ad')
   mapped=raw[:,1:].copy();mapped.obsm['spatial']=np.column_stack([np.arange(n)%10,np.arange(n)//10]);mapped.obsm['q05_cell_abundance_w_sf']=pd.DataFrame(rng.uniform(1,20,(n,2)),index=obs.index,columns=m.safe_label);mapped.write_h5ad(d/'mapped.h5ad')
   ref=pd.DataFrame(np.ones((8,2)),index=mapped.var_names,columns=m.safe_label);ref.to_csv(d/'reference.tsv',sep='\t')
   design={'approved':True,'min_spots':100,'min_positive_spots':20,'interaction_model':'none_additive_section','section_to_person':'mean_section_z_same_person_complete_support','bh_family':'all_target_state_pairs_person_empirical','variance_method':'spatial_block_cluster'}
   sp={'reference_mode':'target_excluded','excluded_targets':['T'],'statistical_design':design};params=effective_parameters(sp);params.update(reference_mode='target_excluded',excluded_targets=['T'],forced_genes=[],only_gses=[])
   cfg={'atlas':{'input_sha256':digest('synthetic atlas')},'external':{'spatial':sp}}
   rec={'schema':'spatial_provenance_v2','atlas_sha256':cfg['atlas']['input_sha256'],'ontology_sha256':digest(o.fillna('').astype(str).sort_values(['level1_annotation','level2_annotation']).to_dict('records')),'parameters_sha256':digest(params),'reference_artifact_sha256':file_digest(d/'reference.tsv'),'posterior_sha256':file_digest(d/'mapped.h5ad'),'reference_genes_sha256':digest(list(mapped.var_names)),'targets_sha256':digest(['T']),'spot_ids_sha256':digest(list(obs.index)),'sample_ids_sha256':digest(obs.sample_id.tolist()),'spatial_units_sha256':digest(obs[['section_id','biological_unit_id']].values.tolist()),'reference_mode':'target_excluded','excluded_targets':['T']}
   (d/'provenance.json').write_text(json.dumps(rec));manifest=pd.DataFrame([{'dataset':'synthetic','mapped_h5ad':str(d/'mapped.h5ad'),'full_raw_h5ad':str(d/'raw.h5ad'),'full_raw_sha256':file_digest(d/'raw.h5ad'),'reference_signatures':str(d/'reference.tsv'),'provenance_json':str(d/'provenance.json')}]);roster=pd.DataFrame([{'gene':'T','cell_state':m.state_id.iloc[0],'annotation_level':'level2_annotation'}])
   entry.calculate(manifest,roster,m,d/'valid',cfg,o)
   result=pd.read_csv(d/'valid/section_count_models.tsv',sep='\t');self.assertEqual(len(result),2);self.assertTrue((result.variance_method=='spatial_block_cluster').all())
   rec['source']='locked71';(d/'provenance.json').write_text(json.dumps(rec))
   with self.assertRaisesRegex(ValueError,'locked71'):entry.calculate(manifest,roster,m,d/'invalid',cfg,o)
 def test_native_small_glm(self):
  rng=np.random.default_rng(4);n=100;context=rng.uniform(1,20,n);y=rng.poisson(4+.2*context);xy=np.column_stack([np.arange(n)%10,np.arange(n)//10])
  z=strict_fit_count_models(y,context,np.ones(n)*100,xy,20,100)
  self.assertEqual(z['variance_method'],'spatial_block_cluster');self.assertGreater(z['n_blocks'],1);self.assertTrue(np.isfinite(z['nb_se']))
if __name__=='__main__':unittest.main()
