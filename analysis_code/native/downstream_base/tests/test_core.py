"""Synthetic tests only. No real atlas path is ever used here."""
import unittest,tempfile,json,copy,os
from pathlib import Path
import numpy as np,pandas as pd
from scipy import sparse
import anndata as ad
from src.core import Blocked,check_annotations,preflight,sha,digest_ids,atomic_json
from src.cohort import cohort_from_obs,annotate_units,apply_unit_review,select_contrast
from src.counts import state_counts,biological_counts
from src.expression import aggregate_sparse,export_pseudobulk,marker_summary

class CoreTests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
  n=24;g=8;rng=np.random.default_rng(10);self.x=sparse.csr_matrix(rng.poisson(2,(n,g)),dtype='int32')
  obs=pd.DataFrame(index=pd.Index([f'c{i}' for i in range(n)],name='cell_id'))
  obs['sample_uid']=[f's{i//4}' for i in range(n)]
  for col,values in {'biological_sample_uid':[f'b{i//4}' for i in range(n)],'donor_uid':[f'd{i//4}' for i in range(n)],'study':['A']*12+['B']*12,'atlas_source_id':['A']*12+['B']*12,'tissue':['aorta']*n,'vascular_bed':['aortic']*n,'vascular_bed_group':['arterial']*n,'comparison_group':['reference']*12+['case']*12,'cohort_disease_context':['disease_context']*n,'assay_type':['scRNA-seq']*n,'age':[40]*n,'sex':['female']*n}.items():obs[col]=values
  obs['level1_annotation']=['Fibroblast']*8+['Macrophage']*8+['Unknown']*4+['Schwann cell']*4
  obs['level2_annotation']=['COL1 Fib']*8+['C1QC Mac']*4+['Unknown']*8+['Not_assessed']*4
  self.obs=obs;self.ont=obs[['level1_annotation','level2_annotation']].drop_duplicates().copy();self.ont['approved']=True
  self.ont['n_cells']=[((obs.level1_annotation==p)&(obs.level2_annotation==s)).sum() for p,s in self.ont[['level1_annotation','level2_annotation']].itertuples(index=False,name=None)]
  a=ad.AnnData(self.x.copy(),obs=obs.copy(),var=pd.DataFrame({'gene_symbol':[f'G{i}' for i in range(g)]},index=pd.Index([f'G{i}' for i in range(g)],name='gene_id')))
  z=self.x.astype(float).multiply(10000/self.x.sum(axis=1).A1[:,None]).tocsr();z.data=np.log1p(z.data);a.layers['log1p']=z
  self.path=self.root/'synthetic.h5ad';a.write_h5ad(self.path)
 def tearDown(self):self.tmp.cleanup()
 def test_annotation_hierarchy(self):self.assertEqual(len(check_annotations(self.obs,self.ont)),24)
 def test_pending_rejected(self):
  o=self.obs.copy();o.iloc[0,o.columns.get_loc('level2_annotation')]='Pending'
  with self.assertRaises(Blocked):check_annotations(o,self.ont)
 def test_alias_conflict(self):
  o=self.obs.copy();o['L2_annotation']='x'
  with self.assertRaises(Blocked):check_annotations(o,self.ont)
 def test_wrong_parent(self):
  o=self.obs.copy();o.iloc[0,o.columns.get_loc('level1_annotation')]='Macrophage'
  with self.assertRaises(Blocked):check_annotations(o,self.ont)
 def test_denominator_unknown_and_missing_parent(self):
  t=state_counts(self.obs)
  x=t[(t.sample_uid=='s0')&(t.level1_annotation=='Macrophage')]
  self.assertTrue(x.prop_within_parent.isna().all());self.assertTrue(x.n_cells.eq(0).all())
  for sample,g in t.groupby('sample_uid'):self.assertEqual(g.n_cells.sum(),4);self.assertAlmostEqual(g.prop_all_qc.sum(),1.)
  y=t[(t.sample_uid=='s3')&(t.level1_annotation=='Macrophage')];self.assertTrue(y.n_parent.eq(4).all())
 def test_counts_stream_vs_dense(self):
  sums,u,v=aggregate_sparse(self.path,self.obs,['sample_uid','level2_annotation'],integer=True,block_rows=3)
  self.assertEqual(sums.sum(),self.x.sum())
  for i,row in u.iterrows():
   ids=(self.obs.sample_uid==row.sample_uid)&(self.obs.level2_annotation==row.level2_annotation)
   np.testing.assert_array_equal(sums[i].toarray(),self.x[ids].sum(axis=0))
 def test_detection_matches_raw(self):
  sums,u,v=aggregate_sparse(self.path,self.obs,['level2_annotation'],detection=True,block_rows=5)
  for i,row in u.iterrows():
   ids=self.obs.level2_annotation==row.level2_annotation
   np.testing.assert_array_equal(sums[i].toarray(),(self.x[ids]>0).sum(axis=0))
 def test_fractional_raw_rejected(self):
  a=ad.read_h5ad(self.path);a.X=a.X.astype(float);a.X.data[0]=.5;a.write_h5ad(self.root/'bad.h5ad')
  with self.assertRaises(Blocked):aggregate_sparse(self.root/'bad.h5ad',self.obs,['sample_uid'],integer=True)
 def test_cohort_conflict_rejected(self):
  o=self.obs.copy();o.loc['c0','donor_uid']='other'
  with self.assertRaises(Blocked):cohort_from_obs(o)
 def test_technical_replicates_require_review(self):
  c=cohort_from_obs(self.obs);c.loc[c.sample_uid=='s1','biological_sample_uid']='b0';c=annotate_units(c)
  self.assertFalse(c.loc[c.sample_uid.isin(['s0','s1']),'aggregate_approved'].any())
 def test_missing_donor_never_new_person(self):
  c=cohort_from_obs(self.obs);c.loc[0,'donor_uid']='NR';c=annotate_units(c)
  self.assertFalse(c.iloc[0].inference_eligible);self.assertEqual(c.iloc[0].donor_uid,'NR')
 def test_missing_covariates_do_not_global_drop(self):
  c=cohort_from_obs(self.obs);c.loc[0,'age']=np.nan
  self.assertEqual(len(c),6);self.assertTrue(c.inference_eligible.all())
 def test_pseudobulk_export_preserves_input(self):
  before=sha(self.path);c=cohort_from_obs(self.obs)
  cfg={'atlas':{'path':str(self.path),'level1_key':'level1_annotation','level2_key':'level2_annotation'},'execution':{'block_rows':3},'pseudobulk':{'min_cells':2}}
  export_pseudobulk(cfg,self.obs,c,self.root/'pb');self.assertEqual(before,sha(self.path))
  marker_summary(cfg,self.obs,self.root/'markers');self.assertEqual(before,sha(self.path))
 def test_unapproved_contrast(self):
  with self.assertRaises(Blocked):select_contrast(cohort_from_obs(self.obs),{'approved':False})
 def test_unpaired_repeated_donor_rejected(self):
  m=cohort_from_obs(self.obs);m.loc[m.sample_uid=='s1','donor_uid']='d0'
  with self.assertRaises(Blocked):select_contrast(m,{'approved':True,'reference':'reference','case':'case','min_independent_per_arm':1,'paired':False,'covariates':''})
 def test_unknown_retained(self):self.assertEqual(state_counts(self.obs).query('level2_annotation == "Unknown"').n_cells.sum(),8)
 def test_shuffled_expression_metadata_rejected(self):
  with self.assertRaisesRegex(Blocked,'order'):aggregate_sparse(self.path,self.obs.iloc[::-1],['sample_uid'],integer=True)
 def test_fractional_detection_rejected(self):
  a=ad.read_h5ad(self.path);a.X=a.X.astype(float);a.X.data[0]=.2;a.write_h5ad(self.root/'bad_detection.h5ad')
  with self.assertRaises(Blocked):aggregate_sparse(self.root/'bad_detection.h5ad',self.obs,['sample_uid'],detection=True)
 def test_provisional_donor_not_independent(self):
  c=cohort_from_obs(self.obs);c.loc[0,'donor_uid']='PROVISIONAL::s0';c=annotate_units(c)
  self.assertFalse(c.iloc[0].inference_eligible)
 def make_release(self):
  cfg={'atlas':{'path':str(self.path),'manifest':str(self.root/'manifest.json'),'checksum':str(self.root/'atlas.sha256'),'schema_contract':str(self.root/'contract.json'),'expected_shape':[24,8],'expected_release_version':'L1L2_final_v1','sample_key':'sample_uid','level1_key':'level1_annotation','level2_key':'level2_annotation'},'outputs':{'release_support':str(self.root/'support')}}
  s=self.root/'support';s.mkdir();self.ont.to_csv(s/'level2_ontology.tsv',sep='\t',index=False)
  self.obs[['level1_annotation','level2_annotation']].to_csv(s/'cell_annotations.tsv.gz',sep='\t')
  self.obs.drop(columns=['level1_annotation','level2_annotation']).drop_duplicates('sample_uid').to_csv(s/'sample_metadata.csv',index=False)
  contract={'expected_identity':{'n_sample_uid':6,'expected_level1_counts':self.obs.level1_annotation.value_counts().to_dict()},'data_layout':{'metadata_columns_expected_from_v4':[c for c in self.obs if not c.startswith('level')]},'release_sidecars':{'level2_ontology.tsv':[],'cell_annotations.tsv.gz':[],'sample_metadata.csv':[]}}
  atomic_json(self.root/'contract.json',contract)
  a=ad.read_h5ad(self.path);a.uns.update(release_status='FINAL_APPROVED',level1_frozen=True,level2_frozen=True,atlas_release_version='L1L2_final_v1',annotation_keys={'level1':'level1_annotation','level2':'level2_annotation'},neighbors_cpu_n30={'connectivities_key':'C','distances_key':'D'})
  for k in ['P','H','U']:a.obsm[k]=np.zeros((24,2))
  a.obsp['C']=sparse.eye(24,format='csr');a.obsp['D']=sparse.eye(24,format='csr');a.write_h5ad(self.path)
  h=sha(self.path);cfg['atlas']['input_sha256']=h;(self.root/'atlas.sha256').write_text(h+'  synthetic.h5ad\n')
  atomic_json(self.root/'manifest.json',{'h5ad_path':str(self.path),'h5ad_sha256':h,'release_status':'FINAL_APPROVED','release_version':'L1L2_final_v1','ready_for_downstream':True,'level1_frozen':True,'level2_frozen':True,'cell_id_order_sha256':digest_ids(self.obs.index),'gene_id_order_sha256':digest_ids(a.var_names),'global_key_map':{'pca':'P','harmony':'H','umap':'U','neighbors_uns':'neighbors_cpu_n30'},'release_sidecar_sha256':{p.name:sha(p) for p in s.iterdir()}})
  return cfg
 def test_complete_synthetic_release_preflight(self):
  cfg=self.make_release();o,_,r=preflight(cfg,self.root/'validation.json');self.assertEqual(len(o),24);self.assertEqual(r['status'],'PASS_RELEASE_PREFLIGHT')
 def test_release_without_approval_rejected(self):
  cfg=self.make_release();p=Path(cfg['atlas']['manifest']);m=json.loads(p.read_text());m['level2_frozen']=False;atomic_json(p,m)
  with self.assertRaises(Blocked):preflight(cfg,self.root/'validation.json')
 def test_checksum_tamper_rejected(self):
  cfg=self.make_release();cfg['atlas']['input_sha256']='0'*64
  with self.assertRaises(Blocked):preflight(cfg,self.root/'validation.json')
 def test_no_fallback_input(self):
  cfg=self.make_release();cfg['atlas']['path']=str(self.root/'not_released.h5ad')
  with self.assertRaisesRegex(Blocked,'WAITING_FOR_FINAL_L2_ATLAS'):preflight(cfg,self.root/'validation.json')

if __name__=='__main__':unittest.main(verbosity=2)
