"""Native tiny AnnData and typed identity regressions; no production analysis."""
import unittest,tempfile,json,copy,hashlib
from pathlib import Path
import pandas as pd,numpy as np,anndata as ad
from scipy import sparse
from src.core import Blocked,sha,read_frame,digest_ids
from src.cohort import identical_identity,annotate_units,apply_unit_review,cohort_from_obs,read_unit_review,FIELDS
from src.counts import state_counts,biological_counts,coverage_table
from src.expression import marker_summary
from src.figures import export_figures,dotplot
from src.release_compat import encoded_digest,LENGTH_ENCODING,NEWLINE_ENCODING,release_digest,bound_file,validate_metadata_sources,validate_registry,validate_manifest_binding

class CompatTests(unittest.TestCase):
 def setUp(self):self.tmp=tempfile.TemporaryDirectory();self.p=Path(self.tmp.name)
 def tearDown(self):self.tmp.cleanup()
 def typed(self,dtype,missing=False):return pd.Series(['001','供者甲',None if missing else 'd3'],dtype=dtype)
 def test_category_no_missing(self):self.assertTrue(identical_identity(self.typed('category'),self.typed('object')))
 def test_category_missing(self):self.assertTrue(identical_identity(self.typed('category',True),self.typed('string',True)))
 def test_object_none_pdNA(self):self.assertTrue(identical_identity(pd.Series(['001',None]),pd.Series(['001',pd.NA],dtype='string')))
 def test_nullable_string(self):self.assertTrue(identical_identity(self.typed('string'),self.typed('string')))
 def test_donor_changed(self):self.assertFalse(identical_identity(pd.Series(['001']),pd.Series(['002'])))
 def test_missing_position(self):self.assertFalse(identical_identity(pd.Series(['001',None]),pd.Series([None,'001'])))
 def test_missing_not_empty(self):self.assertFalse(identical_identity(pd.Series([None]),pd.Series([''])))
 def test_literal_nan_not_missing(self):self.assertFalse(identical_identity(pd.Series([None]),pd.Series(['nan'])))
 def test_leading_zero_changed(self):self.assertFalse(identical_identity(pd.Series(['001']),pd.Series(['1'])))
 def test_exact_no_trim_casefold(self):self.assertFalse(identical_identity(pd.Series(['A ']),pd.Series(['a'])))
 def meta(self):
  n=12;o=pd.DataFrame(index=pd.Index(['c'+str(i) for i in range(n)],name='cell_id'))
  for k in FIELDS:o[k]=['x']*n
  o['sample_uid']=np.repeat(['001','002','003'],4);o['donor_uid']=np.repeat(['0001','供者甲',None],4);o['biological_sample_uid']=np.repeat(['b01','b02','b03'],4);o['age']=[40.]*n;o['sex']=['female']*8+[None]*4
  o['level1_annotation']=['Fibroblast']*6+['Unknown']*2+['Schwann cell']*4
  o['level2_annotation']=['G1 Fib']*4+['G2 Fib']*2+['Unknown']*2+['Not_assessed']*4
  for k in o.select_dtypes('object'):o[k]=pd.Categorical(o[k])
  return o
 def review(self,m):
  r=m.copy();r['aggregate_approved']=True;r['inference_eligible']=False;r['review_evidence']='synthetic single sample unit';r['design_status']='APPROVED';return r
 def test_key_alignment(self):
  m=cohort_from_obs(self.meta());r=self.review(m).iloc[::-1];out=apply_unit_review(m,r);self.assertEqual(list(out.sample_uid),list(m.sample_uid))
 def test_true_donor_mutation_blocked(self):
  m=cohort_from_obs(self.meta());r=self.review(m);r['donor_uid']=r.donor_uid.astype('string');r.loc[0,'donor_uid']='0002'
  with self.assertRaises(Blocked):apply_unit_review(m,r)
 def test_duplicate_sample_blocked(self):
  m=cohort_from_obs(self.meta());r=self.review(m);r['sample_uid']=r.sample_uid.astype('string');r.loc[0,'sample_uid']=r.loc[1,'sample_uid']
  with self.assertRaises(Blocked):apply_unit_review(m,r)
 def test_native_csv_leading_zero(self):
  self.p.joinpath('review.tsv').write_text('sample_uid\tdonor_uid\tbiological_sample_uid\n001\t0001\tb01\n002\t\tb02\n');r=read_unit_review(self.p/'review.tsv');self.assertEqual(r.loc[0,'sample_uid'],'001');self.assertTrue(pd.isna(r.loc[1,'donor_uid']))
 def test_hash_encodings_distinct(self):self.assertNotEqual(encoded_digest(['a','甲'],LENGTH_ENCODING),encoded_digest(['a','甲'],NEWLINE_ENCODING))
 def test_hash_unicode_byte_length(self):
  h=hashlib.sha256();b='甲'.encode();h.update((3).to_bytes(8,'little'));h.update(b);self.assertEqual(h.hexdigest(),encoded_digest(['甲'],LENGTH_ENCODING))
 def test_hash_generic_unchanged(self):self.assertEqual(digest_ids(['a']),hashlib.sha256(b'a\n').hexdigest())
 def test_hash_shuffle_delete_add(self):
  a=encoded_digest(['a','b'],LENGTH_ENCODING)
  for v in [['b','a'],['a'],['a','b','c']]:self.assertNotEqual(a,encoded_digest(v,LENGTH_ENCODING))
 def test_unknown_encoding_rejected(self):
  with self.assertRaises(Blocked):encoded_digest(['a'],'guess')
 def test_missing_encoding_rejected(self):
  with self.assertRaises(Blocked):release_digest(['a'],{'support_revision':'acceptance_r1'},'cell_id_order_sha256',True)
 def test_compat_manifest_required(self):
  with self.assertRaises(Blocked):release_digest(['a'],{},'cell_id_order_sha256',True)
 def test_bound_file_missing(self):
  with self.assertRaises(Blocked):bound_file({'path':str(self.p/'missing'),'sha256':'0'*64})
 def test_bound_file_mutation(self):
  p=self.p/'source';p.write_text('v1');r={'path':str(p),'sha256':sha(p)};p.write_text('v2')
  with self.assertRaises(Blocked):bound_file(r)
 def metadata_source(self):
  p=self.p/'patch.csv';p.write_text('sample_uid,donor_uid,age\ns1,001,40\n');s=self.p/'sample.csv';s.write_bytes(p.read_bytes());b=self.p/'checks.json';b.write_text(json.dumps({'inputs/patch.csv':sha(p)}));ev=self.p/'evidence.txt';ev.write_text('synthetic recorded import, not historical project evidence')
  files=[]
  for role in ['imported_metadata_csv','import_code','import_run','bundle_checksums','sample_metadata_snapshot','metadata_provenance','user_confirmation','audit_A','audit_B','merge_record','public_GSA_evidence']:
   x={'imported_metadata_csv':p,'sample_metadata_snapshot':s,'bundle_checksums':b}.get(role,ev);files.append({'role':role,'path':str(x),'sha256':sha(x),'scope':'synthetic fixture','version':'test','evidence_type':'synthetic'})
  m=self.p/'metadata.json';m.write_text(json.dumps({'schema':'metadata_review_sources_acceptance_r1','files':files,'import_link_validation':{'status':'PASS','bundle_member':'inputs/patch.csv'}}));return m,ev
 def test_metadata_sources_native_good(self):m,_=self.metadata_source();self.assertEqual(validate_metadata_sources(m)['sample_count'],1)
 def test_metadata_sources_missing(self):
  m,e=self.metadata_source();e.unlink()
  with self.assertRaises(Blocked):validate_metadata_sources(m)
 def test_metadata_sources_tamper(self):
  m,e=self.metadata_source();e.write_text('changed')
  with self.assertRaises(Blocked):validate_metadata_sources(m)
 def registry(self):
  a=ad.AnnData(sparse.eye(3,format='csr'),obs=pd.DataFrame({'cl':pd.Categorical(['C0','C1','C0'])},index=['001','002','甲']),var=pd.DataFrame(index=['g1','g2','g3']));a.obsm['U']=np.zeros((3,2));a.obsp['C']=sparse.eye(3,format='csr');a.obsp['D']=sparse.eye(3,format='csr');a.uns['N']={'connectivities_key':'C','distances_key':'D'};p=self.p/'local.h5ad';a.write_h5ad(p)
  mp=self.p/'map.tsv';a.obs.to_csv(mp,sep='\t',index_label='cell_id');v=self.p/'task.json';v.write_text('{}')
  r=pd.DataFrame([dict(role='lineage_analysis',branch='Test',path=str(p),sha256=sha(p),version=str(self.p),version_source=str(v),version_source_sha256=sha(v),n_obs=3,n_vars=3,label_key='cl',graph_key='N',embedding_key='U',mapping_path=str(mp),mapping_sha256=sha(mp))]);rp=self.p/'reg.tsv';r.to_csv(rp,sep='\t',index=False);return rp,r
 def test_registry_native_good(self):p,_=self.registry();self.assertEqual(validate_registry(p)['rows'],1)
 def test_registry_invalid_graph(self):
  p,r=self.registry();r['graph_key']='global_fake';r.to_csv(p,sep='\t',index=False)
  with self.assertRaises(Blocked):validate_registry(p)
 def test_registry_invalid_embedding(self):
  p,r=self.registry();r['embedding_key']='fake';r.to_csv(p,sep='\t',index=False)
  with self.assertRaises(Blocked):validate_registry(p)
 def test_registry_changed_mapping(self):
  p,r=self.registry();q=Path(r.loc[0,'mapping_path']);q.write_text(q.read_text().replace('C1','C9'));r['mapping_sha256']=sha(q);r.to_csv(p,sep='\t',index=False)
  with self.assertRaises(Blocked):validate_registry(p)
 def test_native_five_module_smoke_with_missing_category(self):
  o=self.meta();x=sparse.csr_matrix(np.arange(36).reshape(12,3)%5,dtype='int32');a=ad.AnnData(x,obs=o,var=pd.DataFrame({'gene_symbol':['G1','G2','G3']},index=['G1','G2','G3']));z=x.astype(float).multiply((10000/np.maximum(x.sum(axis=1).A1,1))[:,None]).tocsr();z.data=np.log1p(z.data);a.layers['log1p']=z;a.obsm['U']=np.arange(24).reshape(12,2).astype(float);p=self.p/'atlas.h5ad';a.write_h5ad(p);h=sha(p)
  obs=read_frame(p);self.assertEqual(str(obs.donor_uid.dtype),'category');cohort=cohort_from_obs(obs);cohort=apply_unit_review(cohort,self.review(cohort));self.assertTrue(cohort.donor_uid.isna().any())
  counts=state_counts(obs);self.assertEqual(counts.n_cells.sum(),12);self.assertTrue(counts.loc[counts.n_parent.eq(0),'prop_within_parent'].isna().all());self.assertTrue(counts.loc[counts.n_parent.gt(0)&counts.n_cells.eq(0),'prop_within_parent'].eq(0).all());self.assertEqual(biological_counts(obs,cohort).n_cells.sum(),12);coverage_table(obs,cohort)
  cfg={'atlas':{'path':str(p),'level1_key':'level1_annotation','level2_key':'level2_annotation'},'execution':{'block_rows':3,'seed':0},'figures':{},'outputs':{'release_support':str(self.p/'support')}}
  marker_summary(cfg,obs,self.p/'markers');(self.p/'preflight.json').write_text(json.dumps({'global_key_map':{'umap':'U'},'sha256':h}));export_figures(cfg,obs,self.p)
  panel=self.p/'panel.tsv';pd.DataFrame({'gene_id':['G1','G2'],'panel':['Fib','Fib']}).to_csv(panel,sep='\t',index=False);cfg['figures']['marker_panel']=str(panel)
  support=self.p/'support';support.mkdir();pd.DataFrame({'level1_annotation':['Fibroblast']*2,'level2_annotation':['G1 Fib','G2 Fib'],'marker_anchor':['G1','G2'],'approved':[True,True]}).to_csv(support/'level2_ontology.tsv',sep='\t',index=False)
  dotplot(cfg,self.p,level='L1');dotplot(cfg,self.p,level='L2');self.assertEqual(h,sha(p));self.assertTrue((self.p/'figures/L2_DotPlot.pdf').is_file())
 def compat_fixture(self):
  from tests.test_core import CoreTests
  fixture=CoreTests();fixture.setUp();self.addCleanup(fixture.tearDown);cfg=fixture.make_release();support=Path(cfg['outputs']['release_support']);oldpath=Path(cfg['atlas']['manifest']);old=json.loads(oldpath.read_text())
  pub=self.p/'committed/audit/finalize_level2.py';pub.parent.mkdir(parents=True);pub.write_text('synthetic publisher: length-prefixed UTF8 fixture')
  old['finalization_run']=str(pub.parent.parent)
  for key,values in [('cell_id_order_sha256',fixture.obs.index),('gene_id_order_sha256',read_frame(fixture.path,'var').index),('annotation_sha256',fixture.obs.level2_annotation),('level1_annotation_sha256',fixture.obs.level1_annotation)]:old[key]=encoded_digest(values,LENGTH_ENCODING)
  oldpath.write_text(json.dumps(old));m=copy.deepcopy(old);m.update(support_revision='acceptance_r1',downstream_execution_authorized=False,original_manifest={'path':str(oldpath),'sha256':sha(oldpath)},original_manifest_sha256=sha(oldpath),publisher_script={'path':str(pub),'sha256':sha(pub)},hash_encoding={k:LENGTH_ENCODING for k in ['cell_id_order_sha256','gene_id_order_sha256','annotation_sha256','level1_annotation_sha256']},release_support_dir=str(support))
  meta,_=self.metadata_source();mdata=json.loads(meta.read_text());patch=self.p/'actual_patch.csv';patch.write_bytes((support/'sample_metadata.csv').read_bytes());bundle=self.p/'checks.json';bundle.write_text(json.dumps({'inputs/patch.csv':sha(patch)}))
  for x in mdata['files']:
   if x['role'] in ['imported_metadata_csv','sample_metadata_snapshot','bundle_checksums']:
    q={'imported_metadata_csv':patch,'sample_metadata_snapshot':support/'sample_metadata.csv','bundle_checksums':bundle}[x['role']];x.update(path=str(q),sha256=sha(q))
  (support/'metadata_review_sources.json').write_text(json.dumps(mdata));rp,_=self.registry();(support/'source_registry.tsv').write_bytes(rp.read_bytes());m['release_sidecar_sha256']={x.name:sha(x) for x in support.iterdir()}
  mp=self.p/'compat.json';mp.write_text(json.dumps(m));cfg['atlas'].update(manifest=str(mp),release_compatibility_required=True)
  contract=Path(cfg['atlas']['schema_contract']);ct=json.loads(contract.read_text());ct['release_sidecars'].update({'metadata_review_sources.json':[],'source_registry.tsv':[]});contract.write_text(json.dumps(ct));return cfg,m,mp,fixture
 def test_native_compat_preflight_end_to_end(self):
  from src.core import preflight
  c,m,mp,f=self.compat_fixture();_,_,r=preflight(c,self.p/'pass.json');self.assertEqual(r['status'],'PASS_RELEASE_PREFLIGHT');self.assertIn('compatibility_checks',r)
 def test_bound_encoding_unknown_full_preflight(self):
  from src.core import preflight
  c,m,mp,f=self.compat_fixture();m['hash_encoding']['cell_id_order_sha256']='guess';mp.write_text(json.dumps(m))
  with self.assertRaisesRegex(Blocked,'encoding'):preflight(c,self.p/'no.json')
 def test_bound_source_sha_mismatch_full_preflight(self):
  from src.core import preflight
  c,m,mp,f=self.compat_fixture();Path(m['publisher_script']['path']).write_text('changed')
  with self.assertRaisesRegex(Blocked,'SHA mismatch'):preflight(c,self.p/'no.json')
 def test_axis_reordered_full_preflight(self):
  from src.core import preflight
  c,m,mp,f=self.compat_fixture();a=ad.read_h5ad(f.path);a=a[::-1].copy();a.write_h5ad(f.path);h=sha(f.path)
  # Rebind whole-file bytes to isolate ID-order check; original and compat expected
  # order hashes stay unchanged. Fixture-only, never done to production files.
  op=Path(m['original_manifest']['path']);om=json.loads(op.read_text());om['h5ad_sha256']=h;op.write_text(json.dumps(om));m.update(h5ad_sha256=h,original_manifest={'path':str(op),'sha256':sha(op)},original_manifest_sha256=sha(op));mp.write_text(json.dumps(m));c['atlas']['input_sha256']=h;Path(c['atlas']['checksum']).write_text(h+'  fixture.h5ad')
  with self.assertRaisesRegex(Blocked,'cell_id_order'):preflight(c,self.p/'no.json')
 def test_dependency_manifest_binds_historical_records_as_leaves(self):
  from src.review_dependencies import reviewed_dependencies
  c,m,mp,f=self.compat_fixture();s=Path(c['outputs']['release_support']);p=s/'metadata_review_sources.json';d=json.loads(p.read_text());ev=self.p/'historical.json';ev.write_text(json.dumps({'stage':'/retired_nonexistent/staging.partial.h5ad'}));d['files'].append({'role':'historical_record','path':str(ev),'sha256':sha(ev),'scope':'fixture','version':'test','evidence_type':'synthetic'});p.write_text(json.dumps(d));m['release_sidecar_sha256'][p.name]=sha(p);mp.write_text(json.dumps(m))
  contract=Path(c['atlas']['schema_contract']);v=json.loads(contract.read_text());v['write_policy']={'partial_filename':'never_created.partial.h5ad'};contract.write_text(json.dumps(v))
  cp=self.p/'config.json';cp.write_text(json.dumps(c));dig,proof=reviewed_dependencies(c,cp);self.assertIn(str(ev),proof['files'])
  self.assertIn(str(contract),proof['files']);v['write_policy']['partial_filename']='different.partial.h5ad';contract.write_text(json.dumps(v));second,_=reviewed_dependencies(c,cp);self.assertNotEqual(dig,second)
  ev.write_text('changed historical evidence')
  with self.assertRaisesRegex(Blocked,'SHA differs'):reviewed_dependencies(c,cp)

if __name__=='__main__':unittest.main(verbosity=2)
