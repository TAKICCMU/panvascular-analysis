import sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np,pandas as pd
from scipy import sparse
from external_modules import *
class ExternalTests(unittest.TestCase):
    def setUp(self):
        self.obs=pd.DataFrame({'study':['S1','S1','S2','S2','S3','S3'],
          'sample_uid':['a','a','b','b','c','c'],'biological_sample_uid':['a','a','b','b','c','c'],
          'level1_annotation':['T','T','T','T','T','T'],
          'level2_annotation':['x','y','x','y','x','y']},index=[f'c{i}' for i in range(6)])
        self.X=sparse.csr_matrix([[8,2],[1,9],[7,3],[2,8],[6,4],[3,7]])
    def test_specificity_chunk_equivalence_and_simple_reference(self):
        labels=qualified_states(self.obs)
        args=(self.X,self.obs,['ENSG0001','ENSG0002'],labels)
        a,b,m,e=specificity_from_counts(*args,chunk=1,min_cells=1,min_studies=1,min_study_cells=1)
        c,d,_,_=specificity_from_counts(*args,chunk=100,min_cells=1,min_studies=1,min_study_cells=1)
        np.testing.assert_allclose(a.iloc[:,2:],c.iloc[:,2:]);np.testing.assert_allclose(b.iloc[:,2:],d.iloc[:,2:])
        raw=self.X.toarray();log=np.log1p(raw/raw.sum(axis=1)[:,None]*1e6)
        means=np.stack([log[[0,2,4]].mean(axis=0),log[[1,3,5]].mean(axis=0)])
        np.testing.assert_allclose(a.iloc[:,2:].to_numpy(),(means/means.sum(axis=0)).T)
    def test_parent_qualified_no_label_collision(self):
        self.obs.loc['c0','level1_annotation']='NK'
        self.assertNotEqual(qualified_states(self.obs).loc['c0'],qualified_states(self.obs).loc['c2'])
    def test_unknown_not_estimable_not_false_zero(self):
        self.obs.loc['c0','level1_annotation']='Unknown'
        t=communication_eligibility(self.obs,min_cells=2)
        self.assertFalse(t.loc[t.biological_sample_uid=='a','eligible'].any())
        self.assertFalse((t.sender=='Unknown').any())
    def test_bad_metadata_blocked(self):
        self.obs.loc['c0','study']='wrong'
        with self.assertRaises(ValueError):validate_sample_metadata(self.obs,'sample_uid',['study'])
    def test_duplicate_symbol_counts_aggregate(self):
        gm=pd.DataFrame({'feature_id':['g1','g2'],'symbol':['G','G']})
        x,n=prepare_scdrs_counts(self.X,self.obs,gm)
        np.testing.assert_array_equal(x.toarray(),np.full((6,1),10))
    def test_missing_resource_blocked_and_hash_mismatch(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'f';p.write_text('x')
            r=audit_resources([{'path':str(p),'expected_sha256':'bad'},{'path':str(p)+'missing'}],True)
            self.assertEqual(r[0]['status'],'BLOCKED_HASH_MISMATCH');self.assertEqual(r[1]['status'],'BLOCKED_MISSING_RESOURCE')
    def test_target_missing_evidence_stays_na(self):
        out=target_evidence_join({'de':pd.DataFrame({'gene_id':['G1'],'q':[0.1]}),'coloc':pd.DataFrame({'gene_id':['G2'],'h4':[0.9]})})
        self.assertTrue(pd.isna(out.loc[out.gene_id=='G1','coloc__h4']).all())
    def test_raw_input_invalid_rejected(self):
        X=self.X.astype(float);X.data/=3
        with self.assertRaises(ValueError):specificity_from_counts(X,self.obs,['ENSG0001','ENSG0002'],qualified_states(self.obs),min_cells=1,min_studies=1,min_study_cells=1)
    def test_native_commands_preserve_settings(self):
        cmd=scdrs_commands('scdrs','input','cov','genes','out')[0]
        self.assertEqual(cmd[cmd.index('--n-ctrl')+1],'1000')
        self.assertIn('condition=Average',magma_command('m','g','s','out'))
    def test_unmapped_gene_keeps_full_library_denominator(self):
        X=sparse.hstack([self.X,sparse.csr_matrix([[100],[0],[50],[0],[10],[0]])]).tocsr()
        a,_,_,_=specificity_from_counts(X,self.obs,['ENSG0001','ENSG0002',None],qualified_states(self.obs),min_cells=1,min_studies=1,min_study_cells=1)
        full=X.toarray();log=np.log1p(full/full.sum(axis=1)[:,None]*1e6)
        means=np.stack([log[[0,2,4],:2].mean(axis=0),log[[1,3,5],:2].mean(axis=0)])
        np.testing.assert_allclose(a.iloc[:,2:],(means/means.sum(axis=0)).T)
        self.assertEqual(len(a),2)
        self.assertEqual(int(a.attrs['gene_mapping_audit'].retained_in_specificity_output.sum()),2)
if __name__=='__main__':unittest.main()
