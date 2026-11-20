import sys, unittest, tempfile, json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
from src.covarnet import prepare,write_prepared,run_nmf,execute

def fixture():
    samples=[f's{i}' for i in range(24)]
    m=pd.DataFrame(dict(sample_uid=samples,biological_sample_uid=samples,donor_uid=[f'd{i}' for i in range(24)],donor_resolved=True,study=[f'study{i//8}' for i in range(24)],tissue='artery',condition='group',aggregate_approved=False))
    rows=[]
    for i,s in enumerate(samples):
        for j in range(6): rows.append([s,'Fibroblast' if j<3 else 'Macrophage',f'state{j}',25+(i*j+i+j)%11])
        rows.append([s,'Fibroblast','Unknown',20])
        rows.append([s,'Unknown','Unknown',20])
    c=pd.DataFrame(rows,columns=['sample_uid','level1_annotation','level2_annotation','count'])
    c['n_all_qc']=c.sample_uid.map(c.groupby('sample_uid')['count'].sum())
    ont=c[['level1_annotation','level2_annotation']].drop_duplicates().copy(); ont['covarnet_eligible']=~ont.level2_annotation.eq('Unknown')
    return c,m,ont

class TestCoVarNet(unittest.TestCase):
    def test_all_denominator(self):
        c,m,o=fixture(); r=prepare(c,m,o)
        self.assertEqual(r['manifest']['status'],'PREPARED_NOT_EXECUTED')
        self.assertEqual(r['manifest']['n_states'],6)
        self.assertTrue((r['prop'].sum()<1).all())
        self.assertEqual(r['parents'].loc['Fibroblast'].iloc[0],sum(c[(c.sample_uid=='s0')&(c.level1_annotation=='Fibroblast')]['count']))
    def test_missing_eligibility_review(self):
        c,m,o=fixture();r=prepare(c,m,o.drop(columns='covarnet_eligible'))
        self.assertIn('REVIEW_COVARNET_NODE_ELIGIBILITY',r['manifest']['blocking_reasons'])
        self.assertEqual(r['manifest']['n_states'],6)
    def test_named_rule_approval(self):
        c,m,o=fixture();r=prepare(c,m,o.drop(columns='covarnet_eligible'),{'approve_named_state_frequency_rule':True})
        self.assertEqual(r['manifest']['status'],'PREPARED_NOT_EXECUTED')
    def test_unknown_donor_not_pooled(self):
        c,m,o=fixture();m['donor_uid']='NR';m['donor_resolved']=False
        self.assertEqual(len(prepare(c,m,o)['profiles']),24)
    def test_unapproved_merge_blocked(self):
        c,m,o=fixture();m.loc[1,'donor_uid']='d0'
        with self.assertRaisesRegex(ValueError,'profile_aggregate_approved'): prepare(c,m,o)
    def test_merge_approved_preserves_cells(self):
        c,m,o=fixture();m.loc[1,'donor_uid']='d0';m['profile_aggregate_approved']=False;m.loc[:1,'profile_aggregate_approved']=True
        r=prepare(c,m,o);self.assertEqual(len(r['profiles']),23); self.assertEqual(r['all_profiles'].n_total_atlas_cells.sum(),c['count'].sum())
    def test_bad_count_rejected(self):
        c,m,o=fixture();c['count']=c['count'].astype(float);c.loc[0,'count']=1.2
        with self.assertRaisesRegex(ValueError,'integers'): prepare(c,m,o)
    def test_duplicate_sample_rejected(self):
        c,m,o=fixture()
        with self.assertRaisesRegex(ValueError,'Nonunique'): prepare(c,pd.concat([m,m.iloc[:1]]),o)
    def test_absent_parent_is_NA(self):
        c,m,o=fixture();c.loc[(c.sample_uid=='s0')&(c.level1_annotation=='Fibroblast'),'count']=0
        c['n_all_qc']=c.sample_uid.map(c.groupby('sample_uid')['count'].sum())
        r=prepare(c,m,o);pid=r['sample_profile_map'].set_index('sample_uid').loc['s0','covarnet_sample_uid']
        self.assertTrue(np.isnan(r['within_parent'].loc['state0',pid]))
    def test_missing_state_rejected(self):
        c,m,o=fixture()
        with self.assertRaisesRegex(ValueError,'missing from ontology'): prepare(c,m,o.iloc[1:])
    def test_hash_and_overwrite_guard(self):
        c,m,o=fixture();r=prepare(c,m,o)
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'prepared';write_prepared(r,p,{'fixture':True})
            with self.assertRaises(FileExistsError):write_prepared(r,p,{})
            with (p/'CoVarNet_node_counts.tsv').open('a') as f:f.write('\n')
            with self.assertRaisesRegex(ValueError,'hash mismatch'):run_nmf(p)
    def test_technical_approval_cannot_approve_distinct_bio(self):
        c,m,o=fixture();m.loc[1,'donor_uid']='d0';m['aggregate_approved']=True
        with self.assertRaisesRegex(ValueError,'profile_aggregate_approved'): prepare(c,m,o)
    def test_same_bio_requires_technical_approval(self):
        c,m,o=fixture();m.loc[1,['donor_uid','biological_sample_uid']]=['d0','s0']
        with self.assertRaisesRegex(ValueError,'technical partitions'): prepare(c,m,o)
        m.loc[:1,'aggregate_approved']=True
        self.assertEqual(len(prepare(c,m,o)['profiles']),23)
    def test_same_bio_conflicting_metadata_rejected(self):
        c,m,o=fixture();m.loc[1,'biological_sample_uid']='s0';m['aggregate_approved']=True
        with self.assertRaisesRegex(ValueError,'identity conflicts'): prepare(c,m,o)
    def test_cross_study_same_known_donor_one_bootstrap_block(self):
        c,m,o=fixture();m.loc[8,'donor_uid']='d0'
        r=prepare(c,m,o);x=r['sample_profile_map'].set_index('sample_uid')
        self.assertEqual(x.loc['s0','covarnet_donor_key'],x.loc['s8','covarnet_donor_key'])
        self.assertNotEqual(x.loc['s0','covarnet_sample_uid'],x.loc['s8','covarnet_sample_uid'])
    def test_execute_null_path_and_comparison_group_alias(self):
        c,m,o=fixture();m=m.rename(columns={'condition':'comparison_group'})
        with tempfile.TemporaryDirectory() as td:
            r=Path(td);(r/'counts').mkdir();(r/'cohort').mkdir()
            c.to_csv(r/'counts/sample_state_counts.tsv',sep='\t',index=False)
            m.to_csv(r/'cohort/cohort_design.tsv',sep='\t',index=False)
            o.to_csv(r/'counts/ontology.tsv',sep='\t',index=False)
            result=execute({'covarnet':{'ontology_path':None}},r)
            self.assertEqual(result['status'],'PREPARED_NOT_EXECUTED')
            self.assertEqual(result['source_identity']['cohort_column_aliases']['condition'],'comparison_group')
    def test_execute_blocked_raises_with_audit_preserved(self):
        c,m,o=fixture();o=o.drop(columns='covarnet_eligible')
        with tempfile.TemporaryDirectory() as td:
            r=Path(td);(r/'counts').mkdir();(r/'cohort').mkdir()
            c.to_csv(r/'counts/sample_state_counts.tsv',sep='\t',index=False)
            m.to_csv(r/'cohort/cohort_design.tsv',sep='\t',index=False)
            o.to_csv(r/'counts/ontology.tsv',sep='\t',index=False)
            with self.assertRaisesRegex(RuntimeError,'BLOCKED'): execute({},r)
            self.assertTrue((r/'covarnet/INPUT_MANIFEST.json').exists())
if __name__=='__main__':unittest.main()
