"""Revision-specific synthetic assertions; no previous log is reused as evidence."""
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
import numpy as np
import pandas as pd
from src.covarnet import prepare, write_prepared
from src.counts import state_counts
from src.state_identity import build_state_identity_map, map_state_axis
from tests.test_covarnet import fixture

ROOT=Path(__file__).resolve().parents[1]

class CoVarNetRevisionTests(unittest.TestCase):
    def test_actual_cell_count_ledger_to_prepare_preserves_special_denominators(self):
        c,m,o=fixture()
        c.loc[c.level1_annotation.eq('Unknown'),'level2_annotation']='Not_assessed'
        # Construct cell rows and invoke the real cell-to-count algorithm, not a
        # pre-compatible manually generated matrix. Dispatcher ontology fill is
        # separately tested by the release->counts integration regression.
        obs=c.loc[c.index.repeat(c['count'])].reset_index(drop=True)
        counts=state_counts(obs)
        o=obs[['level1_annotation','level2_annotation']].drop_duplicates().copy()
        o['covarnet_eligible']=~o.level2_annotation.isin(['Unknown','Not_assessed'])
        r=prepare(counts,m,o)
        self.assertEqual(r['manifest']['status'],'PREPARED_NOT_EXECUTED')
        self.assertEqual(r['manifest']['phase'],'INPUTS_READY')
        self.assertEqual(r['manifest']['n_cells'],len(obs))
        self.assertNotIn('Unknown',r['counts'].index)
        self.assertNotIn('Not_assessed',r['counts'].index)
        self.assertEqual(r['parents'].sum().sum(),len(obs))
        self.assertTrue((r['prop'].sum(axis=0)<1).all())

    def test_actual_prepare_state_axis_roundtrip_with_punctuation(self):
        c,m,o=fixture()
        names={'state0':'C1QC / HLA-II + resident','state1':'APC: state (1)'}
        c['level2_annotation']=c.level2_annotation.replace(names)
        o['level2_annotation']=o.level2_annotation.replace(names)
        r=prepare(c,m,o)
        sm=r['state_identity_map']
        spatial=map_state_axis(r['counts'].index,sm,'state_id','safe_label')
        restored=map_state_axis(spatial,sm,'safe_label','state_id')
        self.assertEqual(list(restored),list(r['counts'].index))
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'prepared';man=write_prepared(r,p,{'fixture':True})
            self.assertIn('CoVarNet_state_identity_map.tsv',man['files'])
        with self.assertRaisesRegex(ValueError,'missing from final ontology'):
            map_state_axis(['misspelled state'],sm,'state_id')

    def test_explicit_ontology_ids_are_weights_axis(self):
        c,m,o=fixture()
        o['state_id']=['ID_'+str(i) for i in range(len(o))]
        c.loc[c.level2_annotation.isin(['state0','state3']),'level2_annotation']='shared text'
        o.loc[o.level2_annotation.isin(['state0','state3']),'level2_annotation']='shared text'
        r=prepare(c,m,o)
        ids=set(o.loc[o.covarnet_eligible,'state_id'])
        self.assertEqual(set(r['counts'].index),ids)
        self.assertEqual(set(r['mapping'].cell_state),ids)
        self.assertIn('level2_annotation',r['mapping'])
        self.assertNotIn('shared text',r['counts'].index)

    def test_counts_added_special_false_does_not_destroy_optional_eligibility_rule(self):
        c,m,o=fixture()
        # Exact output shape after counts concatenates an unqualified named
        # ontology with explicit False special rows. Blank represents TSV reload.
        for blank in [np.nan,'']:
            outont=o.copy();outont['covarnet_eligible']=outont.covarnet_eligible.astype(object)
            outont.loc[~outont.level2_annotation.eq('Unknown'),'covarnet_eligible']=blank
            for approved in [False,True]:
                with self.subTest(blank=repr(blank),approved=approved):
                    result=prepare(c,m,outont,{'approve_named_state_frequency_rule':approved})
                    self.assertEqual(result['manifest']['status'],'PREPARED_NOT_EXECUTED' if approved else 'BLOCKED')
                    self.assertEqual(result['manifest']['eligibility_basis'],'named_state_frequency_rule_APPROVED' if approved else 'named_state_frequency_rule_PROPOSED')
                    self.assertEqual(result['manifest']['n_states'],6)
                    self.assertEqual(result['manifest']['n_cells'],int(c['count'].sum()))
                    self.assertNotIn('Unknown',result['counts'].index)

    def test_partial_named_qualification_missing_stays_strict(self):
        c,m,o=fixture()
        for value in [None,np.nan,'']:
            outont=o.copy();outont['covarnet_eligible']=outont.covarnet_eligible.astype(object)
            outont.loc[outont.level2_annotation.eq('state0'),'covarnet_eligible']=value
            for approved in [False,True]:
                with self.subTest(value=repr(value),approved=approved),self.assertRaisesRegex(ValueError,'Partial named'):
                    prepare(c,m,outont,{'approve_named_state_frequency_rule':approved})

    def test_invalid_named_eligibility_never_coerces(self):
        c,m,o=fixture()
        for value in ['yes','False ','unknown','',None,np.nan,2]:
            bad=o.copy();bad['covarnet_eligible']=bad.covarnet_eligible.astype(object)
            bad.loc[bad.level2_annotation.eq('state0'),'covarnet_eligible']=value
            with self.subTest(value=value),self.assertRaisesRegex(ValueError,'[Ee]xplicit boolean'):
                prepare(c,m,bad)

    @unittest.skipUnless(shutil.which('Rscript'),'Native R unavailable')
    def test_native_r_descriptive_and_matching(self):
        script=r'''
args<-commandArgs(TRUE);source(args[1])
# Three donors each contribute two tissue profiles. n_profiles is descriptive,
# not an asserted independent-person sample size. No ordinary P/FDR survives.
m<-rbind(a=c(.1,.2,.4,.3,.8,.7),b=c(.3,.4,.2,.1,.9,.8),c=c(.9,.8,.7,.5,.3,.2))
r<-descriptive_profile_spearman(m,"raw_proportion",rep(c("d1","d2","d3"),each=2))
l<-descriptive_profile_spearman(log((m+1e-5)/(1-m+1e-5)),"logit_proportion")
stopifnot(all(r$n_profiles==6),all(r$n_resampling_identity_blocks==3L),!any(c("P","FDR","pvalue")%in%names(r)),
          all(r$inference=="DESCRIPTIVE_NO_P_OR_FDR"),max(abs(r$rho-l$rho))<1e-14)
ref<-cbind(R1=c(1,.1,0),R2=c(1,.2,0));rownames(ref)<-c("s1","s2","s3")
cand<-cbind(C1=c(1,.15,0),C2=c(0,0,1));rownames(cand)<-rownames(ref)
a<-compare_basis_review(ref,cand,"synthetic_repeat_donor",top_nodes=2L)
many<-a[a$matching_mode=="inherited_many_to_one_primary",]
one<-a[a$matching_mode=="one_to_one_hungarian_sensitivity",]
stopifnot(all(many$matched_module==1L),all(many$duplicate_match),
 length(unique(one$matched_module))==2L,!any(one$duplicate_match))
write.table(a,args[2],sep="\t",row.names=FALSE,quote=FALSE)
cat("SYNTHETIC_NATIVE_R_DESCRIPTIVE_MATCHING_PASS\n")
'''
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);(p/'test.R').write_text(script)
            env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
            done=subprocess.run(['Rscript',str(p/'test.R'),str(ROOT/'src/covarnet_review_helpers.R'),str(p/'out.tsv')],env=env,capture_output=True,text=True)
            self.assertEqual(done.returncode,0,done.stdout+done.stderr)
            self.assertIn('PASS',done.stdout)
            a=pd.read_csv(p/'out.tsv',sep='\t')
            self.assertEqual(set(a.matching_mode),{'inherited_many_to_one_primary','one_to_one_hungarian_sensitivity'})

    @unittest.skipUnless(shutil.which('Rscript'),'Native R unavailable')
    def test_hungarian_native_matches_exhaustive_optimum(self):
        rng=np.random.default_rng(7)
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)
            for i,(n,m) in enumerate([(2,2),(3,3),(4,6),(6,6)]):
                mat=rng.uniform(-1,1,(n,m));np.savetxt(p/f'{i}.tsv',mat,delimiter='\t')
            script='source(commandArgs(TRUE)[1]);for(i in 0:3){m<-as.matrix(read.table(paste0(i,".tsv"),sep="\\t"));write(maximum_one_to_one(m),paste0(i,".out"),ncolumns=nrow(m))}'
            (p/'test.R').write_text(script)
            done=subprocess.run(['Rscript','test.R',str(ROOT/'src/covarnet_review_helpers.R')],cwd=p,capture_output=True,text=True)
            self.assertEqual(done.returncode,0,done.stderr)
            for i in range(4):
                mat=np.loadtxt(p/f'{i}.tsv');assignment=np.loadtxt(p/f'{i}.out',dtype=int)-1
                optimum=max(mat[np.arange(mat.shape[0]),list(cols)].sum() for cols in itertools.permutations(range(mat.shape[1]),mat.shape[0]))
                self.assertAlmostEqual(mat[np.arange(len(assignment)),assignment].sum(),optimum,places=12)

    def test_active_r_relative_coefficient_name_and_no_naive_test(self):
        s=(ROOT/'src/covarnet_nmf.R').read_text()
        self.assertIn('sample_module_relative_coefficient_K',s)
        self.assertNotIn('sample_module_abundance_K',s)
        self.assertNotIn('cor.test(',s)
        self.assertIn('module_matching_collapse_audit',s)
        self.assertIn('matching_mode, reference_module',s)

if __name__=='__main__':unittest.main()
