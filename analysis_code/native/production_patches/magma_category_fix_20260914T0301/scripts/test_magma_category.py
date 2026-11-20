import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='2'
import sys,fcntl,json,time,hashlib
from pathlib import Path
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908');W=R/'production_launches/magma_category_fix_20260914T0301';C=W/'magma_adapter_v1';lock=open(R/'production_launches/execution_control_20260914/.auxiliary_compute.lock','a');print('Waiting auxiliary lock',flush=True);fcntl.flock(lock,fcntl.LOCK_EX)
sys.path.insert(0,str(C));import numpy as np,pandas as pd,h5py
from scipy import sparse
from anndata.io import read_elem
from src.external_modules import specificity_from_counts,identity_codes
passed=[]
v=pd.Series(pd.Categorical(['001','01','1',None,'unused_observed'],categories=['unused_level','unused_observed','1','01','001']))
orig=v.copy();mapping={'001':0,'01':1,'1':2};assert np.array_equal(identity_codes(v,mapping),[0,1,2,-1,-1]);pd.testing.assert_series_equal(v,orig);passed.append('category_missing_unused_and_leading_zero_exact_codes_no_mutation')
rng=np.random.default_rng(7);X=sparse.csr_matrix(rng.poisson(3,(13,4)));genes=pd.Index(['ENSG00001','ENSG00002','ENSG00003','UNMAPPED'])
obs=pd.DataFrame({'study':pd.Categorical(['001','01','1']*4+[None],categories=['unused','1','01','001'])});labels=pd.Series(pd.Categorical(['A']*6+['B']*6+['Unknown'],categories=['A','B','Unknown','unused_label']))
a=specificity_from_counts(X,obs,genes,labels,min_cells=2,min_studies=3,min_study_cells=1);oo=obs.copy();oo['study']=oo.study.astype(object);b=specificity_from_counts(X,oo,genes,labels.astype(object),min_cells=2,min_studies=3,min_study_cells=1)
for x,y in zip(a,b):
 x.attrs={};y.attrs={};pd.testing.assert_frame_equal(x,y)
passed.append('actual_specificity_function_category_vs_object_all_four_outputs_identical')
try:
 bad=obs.copy();bad.loc[0,'study']=None;specificity_from_counts(X,bad,genes,labels,min_cells=2,min_studies=3,min_study_cells=1)
 raise AssertionError('valid biological state missing study unexpectedly accepted')
except ValueError as e:assert 'Missing study identity' in str(e)
passed.append('missing_study_in_included_cells_still_blocks')
try:
 bad=X.copy().astype(float);bad.data[0]=-.5;specificity_from_counts(bad,obs,genes,labels,min_cells=2,min_studies=3,min_study_cells=1)
 raise AssertionError('invalid raw counts accepted')
except ValueError:pass
passed.append('actual_function_raw_count_validation_preserved')
with h5py.File(R/'panvascular_CPU_theta2_L1L2_final_v1.h5ad','r') as f:
 studies=read_elem(f['obs/study']);l1=read_elem(f['obs/level1_annotation']);n=len(studies)
ss=pd.Series(studies);ll=pd.Series(l1);before_s=ss.copy();before_l=ll.copy()
for series in [ss,ll]:
 keys=sorted(series.dropna().unique().tolist());mapr={k:i for i,k in enumerate(keys)};got=identity_codes(series,mapr);expected=series.astype(object).map(mapr).fillna(-1).to_numpy(int);assert np.array_equal(got,expected);assert np.all(got>=0)
pd.testing.assert_series_equal(ss,before_s);pd.testing.assert_series_equal(ll,before_l);passed.append('real_1175115_native_categorical_study_and_L1_encoding_entry_pass_no_mutation')
from src.dispatcher import code_digest
base=R/'downstream_CPU_L1L2_v2_1_release_compat_20260914';assert code_digest(base)=='7ecd30a62e3528a4c3c6ec8329124c8f2754ff90f65c8b28e75fcfc0f7e695df';changed=[]
for d in ['src','scripts','external_legacy']:
 for p in (C/d).rglob('*'):
  if p.is_file() and p.suffix in ['.py','.R','.sh'] and p.read_bytes()!=(base/p.relative_to(C)).read_bytes():changed.append(str(p.relative_to(C)))
assert changed==['src/external_modules.py'];passed.append('only_one_cloned_source_file_changed_original_code_digest_unchanged')
out={'status':'PASS','tests_passed':len(passed),'tests':passed,'n_real_obs':n,'no_real_expression_analysis':True,'new_code_digest':code_digest(C),'baseline_code_digest':code_digest(base),'changed_files':changed};(W/'evidence/CATEGORY_TEST_RESULTS.json').write_text(json.dumps(out,indent=2));print(json.dumps(out),flush=True);fcntl.flock(lock,fcntl.LOCK_UN)
