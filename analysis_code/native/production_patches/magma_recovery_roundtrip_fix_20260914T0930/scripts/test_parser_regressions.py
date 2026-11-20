import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='2'
from pathlib import Path
import sys,json,importlib.util
import pandas as pd,numpy as np
from scipy.stats import false_discovery_control
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908');W=R/'production_launches/magma_recovery_roundtrip_fix_20260914T0930';C=W/'magma_adapter_v3';S=R/'production_launches/magma_category_fix_20260914T0301/outputs/runs/20260913T190814_808850_7ccf71427e40/magma';sys.path.insert(0,str(C/'src'))
from magma_recovery import parse_conditional_results,_read_and_validate
m=pd.read_csv(S/'L1_variables.tsv',sep='\t');r=pd.read_csv(S/'primary_r0.1_L1_cell_weighted.gsa.out',sep=r'\s+',comment='#');passed=[]
t,n,models=_read_and_validate(S/'primary_r0.1_L1_cell_weighted.gsa.out',m,'actual_smoke_family');assert len(t)==len(n)==len(models)==15;assert not t.VARIABLE.eq('Average').any();assert n.VARIABLE.eq('Average').all();assert set(t.group)==set(m.group);assert np.allclose(t.FDR_BH,false_discovery_control(t.P.to_numpy(),method='bh'),rtol=1e-14,atol=1e-14);assert np.array_equal(t.P.to_numpy(),r.loc[r.VARIABLE.ne('Average'),'P'].to_numpy());passed.append('real_gsa_15_targets_15_condition_rows_preserved_one_to_one_BH_target_only')
def rejects(name,frame,msg):
 try:parse_conditional_results(frame,m,'synthetic_mutation')
 except ValueError as e:
  assert msg in str(e),(name,str(e));passed.append(name);return
 raise AssertionError(name+' unexpectedly accepted')
rejects('duplicate_target_variable_model_rejected',pd.concat([r,r.loc[r.VARIABLE.eq('S0000')]]),'Duplicate variable-model')
extra=r.loc[r.VARIABLE.eq('S0000')].copy();extra.MODEL=99;rejects('same_target_in_second_model_rejected',pd.concat([r,extra]),'Target appears in repeated models')
bad=r.copy();bad.loc[bad.VARIABLE.eq('S0001'),'MODEL']=1;rejects('two_targets_in_same_model_rejected',bad,'Multiple target variables')
bad=r.copy();bad.loc[bad.VARIABLE.eq('S0000'),'VARIABLE']='UNEXPECTED_GENESET';rejects('unknown_variable_rejected',bad,'Unexpected MAGMA output variable')
rejects('duplicate_nuisance_model_rejected',pd.concat([r,r.loc[(r.VARIABLE.eq('Average'))&(r.MODEL.eq(1))]]),'Duplicate variable-model')
rejects('missing_target_rejected',r.loc[r.VARIABLE.ne('S0000')],'Missing mapped target')
rejects('missing_condition_rejected',r.loc[~(r.VARIABLE.eq('Average')&r.MODEL.eq(1))],'Every target model requires')
bad=r.copy();bad.loc[0,'P']=np.nan;rejects('nonfinite_P_rejected',bad,'Nonfinite MAGMA output')
t.to_csv(W/'evidence/REAL_GSA_PARSED_TARGETS_TEST.tsv',sep='\t',index=False);n.to_csv(W/'evidence/REAL_GSA_CONDITIONS_TEST.tsv',sep='\t',index=False)
# The unchanged category/specificity implementation must still pass the previous tests.
old=R/'production_launches/magma_category_fix_20260914T0301/magma_adapter_v1/src/external_modules.py';assert (C/'src/external_modules.py').read_bytes()==old.read_bytes();passed.append('specificity_and_category_fix_code_exactly_unchanged')
sys.path.insert(0,str(C));from src.dispatcher import code_digest
result={'status':'PASS','passed':len(passed),'tests':passed,'new_code_digest':code_digest(C),'source_raw_gsa':str(S/'primary_r0.1_L1_cell_weighted.gsa.out'),'full_scientific_run_completed':False};(W/'evidence/PARSER_TEST_RESULTS.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
