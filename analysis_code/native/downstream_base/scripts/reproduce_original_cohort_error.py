from pathlib import Path
import sys,json,traceback
sys.dont_write_bytecode=True
B=Path(__file__).resolve().parents[1];R=B.parent;OLD=R/'downstream_CPU_L1L2_v2_1_revision_20260913';sys.path.insert(0,str(OLD))
import pandas as pd
from src.core import read_frame,sha
from src.cohort import cohort_from_obs,apply_unit_review,FIELDS
p=R/'panvascular_CPU_theta2_L1L2_final_v1.h5ad'
r={'original_code':str(OLD/'src/cohort.py'),'original_code_sha256':sha(OLD/'src/cohort.py'),'input':str(p),'matrix_read':False}
try:
 obs=read_frame(p,columns=FIELDS);m=cohort_from_obs(obs);review=pd.read_csv(OLD/'configs/cohort_design.tsv',sep='\t');apply_unit_review(m,review);r['status']='UNEXPECTED_PASS'
except Exception as e:r.update(status='ORIGINAL_FAILURE_REPRODUCED',error=str(e),traceback=traceback.format_exc())
(B/'evidence/ORIGINAL_COHORT_FAILURE_REPRODUCED.json').write_text(json.dumps(r,indent=2));print(json.dumps(r),flush=True)
assert r['status']=='ORIGINAL_FAILURE_REPRODUCED' and 'Categorical' in r['error']
