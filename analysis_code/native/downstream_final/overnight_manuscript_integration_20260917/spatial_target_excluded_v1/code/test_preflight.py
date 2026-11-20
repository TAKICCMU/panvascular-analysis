import json,sys,hashlib,ast
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad,torch
from exclusion import token,excluded_mask,assert_excluded
J=Path(__file__).resolve().parents[1];c=json.loads((J/'configs/config.json').read_text());exc=set(c['excluded_feature_tokens'])
for p in (J/'code').glob('*.py'):ast.parse(p.read_text())
x=pd.DataFrame({'stable_ensembl_id':['ENSG00000186665.8','ENSG_NOT_TARGET']},index=['C17orf58','CONTROL'])
assert excluded_mask(x,exc).tolist()==[True,False]
try:assert_excluded(x,exc)
except AssertionError:pass
else:raise AssertionError('Residual identity guard did not fail')
a=ad.read_h5ad(c['reference_counts'],backed='r');mask=excluded_mask(a.var,exc);retained=a.var.loc[~mask];assert_excluded(retained,exc)
assert a.obs.reference_state.nunique()==128
assert int(mask.sum())>0
nobs,nvars=a.shape;nremoved=int(mask.sum());a.file.close()
# No optimization/training: a single CUDA tensor validates the explicit device restoration.
torch.set_num_threads(4);assert torch.cuda.is_available();t=torch.arange(128,device='cuda').float();assert float((t+1).sum().cpu())==8256
result={'status':'PASS','source_reference_shape':[nobs,nvars],'target_feature_rows_removed':nremoved,'reference_retained_genes':len(retained),'target_id_count':len(c['target_ids']),'fixed_factors':128,'synthetic_versioned_ensembl_alias_test':'PASS','residual_alias_guard':'PASS','cuda_tensor_device':str(t.device),'cuda_device_name':torch.cuda.get_device_name(0),'trained_models':False,'full_gpu_smoke':'NOT_RUN; wrapper will run before production','biological_independence_boundary':'conditional on frozen upstream atlas labels; targets were not excluded from historical atlas annotation'}
(J/'evidence/PREFLIGHT_TEST.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
