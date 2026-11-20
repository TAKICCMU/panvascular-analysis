import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='2'
import sys,json
from pathlib import Path
W=Path(__file__).resolve().parents[1];R=W.parents[1];sys.path.insert(0,str(W/'code'))
from scdrs_adapter import native_sequence,sha
out=W/'tests/actual_obs_smoke_v1';result=json.loads((out/'SMOKE_RESULT.json').read_text());assert result['status']=='PASS_NATIVE_SHARED_ADAPTER_ACTUAL_OBS'
m=json.loads((R/'enrichment/scdrs_label_independent/genesets/GENESETS_FROZEN.json').read_text());s=[{'name':v['scenario'],'geneset':v['gs_file'],'sha256':v['gs_sha256']} for v in m['scenarios']]
native_sequence(out/'derived.h5ad',out/'native_covariates.cov.gz',s,out/'native','/home/anzhen_vhd/miniconda3/envs/scdrs_env/bin/scdrs',result['group_keys'],2,validate_existing=True)
report={'status':'PASS','shared_adapter_sha256':sha(W/'code/scdrs_adapter.py'),'validation':'Exact same derived-builder and native commands passed; added validators applied to all actual native outputs without redundant scientific recomputation','traits':[json.loads((out/'native'/v['name']/'SCHEMA_VALIDATION.json').read_text()) for v in s]}
(W/'tests/FINAL_SMOKE_VALIDATION.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
