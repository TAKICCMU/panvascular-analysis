import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import ast,json,tempfile,gzip
import numpy as np,pandas as pd
import validate_native as v
W=Path(__file__).resolve().parents[1];O=W.parents[1]
old=ast.parse((O/'code/scdrs_adapter.py').read_text());new=ast.parse((W/'code/coordinator.py').read_text())
def cmds(tree):
 out=[]
 for x in ast.walk(tree):
  if isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='cmd' for t in x.targets) and isinstance(x.value,ast.List):
   vals=[z.value for z in x.value.elts if isinstance(z,ast.Constant)]
   if 'compute-score' in vals or 'perform-downstream' in vals:out.append(vals)
 return sorted(out,key=str)
assert cmds(old)==cmds(new)
bh_old=next(x for x in old.body if isinstance(x,ast.FunctionDef) and x.name=='bh');bh_new=next(x for x in ast.parse(Path(v.__file__).read_text()).body if isinstance(x,ast.FunctionDef) and x.name=='bh')
assert ast.dump(bh_old)==ast.dump(bh_new)
with tempfile.TemporaryDirectory(prefix='scdrs_validator_fixture_') as tmp:
 run=Path(tmp);(run/'scores').mkdir();pd.DataFrame(np.zeros((4,1000)),index=['a','b','c','d'],columns=['ctrl_norm_score_'+str(i) for i in range(1000)]).to_csv(run/'scores/test.full_score.gz',sep='\t',compression='gzip')
 (run/'compute-score.log').write_text('--h5ad-file loaded: n_cell=4, n_gene=3000\nTrait=test, n_gene=25\n')
 s=v.score(run);assert s['n_scored_cells']==4 and s['n_control_columns']==1000
 (run/'L1').mkdir();pd.DataFrame({'group':['Mast','Fib','Excluded_Not_eligible'],'assoc_mcp':[.01,.04,.8],'hetero_mcp':[.3,.8,.8]}).to_csv(run/'L1/test.scdrs_group.fixture',sep='\t',index=False)
 g=v.group(run,'L1');assert g['n_eligible_groups']==2
 bad=pd.DataFrame(np.zeros((2,999)),index=['a','b'],columns=['ctrl_norm_score_'+str(i) for i in range(999)])
 bad.to_csv(run/'scores/test.full_score.gz',sep='\t',compression='gzip')
 try:v.score(run)
 except ValueError:pass
 else:raise AssertionError('999 controls accepted')
print(json.dumps(dict(status='PASS',tests=['exact_original_CLI_literal_options','exact_original_BH_AST','1000_control_complete_score_validation','original_eligible_group_scope','reject_999_control_incomplete_score'],first_live_score_recomputed=False),indent=2))
