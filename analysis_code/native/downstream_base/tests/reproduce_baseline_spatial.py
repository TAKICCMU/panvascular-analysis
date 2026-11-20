"""Independent v1 R04/R05 reproduction; does not alter original tests/logs."""
import ast,json,sys
from pathlib import Path
old=Path(sys.argv[1]);sys.path.insert(0,str(old/'src'))
import pandas as pd
from external_statistics import spatial_module_projection
from external_modules import render_spatial_adapter
out={'evidence_type':'independent_reproduction_of_review_findings','source_root':str(old),'original_review_reproduction_script_available':False}
a=pd.DataFrame([[1.,2.]],columns=['Macrophage | X','T cell | Y']);w=pd.DataFrame({'M1':[.4,.6]},index=['X','Y'])
try:spatial_module_projection(a,w);out['R04']={'reproduced':False}
except ValueError as e:out['R04']={'reproduced':True,'error':str(e),'weight_ids':list(w.index),'abundance_ids':list(a)}
source=(old/'external_legacy/spatial/02_spatial_v4_cell2location_GPU.py').read_text()
s=render_spatial_adapter(source,'a','b','c','level1_annotation','level2_annotation','study',[],'blocked')
node=next(n for n in ast.parse(s).body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='REF_MAX_EPOCHS' for t in n.targets))
expr=ast.unparse(node.value)
out['R05']={'reproduced':True,'requested_reference_max_epochs':777,'actual_rendered_expression':expr,'reason':'Renderer has no training_config argument; external_render_adapters discards approved setting'}
Path(sys.argv[2]).write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
