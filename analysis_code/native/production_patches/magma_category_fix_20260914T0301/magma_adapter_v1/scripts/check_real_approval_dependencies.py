from pathlib import Path
import sys,json,time
sys.dont_write_bytecode=True
B=Path(__file__).resolve().parents[1];sys.path.insert(0,str(B))
from src.review_dependencies import reviewed_dependencies
from src.core import atomic_json
p=B/'configs/phase1_config_PROPOSED.yaml';c=json.loads(p.read_text());t=time.time()
try:
 d,m=reviewed_dependencies(c,p);atomic_json(B/'evidence/REVIEW_DEPENDENCIES_PROPOSED.json',m);r={'status':'PASS_CONTENT_BINDING_ONLY','digest':d,'files':len(m['files']),'seconds':time.time()-t,'production_approval_receipt':False}
except Exception as e:
 import traceback
 r={'status':'BLOCKED','error':str(e),'traceback':traceback.format_exc(),'seconds':time.time()-t};print(r['traceback'],flush=True)
atomic_json(B/'evidence/REVIEW_DEPENDENCY_CHECK.json',r);print(json.dumps(r),flush=True)
if r['status']=='BLOCKED':sys.exit(2)
