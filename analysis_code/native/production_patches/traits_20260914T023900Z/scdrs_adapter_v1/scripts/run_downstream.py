#!/usr/bin/env python3
"""Default is status only. --execute requires a later user approval receipt."""
import argparse,os,sys,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--execute',action='store_true');p.add_argument('--worker');p.add_argument('--run-dir');p.add_argument('--modules',nargs='+');args=p.parse_args()
import yaml
cfgpath=Path(args.config).resolve();c=yaml.safe_load(cfgpath.read_text())
# Set limits before scientific imports.
threads=str(c['execution']['per_worker_threads'])
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]=threads
os.environ['CUDA_VISIBLE_DEVICES']='';os.environ['DOWNSTREAM_CONFIG_PATH']=str(cfgpath)
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.core import Blocked,atomic_json
from src.dispatcher import authorization,run,worker
try:
    if args.worker:
        authorization(c,cfgpath)
        if not args.run_dir:raise Blocked('Missing run directory')
        result=worker(c,args.worker,args.run_dir);atomic_json(Path(args.run_dir)/(args.worker+'.receipt.json'),{'status':result.get('status','FAILED') if isinstance(result,dict) else 'FAILED','result':result});print(json.dumps(result,default=str))
    elif args.execute:
        rd,status=run(c,cfgpath,args.modules);print(json.dumps({'run_dir':str(rd),'status':status},indent=2))
        if any(v['status']!='COMPLETE' for v in status.values()):sys.exit(3)
    else:
        a=c['atlas'];exists=Path(a['path']).exists() and Path(a['manifest']).exists()
        m=json.load(open(a['manifest'])) if exists else {}
        ready=m.get('release_status')=='FINAL_APPROVED' and m.get('level2_frozen') is True and m.get('ready_for_downstream') is True
        print(json.dumps({'status':'CODE_READY_FOR_REVIEW_NOT_EXECUTED' if ready else 'WAITING_FOR_FINAL_L2_ATLAS','input_exists':exists,'executed':False,'policy':'No automatic fallback to working/partial/legacy atlas'},indent=2))
except Blocked as e:
    print(json.dumps({'status':'BLOCKED','reason':str(e)},ensure_ascii=False));sys.exit(3)
