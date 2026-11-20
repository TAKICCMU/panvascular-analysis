#!/usr/bin/env python3
"""One-thread control process; unchanged reviewed dispatcher allocates its one 12-thread worker."""
import os,sys,json,argparse
from pathlib import Path
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
os.environ['CUDA_VISIBLE_DEVICES']=''
p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--execute',action='store_true',required=True);p.add_argument('--modules',nargs='+',required=True);args=p.parse_args()
if args.modules!=['cohort','magma']:raise SystemExit('Only cohort dependency and MAGMA are authorized in this launcher')
import yaml
cf=Path(args.config).resolve();c=yaml.safe_load(cf.read_text());os.environ['DOWNSTREAM_CONFIG_PATH']=str(cf);C=Path(__file__).resolve().parents[1]/'magma_adapter_v3';sys.path.insert(0,str(C));from src.dispatcher import run
rd,status=run(c,cf,args.modules);print(json.dumps({'run_dir':str(rd),'status':status},indent=2));sys.exit(0 if all(x['status']=='COMPLETE' for x in status.values()) else 3)
