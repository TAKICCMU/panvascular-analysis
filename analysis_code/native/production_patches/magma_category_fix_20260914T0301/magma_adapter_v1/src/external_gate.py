"""Single external authorization bridge to the package's reviewed release contract.
Only standard-library imports occur before configure_cpu; no scientific code in gates.
"""
import argparse,json,os,sys,tempfile
from pathlib import Path
THREAD_KEYS=['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS','VECLIB_MAXIMUM_THREADS']

def load_config(path):
    import yaml
    p=Path(path).resolve();c=yaml.safe_load(p.read_text())
    n=int(c['execution']['per_worker_threads'])
    if n<1 or n>20:raise RuntimeError('External worker CPU thread cap is20')
    for key in THREAD_KEYS:os.environ[key]=str(n)
    os.environ['CUDA_VISIBLE_DEVICES']='';os.environ['DOWNSTREAM_CONFIG_PATH']=str(p)
    return c,p

def bootstrap_cli():
    p=argparse.ArgumentParser(add_help=False);p.add_argument('--config',required=True)
    a,_=p.parse_known_args();return load_config(a.config)

def imports():
    root=Path(__file__).resolve().parents[1]
    if str(root) not in sys.path:sys.path.insert(0,str(root))
    from src.dispatcher import authorization,code_digest
    from src.core import preflight,sha,Blocked,atomic_json
    return root,authorization,code_digest,preflight,sha,Blocked,atomic_json

def authorize(config_path):
    c,p=load_config(config_path);_,authorization,*_=imports();authorization(c,p);return c,p

def approved_preflight(config_path,receipt_dir):
    c,p=load_config(config_path)
    root,auth,digest,preflight,sha,Blocked,atomic_json=imports()
    auth(c,p)  # exact final manifest, atlas/config/code SHA, receipt, CPU budget
    output=Path(receipt_dir)
    if output.exists():raise Blocked('Non-overwrite: external run directory exists')
    output.mkdir(parents=True,exist_ok=False)
    obs,ontology,info=preflight(c,output/'preflight.json')
    record={'config_sha256':sha(p),'code_sha256':digest(root),'atlas_sha256':c['atlas']['input_sha256'],'preflight_sha256':sha(output/'preflight.json')}
    atomic_json(output/'EXTERNAL_AUTHENTICATED.json',record)
    os.environ['DOWNSTREAM_AUTH_CONFIG']=str(p)
    os.environ['DOWNSTREAM_AUTH_RECEIPT']=str(output/'EXTERNAL_AUTHENTICATED.json')
    os.environ['DOWNSTREAM_AUTH_PYTHON']=sys.executable
    os.environ['DOWNSTREAM_AUTH_GATE']=str(root/'scripts/external_gate_cli.py')
    os.environ['DOWNSTREAM_AUTH_THREADS']=str(c['execution']['per_worker_threads'])
    return c,obs,ontology,info

def check_receipt(config_path,receipt_path):
    c,p=load_config(config_path);root,auth,digest,preflight,sha,Blocked,_=imports();auth(c,p)
    rp=Path(receipt_path);r=json.loads(rp.read_text());pp=rp.parent/'preflight.json'
    if r!={'config_sha256':sha(p),'code_sha256':digest(root),'atlas_sha256':c['atlas']['input_sha256'],'preflight_sha256':sha(pp)}:
        raise Blocked('External preflight receipt changed or belongs to other code/config')
    q=json.loads(pp.read_text())
    if q.get('status')!='PASS_RELEASE_PREFLIGHT' or q.get('sha256')!=c['atlas']['input_sha256']:raise Blocked('Invalid release preflight')
    return c

def reviewed_cohort(config,obs):
    import pandas as pd
    from src.cohort import cohort_from_obs,apply_unit_review
    review=pd.read_csv(config['statistical_units']['cohort_design_file'],sep='\t',keep_default_na=False)
    return apply_unit_review(cohort_from_obs(obs),review)
