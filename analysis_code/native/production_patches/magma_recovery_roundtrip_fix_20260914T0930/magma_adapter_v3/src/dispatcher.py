"""Dependency-aware CPU dispatcher; full execution requires a separately reviewed receipt."""
from pathlib import Path
import json,os,subprocess,time,hashlib,sys
from concurrent.futures import ThreadPoolExecutor,wait,FIRST_COMPLETED
from .core import Blocked,sha,atomic_json
from .review_dependencies import reviewed_dependencies
from .status import receipt_state
DAG={'cohort':[],'counts':['cohort'],'markers':['cohort'],'pseudobulk':['cohort'],'de':['pseudobulk'],'meta':['de'],'covarnet':['counts'],'composition':['counts'],'figures':['counts'],'dotplot':['markers','figures'],'external_audit':['cohort'],'magma':['cohort'],'scdrs':['cohort'],'communication':['cohort'],'targets':['meta','magma'],'sensitivity':['de','composition']}

def code_digest(root):
    files=sorted(p for folder in ['src','scripts','external_legacy'] for p in (Path(root)/folder).glob('**/*') if p.is_file() and p.suffix in ['.py','.R','.sh'])
    return hashlib.sha256('\n'.join(str(p.relative_to(root))+':'+sha(p) for p in files).encode()).hexdigest()

def authorization(config,config_path):
    e=config['execution'];a=config['atlas']
    if not Path(a['path']).exists() or not Path(a['manifest']).exists():raise Blocked('WAITING_FOR_FINAL_L2_ATLAS')
    m=json.load(open(a['manifest']))
    if m.get('release_status')!='FINAL_APPROVED' or not m.get('level2_frozen') or not m.get('ready_for_downstream'):raise Blocked('WAITING_FOR_FINAL_L2_ATLAS')
    if e.get('mode')!='approved_execution' or e.get('run_full_analysis') is not True:raise Blocked('CODE_READY_FOR_REVIEW_NOT_EXECUTED')
    receipt=e.get('approval_receipt')
    if not receipt or not Path(receipt).exists():raise Blocked('Missing subsequent user approval receipt')
    r=json.load(open(receipt))
    if r.get('approved') is not True or r.get('atlas_sha256')!=a.get('input_sha256') or r.get('config_sha256')!=sha(config_path):raise Blocked('Approval receipt does not match exact reviewed atlas/config')
    dependency_digest,dependency_manifest=reviewed_dependencies(config,config_path)
    if r.get('reviewed_dependency_digest')!=dependency_digest:raise Blocked('Approval receipt reviewed dependency digest differs')
    if r.get('code_sha256')!=code_digest(Path(__file__).resolve().parents[1]):raise Blocked('Approval receipt code hash differs')
    if e.get('cpu_only') is not True or int(e['max_workers'])*int(e['per_worker_threads'])>40 or int(e['max_workers'])>2:raise Blocked('CPU resource budget exceeds reviewed 2x20/40')
    cv=config.get('covarnet',{})
    if int(cv.get('workers',1))*int(cv.get('threads_per_worker',20))>int(e['per_worker_threads']):raise Blocked('CoVarNet nested parallelism exceeds module thread allocation')
    return r

def cache_identity(config,config_path):
    paths=[config_path,config['atlas']['manifest'],config['atlas']['schema_contract']]
    for name in ['cohort_design_file','contrasts_file']:
        p=config['statistical_units'].get(name)
        if p and Path(p).exists():paths.append(p)
    for p in Path(config_path).parent.glob('*'):
        if p.is_file():paths.append(str(p))
    payload={'atlas_sha256':config['atlas']['input_sha256'],'code':code_digest(Path(__file__).resolve().parents[1]),'files':{str(p):sha(p) for p in sorted(set(map(str,paths)))},'external_resources':config.get('external',{})}
    return hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest(),payload

def run(config,config_path,selected=None):
    authorization(config,config_path)
    dag={k:list(v) for k,v in DAG.items()}
    if 'sensitivity' in dag:
        requested=config.get('sensitivity',{}).get('model_modules',['de','composition'])
        if not set(requested)<= {'de','composition','covarnet'}:raise Blocked('Unsupported sensitivity model dependency')
        dag['sensitivity']=['counts']+list(requested)
    from .core import preflight
    import datetime
    key,payload=cache_identity(config,config_path)
    parent=Path(config['outputs']['root'])/'runs';parent.mkdir(parents=True,exist_ok=True)
    run=parent/(datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%f')+'_'+key[:12]);run.mkdir(exist_ok=False)
    atomic_json(run/'RUN_IDENTITY.json',payload)
    obs,ont,info=preflight(config,run/'preflight.json');del obs,ont
    modules=selected or config['execution']['enabled_modules']
    unknown=set(modules)-set(dag)
    if unknown:raise Blocked(f'Unknown modules: {unknown}')
    wanted=set(modules)
    def add_deps(m):
        for d in dag[m]:wanted.add(d);add_deps(d)
    for m in modules:add_deps(m)
    env=os.environ.copy();threads=str(config['execution']['per_worker_threads'])
    for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:env[k]=threads
    env['CUDA_VISIBLE_DEVICES']='';env['DOWNSTREAM_RUN_CONFIG_SHA256']=sha(config_path)
    status={m:{'status':'PENDING'} for m in sorted(wanted)}
    def launch(m):
        log=run/(m+'.log');start=time.time();cmd=[sys.executable,str(Path(__file__).resolve().parents[1]/'scripts/run_downstream.py'),'--config',str(config_path),'--worker',m,'--run-dir',str(run)]
        with log.open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env,check=False)
        receipt=run/(m+'.receipt.json')
        state=receipt_state(r.returncode,json.loads(receipt.read_text()) if receipt.exists() else None)
        return {'status':state,'exit_code':r.returncode,'seconds':time.time()-start,'command':cmd,'log':str(log)}
    with ThreadPoolExecutor(max_workers=config['execution']['max_workers']) as pool:
        futures={}
        while any(v['status'] in ['PENDING','RUNNING'] for v in status.values()):
            for m in status:
                if status[m]['status']!='PENDING':continue
                deps=[status[d]['status'] for d in dag[m]]
                if any(d in ['FAILED','BLOCKED','BLOCKED_DEPENDENCY','PARTIAL','PREPARED_NOT_EXECUTED'] for d in deps):status[m]={'status':'BLOCKED_DEPENDENCY','dependencies':dag[m]}
                elif all(d=='COMPLETE' for d in deps) and len(futures)<config['execution']['max_workers']:
                    status[m]={'status':'RUNNING'};futures[pool.submit(launch,m)]=m
            atomic_json(run/'RUN_STATUS.json',status)
            if futures:
                done,_=wait(futures,timeout=10,return_when=FIRST_COMPLETED)
                for f in done:
                    m=futures.pop(f)
                    try:status[m]=f.result()
                    except Exception as exc:status[m]={'status':'FAILED','reason':str(exc)}
            elif any(v['status']=='PENDING' for v in status.values()):raise Blocked('DAG stalled')
    atomic_json(run/'RUN_STATUS.json',status)
    return run,status

def worker(config,module,run):
    import pandas as pd
    from .core import read_frame
    run=Path(run)
    # An actual run receipt is necessary even for an individually invoked worker.
    if not (run/'preflight.json').exists() or not (run/'RUN_IDENTITY.json').exists():raise Blocked('Worker requires authenticated release/run receipts')
    key,payload=cache_identity(config,os.environ.get('DOWNSTREAM_CONFIG_PATH',''))
    if payload!=json.load(open(run/'RUN_IDENTITY.json')):raise Blocked('Worker run identity differs')
    a=config['atlas'];cols=[a['sample_key'],a['level1_key'],a['level2_key']]
    if module=='cohort':
        from .cohort import cohort_from_obs,apply_unit_review,read_unit_review,FIELDS
        obs=read_frame(a['path'],columns=FIELDS+cols[1:]);m=cohort_from_obs(obs)
        review=read_unit_review(config['statistical_units']['cohort_design_file']);m=apply_unit_review(m,review)
        (run/'cohort').mkdir();m.to_csv(run/'cohort/cohort_design.tsv',sep='\t',index=False);return {'status':'COMPLETE'}
    if module=='counts':
        from .counts import state_counts,biological_counts,coverage_table
        from .cohort import read_unit_review
        obs=read_frame(a['path'],columns=cols);c=read_unit_review(run/'cohort/cohort_design.tsv');out=run/'counts';out.mkdir()
        state_counts(obs).to_csv(out/'sample_state_counts.tsv',sep='\t',index=False)
        biological_counts(obs,c).to_csv(out/'biological_state_counts.tsv',sep='\t',index=False)
        coverage_table(obs,c).to_csv(out/'state_coverage.tsv',sep='\t',index=False)
        ont=pd.read_csv(Path(config['outputs']['release_support'])/'level2_ontology.tsv',sep='\t');special=obs.loc[obs[a['level2_key']].isin(['Unknown','Not_assessed']),[a['level1_key'],a['level2_key']]].drop_duplicates();special['approved']=True;special['covarnet_eligible']=False;ont=pd.concat([ont,special],ignore_index=True).drop_duplicates([a['level1_key'],a['level2_key']]);ont.loc[ont[a['level2_key']].isin(['Unknown','Not_assessed']),'covarnet_eligible']=False;ont.to_csv(out/'ontology.tsv',sep='\t',index=False);return {'status':'COMPLETE'}
    if module in ['markers','pseudobulk']:
        from .expression import marker_summary,export_pseudobulk
        obs=read_frame(a['path'],columns=cols)
        return marker_summary(config,obs,run/'markers') if module=='markers' else export_pseudobulk(config,obs,pd.read_csv(run/'cohort/cohort_design.tsv',sep='\t'),run/'pseudobulk')
    if module=='de':
        from .inference import execute_de
        return execute_de(config,run)
    if module=='meta':
        from .inference import meta_input
        return meta_input(config,run)
    if module=='covarnet':
        from .covarnet import execute
        return execute(config,run)
    if module=='composition':
        from .composition import execute
        return execute(config,run)
    if module=='figures':
        from .figures import export_figures
        return export_figures(config,read_frame(a['path'],columns=cols+['study','tissue','comparison_group']),run)
    if module=='dotplot':
        from .figures import dotplot
        return {'status':'COMPLETE','panels':{level:dotplot(config,run,level) for level in ['L1','L2']}}
    if module=='sensitivity':
        from .sensitivity import execute
        return execute(config,run)
    if module=='external_audit':
        from .external_modules import execute
        result=execute(config,run)
        return {'status':'PREPARED_NOT_EXECUTED','phase':'INPUTS_AUDITED_NOT_MODELS','result':result}
    if module in ['magma','scdrs','communication','targets']:
        cmd=[sys.executable,str(Path(__file__).resolve().parents[1]/'scripts/external_run_module.py'),'--config',os.environ['DOWNSTREAM_CONFIG_PATH'],'--module',module,'--out',str(run/module)]
        result=subprocess.run(cmd,check=False)
        if result.returncode:raise Blocked(f'External {module} failed/blocked; exit={result.returncode}; see module log')
        return {'status':'COMPLETE','command':cmd}
    raise Blocked(f'No module {module}')
