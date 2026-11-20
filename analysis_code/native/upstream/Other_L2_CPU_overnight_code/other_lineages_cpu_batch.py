# -*- coding: utf-8 -*-
"""Independent CPU-only batch for 13 non-macrophage L1 lineages.

Imported by a notebook without importing Scanpy/OmicVerse or changing any environment.
The launcher starts a separate session using the same Python interpreter. Workers never
receive the notebook's adata; only the frozen h5ad is read. No final L2 assignment.
"""
from __future__ import annotations
import argparse, copy, fcntl, hashlib, html, importlib.util, json, os
from pathlib import Path
import shutil, signal, subprocess, sys, time, traceback, zipfile
from datetime import datetime

ROOT=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908')
SOURCE=ROOT/'panvascular_CPU_theta2_level1_frozen_v1.h5ad'
PATCH_DIR=ROOT/'annotation/HVG_clip_review_20260910'
PATCH_SHA256='41db6ac24969e95d6241770ddb88a4530fe6be662e0d3f95bed47b1ba195254b'
FORBIDDEN={'Macrophage','Schwann cell','Unknown'}
ORDER=['Fibroblast','Smooth muscle cell','Endothelial cell','T cell','Fibromyocyte',
       'Monocyte','Natural killer cell','Dendritic cell','B cell','Plasma cell',
       'Mast cell','Neutrophil','Pericyte']
VERSION='Other_L2_CPU_Round1_1.0'


def require(ok,msg):
    if not bool(ok):raise RuntimeError(msg)


def stamp():return datetime.now().strftime('%Y%m%d_%H%M%S_%f')


def sha_file(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(16*1024*1024),b''):h.update(b)
    return h.hexdigest()


def read_json(path):return json.loads(Path(path).read_text(encoding='utf-8'))


def write_json(path,obj):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    temp=p.with_name(p.name+f'.{os.getpid()}.tmp')
    temp.write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=str),encoding='utf-8');temp.replace(p)


def stat_key(path):
    s=Path(path).stat()
    return {'size':s.st_size,'mtime_ns':s.st_mtime_ns,'inode':s.st_ino,'device':s.st_dev}


def cpu_env(threads):
    env=os.environ.copy()
    env.update(CUDA_VISIBLE_DEVICES='',MPLBACKEND='Agg',PYTHONHASHSEED='0',PYTHONUNBUFFERED='1')
    for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'):
        env[k]=str(threads)
    return env


def load_core(code_dir):
    p=Path(code_dir)/'lineage_cpu_core.py'
    spec=importlib.util.spec_from_file_location('_other_lineage_core',p)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def _live_pid(info):
    try:
        import psutil
        p=psutil.Process(int(info['pid']))
        return abs(p.create_time()-float(info['create_time']))<1 and p.status()!=psutil.STATUS_ZOMBIE
    except Exception:return False


def _write_active(run,pid):
    import psutil
    write_json(Path(run)/'scheduler_process.json',{'pid':pid,'create_time':psutil.Process(pid).create_time(),
                                                  'run':str(run),'time':datetime.now().isoformat()})


def launch_other_lineages(source=SOURCE, max_workers=2, total_threads=40,
                           resume_run=None, macrophage_run=None,
                           min_available_memory_gb=48, force_new=False, wait_for_start=True, startup_timeout=180):
    """Start/reattach detached server batch; returns run directory and PID immediately.

    A closed browser is OK while server processes stay alive. Server shutdown/job limits
    still stop work. Recalling with resume_run reuses verified checkpoints.
    """
    source=Path(source).expanduser().resolve();manifest=source.with_suffix('.manifest.json')
    require(source.is_file() and manifest.is_file(),f'冻结h5ad或manifest缺失：{source}')
    require(1<=int(total_threads)<=40 and 1<=int(max_workers)<=4,'线程上限40，分支并发数1–4。')
    require(total_threads>=max_workers,'总线程数小于并发数。')
    code=Path(__file__).resolve().parent
    for f in ('lineage_cpu_core.py','lineage_plans.json','run_defaults.json','sample_metadata_export.py'):
        require((code/f).is_file(),f'请完整解压代码包，缺少：{code/f}')
    # Import only small utility dependencies here, never reload/modify notebook data.
    import psutil
    home=source.parent/'L2_CPU_others';home.mkdir(parents=True,exist_ok=True)
    pointer=home/'ACTIVE_BATCH.json'
    with (home/'launch.lock').open('a+') as lock:
        fcntl.flock(lock.fileno(),fcntl.LOCK_EX)
        if resume_run is None and not force_new and pointer.exists():
            info=read_json(pointer);p=Path(info['run'])
            if (p/'batch_task.json').is_file():
                prior=read_json(p/'batch_task.json')
                if prior['source']==str(source) and prior['code_version']==VERSION:
                    resume_run=p
        if resume_run is None:
            run=home/('Other_lineages_round1_'+stamp());run.mkdir()
            shutil.copytree(code,run/'code',ignore=shutil.ignore_patterns('__pycache__','*.pyc','*.ipynb','TEST*'))
            plans=read_json(run/'code/lineage_plans.json')
            require(set(plans)==set(ORDER) and not FORBIDDEN.intersection(plans),'谱系列表异常。')
            patch=PATCH_DIR/'code/patched_function.py';tests=PATCH_DIR/'audit/TESTS.json'
            require(patch.is_file() and tests.is_file(),'缺少已验证HVG补丁或测试记录。')
            require(sha_file(patch)==PATCH_SHA256 and read_json(tests).get('passed') is True,
                    'HVG补丁或既往测试状态不符；不自动更换算法。')
            (run/'audit').mkdir();shutil.copy2(patch,run/'audit/patched_function.py');shutil.copy2(tests,run/'audit/HVG_TESTS.json')
            protected=[]
            if macrophage_run is not None:
                mp=Path(macrophage_run).expanduser().resolve()
                # Only read fingerprints of selected pre-existing files; no macrophage writes.
                for rel in ('task.json','RESULTS.json','tables/cell_cluster_mapping.tsv.gz',
                            'tables/all_resolution_review_candidates.tsv'):
                    p=mp/rel
                    if p.is_file():protected.append({'path':str(p),'sha256':sha_file(p)})
            task={'code_version':VERSION,'created_at':datetime.now().isoformat(),'source':str(source),
                  'source_manifest':str(manifest),'threads_per_worker':int(total_threads)//int(max_workers),
                  'total_threads':int(total_threads),'max_workers':int(max_workers),
                  'min_available_memory_gb':float(min_available_memory_gb),'min_free_disk_gb':12,
                  'lineages':ORDER,'n_expected':903398,'excluded':sorted(FORBIDDEN),
                  'protected_macrophage_files':protected,'python':sys.executable,
                  'core_sha256':sha_file(run/'code/lineage_cpu_core.py'),
                  'plans_sha256':sha_file(run/'code/lineage_plans.json'),
                  'defaults_sha256':sha_file(run/'code/run_defaults.json')}
            write_json(run/'batch_task.json',task)
            write_json(run/'BATCH_STATUS.json',{'status':'READY_TO_START','run':str(run)})
        else:
            run=Path(resume_run).expanduser().resolve();require((run/'batch_task.json').is_file(),'resume_run不是本批处理目录。')
            task=read_json(run/'batch_task.json')
            require(task['source']==str(source),'继续运行的源h5ad与原任务不同。')
            procfile=run/'scheduler_process.json'
            if procfile.exists() and _live_pid(read_json(procfile)):
                pid=read_json(procfile)['pid'];print(f'已在运行，不重复启动。PID={pid}\n{run}',flush=True);return run,pid
            statusfile=run/'BATCH_STATUS.json'
            if statusfile.exists() and read_json(statusfile).get('status')=='PAUSED_FOR_OTHER_L2_REVIEW':
                print(f'该批13个分支已完成，直接复用：{run}',flush=True);return run,None
        write_json(pointer,{'run':str(run),'source':str(source),'code_version':VERSION})
        log=(run/'scheduler.log').open('a',encoding='utf-8')
        try:
            proc=subprocess.Popen([sys.executable,'-u',str(run/'code/other_lineages_cpu_batch.py'),'--scheduler',str(run)],
                                  env=cpu_env(1),stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,
                                  start_new_session=True,close_fds=True,cwd=str(run))
            _write_active(run,proc.pid)
        finally:log.close()
        time.sleep(0.7)
        if proc.poll() is not None and proc.returncode!=0:
            raise RuntimeError(f'调度器启动失败，请查看：{run/"scheduler.log"}')
    print(f'其他13谱系CPU任务已提交（不是已完成）。PID={proc.pid}\n目录：{run}\n日志：{run/"scheduler.log"}',flush=True)
    if wait_for_start:
        deadline=time.time()+float(startup_timeout); last=''
        while time.time()<deadline:
            st=read_json(run/'BATCH_STATUS.json'); phase=st.get('status','STARTING')
            if phase!=last:print('启动核查：',phase,flush=True);last=phase
            if phase in ('BATCH_STOPPED_WITH_ERROR','PARTIAL_FAILURE_REVIEW_AVAILABLE'):
                raise RuntimeError(f'启动/运行发生错误，请查看{run/"scheduler.log"}\n'+st.get('error',''))
            if phase in ('RUNNING_OTHER_L2','WAITING_FOR_RESOURCES','PAUSED_FOR_OTHER_L2_REVIEW'):break
            if proc.poll() is not None and proc.returncode!=0:
                raise RuntimeError(f'调度器退出，请查看：{run/"scheduler.log"}')
            time.sleep(2)
        else:print('启动核查尚未完成；子进程仍继续，稍后运行batch_status查看。',flush=True)
    return run,proc.pid


def batch_status(run):
    """Read current status only. Does not start/restart jobs."""
    import pandas as pd
    root=Path(run)
    state=read_json(root/'BATCH_STATUS.json') if (root/'BATCH_STATUS.json').exists() else {}
    print('批处理状态：',state.get('status','UNKNOWN'))
    p=root/'batch_progress.csv'
    if p.exists():return pd.read_csv(p,keep_default_na=False)
    return pd.DataFrame([state])


def preflight(run,task):
    import numpy as np
    import h5py
    import scanpy as sc
    import anndata as ad
    import torch, numba, igraph, pynndescent, psutil
    import importlib, inspect
    try:from anndata.io import read_elem
    except ImportError:from anndata._io.specs import read_elem
    require(sc.__version__=='1.11.5','需要先前成功运行的Scanpy 1.11.5环境；不自动升级/安装。')
    require(not torch.cuda.is_available(),'CPU调度器不应看到CUDA。')
    require('mask_var' in inspect.signature(sc.tl.rank_genes_groups).parameters,'Scanpy API缺少mask_var。')
    hm=importlib.import_module('omicverse.external.harmony.harmony')
    require('device' in inspect.signature(hm.run_harmony).parameters,'Harmony未暴露device参数。')
    require(sha_file(run/'code/lineage_cpu_core.py')==task['core_sha256'],'Worker snapshot被修改。')
    require(sha_file(run/'code/lineage_plans.json')==task['plans_sha256'],'Marker计划被修改。')
    require(sha_file(run/'code/run_defaults.json')==task['defaults_sha256'],'参数被修改。')
    source=Path(task['source']);before=stat_key(source)
    print('一次性校验冻结母文件SHA256；所有worker只读复用该校验票据。',flush=True)
    digest=sha_file(source);manifest=read_json(task['source_manifest'])
    require(digest==manifest.get('sha256'),'冻结母图谱SHA256不符。')
    require(manifest.get('level1_frozen') is True and manifest.get('annotation_key')=='level1_annotation','输入未冻结。')
    core=load_core(run/'code')
    with h5py.File(source,'r') as f:
        obs=read_elem(f['obs']);var=read_elem(f['var'])
        require(obs.index.is_unique and var.index.is_unique,'cell/gene ID重复。')
        require((len(obs),len(var))==core.EXPECTED_SHAPE,'输入不是当前全量冻结图谱。')
        require(obs['level1_annotation'].notna().all() and obs['level1_annotation'].astype(str).value_counts().to_dict()==core.EXPECTED_L1,'冻结L1构成不符。')
        require(obs['sample_uid'].notna().all() and obs['sample_uid'].astype(str).nunique()==247,'sample_uid异常。')
        require(bool(read_elem(f['uns/level1_frozen'])),'h5ad冻结标记缺失。')
        shape=tuple(int(x) for x in f['X'].attrs['shape']);require(shape==core.EXPECTED_SHAPE,'raw X维度不符。')
        raw_enc=core.scalar(f['X'].attrs['encoding-type']);require(raw_enc in ('csr_matrix','csc_matrix'),'raw X必须稀疏。')
    require(stat_key(source)==before,'校验期间源文件改变。')
    free=shutil.disk_usage(run).free/1024**3
    require(free>task['min_free_disk_gb'],'工作磁盘可用空间不足12GiB。')
    estimate=3*before['size']/1024**3
    if free<estimate:print(f'注意：当前空闲{free:.1f}GiB，全部谱系raw/checkpoint/结果可能需约{estimate:.1f}GiB或更多。',flush=True)
    for item in task['protected_macrophage_files']:
        require(sha_file(item['path'])==item['sha256'],'预先记录的巨噬细胞文件发生改变。')
    ticket={'source':str(source),'source_sha256':digest,'source_stat':before,
            'manifest_sha256':sha_file(task['source_manifest']),'n_cells':len(obs),'n_genes':len(var),
            'n_lineages':13,'n_expected_other_cells':903398,'verified_at':datetime.now().isoformat()}
    write_json(run/'audit/verified_frozen_source.json',ticket)
    return ticket


def setup_tasks(run,task,ticket):
    plans=read_json(run/'code/lineage_plans.json');defaults=read_json(run/'code/run_defaults.json')
    jobs=[]
    for lineage in task['lineages']:
        require(lineage not in FORBIDDEN,'禁止重复运行巨噬/Schwann/Unknown。')
        plan=plans[lineage];folder=run/'lineages'/plan['slug']
        for sub in ('data','tables','figures','audit'):(folder/sub).mkdir(parents=True,exist_ok=True)
        for src,dst in [(run/'audit/patched_function.py',folder/'audit/patched_function.py'),
                        (run/'audit/HVG_TESTS.json',folder/'audit/HVG_TESTS.json')]:
            if not dst.exists():shutil.copy2(src,dst)
            else:require(sha_file(dst)==sha_file(src),'分支HVG补丁不符。')
        cfg=copy.deepcopy(defaults);cfg['core_panel_key']=plan['core_panel_key']
        spec={'source':task['source'],'source_manifest':task['source_manifest'],'lineage':lineage,
              'expected_n_cells':plan['expected_n_cells'],'uid_prefix':plan['uid_prefix'],
              'state_panel':plan['state'],'boundary_panel':plan['boundary'],'core_panel_key':plan['core_panel_key'],
              'config':cfg,'threads':task['threads_per_worker'],'core_sha256':task['core_sha256'],
              'verified_source_sha256':ticket['source_sha256'],'source_stat':ticket['source_stat'],
              'manifest_sha256':ticket['manifest_sha256'],'batch_verification':str(run/'audit/verified_frozen_source.json'),
              'batch_verification_sha256':sha_file(run/'audit/verified_frozen_source.json')}
        p=folder/'task.json'
        if p.exists():
            old=read_json(p)
            # Verification timestamp may change on resume, but analysis inputs/config/panels must not.
            omit={'batch_verification_sha256'}
            require({k:v for k,v in old.items() if k not in omit}=={k:v for k,v in spec.items() if k not in omit},
                    f'{lineage}断点参数/输入与历史不符。')
        write_json(p,spec)
        write_json(folder/'marker_panels.json',{'state':plan['state'],'boundary':plan['boundary'],'use':'REVIEW_ONLY_NO_AUTOMATIC_LABELS'})
        jobs.append({'lineage':lineage,'slug':plan['slug'],'n_cells':plan['expected_n_cells'],'folder':str(folder)})
    return jobs


def write_progress(run,jobs,phase):
    import pandas as pd
    rows=[]
    for j in jobs:
        p=Path(j['folder'])/'STATUS.json'
        s=read_json(p) if p.exists() else {'status':'QUEUED'}
        rows.append({'lineage':j['lineage'],'n_cells':j['n_cells'],'status':s.get('status','UNKNOWN'),
                     'error':s.get('error',''),'folder':j['folder'],'worker_log':str(Path(j['folder'])/'worker.log')})
    df=pd.DataFrame(rows);temp=run/'batch_progress.csv.tmp';df.to_csv(temp,index=False,encoding='utf-8-sig');temp.replace(run/'batch_progress.csv')
    write_json(run/'BATCH_STATUS.json',{'status':phase,'updated_at':datetime.now().isoformat(),'run':str(run),
               'completed':int(df.status.eq('PAUSED_FOR_LINEAGE_L2_ROUND1_REVIEW').sum()),'n_requested':len(jobs),
               'failed':int(df.status.str.startswith('FAILED').sum()),'L1_modified':False,'Macrophage_modified':False,
               'level2_annotation_applied':False})
    return df


def completed_and_valid(folder,task):
    p=Path(folder)/'RESULTS.json';status=Path(folder)/'STATUS.json'
    if not p.exists() or not status.exists():return False
    r=read_json(p);s=read_json(status)
    if s.get('status')!='PAUSED_FOR_LINEAGE_L2_ROUND1_REVIEW':return False
    require(r.get('source_sha256')==read_json(task['source_manifest'])['sha256'],'已完成分支源文件不符。')
    require(sha_file(r['analysis_file'])==r['analysis_sha256'],'已完成分支h5ad有变动。')
    return True


def finalize_batch(run,task,jobs):
    import pandas as pd
    tables=[];sumrows=[];htmlrows=[];packfiles=[]
    for j in jobs:
        folder=Path(j['folder']);result=folder/'RESULTS.json';s=folder/'STATUS.json'
        state=read_json(s) if s.exists() else {'status':'NOT_STARTED'}
        good=result.exists() and state.get('status')=='PAUSED_FOR_LINEAGE_L2_ROUND1_REVIEW'
        rel=str(folder.relative_to(run))
        if good:
            r=read_json(result);t=pd.read_csv(folder/'tables/all_resolution_review_candidates.tsv',sep='\t')
            t.insert(0,'level1_annotation',j['lineage']);tables.append(t)
            sumrows.append({'lineage':j['lineage'],'n_cells':j['n_cells'],'status':state['status'],
                            'partition_counts':json.dumps(r['partition_counts']),'analysis_file':r['analysis_file']})
            htmlrows.append(f"<li><a href='{rel}/LINEAGE_L2_ROUND1_OFFLINE.html'>{html.escape(j['lineage'])} ({j['n_cells']:,})</a></li>")
            packfiles.append(folder/'lineage_L2_round1_review_materials.zip')
        else:
            sumrows.append({'lineage':j['lineage'],'n_cells':j['n_cells'],'status':state.get('status'),'partition_counts':'','analysis_file':''})
            htmlrows.append(f"<li>{html.escape(j['lineage'])}: {html.escape(state.get('status','NOT_STARTED'))}</li>")
    if tables:pd.concat(tables,ignore_index=True).to_csv(run/'all_other_L2_cluster_candidates.tsv',sep='\t',index=False)
    pd.DataFrame(sumrows).to_csv(run/'all_other_L2_branch_summary.csv',index=False,encoding='utf-8-sig')
    doc="<!doctype html><meta charset='utf-8'><h1>其他13谱系L2首轮审核</h1><p>无正式L2注释；巨噬细胞未重复；Schwann与L1 Unknown未纳入。</p><ul>"+''.join(htmlrows)+"</ul><p>每个分支HTML图片内嵌；请从服务器打开本入口或完整保留相对目录。单独上传每个分支审核zip亦可。</p>"
    (run/'OTHER_L2_REVIEW_INDEX.html').write_text(doc,encoding='utf-8')
    # Per-lineage ZIPs are already compressed. Outer ZIP uses STORE to avoid expensive recompression.
    outer=run/'Other_L2_all_review_materials.zip';temp=outer.with_suffix('.partial.zip')
    with zipfile.ZipFile(temp,'w',compression=zipfile.ZIP_STORED,allowZip64=True) as z:
        for p in packfiles:z.write(p,Path('lineage_review_zips')/(p.parent.name+'__'+p.name))
        for name in ('all_other_L2_cluster_candidates.tsv','all_other_L2_branch_summary.csv','batch_task.json','batch_progress.csv'):
            if (run/name).exists():z.write(run/name,name)
        for code_file in (run/'code').glob('*'):
            if code_file.is_file() and code_file.suffix in ('.py','.json'):
                z.write(code_file,Path('code')/code_file.name)
        z.writestr('READ_ME.txt','每个lineage_review_zips内ZIP为自带图片的独立图集。解压所需谱系ZIP并打开LINEAGE_L2_ROUND1_OFFLINE.html。表达h5ad不包含在审核包中。')
    with zipfile.ZipFile(temp) as z:require(z.testzip() is None,'统一ZIP校验失败。')
    temp.replace(outer)
    return all(x['status']=='PAUSED_FOR_LINEAGE_L2_ROUND1_REVIEW' for x in sumrows)


def scheduler(run):
    import psutil
    run=Path(run).resolve();task=read_json(run/'batch_task.json');children={}
    with (run/'RUN.lock').open('a+') as lock:
        try:fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise RuntimeError('同一批次已经有调度器运行，拒绝重复。')
        _write_active(run,os.getpid())
        write_json(run/'BATCH_STATUS.json',{'status':'PREFLIGHT_RUNNING','run':str(run)})
        # Do not start another worker on top of an orphan from a prior scheduler.
        for p in (run/'lineages').glob('*/worker_process.json') if (run/'lineages').exists() else []:
            if _live_pid(read_json(p)):raise RuntimeError(f'上次子进程仍在运行：{p}。等待结束后再resume。')
        ticket=preflight(run,task);jobs=setup_tasks(run,task,ticket)
        queue=[]
        for j in jobs:
            if completed_and_valid(j['folder'],task):print('已完成，跳过：',j['lineage'],flush=True)
            else:queue.append(j)
        write_progress(run,jobs,'RUNNING_OTHER_L2')
        def stop(signum,frame):raise KeyboardInterrupt(f'Scheduler received signal {signum}')
        oldterm=signal.signal(signal.SIGTERM,stop)
        try:
            while queue or children:
                # Collect completed children. One failed lineage does not cancel the remaining ones.
                for pid,item in list(children.items()):
                    p,j,log=item
                    if p.poll() is None:continue
                    log.close();code=p.returncode
                    if code!=0:
                        st=Path(j['folder'])/'STATUS.json';previous=read_json(st) if st.exists() else {}
                        if not str(previous.get('status','')).startswith('FAILED'):
                            write_json(st,{'status':'FAILED_NOT_FROZEN','error':f'worker exit code {code}',
                                           'previous_stage':previous.get('status')})
                        print('分支失败，保留断点并继续其他谱系：',j['lineage'],flush=True)
                    else:print('分支完成：',j['lineage'],flush=True)
                    del children[pid]
                while queue and len(children)<task['max_workers']:
                    vm=psutil.virtual_memory().available/1024**3
                    disk=shutil.disk_usage(run).free/1024**3
                    if vm<task['min_available_memory_gb'] or disk<task['min_free_disk_gb']:
                        print(f'暂缓启动下一分支：可用内存{vm:.1f}GiB；磁盘{disk:.1f}GiB',flush=True);break
                    j=queue.pop(0);folder=Path(j['folder']);log=(folder/'worker.log').open('a',encoding='utf-8')
                    p=subprocess.Popen([task['python'],'-u',str(run/'code/lineage_cpu_core.py'),'--worker',str(folder)],
                        env=cpu_env(task['threads_per_worker']),stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,
                        close_fds=True,cwd=str(folder))
                    children[p.pid]=(p,j,log)
                    write_json(folder/'worker_process.json',{'pid':p.pid,'create_time':psutil.Process(p.pid).create_time(),'lineage':j['lineage']})
                    print(f"启动 {j['lineage']} {j['n_cells']:,} cells | PID={p.pid} | threads={task['threads_per_worker']}",flush=True)
                write_progress(run,jobs,'RUNNING_OTHER_L2' if children else 'WAITING_FOR_RESOURCES')
                if queue or children:time.sleep(15)
        finally:
            signal.signal(signal.SIGTERM,oldterm)
            for p,j,log in children.values():
                if p.poll() is None:
                    p.terminate()
                    try:p.wait(timeout=20)
                    except subprocess.TimeoutExpired:p.kill();p.wait()
                log.close()
        print('所有已排队分支均已退出，复核冻结母文件……',flush=True)
        require(stat_key(task['source'])==ticket['source_stat'] and sha_file(task['source'])==ticket['source_sha256'],
                '最终冻结母文件校验不符，结果不可当成同源完结版。')
        for x in task['protected_macrophage_files']:require(sha_file(x['path'])==x['sha256'],'记录的巨噬细胞文件改变。')
        write_progress(run,jobs,'PACKAGING_OTHER_L2_REVIEW')
        all_ok=finalize_batch(run,task,jobs)
        final_state='PAUSED_FOR_OTHER_L2_REVIEW' if all_ok else 'PARTIAL_FAILURE_REVIEW_AVAILABLE'
        write_progress(run,jobs,final_state)
        write_json(run/'audit/final_source_verification.json',{'sha256':ticket['source_sha256'],'unchanged':True,
                        'macrophage_recorded_files_unchanged':True,'completed_at':datetime.now().isoformat()})
        (run/'PROJECT_STATE_OTHER_L2_ROUND1.md').write_text('# Other L2 CPU first-round branch\n\n'+
            json.dumps({'status':final_state,'source':task['source'],'L1_unchanged':True,'Macrophage_unchanged':True,
                        'lineages':task['lineages'],'L2_frozen':False,'annotation_applied':False,
                        'review_index':str(run/'OTHER_L2_REVIEW_INDEX.html')},ensure_ascii=False,indent=2),encoding='utf-8')
        print(final_state, '\n审核包：',run/'Other_L2_all_review_materials.zip',flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--scheduler',required=True);args=ap.parse_args()
    try:scheduler(args.scheduler)
    except BaseException as exc:
        write_json(Path(args.scheduler)/'BATCH_STATUS.json',{'status':'BATCH_STOPPED_WITH_ERROR','error':repr(exc),
                    'traceback':traceback.format_exc(),'time':datetime.now().isoformat(),
                    'L1_modified':False,'Macrophage_modified':False})
        traceback.print_exc();sys.exit(1)
