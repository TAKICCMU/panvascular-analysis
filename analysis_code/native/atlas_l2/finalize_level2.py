# -*- coding: utf-8 -*-
"""Publish the approved L2 annotation snapshot; no expression analysis.

Read full current h5ad from an EXPLICIT path. Verify exact reviewed per-cell
membership and previous labels; apply a fixed 75-row decision table and one
approved rename. Never classify by scores, infer clinical data, or alter L1.
Run from the existing notebook via run_from_notebook(); work happens in a CPU
child so existing notebook objects are neither cleared nor silently refreshed.
"""
from __future__ import annotations
import argparse
import contextlib
from datetime import datetime
import gc
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import traceback

BUNDLE = Path(__file__).resolve().parent
ROOT = Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908')
FINAL_NAME = 'panvascular_CPU_theta2_L1L2_final_v1.h5ad'
L1, L2 = 'level1_annotation', 'level2_annotation'
VERSION = 'L1L2_final_v1'
SCHEMA = '1.0'
EXPECTED_SHAPE = (1175115, 37336)
EXPECTED_L1 = {'Endothelial cell':104415,'Fibroblast':251371,'Fibromyocyte':29347,
 'Smooth muscle cell':243346,'Pericyte':4513,'Macrophage':210967,'Monocyte':31046,
 'Neutrophil':18509,'T cell':156349,'Natural killer cell':25719,'B cell':10496,
 'Plasma cell':4802,'Mast cell':12929,'Dendritic cell':10556,'Schwann cell':2832,'Unknown':57918}
META_KEYS = ['sample_uid','analysis_sample_uid','study','study_uid','atlas_source_id','sample_id',
 'sample_name','source_series_accession','source_sample_accession','sra_experiment_accession',
 'biosample_accession','biological_sample_uid','donor_uid','donor_id','age','sex','tissue',
 'tissue_detail','vascular_bed','vascular_bed_group','sample_compartment','condition','disease',
 'comparison_group','disease_context','cohort_disease_context','source_diagnosis',
 'donor_disease_context','tissue_disease_status','assay_type']
SPECIAL = {'Unknown', 'Not_assessed'}
FORBIDDEN = {'Pending','PENDING','Unresolved','NA','NR','nan','None',''}


def require(ok, message):
    if not bool(ok):
        raise RuntimeError(message)


def stamp():
    return datetime.now().strftime('%Y%m%d_%H%M%S_%f')


def load_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def dump_json(path, obj):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=2, default=str), encoding='utf-8')
    os.replace(tmp, path)


def sha_file(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(16 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def sha_index(values):
    h = hashlib.sha256()
    for value in values:
        b = str(value).encode('utf-8')
        h.update(len(b).to_bytes(8, 'little')); h.update(b)
    return h.hexdigest()


def read_table(path):
    import pandas as pd
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def verify_bundle():
    for name, digest in load_json(BUNDLE/'CHECKSUMS.json').items():
        require(sha_file(BUNDLE/name) == digest, f'代码/决定表已变化，请不要绕过校验：{name}')


def io_funcs():
    try:
        from anndata.io import read_elem, write_elem
    except ImportError:
        from anndata._io.specs import read_elem, write_elem
    return read_elem, write_elem


def _json_value(x):
    import numpy as np
    if isinstance(x, bytes): return x.decode('utf-8')
    if isinstance(x, np.ndarray):
        return _json_value(x.item()) if x.ndim == 0 else [_json_value(v) for v in x.tolist()]
    if isinstance(x, np.generic): return x.item()
    return x


def hash_h5_tree(node):
    """Same logical dataset hashing as the preceding verified working update.
    Streams large datasets; rejects external links and virtual/external storage.
    """
    import h5py
    import numpy as np
    h = hashlib.sha256()
    def walk(n, relative):
        h.update(relative.encode())
        h.update(json.dumps({k:_json_value(v) for k,v in n.attrs.items()}, sort_keys=True, default=str).encode())
        if isinstance(n, h5py.Dataset):
            require(not n.is_virtual and not n.external, '表达必须内嵌，禁止virtual/external dataset。')
            h.update(str(n.dtype).encode()); h.update(str(n.shape).encode())
            if n.shape == (): chunks = [n[()]]
            else:
                step = max(1, int(16*1024**2 / max(1, n.dtype.itemsize * int(np.prod(n.shape[1:])))))
                chunks = (n[i:i+step] for i in range(0, len(n), step))
            for x in chunks:
                a = np.asarray(x)
                if a.dtype.kind in 'OUS':
                    for v in a.ravel():
                        b = str(_json_value(v)).encode('utf-8')
                        h.update(len(b).to_bytes(8, 'little')); h.update(b)
                else: h.update(memoryview(np.ascontiguousarray(a)).cast('B'))
        else:
            for k in sorted(n):
                require(isinstance(n.get(k, getlink=True), h5py.HardLink), 'HDF5中存在外部/软链接。')
                walk(n[k], relative+'/'+k)
    walk(node, '')
    return h.hexdigest()


def protected_h5(path):
    import h5py
    with h5py.File(path, 'r') as f:
        return {k:hash_h5_tree(f[k]) for k in ('X','layers','var','varm','varp','obsm','obsp','raw') if k in f}


def verify_raw_storage(path, shape):
    import h5py
    import numpy as np
    with h5py.File(path, 'r') as f:
        enc = _json_value(f['X'].attrs.get('encoding-type'))
        require(enc in ('csr_matrix','csc_matrix'), 'X必须为完整CSR/CSC raw counts。')
        require(tuple(f['X'].attrs['shape']) == tuple(shape), 'raw counts形状不符。')
        require('log1p' in f['layers'], '缺少全基因layers/log1p。')
        require(tuple(f['layers/log1p'].attrs['shape']) == tuple(shape), 'log1p形状不符。')
        for slot in ('X', 'layers/log1p'):
            ds=f[slot]['data']
            for i in range(0,len(ds),4_000_000):
                x=ds[i:i+4_000_000]
                require(np.isfinite(x).all() and (x>=0).all(), f'{slot}含非有限或负值。')
                if slot=='X' and x.dtype.kind not in 'iu':
                    require(np.equal(x,np.floor(x)).all(), 'X不再是整数值raw counts。')


def file_stat(path):
    s=Path(path).stat();return (s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns)


def sync_file(path):
    with Path(path).open('rb') as f: os.fsync(f.fileno())


def sync_dir(path):
    fd=os.open(path,os.O_RDONLY)
    try: os.fsync(fd)
    finally: os.close(fd)


@contextlib.contextmanager
def project_lock(root):
    import fcntl
    lock=Path(root)/'.L1L2_finalization.lock'
    with lock.open('a+') as handle:
        try: fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError: raise RuntimeError('另一个最终化进程正在运行；不要重复启动。')
        handle.seek(0);handle.truncate();handle.write(str(os.getpid()));handle.flush()
        try: yield
        finally: fcntl.flock(handle,fcntl.LOCK_UN)


def labels_from_obs(obs):
    import numpy as np
    keys=[k for k in (L2,'L2_annotation') if k in obs]
    require(keys, '当前文件没有L2工作标签；禁止对L1-only文件直接冻结。')
    x=obs[keys[0]]
    require(x.notna().all(), '当前L2存在真正缺失。')
    x=x.astype(str).rename(L2)
    for key in keys[1:]:
        require(np.array_equal(x.to_numpy(),obs[key].astype(str).to_numpy()), '两个L2注释字段存在冲突。')
    return x,keys


def require_idle_review(review):
    review=Path(review)
    import fcntl
    scheduler_lock=review/'focus_scheduler.lock'
    if scheduler_lock.exists():
        with scheduler_lock.open('r') as handle:
            try: fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
            except BlockingIOError: raise RuntimeError('上一轮定点分析调度器仍在运行；请等其结束。')
            else: fcntl.flock(handle,fcntl.LOCK_UN)
    if (review/'FOCUS_STATUS.json').exists():
        require(load_json(review/'FOCUS_STATUS.json').get('status')=='PAUSED_FOR_L2_COMPLETION_REVIEW',
                '上一轮定点任务未全部完成，请先检查FOCUS_STATUS.json。')
    for lock in ('FOCUS.lock','SCHEDULER.lock'):
        p=review/lock
        if p.is_file():
            try: pid=int(p.read_text().strip());os.kill(pid,0)
            except (ValueError,ProcessLookupError): continue
            except PermissionError: raise RuntimeError(f'无法判断活动任务：{p}')
            else: raise RuntimeError(f'仍有可能读写输入的任务PID {pid}，请等待其结束。')
    contracts=load_json(BUNDLE/'inputs/source_contracts.json')
    hashes=load_json(BUNDLE/'inputs/review_summary_checksums.json')
    for branch in contracts:
        d=review/'branches'/branch
        st=load_json(d/'STATUS.json')
        require(st.get('status')=='PAUSED_FOR_L2_COMPLETION_REVIEW', f'{branch}定点审核未正常完成。')
        require(sha_file(d/'cluster_sources_updated_metadata.tsv')==hashes[branch],
                f'{branch}不是本次上传审核的汇总快照。')


def load_maps(review):
    contracts=load_json(BUNDLE/'inputs/source_contracts.json')
    paths=load_json(Path(review)/'inputs/resolved_sources.json')
    sigs=load_json(BUNDLE/'inputs/reviewed_membership_signatures.json')
    maps={}
    for branch,c in contracts.items():
        require(branch in paths, f'缺少{branch}来源路径。')
        p=Path(paths[branch]['mapping'])
        require(p.is_file(), f'原始成员表不存在：{p}')
        require(sha_file(p)==c['mapping_sha256'], f'{branch}成员文件不属于已审核版本。')
        m=read_table(p)
        require(m.cell_id.is_unique and len(m)==c['n_cells'], f'{branch}细胞ID重复/数量变化。')
        require(sha_index(m.cell_id)==sigs[branch]['cell_order_sha256'], f'{branch}细胞顺序变化。')
        require(sha_index(m.sample_uid)==sigs[branch]['sample_order_sha256'], f'{branch}样本对应变化。')
        for k,h in sigs[branch]['group_hashes'].items():
            require(sha_index(m[k])==h, f'{branch}分区成员变化：{k}')
        maps[branch]=m
    return maps,contracts,paths


def derive_review_units(obs,maps,contracts):
    """Exact previous display units: known labels OR Pending::<source UID>.
    The already-reviewed Fib r1.2_16 residual takes priority over its r0.8 parents.
    No new clustering or guessed numeric suffixes are used.
    """
    import pandas as pd
    import numpy as np
    current,_=labels_from_obs(obs)
    units=pd.Series('',index=obs.index,dtype=object)
    branches=pd.Series('',index=obs.index,dtype=object)
    seen=np.zeros(len(obs),dtype=bool)
    for b,c in contracts.items():
        m=maps[b]; ix=obs.index.get_indexer(m.cell_id)
        require((ix>=0).all() and not seen[ix].any(),f'{b}含缺失/跨谱系重叠细胞。')
        require(obs[L1].iloc[ix].astype(str).eq(c['lineage']).all(),f'{b}父L1与原映射不符。')
        require(np.array_equal(obs.sample_uid.iloc[ix].astype(str).to_numpy(),m.sample_uid.to_numpy()), f'{b}sample_uid错配。')
        seen[ix]=True;branches.iloc[ix]=b
        labels=current.iloc[ix].to_numpy(); u=labels.astype(object,copy=True)
        pending=labels=='Pending'
        u[pending]='Pending::'+m.loc[pending,c['selected_key']].astype(str).to_numpy()
        if b=='Fibroblast':
            special=pending & m['lineage_L2_r1p2'].eq('FibL2_R1_r1p2_16').to_numpy()
            u[special]='Pending::FibL2_R1_r1p2_16'
        units.iloc[ix]=u
    require(seen.sum()==sum(c['n_cells'] for c in contracts.values()),'纳入谱系成员数量错误。')
    require(obs[L1].iloc[np.flatnonzero(~seen)].astype(str).isin(['Unknown','Schwann cell']).all(), '未覆盖某个L2父类。')
    return units,branches


def make_final_plan(obs,maps,contracts,decisions,rename,expected_units):
    import pandas as pd
    import numpy as np
    current,oldkeys=labels_from_obs(obs)
    units,branches=derive_review_units(obs,maps,contracts)
    counts=pd.DataFrame({'branch':branches,'review_unit':units})
    got=counts[counts.branch.ne('')].value_counts().rename('n_cells').sort_index()
    want=expected_units.assign(n_cells=expected_units.n_cells.astype(int)).set_index(['branch','review_unit']).n_cells.sort_index()
    require(got.to_dict()==want.to_dict(), '当前逐细胞标签/分区与已审核170个单元不一致。')
    require(decisions[['branch','review_unit']].duplicated().sum()==0,'决定表含重复单元。')
    expected_pending=set(map(tuple,decisions[['branch','review_unit']].to_numpy()))
    observed_pending=set(map(tuple,counts.loc[current.eq('Pending'),['branch','review_unit']].drop_duplicates().to_numpy()))
    require(expected_pending==observed_pending,'75个待定单元不完整或存在新增待定群。')
    after=current.copy();why=pd.Series('',index=obs.index,dtype=object);handled=np.zeros(len(obs),bool)
    for row in decisions.itertuples(index=False):
        select=(branches.eq(row.branch)&units.eq(row.review_unit)&current.eq('Pending')).to_numpy()
        require(int(select.sum())==int(row.n_cells),f'仅当前Pending残余数量不符：{row.review_unit}')
        require(not handled[select].any(),'决定成员重叠。');handled[select]=True
        require(row.level2_annotation not in FORBIDDEN and 'unspecified' not in row.level2_annotation.lower(), '决定表含不允许最终标签。')
        after.iloc[np.flatnonzero(select)]=row.level2_annotation
        why.iloc[np.flatnonzero(select)]=row.review_unit
    require(np.array_equal(handled,current.eq('Pending').to_numpy()),'尚有Pending未覆盖，不能冻结。')
    for row in rename.itertuples(index=False):
        m=maps[row.branch];ids=m.loc[m[row.cluster_key].eq(row.original_cluster_uid),'cell_id']
        require(len(ids)==int(row.n_cells),'更名单元数量不符。')
        require(current.loc[ids].eq(row.old_name).all(),'待更名的旧标签已发生变化。')
        after.loc[ids]=row.level2_annotation;why.loc[ids]='RENAME::'+row.original_cluster_uid
    require(after.loc[current.eq('Unknown')].eq('Unknown').all(),'原Unknown被错误回收。')
    require(after.loc[obs[L1].astype(str).eq('Schwann cell')].eq('Not_assessed').all(),'Schwann规则改变。')
    require(after.loc[obs[L1].astype(str).eq('Unknown')].eq('Unknown').all(),'L1 Unknown规则改变。')
    require(not after.isin(FORBIDDEN).any() and not after.str.contains('unspecified',case=False).any(),'存在未完成或禁止的标签。')
    changed=current.ne(after)
    audit=pd.DataFrame({'cell_id':obs.index[changed],L1:obs.loc[changed,L1].astype(str).to_numpy(),
        'old_level2_annotation':current[changed].to_numpy(),L2:after[changed].to_numpy(),
        'decision_review_unit':why[changed].to_numpy()})
    return current,after,units,branches,audit


def marker_anchor(label):
    if label in SPECIAL: return ''
    first=label.split()[0]
    for suffix in ('-associated','-expressing','-high'):
        first=first.removesuffix(suffix)
    require(bool(re.fullmatch('[A-Za-z0-9-]+',first)),f'无法严格解析marker锚点：{label}')
    return first


def make_ontology(obs,after,var,decisions,rename):
    import pandas as pd
    frame=pd.DataFrame({L1:obs[L1].astype(str).to_numpy(),L2:after.to_numpy()})
    grouping=frame.groupby([L1,L2],sort=True).size().rename('n_cells').reset_index()
    nonspecial=grouping[~grouping[L2].isin(SPECIAL)]
    require(nonspecial.groupby(L2)[L1].nunique().max()==1, '同一个L2名称对应多个父类。')
    scopes=decisions.groupby('level2_annotation').label_scope.agg(lambda x:';'.join(sorted(set(x)))).to_dict()
    rows=[]
    for r in grouping.itertuples(index=False,name=None):
        parent,name,n=r
        if name in SPECIAL: continue  # special-count file is separate; they are NOT biological states.
        anchor=marker_anchor(name)
        require(anchor in var.index,f'最终锚点不在全基因表中：{anchor}')
        identifier='L2_'+hashlib.sha256((parent+'\0'+name).encode()).hexdigest()[:12]
        rows.append(dict(level2_state_id=identifier,level2_annotation=name,level1_annotation=parent,
            marker_anchor=anchor,state_description=name[len(name.split()[0]):].strip(),n_cells=int(n),approved=True,
            label_scope=scopes.get(name,'previously_approved_expression_state'),
            interpretation='transcriptomic annotation; anchor need not be unique or detected in every cell',
            covarnet_eligibility_review='NOT_YET_ASSESSED'))
    out=pd.DataFrame(rows)
    require(out.level2_state_id.is_unique and out.level2_annotation.is_unique,'ontology键不唯一。')
    return out,grouping


def write_sample_metadata(obs,path):
    """Do not turn per-cell QC into a spurious sample-level conflict."""
    import pandas as pd
    keys=[k for k in META_KEYS if k in obs and k!='sample_uid']
    rows=[]
    for sid,ix in obs.groupby('sample_uid',observed=True,sort=False).indices.items():
        g=obs.iloc[ix];d={'sample_uid':str(sid),'n_cells':len(ix)}
        for k in keys:
            vals=g[k].dropna().unique()
            require(len(vals)<=1,f'临床字段在同一个sample内多值：{sid}/{k}；不自动取第一条。')
            d[k]=vals[0] if len(vals) else None
        rows.append(d)
    result=pd.DataFrame(rows)
    require(result.sample_uid.is_unique,'样本表键重复。')
    result.to_csv(path,index=False,encoding='utf-8-sig')
    return result


def smoke_anndata_io(obs,out):
    """Test the installed writer using ACTUAL obs dtypes before copying a large file."""
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp
    import numpy as np
    read_elem,write_elem=io_funcs()
    import h5py
    p=Path(out)/'anndata_io_smoke.h5ad'
    tiny=ad.AnnData(X=sp.csr_matrix(np.eye(3,dtype=np.int32)),obs=obs.iloc[:3].copy(),var=pd.DataFrame(index=['g1','g2','g3']))
    tiny.layers['log1p']=tiny.X.astype(np.float32);tiny.uns['test']=True
    tiny.write_h5ad(p)
    with h5py.File(p,'r+') as f:
        del f['obs'];write_elem(f,'obs',tiny.obs,dataset_kwargs={'compression':'lzf'})
    check=ad.read_h5ad(p)
    pd.testing.assert_frame_equal(check.obs,tiny.obs,check_dtype=False,check_categorical=False)
    require((check.X!=tiny.X).nnz==0,'安装环境小AnnData读写失败。')
    return {'status':'PASS','test':'native tiny AnnData with actual metadata dtypes'}


def _set_uns(group,key,value,write_elem):
    if key in group: del group[key]
    write_elem(group,key,value)


def stage_h5(source,stage,new_obs,values,protected):
    import h5py
    import pandas as pd
    read_elem,write_elem=io_funcs()
    shutil.copyfile(source,stage)
    with h5py.File(stage,'r+') as f:
        del f['obs'];write_elem(f,'obs',new_obs,dataset_kwargs={'compression':'lzf'})
        u=f['uns']
        for key in ('level2_annotation_colors','L2_annotation_colors'):
            if key in u:del u[key]
        for key,value in values.items():_set_uns(u,key,value,write_elem)
        f.flush()
    require(protected_h5(stage)==protected,'存盘后表达/基因/坐标/图发生变化，拒绝发布。')
    with h5py.File(stage,'r') as f:
        saved=read_elem(f['obs'])
        pd.testing.assert_frame_equal(saved,new_obs,check_dtype=False,check_categorical=False)
        require(read_elem(f['uns']['level2_frozen'])==True,'冻结标志未保存。')
    sync_file(stage)


def transactional_publish(stage,target,canonical,archive,expected_stat):
    """Replace explicit working file with backup. Canonical and working paths become
    hard links to one fully embedded h5ad, NOT two expression copies or HDF5 external links.
    Both names are read-only after commit. Roll back target if canonical install fails.
    """
    stage,target,canonical,archive=map(Path,(stage,target,canonical,archive))
    require(file_stat(target)==tuple(expected_stat),'待覆盖文件被其他任务修改。')
    require(canonical==target or not os.path.lexists(canonical),'正式路径已存在其他文件；禁止无提示覆盖最终release。')
    backup=archive/target.name;require(not backup.exists(),'备份已存在。')
    alias_temp=canonical.with_name(canonical.name+'.link_staging')
    require(not os.path.lexists(alias_temp),'存在未完成的发布链接。')
    linked=False;committed=False;canonical_installed=False
    try:
        # Hard-link preflight BEFORE altering any existing file. Same filesystem required.
        if canonical!=target:os.link(stage,alias_temp);linked=True
        os.link(target,backup)
        sync_dir(archive)
        require(file_stat(target)==tuple(expected_stat),'建立备份期间目标变化。')
        os.replace(stage,target);committed=True
        if linked:
            os.replace(alias_temp,canonical);linked=False;canonical_installed=True
        sync_dir(target.parent)
    except BaseException:
        if canonical_installed and os.path.lexists(canonical): canonical.unlink()
        if committed:
            rollback=target.with_name(target.name+'.rollback_link')
            if os.path.lexists(rollback):rollback.unlink()
            os.link(backup,rollback);os.replace(rollback,target);sync_dir(target.parent)
        if linked and os.path.lexists(alias_temp):alias_temp.unlink()
        raise
    return backup



def finish_publication(manifest, run):
    """Commit/recover sidecars only for a verified full h5ad and approved manifest.
    May be rerun after an interruption between h5ad replacement and sidecar writes.
    """
    run=Path(run);target=Path(manifest['working_path_alias']);canonical=Path(manifest['h5ad_path'])
    root=canonical.parent;release=Path(manifest['release_support_dir']);support=run/'release_staging'
    expected=manifest['h5ad_sha256']
    require(sha_file(target)==expected,'已写入的工作文件与已核验stage摘要不同，拒绝补写发布信息。')
    if canonical!=target:
        if canonical.exists(): require(os.path.samefile(target,canonical),'正式路径已指向其他文件。')
        else:
            alias=canonical.with_name(canonical.name+'.link_staging')
            if os.path.lexists(alias):
                require(os.path.samefile(alias,target),'未完成链接与当前已写文件不同。')
                os.replace(alias,canonical)
            else: os.link(target,canonical)
    backup=Path(manifest['backup_h5ad'])
    require(backup.is_file() and sha_file(backup)==manifest['source_h5ad_sha256'], '旧全量文件备份不符，停止发布。')
    active=release if release.exists() else support
    require(active.is_dir(), '缺少已核验发布附表。')
    for name,digest in manifest['release_sidecar_sha256'].items():
        require(sha_file(active/name)==digest,'发布附表已改变：'+name)
    if active==support:
        release.parent.mkdir(parents=True,exist_ok=True);os.replace(support,release)
    for path in dict.fromkeys((canonical,target)):
        dump_json(path.with_suffix('.manifest.json'),manifest)
        checksum=path.with_suffix('.h5ad.sha256');tmp=checksum.with_name(checksum.name+'.tmp')
        tmp.write_text(expected+'  '+path.name+'\n',encoding='utf8');os.replace(tmp,checksum)
    statepath=root/'PROJECT_STATE_ATLAS_L1L2_FINAL_v1.md'
    prev=run/'backup'/statepath.name
    if statepath.exists() and not prev.exists(): shutil.copy2(statepath,prev)
    special=manifest['level2_special_counts_verified'];named=manifest['shape_verified'][0]-sum(special.values())
    text=f"""# CPU L1/L2 final release v1

L1_frozen = TRUE
L2_frozen = TRUE
Downstream execution = NOT AUTHORIZED

Canonical full h5ad: {canonical}
SHA256: {expected}
Source backup: {backup}

X: unchanged all-gene raw counts. layers[log1p]: unchanged CP10k+log1p.
Clinical metadata, sample/donor identities and global coordinates/graphs: unchanged.
Formal annotation fields: level1_annotation and level2_annotation.
No Pending/unspecified label. Schwann remains Not_assessed; L1 Unknown remains Unknown.

Named cells: {named}
Unknown cells: {special['Unknown']}
Not_assessed cells: {special['Not_assessed']}
Named L2 labels: {manifest['named_level2_states_verified']}

Some labels are marker-parent groupings, not independently resolved functional states.
Do not treat all named labels as automatically eligible CoVarNet nodes. See external ontology.
No new clustering, integration, differential testing or clinical imputation was performed.
Release sidecars: {release}
This release does not authorize executing downstream scientific models.
"""
    temp=statepath.with_name(statepath.name+'.tmp');temp.write_text(text,encoding='utf8');os.replace(temp,statepath)
    journal=load_json(run/'FREEZE_STATUS.json')
    journal.update(status='COMMITTED_FINAL_APPROVED',h5ad=str(canonical),backup_h5ad=str(backup),
                   release_support_dir=str(release),validation=load_json(release/'validation_report.json'))
    dump_json(run/'FREEZE_STATUS.json',journal);sync_dir(root)
    return manifest


def finalize(target,root=ROOT,review_run=None):
    import numpy as np
    import pandas as pd
    import h5py
    verify_bundle()
    root=Path(root).resolve();target=Path(target).resolve();canonical=root/FINAL_NAME
    require(target.is_file() and target.parent==root,'输入必须是指定根目录下明确的全量L2工作h5ad。')
    require(target.name!='panvascular_CPU_theta2_level1_frozen_v1.h5ad' and 'raw_counts' not in target.name.lower(),
            '禁止覆盖raw源文件或冻结L1-only母图谱。')
    read_elem,write_elem=io_funcs()
    decisions=read_table(BUNDLE/'inputs/final_pending_decisions.tsv')
    rename=read_table(BUNDLE/'inputs/approved_existing_name_refinement.tsv')
    spec_sha=hashlib.sha256(json.dumps({p:sha_file(BUNDLE/'inputs'/p) for p in (
        'final_pending_decisions.tsv','approved_existing_name_refinement.tsv','reviewed_current_unit_counts.tsv',
        'source_contracts.json')},sort_keys=True).encode()).hexdigest()
    with project_lock(root):
        with h5py.File(target,'r') as f:
            obs=read_elem(f['obs']);var=read_elem(f['var']);uns=read_elem(f['uns'])
        if bool(uns.get('level2_frozen',False)):
            require(uns.get('L2_final_decision_sha256')==spec_sha,'已有不同冻结版本；不重复覆盖。')
            prior_run=Path(uns['L2_finalization_run'])
            ready=load_json(prior_run/'audit/PREPARED_RELEASE_MANIFEST.json')
            require(ready['decision_spec_sha256']==spec_sha,'发布决定摘要不符。')
            # Also completes a previously interrupted sidecar transaction, never reclassifies cells.
            manifest=finish_publication(ready,prior_run)
            print('本次最终决定已发布并复核；不重复写入表达或标签。',flush=True)
            return prior_run,canonical
        require(not canonical.exists() or canonical==target,'正式v1文件已存在，不能自动覆盖不同发布。')
        require((len(obs),len(var))==EXPECTED_SHAPE and obs.index.is_unique and var.index.is_unique,'图谱维度或唯一ID不符。')
        require(L1 in obs and obs[L1].notna().all() and obs[L1].astype(str).value_counts().to_dict()==EXPECTED_L1,'冻结L1变化。')
        require(bool(uns.get('level1_frozen',uns.get('L1_frozen',False))),'L1尚未冻结。')
        require(obs.sample_uid.notna().all() and obs.sample_uid.nunique()==247,'样本分组基线变化。')
        review=Path(review_run or str(uns.get('L2_last_update_run',''))).resolve()
        require((review/'UPDATE.json').is_file(),'请显式提供上一轮L2_NEXT_RUN目录。')
        rstate=load_json(review/'UPDATE.json')
        require(rstate.get('status')=='COMMITTED','上一次L2工作文件未完成安全保存。')
        require_idle_review(review)
        # This registry binds the reviewed full labels; metadata-only newer values can be retained.
        previous=read_table(review/'tables/all_cells_level2_annotation.tsv.gz')
        current,oldkeys=labels_from_obs(obs)
        require(np.array_equal(previous.cell_id.to_numpy(),obs.index.astype(str).to_numpy()),'母细胞顺序与上轮审批登记不同。')
        for col in ('sample_uid',L1):require(np.array_equal(previous[col].to_numpy(),obs[col].astype(str).to_numpy()), f'上轮登记中的{col}不符。')
        require(np.array_equal(previous[L2].to_numpy(),current.to_numpy()), '输入不是本次审核的完整标签版本。')
        del previous
        maps,contracts,paths=load_maps(review)
        expected_units=read_table(BUNDLE/'inputs/reviewed_current_unit_counts.tsv')
        before,after,units,branches,audit=make_final_plan(obs,maps,contracts,decisions,rename,expected_units)
        expect=load_json(BUNDLE/'inputs/expected_counts.json')
        require(after.eq('Unknown').sum()==expect['final_Unknown_cells'] and after.eq('Not_assessed').sum()==2832,'最终特殊状态数量不符。')
        require((~after.isin(SPECIAL)).sum()==expect['final_named_cells'],'最终命名细胞数量不符。')
        run=root/('L2_finalization_'+stamp());run.mkdir()
        for d in ('backup','audit','tables','release_staging'):(run/d).mkdir()
        shutil.copytree(BUNDLE/'inputs',run/'audit/decision_inputs')
        shutil.copy2(__file__,run/'audit'/Path(__file__).name)
        audit.to_csv(run/'tables/final_label_changes.tsv.gz',sep='\t',index=False)
        before.rename('previous_level2_annotation').to_csv(run/'tables/previous_labels.tsv.gz',sep='\t',index_label='cell_id')
        # Do not broadly regex-delete metadata. Only exact redundant identity keys are migrated.
        drop=list(oldkeys)
        if 'L1' in obs:
            require(obs['L1'].astype(str).equals(obs[L1].astype(str)),'历史L1列与冻结标签不一致。')
            drop.append('L1')
        new_obs=obs.drop(columns=drop).copy();new_obs[L2]=pd.Categorical(after.to_numpy(),categories=sorted(after.unique()))
        pd.testing.assert_frame_equal(new_obs.drop(columns=L2),obs.drop(columns=drop),check_exact=True)
        ontology,counts=make_ontology(new_obs,after,var,decisions,rename)
        require(len(ontology)==expect['named_level2_labels'],'预期最终名称数量不符。')
        support=run/'release_staging';release=root/'atlas_release'/VERSION
        require(not release.exists(),'正式support目录已存在；请不要覆盖历史发布。')
        ontology.to_csv(support/'level2_ontology.tsv',sep='\t',index=False)
        counts.to_csv(support/'level2_counts.tsv',sep='\t',index=False)
        counts[counts[L2].isin(SPECIAL)].to_csv(support/'special_label_counts.tsv',sep='\t',index=False)
        obs[L1].astype(str).value_counts().rename_axis(L1).rename('n_cells').to_csv(support/'level1_counts.tsv',sep='\t')
        new_obs[[L1,L2]].to_csv(support/'cell_annotations.tsv.gz',sep='\t',index_label='cell_id')
        write_sample_metadata(new_obs,support/'sample_metadata.csv')
        var.to_csv(support/'gene_metadata.tsv.gz',sep='\t',index_label='gene_id')
        approvals=pd.concat([decisions,rename],ignore_index=True).fillna('')
        approvals.to_csv(support/'approved_level2_mapping.tsv',sep='\t',index=False)
        pd.DataFrame([{'column':c,'action':'REPLACED_BY_'+L2 if c in oldkeys else 'ARCHIVED_REDUNDANT_L1'} for c in drop]).to_csv(support/'annotation_column_migration.tsv',sep='\t',index=False)
        pd.DataFrame([{'role':'lineage_analysis','branch':b,'path':paths[b]['analysis'],'sha256':contracts[b]['analysis_sha256'],
          'mapping_path':paths[b]['mapping'],'mapping_sha256':contracts[b]['mapping_sha256']} for b in contracts]).to_csv(support/'source_registry.tsv',sep='\t',index=False)
        sample_meta_sha=sha_file(support/'sample_metadata.csv')
        dump_json(run/'audit/native_io_smoke.json',smoke_anndata_io(new_obs,run/'audit'))
        source_sig=file_stat(target);source_sha=sha_file(target)
        print('只读核验全量raw、log1p与原坐标/图。不会运行任何模型……',flush=True)
        verify_raw_storage(target,EXPECTED_SHAPE)
        protected=protected_h5(target)
        previous_protected=load_json(review/'audit/protected_before.json')
        require(protected==previous_protected,'表达/基因/坐标/图与上一轮审核工作版本不符。')
        keymap={'pca':'X_pca_cpu','harmony':'X_harmony_cpu','umap':'X_umap_cpu','neighbors_uns':'neighbors_cpu_n30'}
        with h5py.File(target,'r') as f:
            for role,key in keymap.items():
                slot='uns' if role=='neighbors_uns' else 'obsm'
                require(key in f[slot],f'最终契约缺少{slot}/{key}。')
            neighbors=read_elem(f['uns'][keymap['neighbors_uns']])
            for field in ('distances_key','connectivities_key'):
                require(field in neighbors and neighbors[field] in f['obsp'],'邻居图指针失效。')
            keymap.update({f:neighbors[f] for f in ('distances_key','connectivities_key')})
            inventory={k:list(f[k]) for k in ('layers','obsm','obsp','varm','varp','uns') if k in f}
        metadata_version=str(uns.get('metadata_version',uns.get('source_metadata_version','inherited_as_saved')))
        values={'atlas_schema_version':SCHEMA,'atlas_release_version':VERSION,'release_status':'FINAL_APPROVED',
            'level1_frozen':True,'level2_frozen':True,'L1_frozen':True,'L2_frozen':True,
            'annotation_keys':{'level1':L1,'level2':L2},'L2_annotation_key':L2,'level2_annotation_key':L2,
            'level1_annotation_key':L1,'active_annotation_key':L2,'L2_status':'FINAL_FROZEN',
            'status':'L1L2_FINAL_APPROVED','current_full_h5ad':str(canonical),'frozen_h5ad':str(canonical),
            'source_h5ad':str(target),'source_h5ad_sha256':source_sha,'L2_finalization_run':str(run),
            'L2_final_decision_sha256':spec_sha,'level2_freeze_time':datetime.now().isoformat(),
            'source_metadata_version':metadata_version,'source_metadata_snapshot_sha256':sample_meta_sha,
            'global_embedding_keys':keymap,'level2_ontology_file':str(release/'level2_ontology.tsv'),
            'matrix_semantics':{'X':'unchanged all-gene raw counts','layers/log1p':'all-gene CP10k then natural log1p; not Harmony corrected'},
            'downstream_execution_authorized':False,'L2_pending_cells':0}
        stage=root/('.'+FINAL_NAME+'.'+stamp()+'.partial.h5ad')
        require(shutil.disk_usage(root).free>target.stat().st_size+2*1024**3,'磁盘不足：需要至少一份全量临时h5ad的空间。')
        journal={'status':'PREPARING','source':str(target),'canonical_h5ad':str(canonical),'stage':str(stage),'run':str(run)}
        dump_json(run/'FREEZE_STATUS.json',journal)
        try:
            print('创建完整临时副本并写入level2_annotation……',flush=True)
            stage_h5(target,stage,new_obs,values,protected)
            require(file_stat(target)==source_sig and sha_file(target)==source_sha,'源文件在保存期间被修改。')
            finalsha=sha_file(stage);stage.chmod(stat.S_IMODE(target.stat().st_mode)&~0o222)
            validation={'status':'PASS','shape':list(EXPECTED_SHAPE),'sample_count':247,'L1_unchanged':True,
                'metadata_unchanged':True,'protected_datasets_unchanged':True,'no_Pending':True,'no_unspecified':True,
                'n_named':int((~after.isin(SPECIAL)).sum()),'n_unknown':int(after.eq('Unknown').sum()),'n_not_assessed':2832,
                'n_named_states':len(ontology),'n_changed':len(audit),'new_scientific_analysis_performed':False,
                'new_donor_inference_performed':False,'source_sha256':source_sha,'output_sha256':finalsha}
            dump_json(support/'validation_report.json',validation)
            dump_json(support/'schema_inventory.json',dict(obs_columns=list(new_obs.columns),var_columns=list(var.columns),slots=inventory))
            dump_json(support/'matrix_fingerprints.json',protected)
            sidecars={p.name:sha_file(p) for p in support.iterdir() if p.is_file()}
            manifest={'schema_version':SCHEMA,'release_version':VERSION,'release_status':'FINAL_APPROVED','ready_for_downstream':True,
                'h5ad_path':str(canonical),'h5ad_sha256':finalsha,'source_h5ad_path':str(target),'source_h5ad_sha256':source_sha,
                'source_metadata_version':metadata_version,'source_metadata_sha256':sample_meta_sha,
                'level1_key':L1,'level2_key':L2,'level1_frozen':True,'level2_frozen':True,
                'shape_expected':list(EXPECTED_SHAPE),'shape_verified':list(EXPECTED_SHAPE),'sample_count_verified':247,
                'level1_counts_verified':EXPECTED_L1,'named_level2_states_verified':len(ontology),
                'level2_special_counts_verified':{'Unknown':int(after.eq('Unknown').sum()),'Not_assessed':2832},
                'matrix_fingerprints':protected,'cell_id_order_sha256':sha_index(obs.index),
                'gene_id_order_sha256':sha_index(var.index),'annotation_sha256':sha_index(after),
                'level1_annotation_sha256':sha_index(obs[L1].astype(str)),'global_key_map':keymap,
                'release_sidecar_sha256':sidecars,'release_support_dir':str(release),
                'validation_report_path':str(release/'validation_report.json'),'finalization_run':str(run),
                'decision_table_sha256':sha_file(BUNDLE/'inputs/final_pending_decisions.tsv'),
                'decision_spec_sha256':spec_sha,'downstream_execution_authorized':False,
                'working_path_alias':str(target),'alias_type':'same_filesystem_hardlink' if target!=canonical else 'same_path',
                'annotation_interpretation':'mixed resolution: named states plus marker-parent groups; not all 116 labels are independently validated biological subtypes',
                'covarnet_eligibility':'separate review; not inferred from annotation approval'}
            # Preserve old sidecars before replacing them.
            for path in (target.with_suffix('.manifest.json'),target.with_suffix('.h5ad.sha256')):
                if path.exists():shutil.copy2(path,run/'backup'/path.name)
            manifest['backup_h5ad']=str(run/'backup'/target.name)
            dump_json(run/'audit/PREPARED_RELEASE_MANIFEST.json',manifest)
            journal.update(status='VERIFIED_STAGING',output_sha256=finalsha);dump_json(run/'FREEZE_STATUS.json',journal)
            backup=transactional_publish(stage,target,canonical,run/'backup',source_sig)
            journal.update(status='H5AD_COMMITTED_SIDECARS_PENDING',backup_h5ad=str(backup));dump_json(run/'FREEZE_STATUS.json',journal)
            finish_publication(manifest,run)
        except BaseException:
            journal=load_json(run/'FREEZE_STATUS.json')
            journal['error']=traceback.format_exc();dump_json(run/'FREEZE_STATUS.json',journal)
            raise
        print('\n最终L2保存完成。',flush=True)
        print('正式文件：',canonical,flush=True)
        print('原工作路径已替换（指向同一份完整文件）：',target,flush=True)
        print('原版备份：',backup,flush=True)
        print('labels=',len(ontology),'named=',validation['n_named'],'Unknown=',validation['n_unknown'],'Not_assessed=2832',flush=True)
        print('未运行下游；kernel中旧对象未自动刷新，请勿再用旧对象覆盖新文件。',flush=True)
        return run,canonical


def run_from_notebook(*,target_h5ad,root=ROOT,review_run=None,threads=4):
    """Blocking CPU subprocess; no live notebook AnnData mutation."""
    require(isinstance(threads,int) and 1 <= threads <= 40,'threads必须在1–40之间。')
    root=Path(root).resolve();target=Path(target_h5ad).resolve();verify_bundle()
    launch=root/('L2_final_freeze_launch_'+stamp());launch.mkdir(parents=True,exist_ok=False)
    request={'target':str(target),'root':str(root),'review_run':str(review_run) if review_run else None,
             'result':str(launch/'RESULT.json')}
    dump_json(launch/'REQUEST.json',request)
    env=os.environ.copy();env['CUDA_VISIBLE_DEVICES']=''
    for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'):
        env[k]=str(int(threads))
    log=launch/'run.log'
    command=[sys.executable,str(Path(__file__).resolve()),'--request',str(launch/'REQUEST.json')]
    with log.open('w',encoding='utf8') as handle:
        proc=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1,env=env)
        try:
            for line in proc.stdout:
                print(line,end='',flush=True);handle.write(line);handle.flush()
            rc=proc.wait()
        except BaseException:
            proc.terminate();proc.wait();raise
    require(rc==0,f'发布没有完成。检查日志：{log}；不要把partial或不完整发布当final。')
    res=load_json(launch/'RESULT.json')
    return Path(res['run']),Path(res['final_h5ad'])


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--request',required=True);args=parser.parse_args()
    req=load_json(args.request)
    require(os.environ.get('CUDA_VISIBLE_DEVICES')=='','请使用run_from_notebook启动CPU子进程。')
    out,path=finalize(req['target'],root=req['root'],review_run=req.get('review_run'))
    dump_json(req['result'],{'status':'COMMITTED_FINAL_APPROVED','run':str(out),'final_h5ad':str(path)})
