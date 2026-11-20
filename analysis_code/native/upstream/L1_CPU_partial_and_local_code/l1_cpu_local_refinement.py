# -*- coding: utf-8 -*-
"""CPU theta=2 atlas: approved partial L1 -> independent local raw-count review.

Import into the EXISTING notebook; call prepare_partial_L1(adata, markers, work_parent)
then run_local_batch(run_dir). No global reintegration and no automatic local writeback.
Local jobs run in new Python processes using the same interpreter, with CUDA hidden.
"""
from __future__ import annotations
import argparse
import base64
import copy
import gc
import hashlib
import html
import importlib
import importlib.metadata as metadata
import inspect
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys
import time
import traceback
import zipfile
from datetime import datetime

# Worker environment is configured by the scheduler BEFORE Python/library startup.
# Importing this module in the notebook does NOT clear or overwrite notebook variables.
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/'
            '06_reanalysis_takeover_20260908')
DEFAULT_PATCH_DIR = ROOT / 'annotation/HVG_clip_review_20260910'
GROUP = 'leiden_res1'
BATCH = 'sample_uid'
EXPECTED_SHAPE = (1175115, 37336)
EXPECTED_COUNTS = {
    '0':28042, '1':138194, '2':25610, '3':46641, '4':125250,
    '5':103761, '6':87560, '7':77491, '8':48146, '9':6639,
    '10':15272, '11':76270, '12':34866, '13':193769, '14':1682,
    '15':4099, '16':22978, '17':12432, '18':39694, '19':4945,
    '20':11973, '21':6066, '22':27074, '23':4207, '24':12929,
    '25':11967, '26':5400, '27':2088, '28':70,
}
APPROVED = {
    '13':'Fibroblast', '17':'Fibroblast',
    '6':'Smooth muscle cell', '7':'Smooth muscle cell',
    '11':'Endothelial cell', '22':'Endothelial cell',
    '4':'Macrophage', '12':'Macrophage', '15':'Macrophage', '20':'Macrophage',
    '1':'T cell', '10':'T cell', '24':'Mast cell',
    '21':'Neutrophil', '14':'Schwann cell',
}
LOCAL_TARGETS = ['0','2','3','5','8','9','16','18','19','23','25','26','27']
EXPECTED_L1 = {'Fibroblast':206201, 'Smooth muscle cell':165051,
               'Endothelial cell':103344, 'Macrophage':176188, 'T cell':153466,
               'Mast cell':12929, 'Neutrophil':6066, 'Schwann cell':1682,
               'Unknown':350188}
RESOLUTIONS = (0.3, 0.5, 0.8)
BASELINE_KEY = 'local_leiden_r0p5'
MURAL = ('Fibroblast','Fibromyocyte','Smooth muscle cell','Pericyte')
FOCUS = {
    'C0':'NK/T boundary', 'C2':'B/Plasma boundary',
    'C3':'Macrophage/Monocyte/Neutrophil boundary',
    'C5':'Fibromyocyte/SMC/Fibroblast boundary',
    'C8':'Monocyte/Macrophage/DC boundary', 'C9':'T/Myeloid mixed program',
    'C16':'Fibroblast/Fibromyocyte/SMC boundary', 'C18':'Pericyte/SMC boundary',
    'C19':'T/Stromal mixed program', 'C23':'T/B mixed program',
    'C25':'Macrophage/Fibroblast mixed program',
    'C26':'Endothelial/Mural mixed program',
    'C27':'SPP1 macrophage/Monocyte boundary',
}


def require(condition, message):
    if not bool(condition):
        raise ValueError(message)


def stamp():
    return datetime.now().strftime('%Y%m%d_%H%M%S_%f')


def dump_json(path, obj):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    temp = p.with_name(p.name + '.tmp')
    temp.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=str), encoding='utf-8')
    temp.replace(p)


def file_sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(16 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def index_sha256(index):
    h = hashlib.sha256()
    for value in index:
        b = str(value).encode('utf-8')
        h.update(len(b).to_bytes(8, 'little')); h.update(b)
    return h.hexdigest()


def _sparse_hash(fmt, shape, arrays):
    # Compatible with cpu_atlas_theta2_raw_locked.py's raw_integrity fingerprint.
    h = hashlib.sha256()
    h.update(json.dumps([fmt, [int(v) for v in shape]]).encode())
    for name, arr in arrays:
        h.update(json.dumps([name, np.dtype(arr.dtype).str, int(arr.shape[0])]).encode())
        for start in range(0, arr.shape[0], 4_000_000):
            b = np.ascontiguousarray(arr[start:start + 4_000_000])
            h.update(memoryview(b).cast('B'))
    return h.hexdigest()


def sparse_sha256(x):
    require(sp.issparse(x) and x.format in ('csr','csc'), 'Expected CSR/CSC sparse matrix.')
    return _sparse_hash(x.format, x.shape, [(k, getattr(x, k)) for k in ('data','indices','indptr')])


def h5_sparse_sha256(path, key='X'):
    import h5py
    with h5py.File(path, 'r') as f:
        g = f[key]
        enc = g.attrs['encoding-type']
        if isinstance(enc, bytes): enc = enc.decode()
        require(enc in ('csr_matrix','csc_matrix'), f'{key}: unsupported encoding {enc}')
        return _sparse_hash(enc[:3], g.attrs['shape'], [(k, g[k]) for k in ('data','indices','indptr')])


def dense_sha256(x):
    x = np.asarray(x)
    h = hashlib.sha256(json.dumps([list(x.shape), x.dtype.str]).encode())
    for start in range(0, x.shape[0], 100_000):
        h.update(memoryview(np.ascontiguousarray(x[start:start+100_000])).cast('B'))
    return h.hexdigest()


def protected_hashes(a):
    # Includes all current coordinates/graphs, not just the selected UMAP.
    return {
        'X':sparse_sha256(a.X), 'log1p':sparse_sha256(a.layers['log1p']),
        'cells':index_sha256(a.obs_names), 'genes':index_sha256(a.var_names),
        'obsm':{k:(sparse_sha256(v) if sp.issparse(v) else dense_sha256(v)) for k,v in a.obsm.items()},
        'obsp':{k:(sparse_sha256(v) if sp.issparse(v) else dense_sha256(v)) for k,v in a.obsp.items()},
    }


def lock_sparse(x):
    for key in ('data','indices','indptr'):
        getattr(x,key).flags.writeable = False


def readonly_file(path):
    p = Path(path)
    p.chmod(stat.S_IMODE(p.stat().st_mode) & ~0o222)


def validate_counts(x):
    require(sp.issparse(x) and x.format in ('csr','csc'), 'Raw X must be sparse, not normalized expression.')
    for i in range(0, x.nnz, 4_000_000):
        b = x.data[i:i+4_000_000]
        require(np.isfinite(b).all() and (b >= 0).all(), 'Invalid raw counts.')
        if not np.issubdtype(b.dtype, np.integer):
            require(np.equal(b, np.floor(b)).all(), 'X is not integer-valued raw counts.')


def save_h5ad_atomic(a, path, raw_hash, lock=True):
    p = Path(path)
    require(not p.exists(), f'Refusing to overwrite {p}')
    temp = p.with_name(p.stem + '.partial.h5ad')
    if temp.exists():  # failed prior automatic write only, never a complete approved object
        temp.rename(p.with_name(p.stem + '.incomplete_' + stamp() + '.h5ad'))
    kw = {'compression':'lzf'}
    if 'convert_strings_to_categoricals' in inspect.signature(a.write_h5ad).parameters:
        kw['convert_strings_to_categoricals'] = False
    a.write_h5ad(temp, **kw)
    require(h5_sparse_sha256(temp) == raw_hash, 'Saved raw counts differ from memory.')
    temp.rename(p)
    if lock: readonly_file(p)


def clean_panel(markers, var_names):
    available = set(var_names)
    panel = {str(k):list(dict.fromkeys(g for g in gs if g in available)) for k,gs in markers.items()}
    panel = {k:gs for k,gs in panel.items() if gs}
    missing = {str(k):[g for g in gs if g not in available] for k,gs in markers.items()}
    return panel, {k:v for k,v in missing.items() if v}


def aggregate_markers(a, labels, genes):
    """All cells, raw detection, mean CP10k-log1p; only marker slice becomes dense."""
    labels = pd.Series(np.asarray(labels, dtype=str), index=a.obs_names)
    order = list(dict.fromkeys(labels))
    try: order = sorted(order, key=lambda s:int(s.rsplit('_',1)[-1]))
    except ValueError: order = sorted(order)
    cats = pd.Categorical(labels, categories=order)
    counts = np.bincount(cats.codes, minlength=len(order))
    indexer = a.var_names.get_indexer(genes)
    require((indexer >= 0).all(), 'Unknown marker name.')
    raw = a.X[:, indexer].tocsr(copy=True)
    raw.sum_duplicates(); raw.eliminate_zeros(); raw.sort_indices()
    log = a.layers['log1p'][:, indexer].tocsr(copy=True)
    log.sum_duplicates(); log.eliminate_zeros(); log.sort_indices()
    binary = raw.copy(); binary.data = np.ones(binary.nnz, dtype=np.float64)
    binary_log = log.copy(); binary_log.data = np.ones(binary_log.nnz, dtype=np.float64)
    require((binary != binary_log).nnz == 0, 'Raw/log1p nonzero supports do not agree.')
    membership = sp.csr_matrix((np.ones(a.n_obs, dtype=np.float64),
                                (cats.codes, np.arange(a.n_obs))),shape=(len(order),a.n_obs))
    def means(x):
        return pd.DataFrame((membership @ x).toarray()/counts[:,None],
                            index=pd.Index(order,name='cluster_uid'), columns=genes)
    return means(log), means(binary), means(raw), pd.Series(counts,index=order,name='n_cells')


def save_dotplot(panel, mean, frac, stem, relative=False):
    import anndata as ad
    import scanpy as sc
    import matplotlib.pyplot as plt
    if not panel: return
    cols = [g for gs in panel.values() for g in gs]
    vals = mean.copy()
    if relative:
        vals = (vals-vals.min()).div((vals.max()-vals.min()).replace(0,np.nan)).fillna(0)
    ids = [f'plot_gene_{i:04d}' for i in range(len(cols))]
    plot_panel, pos = {}, 0
    for name, gs in panel.items():
        plot_panel[name] = ids[pos:pos+len(gs)]; pos += len(gs)
    color_df = vals.loc[:,cols].copy(); size_df = frac.loc[:,cols].copy()
    color_df.columns = size_df.columns = ids
    order = mean.index.astype(str).tolist()
    dummy = ad.AnnData(X=color_df.to_numpy(dtype=np.float32),
                      obs=pd.DataFrame({'group':pd.Categorical(order,categories=order)},
                                       index=['summary_'+str(i) for i in range(len(order))]),
                      var=pd.DataFrame(index=ids))
    dp = sc.pl.DotPlot(dummy, var_names=plot_panel, groupby='group', categories_order=order,
                       use_raw=False, standard_scale=None, dot_color_df=color_df,
                       dot_size_df=size_df,
                       figsize=(max(12,0.21*len(cols)), max(5,0.32*len(order)+3)))
    dp.style(dot_min=0,dot_max=1,dot_edge_lw=0)
    dp.legend(colorbar_title=('Gene-wise relative mean (0-1)' if relative else 'Mean log1p(CP10k)'),
              size_title='Cells with raw counts > 0 (%)')
    dp.make_figure()
    ax = dp.get_axes()['mainplot_ax']; ax.set_xticklabels(cols,rotation=90)
    fig = ax.figure
    try:
        for ext in ('png','pdf','svg'):
            fig.savefig(str(stem)+'.'+ext, dpi=250, bbox_inches='tight')
    finally: plt.close(fig)


def save_umap(a, basis, key, stem, point_size=1.0):
    import anndata as ad
    import scanpy as sc
    import matplotlib.pyplot as plt
    from matplotlib.collections import PathCollection
    order = np.random.default_rng(0).permutation(a.n_obs)
    vis = ad.AnnData(obs=a.obs[[key]].iloc[order].copy(),
                     obsm={'X_plot':np.asarray(a.obsm[basis])[order].copy()})
    # Missing metadata are displayed as Missing, never inferred.
    vis.obs[key] = vis.obs[key].astype('string').fillna('Missing').astype('category')
    n = vis.obs[key].nunique()
    palette = sc.pl.palettes.default_102[:n]
    require(len(palette)==n, f'{key}: too many categories for static palette.')
    fig, ax = plt.subplots(figsize=(9,8))
    try:
        sc.pl.embedding(vis,basis='X_plot',color=key,palette=palette,size=point_size,
                        marker='.',alpha=1,frameon=False,title=key,legend_fontsize=8,
                        show=False,save=False,ax=ax)
        ax.set_aspect('equal',adjustable='box')
        for collection in ax.collections:
            if isinstance(collection,PathCollection):
                collection.set_linewidths(0); collection.set_edgecolor('none')
                collection.set_rasterized(True)
        for ext in ('png','pdf','svg'):
            fig.savefig(str(stem)+'.'+ext,dpi=300,bbox_inches='tight',facecolor='white')
    finally: plt.close(fig)


def prepare_partial_L1(adata, markers, work_parent=None, patch_dir=DEFAULT_PATCH_DIR):
    """Apply only approved parent labels in memory; save FULL mother + 13 raw inputs.

    Does not clear notebook variables or change X/log1p/coordinates/graphs.
    Counts per parent must match the exact 29-cluster theta=2 review.
    """
    import anndata as ad
    import scanpy as sc
    require(not adata.isbacked and not adata.is_view, 'Use your existing in-memory full adata, not a backed/view object.')
    require(adata.shape == EXPECTED_SHAPE, f'Unexpected shape {adata.shape}; must use latest 29-cluster atlas.')
    require(adata.obs_names.is_unique and adata.var_names.is_unique, 'Non-unique cell/gene IDs.')
    require(GROUP in adata.obs and not adata.obs[GROUP].isna().any(), 'Missing leiden_res1 labels.')
    require('log1p' in adata.layers and adata.layers['log1p'].shape==adata.shape, 'Full-gene log1p layer required.')
    require(BATCH in adata.obs and not adata.obs[BATCH].isna().any(), 'Missing sample_uid.')
    labels = adata.obs[GROUP].astype(str)
    actual = {str(k):int(v) for k,v in labels.value_counts().items()}
    require(actual == EXPECTED_COUNTS, 'Cluster membership counts differ from reviewed 29-cluster version; STOP, do not map old IDs.')
    require(adata.obs[BATCH].astype(str).nunique()==247,'Expected 247 samples.')
    require('CPU_pipeline_json' in adata.uns, 'Cannot identify the parent CPU theta=2 pipeline.')
    parent_config = json.loads(str(adata.uns['CPU_pipeline_json']))
    hp = dict(parent_config['harmony'])
    require(float(hp['theta'])==2.0 and hp['device']=='cpu', 'Parent pipeline must be CPU Harmony theta=2.')
    require(sc.__version__=='1.11.5', 'Use the parent Scanpy 1.11.5 environment.')
    patch_dir = Path(patch_dir)
    patch, tests = patch_dir/'code/patched_function.py', patch_dir/'audit/TESTS.json'
    require(patch.is_file() and tests.is_file(), 'Missing audited HVG patch/tests.')
    require(json.loads(tests.read_text())['passed'] is True, 'HVG tests not recorded as passed.')
    old_patch_sha = parent_config.get('hvg_patch_sha256')
    if old_patch_sha:
        require(file_sha256(patch)==old_patch_sha, 'HVG patch changed since mother atlas construction.')
    # Exact membership guard: use the original theta=2 run's cell-to-cluster registry,
    # not merely identical cluster sizes. Never choose an older file by fuzzy search.
    source_run = Path(str(adata.uns.get('output_directory','')))
    source_map = source_run / 'cell_cluster_map.tsv.gz'
    require(source_map.is_file(), f'Missing original CPU cell-cluster registry: {source_map}')
    frozen_map = pd.read_csv(source_map, sep='\t', dtype={'cell_id':str, GROUP:str},
                             usecols=['cell_id',GROUP])
    require(len(frozen_map)==adata.n_obs and frozen_map.cell_id.is_unique,
            'Original cell-cluster registry has invalid membership.')
    require(np.array_equal(frozen_map.cell_id.to_numpy(), adata.obs_names.astype(str).to_numpy())
            and np.array_equal(frozen_map[GROUP].to_numpy(), labels.to_numpy()),
            'Current cells/cluster labels differ from the original CPU theta=2 registry.')
    del frozen_map
    validate_counts(adata.X)
    before = protected_hashes(adata)
    ri = adata.uns.get('raw_integrity',{})
    for field,new in [('sparse_X_sha256',before['X']),('cell_order_sha256',before['cells']),('gene_order_sha256',before['genes'])]:
        require(field in ri and str(ri[field])==new, f'Parent raw_integrity mismatch: {field}')
    source = Path(str(ri['source_path']))
    require(source.is_file(), f'Original immutable counts source is missing: {source}')
    require(file_sha256(source)==str(ri['source_sha256']), 'Original counts FILE hash changed.')
    work_parent = Path(work_parent or adata.uns.get('output_directory',ROOT))
    run = work_parent / ('L1_partial_and_local_' + stamp())
    run.mkdir(parents=True,exist_ok=False)
    for sub in ('audit','tables','figures','clusters'):
        (run/sub).mkdir()
    shutil.copy2(patch,run/'audit/patched_function.py')
    shutil.copy2(tests,run/'audit/HVG_TESTS.json')
    shutil.copy2(Path(__file__),run/'audit/l1_cpu_local_refinement.py')
    marker_panel = {str(k):list(v) for k,v in markers.items()}
    dump_json(run/'marker_panel.json',marker_panel)

    # Only annotation metadata are changed; preserve any existing L1 in a new column.
    if 'L1' in adata.obs:
        old_key = 'L1_before_partial_' + stamp()
        adata.obs[old_key] = adata.obs['L1'].copy()
    adata.obs['L1'] = pd.Categorical(labels.map(APPROVED).fillna('Unknown'))
    adata.obs['L1_parent_cluster'] = 'C' + labels
    adata.obs['L1_review_status'] = pd.Categorical(np.where(labels.isin(APPROVED),
        'APPROVED_PARENT',np.where(labels.eq('28'),'APPROVED_UNKNOWN','PENDING_LOCAL_REVIEW')))
    adata.obs['L1_review_note'] = np.where(labels.eq('10'),'Approved T cell; weaker marker detection',
                                         np.where(labels.eq('28'),'Retained Unknown; excluded from local refinement',''))
    adata.uns['active_L1_key'] = 'L1'
    adata.uns['CPU_L1_status'] = 'PARTIALLY_APPLIED_PENDING_LOCAL_REVIEW'
    adata.uns['L1_frozen'] = False
    adata.uns['L2_started'] = False
    adata.uns['L1_working_directory'] = str(run)
    adata.uns['L1_partial_approved_mapping_json'] = json.dumps(APPROVED,ensure_ascii=False)
    lock_sparse(adata.X)
    composition = adata.obs['L1'].value_counts()
    require({str(k):int(v) for k,v in composition.items()}==EXPECTED_L1,'Unexpected preliminary L1 counts.')
    composition.rename_axis('L1').reset_index(name='n_cells').assign(
        fraction=lambda d:d.n_cells/adata.n_obs).to_csv(run/'tables/preliminary_L1_counts.tsv',sep='\t',index=False)
    pd.DataFrame({'cell_id':adata.obs_names,'parent_cluster':labels.to_numpy(),
                  'L1':adata.obs['L1'].astype(str).to_numpy(),
                  'status':adata.obs['L1_review_status'].astype(str).to_numpy()}).to_csv(
                  run/'tables/preliminary_L1_cell_mapping.tsv.gz',sep='\t',index=False)
    basis = 'X_umap_cpu'
    require(basis in adata.obsm,'Expected parent X_umap_cpu.')
    parent_file = run/'CPU_theta2_L1_partial_full.h5ad'
    print('Saving full mother atlas with raw X, log1p layer, graphs, UMAP and preliminary L1 ...',flush=True)
    require(protected_hashes(adata)==before,'A protected matrix/coordinate changed before saving.')
    save_h5ad_atomic(adata,parent_file,before['X'])
    parent_file_sha = file_sha256(parent_file)
    save_umap(adata,basis,'L1',run/'figures/UMAP_preliminary_L1')

    # Approved parent references for expression-only comparison (NOT extra clustering cells).
    panel, missing = clean_panel(marker_panel,adata.var_names)
    genes = list(dict.fromkeys(g for gs in panel.values() for g in gs))
    ml,fr,_,_ = aggregate_markers(adata,adata.obs['L1'].astype(str),genes)
    ml=ml.drop(index='Unknown'); fr=fr.drop(index='Unknown')
    ml.index=['REF '+str(x) for x in ml.index]; fr.index=ml.index.copy()
    ml.to_csv(run/'tables/approved_reference_mean_log1p.tsv',sep='\t')
    fr.to_csv(run/'tables/approved_reference_detection_fraction.tsv',sep='\t')
    dump_json(run/'audit/protected_mother_hashes.json',before)
    dump_json(run/'audit/missing_markers.json',missing)

    config = {
        'seed':0, 'n_hvg':2000, 'n_pcs':30, 'pearson_theta':100,
        'pearson_clip':'sqrt(local n_cells per sample)', 'hvg_min_sample_cells':2,
        'normalize_target_sum':10000, 'scale_max_value':10,
        'harmony':hp, 'n_neighbors':30, 'metric':'euclidean',
        'leiden_resolutions':list(RESOLUTIONS),'baseline_resolution':0.5,
        'leiden_flavor':'igraph','leiden_n_iterations':100,
        'umap_min_dist':0.3,'umap_maxiter':500,
        'findmarkers':'scanpy.tl.rank_genes_groups wilcoxon on log1p, all genes detected in >=3 cells',
        'sample_metadata_use':'descriptive review only',
        'script_sha256':file_sha256(run/'audit/l1_cpu_local_refinement.py'),
        'patch_sha256':file_sha256(run/'audit/patched_function.py'),
    }
    tasks=[]
    for cid in LOCAL_TARGETS:
        name='C'+cid
        folder=run/'clusters'/name
        (folder/'tables').mkdir(parents=True)
        (folder/'figures').mkdir()
        ix=np.flatnonzero(labels.to_numpy()==cid)
        require(len(ix)==EXPECTED_COUNTS[cid],f'{name} member count mismatch.')
        obs=adata.obs.iloc[ix].copy()
        obs['parent_cell_id']=adata.obs_names[ix].astype(str)
        obs['parent_row']=ix.astype(np.int64)
        obs['parent_cluster']=cid
        obs['parent_cluster_uid']=name
        # Preserve gene metadata but label parent preprocessing columns as historical.
        var=adata.var.copy()
        for c in list(var.columns):
            if c.startswith(('highly_variable','residual','cpu_scale')) or c in ('means','variances','variances_norm'):
                var.rename(columns={c:'parent_'+c},inplace=True)
        C = ad.AnnData(X=adata.X[ix,:].copy(),obs=obs,var=var)
        C.obsm['X_umap_parent']=np.asarray(adata.obsm[basis])[ix,:].copy()
        rh=sparse_sha256(C.X)
        C.uns['object_name']=name
        C.uns['raw_integrity_local']={
            'raw_sha256':rh,'cell_order_sha256':index_sha256(C.obs_names),
            'gene_order_sha256':index_sha256(C.var_names),
            'parent_file':str(parent_file),'parent_file_sha256':parent_file_sha,
            'parent_raw_sha256':before['X'], 'source_raw_file':str(source),
        }
        rawfile=folder/f'{name}_raw_counts.h5ad'
        save_h5ad_atomic(C,rawfile,rh)
        task={'name':name,'raw_input':str(rawfile),'input_sha256':file_sha256(rawfile),
              'raw_matrix_sha256':rh,'n_cells':len(ix),'run_dir':str(run),
              'folder':str(folder),'focus':FOCUS[name],'config':config}
        spec=folder/'task.json'; dump_json(spec,task)
        tasks.append(str(spec))
        print(f'Extracted {name}: {len(ix):,} cells x {C.n_vars:,} genes',flush=True)
        del C; gc.collect()
    require(sum(EXPECTED_COUNTS[c] for c in LOCAL_TARGETS)==350118,'Incorrect target size.')
    require(protected_hashes(adata)==before,'Mother matrices changed during extraction.')
    runmeta={'status':'READY_FOR_LOCAL_CPU_ANALYSIS','mother_file':str(parent_file),
             'mother_file_sha256':parent_file_sha,'tasks':tasks,
             'unknown_total':350188,'local_cells':350118,'C28_retained_unknown':70,
             'approved_cells':824927,'all_cells':1175115,
             'approved_mapping':APPROVED,'raw_integrity':dict(ri),
             'config':config,'script':str(run/'audit/l1_cpu_local_refinement.py')}
    dump_json(run/'RUN.json',runmeta)
    print(f'Preliminary L1 saved: {parent_file}\nRun folder: {run}',flush=True)
    return run


def _cpu_setup(threads):
    import matplotlib
    matplotlib.use('Agg')
    import scanpy as sc
    import torch
    import numba
    import omicverse as ov
    ov.settings.cpu_init()
    require(not torch.cuda.is_available() and not torch.cuda.is_initialized(), 'CPU-only worker unexpectedly sees/uses CUDA.')
    require(sc.__version__=='1.11.5','HVG patch is pinned to Scanpy 1.11.5.')
    numba.set_num_threads(threads); torch.set_num_threads(threads)
    sc.settings.n_jobs=threads; sc.settings.autoshow=False; sc.settings.autosave=False
    np.random.seed(0)
    return ov,torch


def local_hvg(work, task, patch, threads):
    import scanpy as sc
    from threadpoolctl import threadpool_limits
    cfg=task['config']; sizes=work.obs[BATCH].astype(str).value_counts()
    # A singleton cannot supply a within-sample variance estimate. Its cell is
    # retained in all other steps, never deleted or reassigned to a fictitious batch.
    eligibility=pd.DataFrame({'sample_uid':sizes.index,'n_local_cells':sizes.to_numpy()})
    eligibility['used_for_HVG_fit']=eligibility['n_local_cells']>=2
    eligibility.to_csv(Path(task['folder'])/'tables/HVG_sample_eligibility.tsv',sep='\t',index=False)
    eligible=sizes[sizes>=2].index
    require(len(eligible)>0,'No sample with >=2 local cells; HVG cannot be estimated.')
    mask=work.obs[BATCH].astype(str).isin(eligible).to_numpy()
    fit=work if mask.all() else work[mask,:].copy()
    fit.obs[BATCH]=fit.obs[BATCH].astype(str)
    fs=fit.obs[BATCH].value_counts()
    clips=[]
    def observe(batch, clip, variances):
        require(np.isfinite(variances).all(),f'Invalid residual variances: {batch}')
        clips.append({'sample_uid':str(batch),'n_local_cells':int(fs[str(batch)]),'clip':float(clip)})
    mod=importlib.import_module('scanpy.experimental.pp._highly_variable_genes')
    ns=dict(mod.__dict__); ns.update(_batch_order=lambda batches:batches,_observe=observe)
    exec(compile(patch.read_text(),str(patch),'exec'),ns)
    with threadpool_limits(limits=threads):
        hvg=ns['_highly_variable_pearson_residuals'](fit,batch_key=BATCH,theta=cfg['pearson_theta'],
            clip=None,n_top_genes=cfg['n_hvg'],chunksize=1000,check_values=False,inplace=False)
    require(hvg.index.equals(work.var_names),'Local HVG gene order mismatch.')
    require(len(clips)==len(fs) and {x['sample_uid'] for x in clips}==set(fs.index),'HVG fitting sample coverage mismatch.')
    require(all(np.isclose(x['clip'],np.sqrt(x['n_local_cells']),rtol=0,atol=1e-12) for x in clips),'Incorrect local per-sample clip.')
    hv=hvg.highly_variable.to_numpy(dtype=bool)
    require(hv.sum()==cfg['n_hvg'],'Unexpected number of local HVGs.')
    hvg.to_csv(Path(task['folder'])/'tables/HVG_statistics.tsv.gz',sep='\t')
    pd.DataFrame(clips).to_csv(Path(task['folder'])/'tables/sample_clips.tsv',sep='\t',index=False)
    return hvg


def create_local_review(C, task):
    import anndata as ad
    import scanpy as sc
    folder=Path(task['folder']); run=Path(task['run_dir']); name=task['name']
    markers=json.loads((run/'marker_panel.json').read_text())
    panel,missing=clean_panel(markers,C.var_names)
    dump_json(folder/'tables/missing_markers.json',missing)
    genes=list(dict.fromkeys(g for gs in panel.values() for g in gs))
    baseline_mean=baseline_frac=None
    for r in task['config']['leiden_resolutions']:
        tag='r'+str(r).replace('.','p'); key='local_leiden_'+tag
        rd=folder/'tables'/tag; rd.mkdir(exist_ok=True)
        mean,frac,rawmean,n=aggregate_markers(C,C.obs[key].astype(str),genes)
        for fn,df in [('marker_mean_log1p.tsv',mean),('marker_detection_fraction.tsv',frac),('marker_mean_raw_counts.tsv',rawmean)]:
            df.to_csv(rd/fn,sep='\t')
        long=pd.concat([df.rename_axis(index='cluster_uid',columns='gene').stack().rename(v)
                        for df,v in [(mean,'mean_log1p_CP10k'),(frac,'detection_fraction'),(rawmean,'mean_raw_counts')]],axis=1).reset_index()
        long.to_csv(rd/'marker_evidence_long.tsv',sep='\t',index=False)
        n.rename_axis('cluster_uid').reset_index().to_csv(rd/'local_cluster_counts.tsv',sep='\t',index=False)
        save_dotplot(panel,mean,frac,folder/'figures'/f'{name}_{tag}_dotplot_relative',True)
        save_umap(C,'X_umap_local',key,folder/'figures'/f'{name}_{tag}_UMAP')
        if key==BASELINE_KEY: baseline_mean,baseline_frac=mean,frac
    save_dotplot(panel,baseline_mean,baseline_frac,folder/'figures'/f'{name}_baseline_dotplot_mean_log1p')
    mural={k:panel[k] for k in MURAL if k in panel}
    save_dotplot(mural,baseline_mean,baseline_frac,folder/'figures'/f'{name}_baseline_mural_dotplot_relative',True)
    refmean=pd.read_csv(run/'tables/approved_reference_mean_log1p.tsv',sep='\t',index_col=0)
    reffrac=pd.read_csv(run/'tables/approved_reference_detection_fraction.tsv',sep='\t',index_col=0)
    rm=pd.concat([baseline_mean,refmean.loc[:,genes]])
    rf=pd.concat([baseline_frac,reffrac.loc[:,genes]])
    save_dotplot(panel,rm,rf,folder/'figures'/f'{name}_baseline_with_approved_references',True)
    rm.to_csv(folder/'tables/marker_means_with_approved_references.tsv',sep='\t')
    for key in ('study','tissue'):
        if key in C.obs: save_umap(C,'X_umap_local',key,folder/'figures'/f'{name}_UMAP_{key}')
    for r in (0.3,0.8):
        key='local_leiden_r'+str(r).replace('.','p')
        pd.crosstab(C.obs[BASELINE_KEY],C.obs[key],dropna=False).to_csv(
            folder/'tables'/f'baseline_vs_r{str(r).replace(".","p")}_members.tsv',sep='\t')

    # Discovery compares local clusters on ALL expressed genes, not only the L1 panel/HVGs.
    labels=C.obs[BASELINE_KEY].astype(str); counts=labels.value_counts()
    tested_groups=counts[counts>=2].index.tolist()
    keepgene=np.asarray(C.layers['log1p'].getnnz(axis=0)).ravel()>=3
    rankkey='rank_genes_local_baseline'
    if len(counts)>1 and len(tested_groups)>0 and keepgene.sum()>0:
        sc.tl.rank_genes_groups(C,groupby=BASELINE_KEY,groups=tested_groups,reference='rest',
            layer='log1p',use_raw=False,mask_var=keepgene,method='wilcoxon',tie_correct=True,
            pts=True,n_genes=min(100,int(keepgene.sum())),key_added=rankkey,corr_method='benjamini-hochberg')
        rank=sc.get.rank_genes_groups_df(C,group=None,key=rankkey)
        rank.to_csv(folder/'tables/FindMarkers_top100.tsv',sep='\t',index=False)
        positive=rank[rank.logfoldchanges.gt(0) & rank.pvals_adj.lt(0.05)]
        top20=positive.groupby('group',observed=True,sort=False).head(20)
        top20.to_csv(folder/'tables/FindMarkers_top20_positive.tsv',sep='\t',index=False)
        top_map=top20.groupby('group',observed=True)['names'].agg(lambda x:';'.join(x)).to_dict()
    else:
        top_map={}
        (folder/'tables/FindMarkers_not_tested.txt').write_text(
            'No valid within-parent between-subcluster comparison; no forced extra clusters.\n')
    # Per-sample counts and marker evidence: sources are descriptive, not label veto rules.
    rows=[]
    for uid,ix in C.obs.groupby(BASELINE_KEY,observed=True).indices.items():
        block=C.obs.iloc[ix]
        row={'parent_cluster':name,'cluster_uid':str(uid),'n_cells':len(ix),
             'n_samples':block[BATCH].nunique(),'top_markers':top_map.get(str(uid),''),
             'proposed_L1':'','decision':'PENDING','comment':''}
        for k in ('study','tissue','condition','disease_context','donor_uid_v2'):
            if k in block:
                vals=block[k].astype('string').dropna()
                vals=vals[~vals.str.lower().isin(['nr','na','nan','none','unknown','not reported',''])]
                vc=vals.value_counts()
                row['n_'+k]=len(vc)
                row['top_'+k]=str(vc.index[0]) if len(vc) else ''
                row['top_'+k+'_fraction']=float(vc.iloc[0]/len(block)) if len(vc) else np.nan
        row['median_total_raw_counts']=float(np.median(block['local_total_raw_counts']))
        row['median_detected_genes']=float(np.median(block['local_detected_genes']))
        rows.append(row)
    review=pd.DataFrame(rows)
    review.to_csv(folder/'tables/local_cluster_evidence_summary.tsv',sep='\t',index=False)
    template=folder/'tables/local_L1_manual_decisions.tsv'
    if not template.exists(): review.to_csv(template,sep='\t',index=False)
    for k in (BATCH,'study','tissue','condition','disease_context','donor_uid_v2'):
        if k in C.obs:
            pd.crosstab(C.obs[BASELINE_KEY],C.obs[k].astype('string').fillna('Missing')).to_csv(
                folder/'tables'/f'cluster_by_{k}.tsv',sep='\t')
    # Marker means by local cluster x sample, for source-specific verification.
    combo=labels.to_numpy()+' || '+C.obs[BATCH].astype(str).to_numpy()
    mean_s,frac_s,_,n_s=aggregate_markers(C,combo,genes)
    mean_s.to_csv(folder/'tables/cluster_sample_marker_mean_log1p.tsv.gz',sep='\t')
    frac_s.to_csv(folder/'tables/cluster_sample_marker_detection_fraction.tsv.gz',sep='\t')
    n_s.to_csv(folder/'tables/cluster_sample_cell_counts.tsv',sep='\t')
    mapping=C.obs[['parent_cell_id','parent_row','parent_cluster','parent_cluster_uid','L1',
                   'local_leiden_r0p3','local_leiden_r0p5','local_leiden_r0p8']].copy()
    mapping.insert(0,'cell_id',C.obs_names)
    mapping.to_csv(folder/'tables/cell_local_cluster_mapping.tsv.gz',sep='\t',index=False)
    msg=(f'# {name} local L1 review\n\nFocus: {task["focus"]}\n'
         f'Cells retained: {C.n_obs:,}; genes: {C.n_vars:,}. Baseline: r0.5 (review anchor, not a proven optimum).\n'
         'All original raw counts are in X; full-gene CP10k-log1p is in layers[log1p].\n'
         'No local label was assigned automatically; L1 remains Unknown.\n'
         'Subset graph/coordinates were recomputed, not inherited from the global graph.\n'
         'REF rows are approved mother-atlas expression summaries, NOT extra integration cells.\n'
         'Different resolutions have different memberships; equal suffix IDs do NOT imply identical groups.\n'
         'Singleton samples are omitted ONLY from within-sample HVG variance estimation; their cells remain in PCA/Harmony/plots.\n'
         'FindMarkers uses within-parent one-versus-rest Wilcoxon; p-values are descriptive annotation evidence, not donor-level disease inference.\n'
         'No Scrublet/doublet filtering, no automatic rejection by source dominance, no target recovery fraction.\n'
         'No L1 mother writeback and no L2. A mixed branch may remain unresolved after reanalysis.\n')
    (folder/'REVIEW_NOTES.md').write_text(msg,encoding='utf-8')
    document=['<!doctype html><meta charset="utf-8"><style>body{font:15px sans-serif;margin:24px}img{max-width:100%}pre{white-space:pre-wrap}td,th{padding:4px;border:1px solid #ccc}</style>',
              '<pre>'+html.escape(msg)+'</pre>',review.to_html(index=False,escape=True)]
    for path in sorted((folder/'figures').glob('*.png')):
        document += ['<h2>'+html.escape(path.stem)+'</h2>',
                     '<img src="data:image/png;base64,'+base64.b64encode(path.read_bytes()).decode()+'">']
    (folder/f'{name}_OFFLINE_REVIEW.html').write_text('\n'.join(document),encoding='utf-8')
    return review


def worker(spec_path):
    import anndata as ad
    import scanpy as sc
    from threadpoolctl import threadpool_limits, threadpool_info
    from pynndescent import PyNNDescentTransformer
    task=json.loads(Path(spec_path).read_text()); cfg=task['config']
    folder=Path(task['folder']); name=task['name']; run=Path(task['run_dir'])
    threads=int(os.environ['L1_CPU_THREADS'])
    ov,torch=_cpu_setup(threads)
    started=time.time()
    def status(s):
        dump_json(folder/'STATUS.json',{'status':s,'time':datetime.now().isoformat(),'name':name})
        print(f'[{name}] {s}',flush=True)
    status('LOADING_LOCAL_RAW_COUNTS')
    require(file_sha256(task['raw_input'])==task['input_sha256'],'Local raw input file changed.')
    patch=run/'audit/patched_function.py'
    require(file_sha256(patch)==cfg['patch_sha256'],'HVG patch snapshot changed.')
    C=sc.read_h5ad(task['raw_input']); validate_counts(C.X)
    raw_sha=sparse_sha256(C.X)
    require(raw_sha==task['raw_matrix_sha256'],'Local raw matrix hash mismatch.')
    require(C.n_obs==task['n_cells'] and C.obs['L1'].astype(str).eq('Unknown').all(),'Incorrect local membership.')
    lock_sparse(C.X)
    work=ad.AnnData(X=C.X.tocsr(copy=True),obs=C.obs[[BATCH]].copy(),var=pd.DataFrame(index=C.var_names.copy()))
    work.X.sum_duplicates(); work.X.eliminate_zeros(); work.X.sort_indices()
    work.obs[BATCH]=work.obs[BATCH].astype(str)
    totals=np.asarray(work.X.sum(axis=1,dtype=np.float64)).ravel()
    require(np.isfinite(totals).all() and (totals>0).all(),'Zero/invalid library size; stop, no cell deletion.')
    C.obs['local_total_raw_counts']=totals
    C.obs['local_detected_genes']=work.X.getnnz(axis=1)
    status('LOCAL_HVG_ON_RAW_COUNTS')
    hvg=local_hvg(work,task,patch,threads)
    for col in hvg: C.var[col]=hvg[col]
    hv=hvg.highly_variable.to_numpy(dtype=bool)
    status('FULL_GENE_NORMALIZATION_COPY')
    with threadpool_limits(limits=threads):
        sc.pp.normalize_total(work,target_sum=cfg['normalize_target_sum'],exclude_highly_expressed=False)
        sc.pp.log1p(work)
    C.layers['log1p']=work.X
    C.uns['log1p']={'base':None}
    del work;gc.collect()
    status('LOCAL_SCALE_PCA')
    p=ad.AnnData(X=C.layers['log1p'][:,hv].copy(),obs=pd.DataFrame(index=C.obs_names.copy()),
                 var=pd.DataFrame(index=C.var_names[hv].copy()))
    with threadpool_limits(limits=threads):
        sc.pp.scale(p,zero_center=True,max_value=cfg['scale_max_value'])
        sc.pp.pca(p,n_comps=cfg['n_pcs'],svd_solver='covariance_eigh',zero_center=True,
                  mask_var=None,random_state=cfg['seed'],dtype='float32')
    C.obsm['X_pca_local']=np.asarray(p.obsm['X_pca'],dtype=np.float32).copy()
    loadings=np.zeros((C.n_vars,cfg['n_pcs']),dtype=np.float32);loadings[hv]=p.varm['PCs']
    C.varm['PCs_local']=loadings;C.uns['pca_local']=copy.deepcopy(p.uns['pca'])
    for key in ('mean','std'):
        v=np.full(C.n_vars,np.nan);v[hv]=p.var[key].to_numpy();C.var['local_scale_'+key]=v
    del p,loadings;gc.collect()
    status('LOCAL_HARMONY_CPU')
    hm=importlib.import_module('omicverse.external.harmony.harmony')
    require('device' in inspect.signature(hm.run_harmony).parameters,'Harmony must accept explicit CPU device.')
    hp=dict(cfg['harmony']);hp['device']='cpu';hp['theta']=2.0
    if C.obs[BATCH].nunique()>1:
        with threadpool_limits(limits=threads):
            result=hm.run_harmony(C.obsm['X_pca_local'].copy(),
                C.obs[[BATCH]].assign(**{BATCH:C.obs[BATCH].astype(str)}),BATCH,**hp)
        z=result.result()
        if torch.is_tensor(z):z=z.detach().cpu().numpy()
        z=np.asarray(z,dtype=np.float32)
        require(z.shape==(C.n_obs,cfg['n_pcs']) and np.isfinite(z).all(),'Invalid Harmony embedding.')
        C.obsm['X_harmony_local']=z.copy()
        cfg['actual_harmony_class']=type(result).__module__+'.'+type(result).__name__
        try:shutil.copy2(inspect.getfile(type(result)),folder/'harmony_actual_source.py')
        except (OSError,TypeError):pass
        obj=getattr(result,'objective_harmony',None)
        if obj is not None:
            if torch.is_tensor(obj):obj=obj.detach().cpu().numpy()
            pd.DataFrame({'record':np.arange(len(obj)),'objective':np.asarray(obj).ravel()}).to_csv(
                folder/'tables/harmony_objective_history.tsv',sep='\t',index=False)
        del result,z;gc.collect()
    else:
        cfg['actual_harmony_class']='SKIPPED_SINGLE_SAMPLE'
        C.obsm['X_harmony_local']=C.obsm['X_pca_local'].copy()
    status('LOCAL_NEIGHBORS_LEIDEN_UMAP_CPU')
    transformer=PyNNDescentTransformer(n_neighbors=cfg['n_neighbors'],metric=cfg['metric'],
                                      random_state=cfg['seed'],n_jobs=threads)
    with threadpool_limits(limits=1,user_api='blas'):
        sc.pp.neighbors(C,n_neighbors=cfg['n_neighbors'],n_pcs=cfg['n_pcs'],use_rep='X_harmony_local',
                        metric=cfg['metric'],method='umap',transformer=transformer,
                        random_state=cfg['seed'],key_added='neighbors_local')
    for r in cfg['leiden_resolutions']:
        tag='r'+str(r).replace('.','p');key='local_leiden_'+tag
        sc.tl.leiden(C,resolution=r,random_state=cfg['seed'],n_iterations=cfg['leiden_n_iterations'],
                     flavor='igraph',directed=False,use_weights=True,neighbors_key='neighbors_local',key_added=key)
        codes=C.obs[key].astype(str)
        order=sorted(codes.unique(),key=int)
        C.obs[key]=pd.Categorical(name+'_'+tag+'_'+codes,
                                  categories=[name+'_'+tag+'_'+x for x in order],ordered=True)
    with threadpool_limits(limits=threads):
        sc.tl.umap(C,min_dist=cfg['umap_min_dist'],spread=1.0,maxiter=cfg['umap_maxiter'],
                   init_pos='spectral',n_components=2,alpha=1.0,negative_sample_rate=5,
                   random_state=cfg['seed'],neighbors_key='neighbors_local')
    C.obsm['X_umap_local']=C.obsm['X_umap'].copy()
    status('MARKER_DOTPLOTS_AND_FINDMARKERS')
    with threadpool_limits(limits=threads):
        review=create_local_review(C,task)
    cfg['threads']=threads;cfg['threadpools']=threadpool_info()
    cfg['versions']={}
    for pkg in ('omicverse','scanpy','anndata','numpy','scipy','pandas','scikit-learn',
                'pynndescent','numba','igraph','umap-learn','torch'):
        try:cfg['versions'][pkg]=metadata.version(pkg)
        except metadata.PackageNotFoundError:cfg['versions'][pkg]='distribution metadata unavailable'

    C.uns['local_config_json']=json.dumps(cfg,ensure_ascii=False,default=str)
    C.uns['active_L1_key']='L1'
    C.uns['local_L1_status']='PAUSED_FOR_LOCAL_L1_REVIEW'
    C.uns['formal_local_L1_applied']=False
    C.uns['expression_layout']={'X':'unchanged full-gene raw counts',
                                'log1p':'full-gene CP10k + log1p in layers[log1p]',
                                'parent_umap':'display reference only; local graph recomputed'}
    require(sparse_sha256(C.X)==raw_sha,'Local raw counts changed in memory.')
    require(C.obs['L1'].astype(str).eq('Unknown').all(),'Local worker altered a formal L1 label.')
    final=folder/f'{name}_analysis.h5ad'
    # Only worker-generated non-final review objects may be archived on explicit retry.
    if final.exists(): final.rename(folder/f'{name}_analysis_previous_{stamp()}.h5ad')
    save_h5ad_atomic(C,final,raw_sha)
    outmanifest={'status':'PAUSED_FOR_LOCAL_L1_REVIEW','name':name,'n_cells':C.n_obs,
                 'analysis_file':str(final),'analysis_sha256':file_sha256(final),
                 'raw_matrix_sha256':raw_sha,'elapsed_seconds':time.time()-started,
                 'task_sha256':file_sha256(spec_path),'config':cfg}
    dump_json(folder/'COMPLETED.json',outmanifest)
    status('PAUSED_FOR_LOCAL_L1_REVIEW')


def build_index(run_dir):
    run=Path(run_dir);meta=json.loads((run/'RUN.json').read_text());rows=[]
    for spec in meta['tasks']:
        t=json.loads(Path(spec).read_text());folder=Path(t['folder'])
        status=json.loads((folder/'STATUS.json').read_text()) if (folder/'STATUS.json').exists() else {'status':'NOT_STARTED'}
        rows.append({'parent':t['name'],'n_cells':t['n_cells'],'status':status['status'],
                     'focus':t['focus'],'review':str(folder/f'{t["name"]}_OFFLINE_REVIEW.html')})
    df=pd.DataFrame(rows);df.to_csv(run/'tables/local_batch_status.tsv',sep='\t',index=False)
    text=['<!doctype html><meta charset="utf-8"><h1>CPU local L1 review</h1>',
          '<p>Mother L1 is preliminary: 824,927 approved cells; 350,188 Unknown, including C28 (70 cells) excluded from refinement. Local labels have NOT been written back.</p>',
          df.drop(columns='review').to_html(index=False)]
    for row in rows:
        rel=Path(row['review']).relative_to(run)
        text.append('<p><a href="'+html.escape(str(rel))+'">'+html.escape(row['parent']+' - '+row['focus'])+'</a></p>')
    (run/'GLOBAL_LOCAL_L1_REVIEW.html').write_text('\n'.join(text),encoding='utf-8')
    reviewzip=run/'Local_L1_review_materials.zip'
    with zipfile.ZipFile(reviewzip,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=3) as z:
        for p in run.rglob('*'):
            if not p.is_file() or p==reviewzip:continue
            if p.suffix in ('.h5ad','.npy','.npz','.zip','.pdf','.svg'):continue
            if '.incomplete_' in p.name or '_previous_' in p.name:continue
            if p.suffix=='.png':continue # embedded in the branch offline HTML
            z.write(p,arcname=str(p.relative_to(run)))
    return df


def run_local_batch(run_dir, max_workers=2, total_threads=40):
    """Blocking notebook call; launches at most two CPU-only jobs, 20 threads each.
    Completed branches are reused only after task/output hash checks.
    Re-run the SAME call after a failure; successful branches will not be recomputed.
    """
    run=Path(run_dir).resolve();meta=json.loads((run/'RUN.json').read_text())
    available=len(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else (os.cpu_count() or 1)
    require(max_workers>=1 and total_threads>=max_workers,'Invalid worker/thread configuration.')
    require(total_threads<=available,f'Only {available} CPU logical cores are available to this process.')
    per=max(1,total_threads//max_workers)
    script=Path(meta['script'])
    require(file_sha256(script)==meta['config']['script_sha256'],'Saved worker code changed.')
    require(file_sha256(meta['mother_file'])==meta['mother_file_sha256'],'Preliminary mother file changed.')
    waiting=[]
    for spec in meta['tasks']:
        t=json.loads(Path(spec).read_text());done=Path(t['folder'])/'COMPLETED.json'
        if done.exists():
            d=json.loads(done.read_text())
            require(d['task_sha256']==file_sha256(spec),'Completed task configuration changed.')
            require(file_sha256(d['analysis_file'])==d['analysis_sha256'],'Completed local output changed.')
            print('Already complete, kept:',t['name'],flush=True)
        else:waiting.append((spec,t))
    # Largest jobs first to avoid one long final tail; names/biological membership unchanged.
    waiting.sort(key=lambda item:item[1]['n_cells'],reverse=True)
    active=[];failures=[]
    env=os.environ.copy();env['CUDA_VISIBLE_DEVICES']='';env['MPLBACKEND']='Agg'
    env['PYTHONHASHSEED']='0';env['PYTHONUNBUFFERED']='1';env['L1_CPU_THREADS']=str(per)
    for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'):
        env[k]=str(per)
    meta['status']='RUNNING_LOCAL_CPU_ANALYSES';meta['parallel']={'jobs':max_workers,'threads_per_job':per}
    dump_json(run/'RUN.json',meta)
    try:
        while waiting or active:
            while waiting and len(active)<max_workers:
                spec,t=waiting.pop(0)
                log=Path(t['folder'])/('run_'+stamp()+'.log')
                f=log.open('w',encoding='utf-8')
                try:
                    p=subprocess.Popen([sys.executable,'-u',str(script),'--worker',spec],
                                       stdout=f,stderr=subprocess.STDOUT,env=env,cwd=str(run),start_new_session=True)
                except Exception:
                    f.close();raise
                active.append((p,f,t,str(log)))
                print(f'Started {t["name"]}: {per} threads; log: {log}',flush=True)
            for item in active[:]:
                p,f,t,log=item;rc=p.poll()
                if rc is None:continue
                f.close();active.remove(item)
                if rc!=0:
                    failures.append(t['name'])
                    dump_json(Path(t['folder'])/'STATUS.json',{'status':'FAILED','exitcode':rc,'log':log})
                    print(f'FAILED {t["name"]}; see {log}',flush=True)
                else: print(f'Completed {t["name"]}',flush=True)
            if active:time.sleep(3)
    except BaseException:
        import signal
        for p,f,t,log in active:
            try:
                os.killpg(p.pid,signal.SIGTERM);p.wait(timeout=15)
            except (ProcessLookupError,subprocess.TimeoutExpired):
                if p.poll() is None:os.killpg(p.pid,signal.SIGKILL)
            finally:
                f.close()
                dump_json(Path(t['folder'])/'STATUS.json',{'status':'INTERRUPTED','log':log})
        meta['status']='INTERRUPTED';dump_json(run/'RUN.json',meta);build_index(run)
        raise
    meta['status']='PAUSED_FOR_LOCAL_L1_REVIEW' if not failures else 'PARTIAL_FAILURE'
    meta['failed_branches']=failures
    require(file_sha256(meta['mother_file'])==meta['mother_file_sha256'],'Mother file was changed during local analysis.')
    dump_json(run/'RUN.json',meta)
    summary=build_index(run)
    print(summary.to_string(index=False),flush=True)
    print('Offline index:',run/'GLOBAL_LOCAL_L1_REVIEW.html',flush=True)
    print('Review ZIP (no expression objects):',run/'Local_L1_review_materials.zip',flush=True)
    if failures: raise RuntimeError(f'Failed branches: {failures}. Fix the reported error; rerun run_local_batch on this run directory.')
    return summary


def load_local(run_dir, cluster_name, backed=None):
    import scanpy as sc
    name=str(cluster_name)
    if not name.startswith('C'):name='C'+name
    p=Path(run_dir)/'clusters'/name/f'{name}_analysis.h5ad'
    if not p.is_file():raise FileNotFoundError(p)
    return sc.read_h5ad(p,backed=backed)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--worker',required=True)
    args=parser.parse_args()
    try:worker(args.worker)
    except BaseException:
        traceback.print_exc()
        sys.exit(1)
