# NOTEBOOK CELL 0, execution_count=1
import scanpy as sc
import omicverse as ov
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import cellhint
import pertpy as pt
ov.plot_set()
sc.settings.set_figure_params(dpi=600, facecolor="white")
# NOTEBOOK CELL 1, execution_count=2
adata = sc.read_h5ad("/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908/annotation/full_atlas_CPU_Harmony_v1/data/full_atlas_CPU_review.h5ad")
# NOTEBOOK CELL 2, execution_count=None
###前置处理参数备份
import importlib
import scanpy as sc
from pynndescent import PyNNDescentTransformer

# Harmony：直接调用 OmicVerse 内置实现
harmony_module = importlib.import_module(
    "omicverse.external.harmony.harmony"
)

result = harmony_module.run_harmony(
    adata.obsm["X_pca_cpu"],
    adata.obs.assign(sample_uid=adata.obs["sample_uid"].astype(str)),
    "sample_uid",
    theta=1.0,
    lamb=None,
    sigma=0.1,
    nclust=None,
    tau=0,
    block_size=0.05,
    max_iter_harmony=100,
    max_iter_kmeans=20,
    epsilon_cluster=1e-5,
    epsilon_harmony=1e-4,
    alpha=0.2,
    random_state=0,
    device="cpu",
    verbose=True,
)
adata.obsm["X_harmony_cpu"] = result.result()

# CPU neighbors：48 个并行任务
transformer = PyNNDescentTransformer(
    n_neighbors=30,
    metric="euclidean",
    random_state=0,
    n_jobs=48,
)

sc.pp.neighbors(
    adata,
    n_neighbors=30,
    use_rep="X_harmony_cpu",
    metric="euclidean",
    transformer=transformer,
    random_state=0,
    key_added="neighbors_cpu_n30",
)

# CPU Leiden
sc.tl.leiden(
    adata,
    resolution=0.5,
    random_state=0,
    n_iterations=100,
    flavor="igraph",
    directed=False,
    use_weights=True,
    neighbors_key="neighbors_cpu_n30",
    key_added="leiden_cpu_r0p5",
)

# CPU UMAP
sc.tl.umap(
    adata,
    min_dist=0.3,
    spread=1.0,
    maxiter=500,
    init_pos="spectral",
    n_components=2,
    alpha=1.0,
    negative_sample_rate=5,
    random_state=0,
    neighbors_key="neighbors_cpu_n30",
)
adata.obsm["X_umap_cpu"] = adata.obsm.pop("X_umap")
# NOTEBOOK CELL 3, execution_count=14
adata
# NOTEBOOK CELL 4, execution_count=15
ov.settings.cpu_init()

ov.utils.embedding(
    adata,
    basis="X_umap_cpu",
    color=["study"],
    show=True,
    frameon="small",
)
# NOTEBOOK CELL 5, execution_count=3
adata = sc.read_h5ad(
    "/home/anzhen_vhd/scRNAseq/codex/"
    "panvascular_atlas_rebuild_v2_20260907/"
    "06_reanalysis_takeover_20260908/"
    "CPU_raw_counts_with_clinical_v1/"
    "panvascular_all_cells_raw_counts_with_clinical_v1.h5ad")
adata
# NOTEBOOK CELL 6, execution_count=3
import os
ov.settings.cpu_init()
os.environ["OMP_NUM_THREADS"] = "40"
os.environ["OPENBLAS_NUM_THREADS"] = "40"
os.environ["MKL_NUM_THREADS"] = "40"
os.environ["NUMEXPR_NUM_THREADS"] = "40"
os.environ["NUMBA_NUM_THREADS"] = "40"
ov.settings.cpu_init()

adata.obsm["X_pca"] = (
    adata.obsm["X_pca_cpu"]
    .astype("float32")
)

ov.single.batch_correction(
    adata,
    batch_key="sample_uid",
    max_iter_harmony=100,
    use_rep="X_pca",
    device="cpu",
    theta=2.0,
    methods="harmony",
    n_pcs=30
)

ov.pp.neighbors(
    adata,
    n_neighbors=30,
    n_pcs=30,
    use_rep="X_pca_harmony"
)

ov.pp.umap(adata)
# NOTEBOOK CELL 7, execution_count=None
ov.settings.cpu_gpu_mixed_init()
#ov.settings.gpu_init()
ov.utils.cluster(adata, key_added="leiden_res1", resolution=1)
#ov.utils.cluster(adata, key_added="leiden_res2", resolution=2)
ov.utils.embedding(adata,
    basis='X_umap_harmony',
    color=["leiden_res1"],
    title=['leiden_res1'],
    palette=ov.palette()[12:],
    show=False,
    frameon='small',
)
# 查看所有 clusters 的细胞数分布
print("\nCell counts in all clusters:")
print(adata.obs['leiden_res1'].value_counts())
# NOTEBOOK CELL 8, execution_count=None
# -*- coding: utf-8 -*-
"""CPU atlas rebuild + L1 marker review.

Run in a NEW Python/Jupyter kernel in the existing omicverse environment.
Final AnnData.X = unchanged ALL-gene counts; layers['log1p'] = ALL-gene
CP10k + natural log1p. Scaling happens ONLY in a temporary HVG object.
Neither expression nor labels are inferred from tissue/disease metadata.
"""
# %% 1. CPU configuration: must precede scientific-library imports
import os
N_THREADS = 40
os.environ['CUDA_VISIBLE_DEVICES'] = ''
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = str(N_THREADS)

import copy
import gc
import hashlib
import importlib
import importlib.metadata as metadata
import inspect
import json
import shutil
import stat
from datetime import datetime
from pathlib import Path

import anndata as ad
import h5py
import matplotlib.pyplot as plt
import numba
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from matplotlib.collections import PathCollection
from packaging.version import Version
from pynndescent import PyNNDescentTransformer
from threadpoolctl import threadpool_info, threadpool_limits

ROOT = Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/'
            '06_reanalysis_takeover_20260908')
COUNTS_FILE = ROOT / 'CPU_raw_counts_with_clinical_v1/panvascular_all_cells_raw_counts_with_clinical_v1.h5ad'
PATCH_DIR = ROOT / 'annotation/HVG_clip_review_20260910'
BATCH_KEY = 'sample_uid'
GROUP = 'leiden_res1'                  # matches subsequent manual-review code
NEIGHBORS_KEY = 'neighbors_cpu_n30'
N_HVG, N_PCS, SEED = 2000, 30, 0
THETA = 2.0                           # Harmony theta, NOT Pearson-residual theta
EXPECTED_SHAPE = (1_175_115, 37_336)
EXPECTED_SAMPLES = 247
LOCK_SOURCE_FILE = True              # remove write bits; never rewrite source bytes

markers = {
    'Endothelial cell': ['PECAM1','VWF','KDR','CLDN5','EMCN','RGCC','CA4','ACKR1','PLVAP','SPARCL1','EDN1','ECSCR','FABP4','EFNB2','GJA5','SEMA3G','ESM1'],
    'Fibroblast': ['DCN','LUM','COL1A1','COL1A2','COL3A1','COL5A1','COL14A1','PDGFRA','FBLN1','THY1','C3','C7'],
    'Fibromyocyte': ['FN1','AEBP1','POSTN','COL3A1','LUM','TNFRSF11B','ACTA2','TAGLN','MYL9','TCF21'],
    'Smooth muscle cell': ['ACTA2','MYH11','TAGLN','CNN1','MYL9','TPM2','CALD1','LMOD1','MYLK'],
    'Pericyte': ['RGS5','CSPG4','MCAM','PDGFRB','ABCC9','KCNJ8','COX4I2','NDUFA4L2','NOTCH3','DES'],
    'Macrophage': ['LST1','TYROBP','FCER1G','CTSS','CD68','AIF1','ITGAM','CSF1R','C1QA','C1QB','C1QC','LGALS3','APOE','MS4A7','MSR1','TREM2','GPNMB','SPP1'],
    'Monocyte': ['FCN1','S100A8','S100A9','S100A12','VCAN','CCR2','LYZ','CTSS','CTSD','LGALS3'],
    'Neutrophil': ['FCGR3B','CSF3R','CXCR2','S100A8','S100A9','NAMPT','AQP9','SOD2','FPR1','G0S2'],
    'T cell': ['CD3D','CD3E','TRAC','TRBC1','TRBC2','CD2','CD4','CD8A','CD8B','IL7R','CCR7','LTB'],
    'Natural killer cell': ['NKG7','GNLY','KLRD1','XCL1','XCL2','PRF1','GZMB','GZMK','FCGR3A','CTSW'],
    'B cell': ['MS4A1','CD79A','CD79B','CD22','FCER2','CD37','IGKC'],
    'Plasma cell': ['MZB1','JCHAIN','SDC1','IGHM','IGHA1','IGLC2','IGLC3','XBP1'],
    'Mast cell': ['TPSAB1','TPSB2','KIT','HDC','CMA1','MS4A2','FCER1A'],
    'Dendritic cell': ['FCER1A','CD1C','CLEC10A','HLA-DRA','HLA-DRB1','CST3','CD74'],
    'Schwann cell': ['S100B','PLP1','MPZ','SOX10','MBP','FABP7'],
    'Epithelial check': ['EPCAM','KRT8','KRT18','KRT19'],
}

# %% 2. Integrity / I/O helpers

def require(ok, message):
    if not bool(ok):
        raise ValueError(message)


def dump_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2,
                                   default=str), encoding='utf-8')


def file_sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(16 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def index_sha256(index):
    h = hashlib.sha256()
    for value in index:
        b = str(value).encode('utf-8')
        h.update(len(b).to_bytes(8, 'little'))
        h.update(b)
    return h.hexdigest()


def _sparse_hash(fmt, shape, arrays):
    # Same digest for in-memory CSR/CSC and its HDF5 representation.
    h = hashlib.sha256()
    h.update(json.dumps([fmt, [int(x) for x in shape]]).encode())
    for name, arr in arrays:
        h.update(json.dumps([name, np.dtype(arr.dtype).str,
                             int(arr.shape[0])]).encode())
        for start in range(0, arr.shape[0], 4_000_000):
            block = np.ascontiguousarray(arr[start:start + 4_000_000])
            h.update(memoryview(block).cast('B'))
    return h.hexdigest()


def raw_matrix_sha256(x):
    require(sp.issparse(x) and x.format in ('csr', 'csc'),
            'This pipeline requires CSR/CSC counts; do not densify all genes.')
    return _sparse_hash(x.format, x.shape,
                        [(k, getattr(x, k)) for k in ('data', 'indices', 'indptr')])


def h5_raw_sha256(path):
    with h5py.File(path, 'r') as f:
        g = f['X']
        encoding = g.attrs['encoding-type']
        if isinstance(encoding, bytes):
            encoding = encoding.decode()
        require(encoding in ('csr_matrix', 'csc_matrix'), 'Stored X is not CSR/CSC.')
        return _sparse_hash(encoding[:3], g.attrs['shape'],
                            [(k, g[k]) for k in ('data', 'indices', 'indptr')])


def validate_counts(x):
    require(sp.issparse(x) and x.format in ('csr', 'csc'), 'Raw X must be sparse CSR/CSC.')
    for start in range(0, x.nnz, 4_000_000):
        b = x.data[start:start + 4_000_000]
        require(np.isfinite(b).all() and (b >= 0).all(), 'Invalid raw counts.')
        if not np.issubdtype(b.dtype, np.integer):
            require(np.equal(b, np.floor(b)).all(), 'X is not integer-valued raw counts.')


def readonly_counts(x):
    # Guards in-place array writes in this session; does NOT make AnnData immutable.
    for k in ('data', 'indices', 'indptr'):
        getattr(x, k).flags.writeable = False


def make_readonly(path):
    path = Path(path)
    mode = stat.S_IMODE(path.stat().st_mode)
    if mode & 0o222:
        path.chmod(mode & ~0o222)
    require(not (stat.S_IMODE(path.stat().st_mode) & 0o222),
            f'Could not remove write permission: {path}')


# %% 3. Full atlas construction

def build_cpu_atlas():
    import torch
    import omicverse as ov

    ov.settings.cpu_init()
    require(not torch.cuda.is_available() and not torch.cuda.is_initialized(),
            'CUDA is active: restart the kernel and run this file first.')
    require(sc.__version__ == '1.11.5',
            'The existing HVG patch is audited for Scanpy 1.11.5; do not upgrade silently.')
    require(Version(metadata.version('scikit-learn')) >= Version('1.5'),
            'covariance_eigh needs a compatible scikit-learn version (>=1.5).')
    numba.set_num_threads(N_THREADS)
    torch.set_num_threads(N_THREADS)
    sc.settings.n_jobs = N_THREADS
    sc.settings.verbosity = 2
    sc.settings.autoshow = False
    sc.settings.autosave = False
    np.random.seed(SEED)

    require(COUNTS_FILE.is_file(), f'Missing counts: {COUNTS_FILE}')
    patch = PATCH_DIR / 'code/patched_function.py'
    tests_path = PATCH_DIR / 'audit/TESTS.json'
    require(patch.is_file() and tests_path.is_file(), 'Missing previously audited HVG patch/tests.')
    require(json.loads(tests_path.read_text())['passed'] is True,
            'The recorded HVG tests did not pass.')

    out = ROOT / ('CPU_full_atlas_theta2_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    out.mkdir(parents=True, exist_ok=False)
    (out / 'audit').mkdir()
    shutil.copy2(patch, out / 'audit/patched_function.py')
    shutil.copy2(tests_path, out / 'audit/HVG_TESTS.json')
    if LOCK_SOURCE_FILE:
        make_readonly(COUNTS_FILE)
    print('Hashing source file (read-only I/O) ...', flush=True)
    input_sha = file_sha256(COUNTS_FILE)
    dump_json(out / 'audit/raw_source.json', {
        'file': str(COUNTS_FILE), 'sha256': input_sha,
        'permission_mode': oct(stat.S_IMODE(COUNTS_FILE.stat().st_mode)),
        'limitation': 'File permissions are reversible; keep a separate immutable backup.'})

    hm = importlib.import_module('omicverse.external.harmony.harmony')
    require('device' in inspect.signature(hm.run_harmony).parameters,
            'Installed Harmony lacks explicit device; do not fall back automatically.')
    harmony_params = dict(theta=THETA, lamb=None, sigma=0.1, nclust=None,
                          tau=0, block_size=0.05, max_iter_harmony=100,
                          max_iter_kmeans=20, epsilon_cluster=1e-5,
                          epsilon_harmony=1e-4, alpha=0.2,
                          random_state=SEED, device='cpu', verbose=True)
    packages = {}
    for p in ('omicverse','scanpy','anndata','numpy','scipy','pandas',
              'scikit-learn','pynndescent','numba','igraph','umap-learn','torch'):
        try:
            packages[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            packages[p] = 'not found as distribution'
    config = dict(counts_file=str(COUNTS_FILE), input_sha256=input_sha,
                  hvg_patch_sha256=file_sha256(patch), packages=packages,
                  threads=N_THREADS, seed=SEED, n_hvg=N_HVG, n_pcs=N_PCS,
                  batch_key=BATCH_KEY, pearson_theta=100,
                  pearson_clip='sqrt(n_cells_per_sample)', normalize_target_sum=10000,
                  scale_max_value=10, pca_solver='covariance_eigh',
                  harmony=harmony_params, neighbors_key=NEIGHBORS_KEY,
                  neighbors=dict(n_neighbors=30, metric='euclidean', method='umap',
                                 transformer='PyNNDescentTransformer', n_jobs=N_THREADS),
                  leiden=dict(key=GROUP, resolution=1.0, flavor='igraph',
                              n_iterations=100, directed=False, use_weights=True),
                  umap=dict(min_dist=0.3, spread=1.0, maxiter=500, init_pos='spectral',
                            n_components=2, alpha=1.0, negative_sample_rate=5),
                  L1_status='PENDING_MANUAL_REVIEW', L2_started=False)
    dump_json(out / 'pipeline_parameters.json', config)
    for name, mod in [('harmony', hm)]:
        src = Path(inspect.getfile(mod))
        shutil.copy2(src, out / f'audit/{name}_installed_source.py')
    print('Loading all-gene counts ...', flush=True)
    source = sc.read_h5ad(COUNTS_FILE)
    require(source.shape == EXPECTED_SHAPE, f'Unexpected dimensions: {source.shape}')
    require(source.obs_names.is_unique and source.var_names.is_unique, 'Duplicate cell/gene IDs.')
    require(BATCH_KEY in source.obs and not source.obs[BATCH_KEY].isna().any(),
            'Missing sample_uid values.')
    require(source.obs[BATCH_KEY].astype(str).nunique() == EXPECTED_SAMPLES,
            'Unexpected number of samples.')
    validate_counts(source.X)

    # Fresh analysis container: retain ALL raw counts and ALL original obs/var.
    # Do not inherit old PCA/graphs/coordinates as this run's result.
    atlas = ad.AnnData(X=source.X, obs=source.obs.copy(), var=source.var.copy())
    atlas.uns['source_uns'] = copy.deepcopy(source.uns)
    del source
    raw_sha = raw_matrix_sha256(atlas.X)
    obs_sha = index_sha256(atlas.obs_names)
    var_sha = index_sha256(atlas.var_names)
    readonly_counts(atlas.X)

    # An independent working COPY; no preprocessing function receives atlas.X for writing.
    work = ad.AnnData(X=atlas.X.tocsr(copy=True),
                     obs=atlas.obs[[BATCH_KEY]].copy(),
                     var=pd.DataFrame(index=atlas.var_names.copy()))
    work.X.sum_duplicates()
    work.X.eliminate_zeros()
    work.X.sort_indices()
    work.obs[BATCH_KEY] = work.obs[BATCH_KEY].astype(str)
    totals = np.asarray(work.X.sum(axis=1, dtype=np.float64)).ravel()
    require(np.isfinite(totals).all() and (totals > 0).all(),
            'Zero/invalid library size: stop, do not silently remove cells.')
    atlas.obs['library_size_raw_cpu'] = totals

    sizes = work.obs[BATCH_KEY].value_counts()
    clips = []
    def observe(batch, clip, variances):
        require(np.isfinite(variances).all(), f'Non-finite HVG variances in {batch}')
        batch = str(batch)
        clips.append(dict(sample_uid=batch, n_cells=int(sizes[batch]), clip=float(clip)))
        if len(clips) % 25 == 0:
            print(f'HVG: {len(clips)}/{EXPECTED_SAMPLES} samples', flush=True)

    mod = importlib.import_module('scanpy.experimental.pp._highly_variable_genes')
    ns = dict(mod.__dict__)
    ns.update(_batch_order=lambda batches: batches, _observe=observe)
    exec(compile(patch.read_text(), str(patch), 'exec'), ns)
    with threadpool_limits(limits=N_THREADS):
        hvg = ns['_highly_variable_pearson_residuals'](
            work, batch_key=BATCH_KEY, theta=100, clip=None,
            n_top_genes=N_HVG, chunksize=1000, check_values=False, inplace=False)
    require(hvg.index.equals(atlas.var_names), 'HVG gene order mismatch.')
    require(len(clips) == EXPECTED_SAMPLES and
            {r['sample_uid'] for r in clips} == set(sizes.index), 'HVG sample coverage mismatch.')
    require(all(np.isclose(r['clip'], np.sqrt(r['n_cells']), rtol=0, atol=1e-12)
                for r in clips), 'Per-sample clipping mismatch.')
    require(hvg['highly_variable'].notna().all(), 'Missing HVG flags.')
    hv = hvg['highly_variable'].to_numpy(dtype=bool)
    require(int(hv.sum()) == N_HVG, 'HVG count mismatch.')
    for col in hvg:
        atlas.var[col] = hvg[col]
    hvg.to_csv(out / 'HVG_statistics.tsv.gz', sep='\t')
    pd.DataFrame(clips).to_csv(out / 'sample_clips.tsv', sep='\t', index=False)
    pd.Series(atlas.var_names[hv], name='gene').to_csv(out / 'HVG2000.tsv', sep='\t', index=False)

    print('Normalize/log1p ALL genes on a copy ...', flush=True)
    with threadpool_limits(limits=N_THREADS):
        sc.pp.normalize_total(work, target_sum=10_000, exclude_highly_expressed=False)
        sc.pp.log1p(work)
    for start in range(0, work.X.nnz, 4_000_000):
        require(np.isfinite(work.X.data[start:start + 4_000_000]).all(), 'Invalid log1p expression.')
    atlas.layers['log1p'] = work.X      # ALL 37,336 genes; raw atlas.X unchanged
    del work
    gc.collect()

    print('Scale/PCA on a temporary 2,000-HVG object ...', flush=True)
    pca_work = ad.AnnData(X=atlas.layers['log1p'][:, hv].copy(),
                         obs=pd.DataFrame(index=atlas.obs_names.copy()),
                         var=pd.DataFrame(index=atlas.var_names[hv].copy()))
    with threadpool_limits(limits=N_THREADS):
        sc.pp.scale(pca_work, zero_center=True, max_value=10)
        sc.pp.pca(pca_work, n_comps=N_PCS, svd_solver='covariance_eigh',
                  zero_center=True, random_state=SEED, mask_var=None, dtype='float32')
    atlas.obsm['X_pca_cpu'] = np.asarray(pca_work.obsm['X_pca'], dtype=np.float32).copy()
    atlas.obsm['X_pca'] = atlas.obsm['X_pca_cpu'].copy()
    # Expand PCA loadings to ALL genes. Zero for non-HVG, explicitly documented.
    pcs = np.zeros((atlas.n_vars, N_PCS), dtype=np.float32)
    pcs[hv] = pca_work.varm['PCs']
    atlas.varm['PCs_cpu'] = pcs
    atlas.uns['pca_cpu'] = copy.deepcopy(pca_work.uns['pca'])
    for col in ('mean', 'std'):
        vals = np.full(atlas.n_vars, np.nan)
        vals[hv] = pca_work.var[col].to_numpy()
        atlas.var[f'cpu_scale_{col}_HVG_only'] = vals
    require(atlas.obsm['X_pca_cpu'].shape == (atlas.n_obs, N_PCS) and
            np.isfinite(atlas.obsm['X_pca_cpu']).all(), 'Invalid PCA.')
    np.save(out / 'X_pca_cpu_checkpoint.npy', atlas.obsm['X_pca_cpu'])
    del pca_work, pcs
    gc.collect()

    print('CPU Harmony, theta=2.0 ...', flush=True)
    # Copy the relatively small PCA input to guard against in-place changes in a backend.
    with threadpool_limits(limits=N_THREADS):
        result = hm.run_harmony(atlas.obsm['X_pca_cpu'].copy(),
                                atlas.obs[[BATCH_KEY]].assign(**{BATCH_KEY: atlas.obs[BATCH_KEY].astype(str)}),
                                BATCH_KEY, **harmony_params)
    corrected = np.asarray(result.result())
    require(corrected.shape == (atlas.n_obs, N_PCS) and np.isfinite(corrected).all(),
            'Unexpected Harmony result; do not transpose or continue silently.')
    # Save actual class and installed implementation source, not an assumed backend.
    config['actual_harmony_class'] = type(result).__module__ + '.' + type(result).__name__
    try:
        cls_file = Path(inspect.getfile(type(result)))
        shutil.copy2(cls_file, out / 'audit/harmony_actual_class_source.py')
        config['actual_harmony_source_sha256'] = file_sha256(cls_file)
    except (TypeError, OSError):
        config['actual_harmony_source_sha256'] = 'unavailable'
    obj = getattr(result, 'objective_harmony', None)
    if obj is not None:
        if torch.is_tensor(obj):
            obj = obj.detach().cpu().numpy()
        obj = np.asarray(obj, dtype=np.float64).ravel()
        pd.DataFrame({'record_index': np.arange(len(obj)), 'objective': obj}).to_csv(
            out / 'harmony_objective_history.tsv', sep='\t', index=False)
    atlas.obsm['X_harmony_cpu'] = np.asarray(corrected, dtype=np.float32).copy()
    atlas.obsm['X_pca_harmony'] = atlas.obsm['X_harmony_cpu'].copy()
    np.save(out / 'X_harmony_cpu_checkpoint.npy', atlas.obsm['X_harmony_cpu'])
    del result, corrected
    gc.collect()

    print('CPU neighbors -> Leiden -> UMAP ...', flush=True)
    transformer = PyNNDescentTransformer(n_neighbors=30, metric='euclidean',
                                         random_state=SEED, n_jobs=N_THREADS)
    # Avoid nested BLAS oversubscription; PyNNDescent/Numba keep their worker settings.
    with threadpool_limits(limits=1, user_api='blas'):
        sc.pp.neighbors(atlas, n_neighbors=30, n_pcs=N_PCS, use_rep='X_harmony_cpu',
                        metric='euclidean', method='umap', transformer=transformer,
                        random_state=SEED, key_added=NEIGHBORS_KEY)
    sc.tl.leiden(atlas, resolution=1.0, random_state=SEED, n_iterations=100,
                 flavor='igraph', directed=False, use_weights=True,
                 neighbors_key=NEIGHBORS_KEY, key_added=GROUP)
    with threadpool_limits(limits=N_THREADS):
        sc.tl.umap(atlas, min_dist=0.3, spread=1.0, maxiter=500, init_pos='spectral',
                   n_components=2, alpha=1.0, negative_sample_rate=5,
                   random_state=SEED, neighbors_key=NEIGHBORS_KEY)
    atlas.obsm['X_umap_cpu'] = atlas.obsm['X_umap'].copy()
    atlas.obsm['X_umap_harmony'] = atlas.obsm['X_umap'].copy()
    atlas.uns['umap_cpu'] = copy.deepcopy(atlas.uns['umap'])
    for k in ('X_pca_cpu', 'X_harmony_cpu', 'X_umap_cpu'):
        require(np.isfinite(atlas.obsm[k]).all(), f'Nonfinite {k}.')
    require(atlas.obs[GROUP].notna().all(), 'Missing Leiden labels.')
    require(atlas.shape == EXPECTED_SHAPE and index_sha256(atlas.obs_names) == obs_sha
            and index_sha256(atlas.var_names) == var_sha, 'Cell/gene set or order changed.')
    print('Checking that every stored raw value and source file is unchanged ...', flush=True)
    require(raw_matrix_sha256(atlas.X) == raw_sha, 'Raw counts changed IN MEMORY.')
    require(file_sha256(COUNTS_FILE) == input_sha, 'Source FILE changed.')
    config['threadpools'] = threadpool_info()
    config['n_leiden_clusters'] = int(atlas.obs[GROUP].nunique())
    atlas.uns['CPU_pipeline_json'] = json.dumps(config, ensure_ascii=False, default=str)
    atlas.uns['expression_layout'] = {
        'X': 'ALL-gene QC/deduplicated integer-valued raw counts; DO NOT normalize X',
        'layers_log1p': 'ALL-gene normalize_total(10000), then natural log1p; NOT batch corrected',
        'PCA_loadings': 'varm[PCs_cpu] has zero rows for non-HVG genes',
        'raw_attribute': 'Not populated; raw counts are explicitly stored in X',
    }
    atlas.uns['raw_integrity'] = dict(source_path=str(COUNTS_FILE), source_sha256=input_sha,
                                    sparse_X_sha256=raw_sha, cell_order_sha256=obs_sha,
                                    gene_order_sha256=var_sha)
    atlas.uns['manual_L1_marker_panel'] = markers
    atlas.uns['CPU_L1_status'] = 'PENDING_MANUAL_REVIEW'
    atlas.uns['output_directory'] = str(out)
    final = out / 'CPU_Harmony_theta2_full_raw_and_log1p.h5ad'
    temp = out / 'CPU_Harmony_theta2_full_raw_and_log1p.partial.h5ad'
    print(f'Saving FULL atlas: {final}', flush=True)
    atlas.write_h5ad(temp, compression='lzf')
    require(h5_raw_sha256(temp) == raw_sha, 'Saved raw counts do NOT match the input.')
    with h5py.File(temp, 'r') as f:
        for key in ('layers/log1p','obsm/X_harmony_cpu','obsm/X_umap_cpu', f'obs/{GROUP}'):
            require(key in f, f'Saved object missing {key}.')
        for key in ('connectivities_key', 'distances_key'):
            graph_key = atlas.uns[NEIGHBORS_KEY][key]
            require(f'obsp/{graph_key}' in f, f'Missing stored graph: {graph_key}')
    temp.rename(final)
    make_readonly(final)  # subsequent annotation must save a NEW file
    final_sha = file_sha256(final)
    dump_json(out / 'manifest.json', dict(final_h5ad=str(final), final_h5ad_sha256=final_sha,
                                        raw_integrity=atlas.uns['raw_integrity'],
                                        status='PAUSED_FOR_CPU_L1_REVIEW'))
    dump_json(out / 'pipeline_parameters.json', config)
    pd.DataFrame({'cell_id': atlas.obs_names,
                  GROUP: atlas.obs[GROUP].astype(str).to_numpy()}).to_csv(
        out / 'cell_cluster_map.tsv.gz', sep='\t', index=False)
    print(f'Completed: {atlas.shape}; {config["n_leiden_clusters"]} clusters\n{final}', flush=True)
    return atlas, out, final



adata, OUT, FINAL_FILE = build_cpu_atlas()

# NOTEBOOK CELL 9, execution_count=None
# %% 4. L1 marker plots and expression evidence (NO label writeback)

def create_L1_review(atlas, out, marker_panel=None):
    """All summaries use all cells; no downsampling, no raw X modification."""
    marker_panel = markers if marker_panel is None else marker_panel
    out = Path(out) / 'L1_marker_review'
    out.mkdir(parents=True, exist_ok=True)
    require(GROUP in atlas.obs and 'log1p' in atlas.layers, 'Missing Leiden/log1p layer.')
    require(atlas.obs_names.is_unique and atlas.var_names.is_unique, 'Duplicate IDs.')
    require(not atlas.obs[GROUP].isna().any(), 'Missing cluster labels.')
    # Reapply the in-memory write guard after loading an h5ad.
    readonly_counts(atlas.X)
    labels = atlas.obs[GROUP].astype(str)
    order = sorted(labels.unique(), key=lambda x: int(x))
    cat = pd.Categorical(labels.to_numpy(), categories=order, ordered=True)
    n = np.bincount(cat.codes, minlength=len(order))
    count_table = pd.DataFrame({'cluster': order, 'n_cells': n, 'fraction': n / atlas.n_obs})
    count_table.to_csv(out / 'cluster_cell_counts.tsv', sep='\t', index=False)

    panel = {k: list(dict.fromkeys(g for g in gs if g in atlas.var_names))
             for k, gs in marker_panel.items()}
    panel = {k: gs for k, gs in panel.items() if gs}
    missing = {k: [g for g in gs if g not in atlas.var_names]
               for k, gs in marker_panel.items()}
    dump_json(out / 'missing_markers.json', {k: gs for k, gs in missing.items() if gs})
    dump_json(out / 'marker_panel.json', marker_panel)
    genes = list(dict.fromkeys(g for gs in panel.values() for g in gs))
    require(len(genes) > 0, 'No marker genes found.')
    idx = atlas.var_names.get_indexer(genes)

    # Select only requested genes, NOT the entire raw matrix as dense.
    raw = atlas.X[:, idx].tocsr(copy=True)
    raw.sum_duplicates()
    raw.eliminate_zeros()
    log = atlas.layers['log1p'][:, idx].tocsr(copy=True)
    log.sum_duplicates()
    log.eliminate_zeros()
    positive = raw.copy()
    positive.data = np.ones(positive.nnz, dtype=np.float32)
    positive_log = log.copy()
    positive_log.data = np.ones(positive_log.nnz, dtype=np.float32)
    require((positive != positive_log).nnz == 0, 'Raw and log1p detection supports differ.')
    del positive_log
    # Sparse cluster x cell indicator, then exact all-cell group sums.
    membership = sp.csr_matrix((np.ones(atlas.n_obs, dtype=np.float64),
                                (cat.codes, np.arange(atlas.n_obs))),
                               shape=(len(order), atlas.n_obs))
    def group_mean(x):
        return pd.DataFrame((membership @ x).toarray() / n[:, None],
                            index=pd.Index(order, name='cluster'), columns=genes)
    mean_raw = group_mean(raw)
    mean_log = group_mean(log)
    frac = group_mean(positive)
    for name, df in [('marker_mean_raw_counts.tsv', mean_raw),
                     ('marker_mean_log1p.tsv', mean_log),
                     ('marker_detection_fraction.tsv', frac)]:
        df.to_csv(out / name, sep='\t')
    long = pd.concat([df.rename_axis(index='cluster', columns='gene').stack().rename(name)
                      for name, df in [('mean_log1p', mean_log), ('mean_raw_count', mean_raw),
                                       ('detection_fraction', frac)]], axis=1).reset_index()
    long.to_csv(out / 'marker_evidence_long.tsv', sep='\t', index=False)
    del membership, positive, raw, log
    gc.collect()

    # Tiny plotting-only object; displayed means/fractions are PRECOMPUTED from ALL cells.
    # No million-cell dense plotting DataFrame is needed.
    relative = (mean_log - mean_log.min(axis=0)).div(
        (mean_log.max(axis=0) - mean_log.min(axis=0)).replace(0, np.nan), axis=1).fillna(0)
    mural = {k: panel[k] for k in ('Fibroblast','Fibromyocyte','Smooth muscle cell','Pericyte')
             if k in panel}
    plt.ioff()
    def save_dot(name, subpanel, values, title):
        if not subpanel:
            return
        cols = [g for gs in subpanel.values() for g in gs]
        # Unique plotting-only IDs prevent duplicate-column expansion when a gene
        # occurs in multiple marker groups. True gene names are restored on the axis.
        ids = [f'm{i:04d}' for i in range(len(cols))]
        plot_panel, pos = {}, 0
        for group, gs in subpanel.items():
            plot_panel[group] = ids[pos:pos + len(gs)]
            pos += len(gs)
        colors = values.loc[order, cols].copy()
        sizes = frac.loc[order, cols].copy()
        colors.columns = sizes.columns = ids
        summary = ad.AnnData(X=colors.to_numpy(dtype=np.float32),
                            obs=pd.DataFrame({GROUP: pd.Categorical(order, categories=order)},
                                             index=['cluster_' + c for c in order]),
                            var=pd.DataFrame(index=ids))
        dp = sc.pl.DotPlot(summary, var_names=plot_panel, groupby=GROUP,
                           categories_order=order, use_raw=False, standard_scale=None,
                           dot_color_df=colors, dot_size_df=sizes,
                           figsize=(max(13, 0.22 * len(cols)), max(7, 0.30 * len(order) + 3)))
        dp.style(dot_min=0, dot_max=1, dot_edge_lw=0)
        dp.legend(colorbar_title=title, size_title='Cells with raw counts > 0 (%)')
        dp.make_figure()
        main_ax = dp.get_axes()['mainplot_ax']
        main_ax.set_xticklabels(cols, rotation=90)
        fig = main_ax.figure
        try:
            for ext in ('png', 'pdf', 'svg'):
                fig.savefig(out / f'{name}.{ext}', dpi=300, bbox_inches='tight')
        finally:
            plt.close(fig)
    save_dot('L1_dotplot_mean_log1p', panel, mean_log, 'Mean log1p(CP10k)')
    save_dot('L1_dotplot_relative', panel, relative, 'Gene-wise relative mean (0-1)')
    save_dot('L1_mural_dotplot_relative', mural, relative, 'Gene-wise relative mean (0-1)')
    save_dot('L1_dotplot_raw_counts_reference', panel, mean_raw, 'Mean raw counts (depth-dependent)')

    # Plot-only shuffle; original cell order and coordinates are untouched.
    keys = [k for k in ('study', GROUP) if k in atlas.obs]
    shuffle = np.random.default_rng(SEED).permutation(atlas.n_obs)
    vis = ad.AnnData(obs=atlas.obs[keys].iloc[shuffle].copy(),
                     obsm={'X_umap_cpu': atlas.obsm['X_umap_cpu'][shuffle].copy()})
    for key in keys:
        vis.obs[key] = vis.obs[key].astype('category').cat.remove_unused_categories()
        palette = sc.pl.palettes.default_102[:len(vis.obs[key].cat.categories)]
        require(len(palette) == len(vis.obs[key].cat.categories), 'Not enough unique colors.')
        fig, ax = plt.subplots(figsize=(10, 9))
        try:
            sc.pl.embedding(vis, basis='X_umap_cpu', color=key, palette=palette,
                            size=1.0, marker='.', alpha=1.0, frameon=False,
                            legend_fontsize=8, show=False, save=False, ax=ax)
            ax.set_aspect('equal', adjustable='box')
            for collection in ax.collections:
                if isinstance(collection, PathCollection):
                    collection.set_linewidths(0)
                    collection.set_edgecolor('none')
                    collection.set_rasterized(True)
            for ext in ('png', 'pdf', 'svg'):
                fig.savefig(out / f'UMAP_{key}.{ext}', dpi=400, bbox_inches='tight')
        finally:
            plt.close(fig)
    review = count_table.assign(proposed_L1='', decision='PENDING', comment='')
    # Do not overwrite manually edited decisions if replotting later.
    template = out / 'cluster_L1_review_template.tsv'
    if not template.exists():
        review.to_csv(template, sep='\t', index=False)
    (out / 'README.txt').write_text(
        'Main DotPlot color: all-cell mean log1p(CP10k), NOT corrected expression.\n'
        'Point size: fraction of cells with original counts > 0; range 0-100%.\n'
        'Relative DotPlot is gene-wise scaled; raw-count figure is depth-dependent reference only.\n'
        'UMAP scatter rasterized; text is vector in PDF/SVG. No coordinates changed.\n'
        'No automatic L1 labels; review the cluster table and marker evidence first.\n', encoding='utf-8')
    print(f'L1 review saved: {out}\nNo L1 labels applied.', flush=True)
    return out



REVIEW_OUT = create_L1_review(adata, OUT)

# NOTEBOOK CELL 10, execution_count=None
adata.obsm["X_umap_harmony"] = adata.obsm["X_umap"].copy()
ov.utils.embedding(adata,
    basis='X_umap_harmony',
    color=["study"],  
    palette=ov.palette()[12:],
    show=False,
    frameon='small',)
# NOTEBOOK CELL 11, execution_count=18
markers = {
    "Endothelial cell": [
        "PECAM1", "VWF", "KDR", "CLDN5", "EMCN",
        "RGCC", "CA4", "ACKR1", "PLVAP",
        "SPARCL1", "EDN1", "ECSCR",
        "FABP4", "EFNB2", "GJA5", "SEMA3G", "ESM1"
    ],

    "Fibroblast": [
        "DCN", "LUM", "COL1A1", "COL1A2",
        "COL3A1", "COL5A1", "COL14A1",
        "PDGFRA", "FBLN1", "THY1", "C3", "C7"
    ],

    "Fibromyocyte": [
        "FN1", "AEBP1", "POSTN", "COL3A1",
        "LUM", "TNFRSF11B", "ACTA2",
        "TAGLN", "MYL9", "TCF21"
    ],

    "Smooth muscle cell": [
        "ACTA2", "MYH11", "TAGLN", "CNN1",
        "MYL9", "TPM2", "CALD1",
        "LMOD1", "MYLK"
    ],

    "Pericyte": [
        "RGS5", "CSPG4", "MCAM", "PDGFRB",
        "ABCC9", "KCNJ8", "COX4I2",
        "NDUFA4L2", "NOTCH3", "DES"
    ],

    "Macrophage": [
        "LST1", "TYROBP", "FCER1G", "CTSS",
        "CD68", "AIF1", "ITGAM", "CSF1R",
        "C1QA", "C1QB", "C1QC",
        "LGALS3", "APOE", "MS4A7",
        "MSR1", "TREM2", "GPNMB", "SPP1"
    ],

    "Monocyte": [
        "FCN1", "S100A8", "S100A9", "S100A12",
        "VCAN", "CCR2", "LYZ", "CTSS",
        "CTSD", "LGALS3"
    ],

    "Neutrophil": [
        "FCGR3B", "CSF3R", "CXCR2",
        "S100A8", "S100A9", "NAMPT",
        "AQP9", "SOD2", "FPR1", "G0S2"
    ],

    "T cell": [
        "CD3D", "CD3E", "TRAC", "TRBC1",
        "TRBC2", "CD2", "CD4", "CD8A",
        "CD8B", "IL7R", "CCR7", "LTB"
    ],

    "Natural killer cell": [
        "NKG7", "GNLY", "KLRD1", "XCL1",
        "XCL2", "PRF1", "GZMB", "GZMK",
        "FCGR3A", "CTSW"
    ],

    "B cell": [
        "MS4A1", "CD79A", "CD79B",
        "CD22", "FCER2", "CD37", "IGKC"
    ],

    "Plasma cell": [
        "MZB1", "JCHAIN", "SDC1",
        "IGHM", "IGHA1", "IGLC2",
        "IGLC3", "XBP1"
    ],

    "Mast cell": [
        "TPSAB1", "TPSB2", "KIT",
        "HDC", "CMA1", "MS4A2",
        "FCER1A"
    ],

    "Dendritic cell": [
        "FCER1A", "CD1C", "CLEC10A",
        "HLA-DRA", "HLA-DRB1",
        "CST3", "CD74"
    ],

    "Schwann cell": [
        "S100B", "PLP1", "MPZ",
        "SOX10", "MBP", "FABP7"
    ],

    "Epithelial check": [
        "EPCAM", "KRT8", "KRT18", "KRT19"
    ],
}
# NOTEBOOK CELL 12, execution_count=4
adata
# NOTEBOOK CELL 13, execution_count=5
from pathlib import Path
import importlib.util

if "adata" not in globals() or "markers" not in globals():
    raise RuntimeError(
        "请在仍有全量 adata 和 markers 的原 Notebook 中执行。"
    )

ROOT = Path(
    "/home/anzhen_vhd/scRNAseq/codex/"
    "panvascular_atlas_rebuild_v2_20260907/"
    "06_reanalysis_takeover_20260908"
)

MODULE_FILE = (
    ROOT / "L1_CPU_partial_and_local_code"
    / "l1_cpu_local_refinement.py"
)

if not MODULE_FILE.is_file():
    raise FileNotFoundError(f"请先上传主程序：{MODULE_FILE}")

spec = importlib.util.spec_from_file_location("l1job", MODULE_FILE)
l1job = importlib.util.module_from_spec(spec)
spec.loader.exec_module(l1job)

# 使用母图谱登记的工作文件夹，
# 不依赖可能在之前绘图时被重新赋值的 OUT 变量。
WORK_PARENT = Path(adata.uns["output_directory"])

L1_RUN = l1job.prepare_partial_L1(
    adata,
    markers=markers,
    work_parent=WORK_PARENT,
    patch_dir=ROOT / "annotation/HVG_clip_review_20260910",
)

print("\n当前初步 L1：")
print(adata.obs["L1"].value_counts())

print("\n本次结果目录：")
print(L1_RUN)

print("\n已保存的全量母图谱：")
print(L1_RUN / "CPU_theta2_L1_partial_full.h5ad")
# NOTEBOOK CELL 14, execution_count=None
# 两个独立分支并行，每个分支20线程，总预算40线程。
# 使用当前Notebook的Python环境，不清空当前adata。
adata
LOCAL_STATUS = l1job.run_local_batch(
    L1_RUN,
    max_workers=2,
    total_threads=40,
)

LOCAL_STATUS
# NOTEBOOK CELL 15, execution_count=7
from pathlib import Path

SCRIPT = Path(
    "/home/anzhen_vhd/scRNAseq/codex/"
    "panvascular_atlas_rebuild_v2_20260907/"
    "06_reanalysis_takeover_20260908/"
    "L1_writeback_Unknown_R2_onecell.py"
)

if not SCRIPT.is_file():
    raise FileNotFoundError(f"请先将下载的Python主程序上传到：\n{SCRIPT}")

if "adata" not in globals() or "markers" not in globals():
    raise RuntimeError(
        "请在仍保留全量adata和markers的原Notebook中执行，不要重启kernel。"
    )

exec(
    compile(SCRIPT.read_text(encoding="utf-8"), str(SCRIPT), "exec"),
    globals(),
)
# NOTEBOOK CELL 16, execution_count=8
from pathlib import Path

SCRIPT = Path(
    "/home/anzhen_vhd/scRNAseq/codex/"
    "panvascular_atlas_rebuild_v2_20260907/"
    "06_reanalysis_takeover_20260908/"
    "L1_R2_memory_writeback_Perivascular_R3_onecell.py"
)

if not SCRIPT.is_file():
    raise FileNotFoundError(f"请先上传Python代码文件：\n{SCRIPT}")

if "adata" not in globals() or "markers" not in globals():
    raise RuntimeError(
        "请在仍保留全量adata和markers的原Notebook中运行，不要重启kernel。"
    )

exec(
    compile(SCRIPT.read_text(encoding="utf-8"), str(SCRIPT), "exec"),
    globals(),
)
# NOTEBOOK CELL 17, execution_count=9
from pathlib import Path

SCRIPT = Path(
    "/home/anzhen_vhd/scRNAseq/codex/"
    "panvascular_atlas_rebuild_v2_20260907/"
    "06_reanalysis_takeover_20260908/"
    "L1_R3_writeback_global_review_onecell.py"
)

if not SCRIPT.is_file():
    raise FileNotFoundError(f"请先上传Python主程序到：\n{SCRIPT}")

if "adata" not in globals() or "markers" not in globals():
    raise RuntimeError(
        "请在仍保留全量adata和markers的原Notebook中运行，不要重启kernel。"
    )

exec(
    compile(SCRIPT.read_text(encoding="utf-8"), str(SCRIPT), "exec"),
    globals(),
)
# NOTEBOOK CELL 18, execution_count=10
from pathlib import Path

SCRIPT = Path(
    "/home/anzhen_vhd/scRNAseq/codex/"
    "panvascular_atlas_rebuild_v2_20260907/"
    "06_reanalysis_takeover_20260908/"
    "C0_NK_TCR_audit_onecell.py"
)

if not SCRIPT.is_file():
    raise FileNotFoundError(f"请先上传审核代码：\n{SCRIPT}")

if "adata" not in globals():
    raise RuntimeError(
        "请在仍保留全量adata的原Notebook中运行，不要重启kernel。"
    )

exec(
    compile(SCRIPT.read_text(encoding="utf-8"), str(SCRIPT), "exec"),
    globals(),
)
# NOTEBOOK CELL 19, execution_count=11
from pathlib import Path

SCRIPT = Path(
    "/home/anzhen_vhd/scRNAseq/codex/"
    "panvascular_atlas_rebuild_v2_20260907/"
    "06_reanalysis_takeover_20260908/"
    "freeze_clean_level1_onecell.py"
)

if not SCRIPT.is_file():
    raise FileNotFoundError(f"请先上传Python主程序：\n{SCRIPT}")

if "adata" not in globals():
    raise RuntimeError(
        "请在仍保留当前全量adata的原Notebook中执行，不要重启kernel。"
    )

exec(
    compile(SCRIPT.read_text(encoding="utf-8"), str(SCRIPT), "exec"),
    globals(),
)
# NOTEBOOK CELL 20, execution_count=12
adata
# NOTEBOOK CELL 21, execution_count=None
