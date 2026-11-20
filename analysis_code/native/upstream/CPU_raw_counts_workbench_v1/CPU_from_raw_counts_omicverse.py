# %%
# 1. 线程/CPU设置：必须在科学计算库导入前执行
import os
N_THREADS = min(40, len(os.sched_getaffinity(0)))
os.environ["CUDA_VISIBLE_DEVICES"] = ""
for key in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"]:
    os.environ[key] = str(N_THREADS)

from pathlib import Path
from datetime import datetime
import gc, json, importlib, hashlib
import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc
import omicverse as ov
import torch
from threadpoolctl import threadpool_limits

ov.settings.cpu_init()
sc.settings.n_jobs = N_THREADS
torch.set_num_threads(N_THREADS)
thread_limits = threadpool_limits(limits=N_THREADS)
assert not torch.cuda.is_available()

ROOT = Path("/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907")
ANNOTATION = ROOT / "06_reanalysis_takeover_20260908/annotation"
CPU_BASE = ANNOTATION / "full_atlas_CPU_Harmony_v1"
PATCH_DIR = ANNOTATION / "HVG_clip_review_20260910"
OUT = ANNOTATION / ("CPU_from_counts_manual_" + datetime.now().strftime("%Y%m%d_%H%M%S"))
OUT.mkdir(exist_ok=False)

BATCH_KEY = "sample_uid"
N_HVG = 2000
N_PCS = 30
N_NEIGHBORS = 30
SEED = 0


# %%
# 2. 从审核对象读取 counts 路径；实际加载的是完整 raw counts
review = sc.read_h5ad(CPU_BASE / "data/full_atlas_CPU_review.h5ad", backed="r")
raw_counts_path = Path(str(review.uns["raw_counts_connection"]))
expected_ids = review.obs_names.copy()
review.file.close()
del review

print("Raw counts:", raw_counts_path)
adata = sc.read_h5ad(raw_counts_path)
assert adata.obs_names.equals(expected_ids), "细胞ID或顺序不一致，停止"
assert adata.obs_names.is_unique and adata.var_names.is_unique
assert adata.shape == (1175115, 37336)
assert not adata.obs[BATCH_KEY].isna().any()
assert adata.obs[BATCH_KEY].nunique() == 247
assert sparse.issparse(adata.X)
# 分块检查非负整数 counts，避免产生全矩阵临时副本
for start in range(0, len(adata.X.data), 10_000_000):
    values = adata.X.data[start:start + 10_000_000]
    assert np.isfinite(values).all() and (values >= 0).all()
    assert (values == np.floor(values)).all()


# %%
# 3. sample-aware HVG：保留已验证的逐样本 clip=sqrt(n_sample)
# 不使用未经修正的 Scanpy 1.11.5 batch clip=None 路径。
assert sc.__version__ == "1.11.5", "本补丁限定 Scanpy 1.11.5，版本不符先审核"
assert json.loads((PATCH_DIR / "audit/TESTS.json").read_text())["passed"]
patch_path = PATCH_DIR / "code/patched_function.py"
patch = patch_path.read_text()
archived_patch = CPU_BASE / "audit/HVG_patched_function.py"
assert hashlib.sha256(patch.encode()).hexdigest() == hashlib.sha256(archived_patch.read_bytes()).hexdigest()
(OUT / "HVG_patched_function.py").write_text(patch)

module = importlib.import_module("scanpy.experimental.pp._highly_variable_genes")
namespace = dict(module.__dict__)
namespace["_batch_order"] = lambda batches: batches
sample_sizes = adata.obs[BATCH_KEY].astype(str).value_counts()
clip_records = []
def observe_sample(batch, clip, variances):
    assert np.isfinite(variances).all()
    clip_records.append({"sample_uid": str(batch),
                         "n_cells": int(sample_sizes[str(batch)]),
                         "clip": float(clip)})
    if len(clip_records) % 25 == 0:
        print(f"HVG: {len(clip_records)}/{len(sample_sizes)} samples")
namespace["_observe"] = observe_sample
exec(compile(patch, str(patch_path), "exec"), namespace)

hvg = namespace["_highly_variable_pearson_residuals"](
    adata,
    batch_key=BATCH_KEY,
    theta=100,
    clip=None,  # 在这个项目补丁中，每个样本独立计算 sqrt(n_sample)
    n_top_genes=N_HVG,
    chunksize=1000,
    inplace=False,
)
assert len(clip_records) == 247
assert {x["sample_uid"] for x in clip_records} == set(sample_sizes.index)
assert all(x["clip"] == np.sqrt(x["n_cells"]) for x in clip_records)
assert hvg.index.equals(adata.var_names)
for col in hvg.columns:
    adata.var[col] = hvg[col]
genes = hvg.index[hvg["highly_variable"]].tolist()
assert len(genes) == N_HVG
hvg.to_csv(OUT / "HVG_statistics.tsv.gz", sep="\t")
pd.DataFrame(clip_records).to_csv(OUT / "sample_clips.tsv", sep="\t", index=False)
pd.Series(genes, name="gene").to_csv(OUT / "HVG2000.tsv", sep="\t", index=False)


# %%
# 4. 必须先在全基因范围归一化，再提取 HVG
# Scanpy用于基础normalize/log1p；其余按OmicVerse CPU接口组织。
adata.X.eliminate_zeros()
sc.pp.normalize_total(adata, target_sum=10000)
sc.pp.log1p(adata)

full_log1p_path = OUT / "full_log1p.h5ad"
adata.uns["raw_counts_connection"] = str(raw_counts_path)
adata.write_h5ad(full_log1p_path, compression="lzf")

# 不在内存里重复保留百万细胞的完整 counts 和 log1p 矩阵
adata_hvg = adata[:, genes].copy()
del adata
gc.collect()
adata = adata_hvg
del adata_hvg
adata.uns["expression_connection"] = str(full_log1p_path)


# %%
# 5. OmicVerse CPU scale + PCA
ov.pp.scale(adata, max_value=10, layers_add="scaled", to_sparse=False)
ov.pp.pca(
    adata,
    n_pcs=N_PCS,
    layer="scaled",
    inplace=True,
    random_state=SEED,
    svd_solver="covariance_eigh",
)
assert adata.obsm["X_pca"].shape == (adata.n_obs, N_PCS)
adata.obsm["X_pca_cpu"] = adata.obsm["X_pca"].astype("float32", copy=True)
np.save(OUT / "pca.npy", adata.obsm["X_pca_cpu"])
np.save(OUT / "pca_loadings.npy", adata.varm["PCs"])
# PCA结束后不再需要scaled层，删除以减小Harmony包装器的对象复制开销。
# adata.X 仍为 HVG 的 log1p 表达，不是scaled表达。
if "scaled" in adata.layers:
    del adata.layers["scaled"]
gc.collect()
adata.write_h5ad(OUT / "HVG_PCA_checkpoint.h5ad", compression="lzf")


# %%
# 6. OmicVerse CPU Harmony：显式参数与原全量 CPU 基线一致
HARMONY = dict(
    theta=1.0, lamb=None, sigma=0.1, nclust=None, tau=0,
    block_size=0.05, max_iter_harmony=100, max_iter_kmeans=20,
    epsilon_cluster=1e-5, epsilon_harmony=1e-4, alpha=0.2,
    random_state=SEED, verbose=True, device="cpu",
)
ov.single.batch_correction(
    adata,
    batch_key=BATCH_KEY,
    use_rep="X_pca_cpu",
    methods="harmony",
    n_pcs=N_PCS,
    **HARMONY,
)
assert np.isfinite(adata.obsm["X_pca_harmony"]).all()
np.save(OUT / "harmony.npy", adata.obsm["X_pca_harmony"])


# %%
# 7. CPU neighbors / Leiden / UMAP
# 显式使用CPU transformer，避免不同包装器默认选择不同邻居算法。
from pynndescent import PyNNDescentTransformer
transformer = PyNNDescentTransformer(
    n_neighbors=N_NEIGHBORS, metric="euclidean",
    random_state=SEED, n_jobs=N_THREADS,
)
sc.pp.neighbors(
    adata,
    n_neighbors=N_NEIGHBORS,
    use_rep="X_pca_harmony",
    metric="euclidean",
    transformer=transformer,
    random_state=SEED,
    key_added="neighbors_cpu_n30",
)
sc.tl.leiden(
    adata, resolution=0.5, random_state=SEED,
    n_iterations=100, flavor="igraph", directed=False,
    use_weights=True, neighbors_key="neighbors_cpu_n30",
    key_added="leiden_cpu_r0p5",
)
ov.pp.umap(
    adata,
    neighbors_key="neighbors_cpu_n30",
    min_dist=0.3, spread=1.0, maxiter=500,
    init_pos="spectral", n_components=2,
    alpha=1.0, negative_sample_rate=5, random_state=SEED,
)
adata.obsm["X_umap_harmony"] = adata.obsm["X_umap"].copy()
adata.uns["formal_L1_applied"] = False
adata.write_h5ad(OUT / "CPU_Harmony_manual_review.h5ad", compression="lzf")
import importlib.metadata as im
(OUT / "run_parameters.json").write_text(json.dumps({
    "raw_counts": str(raw_counts_path), "HVG": N_HVG,
    "clip": "per_sample_sqrt_n", "theta_HVG": 100,
    "target_sum": 10000, "scale_max_value": 10,
    "PCA": N_PCS, "PCA_solver": "covariance_eigh",
    "Harmony": HARMONY, "neighbors": N_NEIGHBORS,
    "Leiden_resolution": 0.5, "UMAP_min_dist": 0.3,
    "UMAP_epochs": 500, "seed": SEED, "threads": N_THREADS,
    "versions": {p: im.version(p) for p in ["omicverse", "scanpy", "numpy", "pynndescent", "umap-learn"]},
}, indent=2))
sc.pl.embedding(adata, basis="X_umap_harmony", color=["leiden_cpu_r0p5", "study"])


# %%
# 8. 从全基因log1p取marker，按CPU Leiden绘制DotPlot
panel = {
    "Endothelial": ["PECAM1", "VWF", "KDR"],
    "Fibroblast": ["DCN", "LUM", "COL1A1", "COL1A2"],
    "SMC": ["MYH11", "ACTA2", "TAGLN", "CNN1"],
    "Pericyte program": ["RGS5", "ABCC9", "KCNJ8", "CSPG4"],
    "Myeloid": ["LST1", "TYROBP", "C1QC", "FCN1", "S100A8", "FCGR3B"],
    "T/NK": ["CD3D", "CD3E", "TRAC", "NKG7", "GNLY"],
    "B/Plasma": ["MS4A1", "CD79A", "MZB1", "JCHAIN"],
    "Mast": ["TPSAB1", "KIT"],
    "Schwann": ["S100B", "PLP1", "MPZ"],
}
expr = sc.read_h5ad(adata.uns["expression_connection"], backed="r")
assert expr.obs_names.equals(adata.obs_names)
missing = [g for genes in panel.values() for g in genes if g not in expr.var_names]
print("Missing markers:", missing)
panel = {k: [g for g in v if g in expr.var_names] for k, v in panel.items()}
panel = {k: v for k, v in panel.items() if v}
marker_genes = list(dict.fromkeys(g for genes in panel.values() for g in genes))
try:
    plot_data = expr[:, marker_genes].to_memory()
finally:
    expr.file.close()
plot_data.obs["leiden_cpu_r0p5"] = adata.obs["leiden_cpu_r0p5"].copy()
dp = sc.pl.dotplot(
    plot_data, var_names=panel, groupby="leiden_cpu_r0p5",
    use_raw=False, standard_scale=None, expression_cutoff=0,
    mean_only_expressed=False, dendrogram=False, return_fig=True,
)
dp.savefig(OUT / "L1_marker_dotplot.pdf")
dp.dot_color_df.to_csv(OUT / "marker_mean_log1p.tsv", sep="\t")
dp.dot_size_df.to_csv(OUT / "marker_detection_fraction.tsv", sep="\t")
dp.show()
