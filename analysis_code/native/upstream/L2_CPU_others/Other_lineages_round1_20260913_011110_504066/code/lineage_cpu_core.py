# -*- coding: utf-8 -*-
"""Generic CPU worker for one of the 13 other frozen L1 lineages.
Adapted in a NEW source file from the supplied, already-run Macrophage Round1 pipeline.
Scientific preprocessing is shared; task-specific marker panels and UID prefixes differ.
No original Macrophage code/results or notebook variables are changed.
Called only by the detached scheduler with a read-only source verification ticket.
"""
from __future__ import annotations

import argparse
import base64
import copy
from datetime import datetime
import gc
import hashlib
import html
import importlib
import importlib.metadata as metadata
import inspect
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import time
import traceback
import zipfile

CODE_VERSION = "Other_L2_CPU_Round1_1.0"
ROOT = Path("/home/anzhen_vhd/scRNAseq/codex/"
            "panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908")
SOURCE = ROOT / "panvascular_CPU_theta2_level1_frozen_v1.h5ad"
PATCH_DIR = ROOT / "annotation/HVG_clip_review_20260910"
PATCH_SHA256 = "41db6ac24969e95d6241770ddb88a4530fe6be662e0d3f95bed47b1ba195254b"
L1_KEY = "level1_annotation"
BATCH = "sample_uid"
TARGET = "UNCONFIGURED"
UID_PREFIX = "UNCONFIGURED"
EXPECTED_SHAPE = (1175115, 37336)
EXPECTED_N = 0
EXPECTED_L1 = {
    "Endothelial cell":104415,"Fibroblast":251371,"Fibromyocyte":29347,
    "Smooth muscle cell":243346,"Pericyte":4513,"Macrophage":210967,
    "Monocyte":31046,"Neutrophil":18509,"T cell":156349,"Natural killer cell":25719,
    "B cell":10496,"Plasma cell":4802,"Mast cell":12929,"Dendritic cell":10556,
    "Schwann cell":2832,"Unknown":57918,
}
DEFAULTS = {
    "n_hvg":3000, "n_pcs":50, "n_neighbors":30, "seed":0,
    "pearson_theta":100, "pearson_chunksize":1000, "hvg_min_sample_cells":2,
    "normalize_target_sum":10000, "scale_max_value":10,
    "pca_solver":"covariance_eigh", "metric":"euclidean",
    "leiden_resolutions":[0.4,0.8,1.2], "baseline_resolution":0.8,
    "leiden_flavor":"igraph", "leiden_n_iterations":100,
    "umap_min_dist":0.3, "umap_spread":1.0, "umap_maxiter":500,
    "min_cells_gene_test":3, "findmarkers_method":"wilcoxon", "tie_correct":True,
    "findmarkers_all_three_resolutions":True, "figure_dpi":250,
    "harmony":{
        "theta":2.0, "lamb":None, "sigma":0.1, "nclust":None, "tau":0,
        "block_size":0.05, "max_iter_harmony":100, "max_iter_kmeans":20,
        "epsilon_cluster":1e-5, "epsilon_harmony":1e-4, "alpha":0.2,
        "random_state":0, "device":"cpu", "verbose":True,
    },
}
STATE_PANEL = {}
BOUNDARY_PANEL = {}
MISSING = {"", "nr", "na", "n/a", "nan", "none", "unknown", "not reported", "missing", "<na>"}


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def stamp():
    return datetime.now().strftime("%Y%m%d_%H%M%S_%f")


def sha_file(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def sha_index(values):
    h = hashlib.sha256()
    for value in values:
        b = str(value).encode("utf-8")
        h.update(len(b).to_bytes(8, "little")); h.update(b)
    return h.hexdigest()


def json_write(path, obj):
    p = Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    temp = p.with_name(p.name + ".tmp")
    temp.write_text(json.dumps(obj, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    temp.replace(p)


def readonly(path):
    p = Path(path); p.chmod(stat.S_IMODE(p.stat().st_mode) & ~0o222)


def scalar(value):
    import numpy as np
    if isinstance(value, bytes): return value.decode("utf-8")
    if isinstance(value, np.ndarray) and value.ndim == 0: return scalar(value.item())
    if isinstance(value, np.generic): return value.item()
    return value


def sparse_parts_hash(fmt, shape, arrays):
    import numpy as np
    h = hashlib.sha256(json.dumps([fmt, [int(n) for n in shape]]).encode())
    for key, x in arrays:
        h.update(json.dumps([key, np.dtype(x.dtype).str, int(x.shape[0])]).encode())
        for start in range(0, x.shape[0], 4_000_000):
            h.update(memoryview(np.ascontiguousarray(x[start:start+4_000_000])).cast("B"))
    return h.hexdigest()


def sparse_hash(x):
    import scipy.sparse as sp
    require(sp.issparse(x) and x.format in ("csr","csc"), "Expected sparse expression.")
    return sparse_parts_hash(x.format, x.shape, [(k,getattr(x,k)) for k in ("data","indices","indptr")])


def h5_sparse_hash(group):
    enc = scalar(group.attrs["encoding-type"])
    require(enc in ("csr_matrix","csc_matrix"), "Unexpected HDF5 sparse encoding.")
    return sparse_parts_hash(enc[:3], group.attrs["shape"], [(k,group[k]) for k in ("data","indices","indptr")])


def lock_sparse(x):
    for key in ("data","indices","indptr"):
        getattr(x,key).flags.writeable = False


def validate_raw(x):
    import numpy as np
    import scipy.sparse as sp
    require(sp.issparse(x), "X must contain sparse full-gene raw counts.")
    for start in range(0, x.nnz, 4_000_000):
        b = x.data[start:start+4_000_000]
        require(np.isfinite(b).all() and (b >= 0).all(), "Invalid raw counts.")
        if not np.issubdtype(b.dtype, np.integer):
            require(np.equal(b, np.floor(b)).all(), "X is not integer-valued raw counts.")


def normalized_copy(raw, target=10000):
    """Independent raw working CSR plus full-gene library sizes; normalization occurs later."""
    import numpy as np
    x = raw.tocsr(copy=True).astype(np.float32)
    x.sum_duplicates(); x.eliminate_zeros(); x.sort_indices()
    totals = np.asarray(x.sum(axis=1, dtype=np.float64)).ravel()
    require(np.isfinite(totals).all() and (totals > 0).all(), "Zero/invalid library size; no cells deleted.")
    return x, totals


def csr_rows_from_h5(group, rows, row_block=25000):
    """Read ONLY X from a CSR h5ad, in contiguous row blocks; no full layers loaded."""
    import numpy as np
    import scipy.sparse as sp
    rows = np.asarray(rows, dtype=np.int64)
    require(rows.ndim==1 and len(rows)>0 and np.all(np.diff(rows)>0), "Selected source rows must be unique and ordered.")
    shape = tuple(int(v) for v in group.attrs["shape"])
    require(rows[0]>=0 and rows[-1]<shape[0], "Out-of-range source row.")
    enc = scalar(group.attrs["encoding-type"])
    if enc != "csr_matrix":
        # Official AnnData selection for the uncommon CSC input; no source rewrite.
        require(enc == "csc_matrix", "Source X must be CSR or CSC.")
        from anndata.io import sparse_dataset
        print("Source X is CSC; using AnnData's read-only sparse selection.", flush=True)
        return sp.csr_matrix(sparse_dataset(group)[rows, :])
    ptr = group["indptr"][:]
    blocks=[]
    for lo in range(0, shape[0], row_block):
        hi=min(lo+row_block,shape[0]); a=np.searchsorted(rows,lo); b=np.searchsorted(rows,hi)
        if a==b: continue
        start=int(ptr[lo]); end=int(ptr[hi])
        part=sp.csr_matrix((group["data"][start:end],group["indices"][start:end],
                           ptr[lo:hi+1]-start), shape=(hi-lo,shape[1]))
        blocks.append(part[rows[a:b]-lo,:].copy())
    out=sp.vstack(blocks,format="csr"); del blocks
    require(out.shape==(len(rows),shape[1]), "Raw row extraction changed shape.")
    return out


def select_lineage_metadata(obs, var):
    import numpy as np
    require(obs.index.is_unique and var.index.is_unique, "Nonunique cell/gene IDs.")
    require((len(obs),len(var))==EXPECTED_SHAPE, "Wrong frozen atlas dimensions.")
    require(L1_KEY in obs and obs[L1_KEY].notna().all(), "Missing frozen level1_annotation.")
    require(obs[L1_KEY].astype(str).value_counts().to_dict()==EXPECTED_L1, "Frozen L1 composition differs from the approved version.")
    require(BATCH in obs and obs[BATCH].notna().all() and obs[BATCH].astype(str).nunique()==247,
            "Wrong or missing sample_uid metadata.")
    rows=np.flatnonzero(obs[L1_KEY].astype(str).eq(TARGET).to_numpy())
    require(len(rows)==EXPECTED_N, "Expected frozen lineage cell count mismatch.")
    sub=obs.iloc[rows].copy()
    # Remove unused categories, not observations or clinical information.
    import pandas as pd
    for key in sub:
        if isinstance(sub[key].dtype,pd.CategoricalDtype): sub[key]=sub[key].cat.remove_unused_categories()
    sub["parent_cell_id"]=sub.index.astype(str)
    sub["parent_row"]=rows
    return rows,sub


def read_source(task, run):
    import numpy as np
    import anndata as ad
    import h5py
    try: from anndata.io import read_elem
    except ImportError: from anndata._io.specs import read_elem
    src=Path(task["source"]); manifest=json.loads(Path(task["source_manifest"]).read_text())
    require(manifest.get("level1_frozen") is True and manifest.get("annotation_key")==L1_KEY,
            "Source manifest is not the final frozen level1_annotation version.")
    check_source_snapshot(task)
    source_hash=task["verified_source_sha256"]
    require(source_hash==manifest["sha256"],"Source manifest hash changed.")
    with h5py.File(src,"r") as f:
        obs=read_elem(f["obs"]); var=read_elem(f["var"])
        require(bool(read_elem(f["uns/level1_frozen"])),"h5ad freeze flag missing.")
        rows,sub=select_lineage_metadata(obs,var)
        pipeline=json.loads(str(scalar(read_elem(f["uns/CPU_pipeline_json"]))))
        old_patch=pipeline.get("hvg_patch_sha256")
        require(old_patch is None or old_patch==PATCH_SHA256,"Parent atlas used a different HVG patch.")
        raw=csr_rows_from_h5(f["X"],rows)
        # Reference coordinates are not used for any local neighbor graph or clustering.
        parent_umap=f["obsm/X_umap_cpu"][:][rows,:] if "obsm/X_umap_cpu" in f else None
    # Preserve source gene information with explicit parent prefix for prior HVG/QC columns.
    for key in list(var.columns):
        if key.startswith(("highly_variable","residual_","cpu_scale","local_scale")) or key in ("means","variances"):
            var=var.rename(columns={key:"parent_"+key})
    C=ad.AnnData(X=raw,obs=sub,var=var)
    validate_raw(C.X)
    raw_hash=sparse_hash(C.X); lock_sparse(C.X)
    if parent_umap is not None: C.obsm["X_umap_parent_reference"]=np.asarray(parent_umap).copy()
    C.uns["source_frozen_L1"]={"path":str(src),"sha256":source_hash,
        "level1_key":L1_KEY,"level1_value":TARGET,"parent_pipeline_json":json.dumps(pipeline,ensure_ascii=False)}
    record={"source_h5ad":str(src),"source_sha256":source_hash,"source_n_cells":len(obs),
        "source_n_genes":len(var),"selected_n_cells":C.n_obs,"selected_gene_count":C.n_vars,
        "selected_cell_order_sha256":sha_index(C.obs_names),"gene_order_sha256":sha_index(C.var_names),
        "selected_raw_sha256":raw_hash,"only_frozen_target_lineage":True,"source_read_only":True}
    json_write(run/"audit/input_integrity.json",record)
    sub.to_csv(run/"tables/input_cell_metadata.tsv.gz",sep="\t",index=True,index_label="source_obs_index")
    del obs,sub,raw,parent_umap; gc.collect()
    return C,record


def fit_hvg(work, run, cfg, threads):
    import numpy as np
    import pandas as pd
    from threadpoolctl import threadpool_limits
    sizes=work.obs[BATCH].astype(str).value_counts()
    eligible=sizes[sizes>=cfg["hvg_min_sample_cells"]].index
    require(len(eligible)>0,"No sample with >=2 cells for within-sample HVG variance.")
    sm=pd.DataFrame({"sample_uid":sizes.index,"n_lineage_cells":sizes.values})
    sm["used_for_hvg"]=sm.n_lineage_cells.ge(cfg["hvg_min_sample_cells"])
    sm["cells_retained_in_all_other_steps"]=True
    sm.to_csv(run/"tables/HVG_sample_eligibility.tsv",sep="\t",index=False)
    mask=work.obs[BATCH].astype(str).isin(eligible).to_numpy()
    fit=work if mask.all() else work[mask,:].copy()
    fit.obs[BATCH]=fit.obs[BATCH].astype(str)
    fs=fit.obs[BATCH].value_counts(); clips=[]
    def observe(batch,clip,variances):
        require(np.isfinite(variances).all(),f"Nonfinite HVG variances: {batch}")
        clips.append({"sample_uid":str(batch),"n_cells":int(fs[str(batch)]),"clip":float(clip)})
        if len(clips)%25==0: print(f"HVG: {len(clips)}/{len(fs)} samples",flush=True)
    mod=importlib.import_module("scanpy.experimental.pp._highly_variable_genes")
    ns=dict(mod.__dict__); ns.update(_batch_order=lambda batches:batches,_observe=observe)
    patch=run/"audit/patched_function.py"
    require(sha_file(patch)==PATCH_SHA256,"HVG patch snapshot changed.")
    exec(compile(patch.read_text(),str(patch),"exec"),ns)
    with threadpool_limits(limits=threads):
        hvg=ns["_highly_variable_pearson_residuals"](fit,batch_key=BATCH,theta=cfg["pearson_theta"],
            clip=None,n_top_genes=cfg["n_hvg"],chunksize=cfg["pearson_chunksize"],check_values=False,inplace=False)
    require(hvg.index.equals(work.var_names),"HVG gene order mismatch.")
    require(len(clips)==len(fs) and {x["sample_uid"] for x in clips}==set(fs.index),"HVG batch coverage mismatch.")
    require(all(np.isclose(r["clip"],np.sqrt(r["n_cells"]),rtol=0,atol=1e-12) for r in clips),"Wrong per-sample clip.")
    require(int(hvg.highly_variable.sum())==cfg["n_hvg"],"Wrong number of HVGs.")
    hvg.to_csv(run/"tables/HVG_statistics.tsv.gz",sep="\t")
    pd.DataFrame(clips).to_csv(run/"tables/HVG_sample_clips.tsv",sep="\t",index=False)
    pd.Series(hvg.index[hvg.highly_variable],name="gene").to_csv(run/"tables/HVG3000.tsv",sep="\t",index=False)
    del fit; gc.collect()
    return hvg


def make_embedding(C, run, cfg, threads, progress):
    import numpy as np
    import pandas as pd
    import anndata as ad
    import scanpy as sc
    import torch
    from threadpoolctl import threadpool_limits
    from pynndescent import PyNNDescentTransformer
    raw_before=sparse_hash(C.X)
    x,totals=normalized_copy(C.X,cfg["normalize_target_sum"])
    work=ad.AnnData(X=x,obs=C.obs[[BATCH]].copy(),var=pd.DataFrame(index=C.var_names.copy()))
    work.obs[BATCH]=work.obs[BATCH].astype(str)
    C.obs["L2_total_raw_counts"]=totals
    C.obs["L2_detected_genes"]=work.X.getnnz(axis=1)
    mt=C.var_names.astype(str).str.upper().str.startswith("MT-")
    C.obs["L2_pct_counts_mt"]=100*np.asarray(work.X[:,mt].sum(axis=1,dtype=np.float64)).ravel()/totals
    progress("SAMPLE_AWARE_HVG_ON_COUNTS")
    hvg=fit_hvg(work,run,cfg,threads)
    for col in hvg: C.var[col]=hvg[col].to_numpy()
    progress("NORMALIZE_LOG1P_WORKING_COPY")
    with threadpool_limits(limits=threads):
        sc.pp.normalize_total(work,target_sum=cfg["normalize_target_sum"],exclude_highly_expressed=False)
        sc.pp.log1p(work)
    C.layers["log1p"]=work.X; C.uns["log1p"]={"base":None}
    del work,x,totals; gc.collect()
    progress("LOCAL_CPU_SCALE_PCA")
    hv=C.var.highly_variable.to_numpy(dtype=bool)
    p=ad.AnnData(X=C.layers["log1p"][:,hv].copy(),obs=pd.DataFrame(index=C.obs_names.copy()),
                 var=pd.DataFrame(index=C.var_names[hv]))
    with threadpool_limits(limits=threads):
        sc.pp.scale(p,zero_center=True,max_value=cfg["scale_max_value"])
        sc.pp.pca(p,n_comps=cfg["n_pcs"],svd_solver=cfg["pca_solver"],zero_center=True,
                  mask_var=None,random_state=cfg["seed"],dtype="float32")
    C.obsm["X_pca_lineage_L2"]=np.asarray(p.obsm["X_pca"],dtype=np.float32).copy()
    load=np.zeros((C.n_vars,cfg["n_pcs"]),dtype=np.float32); load[hv]=p.varm["PCs"]
    C.varm["PCs_lineage_L2"]=load; C.uns["pca_lineage_L2"]=copy.deepcopy(p.uns["pca"])
    for key in ("mean","std"):
        v=np.full(C.n_vars,np.nan); v[hv]=p.var[key].to_numpy(); C.var["L2_scale_"+key]=v
    del p,load; gc.collect()
    progress("HARMONY_EXPLICIT_CPU_THETA2")
    hm=importlib.import_module("omicverse.external.harmony.harmony")
    require("device" in inspect.signature(hm.run_harmony).parameters,"Harmony must support explicit device=cpu.")
    hp=copy.deepcopy(cfg["harmony"])
    with threadpool_limits(limits=threads):
        if C.obs[BATCH].nunique()>1:
            result=hm.run_harmony(C.obsm["X_pca_lineage_L2"].copy(),
                C.obs[[BATCH]].assign(**{BATCH:C.obs[BATCH].astype(str)}),BATCH,**hp)
            z=result.result()
            if torch.is_tensor(z):
                require(z.device.type=="cpu","Harmony returned a non-CPU tensor.")
                z=z.detach().cpu().numpy()
            z=np.asarray(z,dtype=np.float32)
            cfg["actual_harmony_class"]=type(result).__module__+"."+type(result).__name__
            require("cuda" not in cfg["actual_harmony_class"].lower() and "gpu" not in cfg["actual_harmony_class"].lower(),
                    "Unexpected GPU Harmony implementation.")
            try: shutil.copy2(inspect.getfile(type(result)),run/"audit/Harmony_actual_source.py")
            except (OSError,TypeError): pass
            hist=getattr(result,"objective_harmony",None)
            if hist is not None:
                if torch.is_tensor(hist): hist=hist.detach().cpu().numpy()
                h=np.asarray(hist).ravel()
                pd.DataFrame({"record":np.arange(len(h)),"objective":h}).to_csv(
                    run/"tables/Harmony_objective_history.tsv",sep="\t",index=False)
                cfg["harmony_objective_records"]=len(h)
            del result
        else:
            z=C.obsm["X_pca_lineage_L2"].copy(); cfg["actual_harmony_class"]="SKIPPED_SINGLE_SAMPLE"
    require(z.shape==(C.n_obs,cfg["n_pcs"]) and np.isfinite(z).all(),"Invalid Harmony output.")
    C.obsm["X_harmony_lineage_L2"]=z.copy(); del z; gc.collect()
    progress("CPU_NEIGHBORS")
    transformer=PyNNDescentTransformer(n_neighbors=cfg["n_neighbors"],metric=cfg["metric"],
                                       random_state=cfg["seed"],n_jobs=threads)
    with threadpool_limits(limits=1,user_api="blas"):
        sc.pp.neighbors(C,n_neighbors=cfg["n_neighbors"],n_pcs=cfg["n_pcs"],use_rep="X_harmony_lineage_L2",
            method="umap",metric=cfg["metric"],transformer=transformer,random_state=cfg["seed"],
            key_added="neighbors_lineage_L2")
    progress("CPU_LEIDEN_THREE_PARTITIONS")
    for resolution in cfg["leiden_resolutions"]:
        tag=res_tag(resolution); key="lineage_L2_"+tag
        sc.tl.leiden(C,resolution=resolution,random_state=cfg["seed"],n_iterations=cfg["leiden_n_iterations"],
            flavor="igraph",directed=False,use_weights=True,neighbors_key="neighbors_lineage_L2",key_added=key)
        old=C.obs[key].astype(str); order=sorted(old.unique(),key=int)
        prefix=UID_PREFIX+"_"+tag+"_"
        C.obs[key]=pd.Categorical(prefix+old,categories=[prefix+k for k in order],ordered=True)
        print(tag,"clusters:",len(order),flush=True)
    progress("CPU_UMAP")
    with threadpool_limits(limits=threads):
        sc.tl.umap(C,min_dist=cfg["umap_min_dist"],spread=cfg["umap_spread"],maxiter=cfg["umap_maxiter"],
            init_pos="spectral",n_components=2,alpha=1.0,negative_sample_rate=5,
            random_state=cfg["seed"],neighbors_key="neighbors_lineage_L2")
    C.obsm["X_umap_lineage_L2"]=np.asarray(C.obsm.pop("X_umap"),dtype=np.float32).copy()
    require(sparse_hash(C.X)==raw_before,"Raw counts changed during preprocessing.")
    C.uns["round1_config_json"]=json.dumps(cfg,ensure_ascii=False,default=str)
    return C


def res_tag(value):
    return "r"+str(value).replace(".","p")


def clean_panel(panel, genes):
    available=set(genes)
    present={k:list(dict.fromkeys(g for g in vals if g in available)) for k,vals in panel.items()}
    missing={k:[g for g in vals if g not in available] for k,vals in panel.items()}
    return {k:v for k,v in present.items() if v},{k:v for k,v in missing.items() if v}


def aggregate_matrix_markers(raw, log, labels, genes, var_names, order=None):
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    labels=np.asarray(labels,dtype=str)
    require(len(labels)==raw.shape[0]==log.shape[0],"Marker aggregation rows mismatch.")
    if order is None: order=list(dict.fromkeys(labels))
    cats=pd.Categorical(labels,categories=order)
    require((cats.codes>=0).all(),"Unrecognized group in marker aggregation.")
    n=np.bincount(cats.codes,minlength=len(order))
    require((n>0).all(),"Empty marker aggregation group.")
    ids=pd.Index(var_names).get_indexer(genes)
    require((ids>=0).all(),"Marker gene is absent.")
    r=raw[:,ids].tocsr(copy=True); r.sum_duplicates(); r.eliminate_zeros()
    l=log[:,ids].tocsr(copy=True); l.sum_duplicates(); l.eliminate_zeros()
    detection=r.copy(); detection.data=np.ones(r.nnz,dtype=np.float64)
    membership=sp.csr_matrix((np.ones(len(labels),dtype=np.float64),(cats.codes,np.arange(len(labels)))),
                              shape=(len(order),len(labels)))
    def mean(x):
        return pd.DataFrame((membership@x).toarray()/n[:,None],index=pd.Index(order,name="cluster_uid"),columns=genes)
    return mean(l),mean(detection),mean(r),pd.Series(n,index=order,name="n_cells")


def dotplot(panel, mean, frac, path, relative=False, dpi=250):
    import numpy as np
    import pandas as pd
    import anndata as ad
    import scanpy as sc
    import matplotlib.pyplot as plt
    if not panel: return
    cols=[g for vals in panel.values() for g in vals]
    vals=mean.copy()
    if relative: vals=(vals-vals.min()).div((vals.max()-vals.min()).replace(0,np.nan)).fillna(0)
    unique_positions=[f"plot_{i:04d}" for i in range(len(cols))]
    pp={}; offset=0
    for key,gs in panel.items():
        pp[key]=unique_positions[offset:offset+len(gs)]; offset+=len(gs)
    colors=vals.loc[:,cols].copy(); sizes=frac.loc[:,cols].copy()
    colors.columns=sizes.columns=unique_positions
    order=list(mean.index.astype(str))
    dummy=ad.AnnData(X=colors.to_numpy(dtype=np.float32),
        obs=pd.DataFrame({"plot_group":pd.Categorical(order,categories=order)},
                         index=["summary_"+str(i) for i in range(len(order))]),
        var=pd.DataFrame(index=unique_positions))
    dp=sc.pl.DotPlot(dummy,var_names=pp,groupby="plot_group",categories_order=order,
        use_raw=False,standard_scale=None,dot_color_df=colors,dot_size_df=sizes,
        figsize=(max(12,.20*len(cols)),max(5,.29*len(order)+3)))
    dp.style(dot_min=0,dot_max=1,dot_edge_lw=0)
    dp.legend(colorbar_title="Gene-wise relative mean (0-1)" if relative else "Mean log1p(CP10k)",
              size_title="Cells with raw counts > 0 (%)")
    dp.make_figure(); ax=dp.get_axes()["mainplot_ax"]; ax.set_xticklabels(cols,rotation=90,fontsize=7)
    fig=ax.figure
    try:
        for ext in ("png","pdf","svg"):
            fig.savefig(str(path)+"."+ext,dpi=dpi,bbox_inches="tight",facecolor="white")
    finally: plt.close(fig)


def umap_plot(C,key,path,cfg):
    import numpy as np
    import pandas as pd
    import anndata as ad
    import scanpy as sc
    import matplotlib.pyplot as plt
    from matplotlib.collections import PathCollection
    ix=np.random.default_rng(cfg["seed"]).permutation(C.n_obs)
    obs=C.obs[[key]].iloc[ix].copy()
    numeric=pd.api.types.is_numeric_dtype(obs[key].dtype)
    palette=None
    if not numeric:
        vals=obs[key].astype("string").fillna("Missing")
        cats=list(C.obs[key].cat.categories) if isinstance(C.obs[key].dtype,pd.CategoricalDtype) else sorted(vals.unique())
        if "Missing" in set(vals) and "Missing" not in cats: cats.append("Missing")
        if len(cats)>60:
            # Metadata only: full values remain in tables; avoid unreadable giant legends.
            keep=set(vals.value_counts().head(25).index)
            vals=vals.where(vals.isin(keep),"Other sources (plot only)")
            cats=sorted(vals.unique())
        obs[key]=pd.Categorical(vals,categories=cats)
        require(len(cats)<=len(sc.pl.palettes.default_102),"Too many plotting categories.")
        palette=sc.pl.palettes.default_102[:len(cats)]
    vis=ad.AnnData(obs=obs,obsm={"X_plot":C.obsm["X_umap_lineage_L2"][ix].copy()})
    fig,ax=plt.subplots(figsize=(10,8))
    try:
        kwargs={"palette":palette} if not numeric else {}
        sc.pl.embedding(vis,basis="X_plot",color=key,size=1.2,marker=".",alpha=1,frameon=False,
            legend_fontsize=7,title=key,show=False,save=False,ax=ax,**kwargs)
        ax.set_aspect("equal",adjustable="box")
        for c in ax.collections:
            if isinstance(c,PathCollection):
                c.set_edgecolor("none");c.set_linewidths(0);c.set_rasterized(True)
        for ext in ("png","pdf","svg"):
            fig.savefig(str(path)+"."+ext,dpi=cfg["figure_dpi"],bbox_inches="tight",facecolor="white")
    finally: plt.close(fig)


def observed_metadata(values):
    s=values.astype("string").fillna("")
    return s[~s.str.strip().str.lower().isin(MISSING)]


def cluster_sources(C,key,folder):
    import numpy as np
    import pandas as pd
    rows=[]
    for uid,ix in C.obs.groupby(key,observed=True,sort=False).indices.items():
        block=C.obs.iloc[ix]
        row={"cluster_uid":str(uid),"n_cells":len(ix),"fraction_within_lineage":len(ix)/C.n_obs}
        for field in (BATCH,"study","tissue","disease","condition","disease_context","donor_uid_v2"):
            if field not in block: continue
            valid=observed_metadata(block[field]); vc=valid.value_counts()
            row["n_"+field]=len(vc); row["missing_"+field]=len(block)-len(valid)
            row["top_"+field]=str(vc.index[0]) if len(vc) else ""
            row["top_"+field+"_fraction"]=float(vc.iloc[0]/len(block)) if len(vc) else np.nan
        if "donor_uid_v2" in block:
            ds=observed_metadata(block["donor_uid_v2"])
            prov=ds.str.contains("PROVISIONAL",case=False,na=False)
            row["n_provisional_donor_keys"]=int(ds[prov].nunique())
            row["n_cells_provisional_donor"]=int(prov.sum())
            row["n_nonmissing_nonprovisional_donor_keys"]=int(ds[~prov].nunique())
            row["donor_interpretation"]="reported keys, not independently verified biological donor count"
        for field in ("L2_total_raw_counts","L2_detected_genes","L2_pct_counts_mt"):
            row["median_"+field]=float(block[field].median())
        rows.append(row)
    out=pd.DataFrame(rows);out.to_csv(folder/"cluster_summary.tsv",sep="\t",index=False)
    for field in (BATCH,"study","tissue","disease","condition","disease_context","donor_uid_v2"):
        if field in C.obs:
            s=C.obs[field].astype("string").fillna("Missing")
            s=s.where(~s.str.strip().str.lower().isin(MISSING),"Missing")
            pd.crosstab(C.obs[key],s,dropna=False).to_csv(folder/f"cluster_by_{field}.tsv",sep="\t")
    return out


def discover_markers(C,key,tag,run,cfg):
    import numpy as np
    import pandas as pd
    import anndata as ad
    import scanpy as sc
    folder=run/"tables"/tag; cache=run/"data"/f"FindMarkers_cache_{tag}.h5ad"
    sizes=C.obs[key].value_counts(); tested=sizes[(sizes>=2)&((C.n_obs-sizes)>=2)].index.tolist()
    skipped=sizes[~sizes.index.isin(tested)]
    skipped.rename_axis("cluster_uid").to_csv(folder/"FindMarkers_untestable_groups.tsv",sep="\t")
    gene_n=C.var["n_cells_detected_lineage"].to_numpy()
    mask=gene_n>=cfg["min_cells_gene_test"]
    rankkey="rank_genes_"+tag
    require(mask.any(),"No expressed genes available for marker discovery.")
    if len(sizes)<2 or not tested:
        pd.DataFrame(columns=["group","names","scores","logfoldchanges","pvals","pvals_adj"]).to_csv(
            folder/"FindMarkers_all_tested.tsv.gz",sep="\t",index=False)
        return pd.DataFrame(),{}
    # Cache stores statistics only; no expression or labels. Verify grouping signature on reuse.
    signature=hashlib.sha256((sha_index(C.obs[key].astype(str))+sha_index(C.var_names[mask])+C.uns["source_frozen_L1"]["sha256"]+json.dumps({k:cfg[k] for k in ("min_cells_gene_test","tie_correct","findmarkers_method")},sort_keys=True)).encode()).hexdigest()
    if cache.exists():
        tiny=sc.read_h5ad(cache)
        require(tiny.uns["grouping_signature"]==signature,"FindMarkers cache membership mismatch.")
        C.uns[rankkey]=tiny.uns["rank_genes_groups"]
        del tiny
    else:
        print(f"FindMarkers {tag}: {len(tested)} clusters, {int(mask.sum())} genes, all {C.n_obs:,} cells",flush=True)
        sc.tl.rank_genes_groups(C,groupby=key,groups=tested,reference="rest",layer="log1p",use_raw=False,
            mask_var=mask,method="wilcoxon",tie_correct=cfg["tie_correct"],pts=True,
            n_genes=int(mask.sum()),rankby_abs=False,corr_method="benjamini-hochberg",key_added=rankkey)
        tiny=ad.AnnData()
        tiny.uns["rank_genes_groups"]=C.uns[rankkey]; tiny.uns["grouping_signature"]=signature
        temp=cache.with_suffix(".partial.h5ad"); tiny.write_h5ad(temp,compression="lzf"); temp.replace(cache)
        del tiny
    rank=sc.get.rank_genes_groups_df(C,group=None,key=rankkey)
    if "pct_nz_group" not in rank:
        pts=C.uns[rankkey]["pts"]; rest=C.uns[rankkey]["pts_rest"]
        rank["pct_nz_group"]=[float(pts.loc[g,c]) for c,g in zip(rank.group,rank.names)]
        rank["pct_nz_reference"]=[float(rest.loc[g,c]) for c,g in zip(rank.group,rank.names)]
    rank["pct_difference"]=rank.pct_nz_group-rank.pct_nz_reference
    rank.to_csv(folder/"FindMarkers_all_tested.tsv.gz",sep="\t",index=False)
    rank.groupby("group",observed=True,sort=False).head(100).to_csv(folder/"FindMarkers_top100.tsv",sep="\t",index=False)
    # Figure selection only: does NOT assign identity or alter clusters.
    noninformative=rank.names.astype(str).str.match(r"^(MT-|RPL|RPS|HBA[12]$|HBB$|MALAT1$)")
    good=rank.loc[(rank.logfoldchanges>0.25)&(rank.pvals_adj<0.05)&(rank.pct_nz_group>=0.10)
                  &(rank.pct_difference>=0.05)&~noninformative].copy()
    top=good.groupby("group",observed=True,sort=False).head(20)
    top.to_csv(folder/"FindMarkers_top20_marker_candidates.tsv",sep="\t",index=False)
    topmap=top.groupby("group",observed=True)["names"].agg(lambda x:";".join(x)).to_dict()
    return top,topmap


def boundary_evidence(C,key,panels,run):
    """Descriptive raw multi-gene detection ONLY; never a classifier or label update."""
    import numpy as np
    import pandas as pd
    records={}; available=set(C.var_names)
    for name,genes in panels.items():
        use=[g for g in genes if g in available]
        if not use: continue
        x=C.X[:,C.var_names.get_indexer(use)].tocsr(copy=True)
        x.sum_duplicates();x.eliminate_zeros()
        records[name+"__detected_genes_ge1"]=np.asarray((x>0).sum(axis=1)).ravel().astype(np.uint8)
        records[name+"__detected_genes_ge2"]=np.asarray((x>=2).sum(axis=1)).ravel().astype(np.uint8)
    percell=pd.DataFrame(records,index=C.obs_names)
    percell.insert(0,key,C.obs[key].astype(str).to_numpy())
    percell.to_csv(run/"tables/boundary_multi_marker_counts_per_cell.tsv.gz",sep="\t",index=True,index_label="cell_id")
    rows=[]
    for uid,ix in C.obs.groupby(key,observed=True).indices.items():
        for name,genes in panels.items():
            for cutoff in (1,2):
                col=name+f"__detected_genes_ge{cutoff}"
                if col not in percell: continue
                v=percell[col].iloc[ix].to_numpy()
                rows.append({"cluster_uid":str(uid),"panel":name,"raw_count_cutoff":cutoff,"n_cells":len(ix),
                    "mean_n_detected_panel_genes":float(v.mean()),
                    "fraction_at_least_2_panel_genes":float((v>=2).mean()),
                    "interpretation":"DESCRIPTIVE_ONLY_NOT_AN_IDENTITY_CLASSIFIER"})
    pd.DataFrame(rows).to_csv(run/"tables/boundary_multi_marker_by_cluster.tsv",sep="\t",index=False)


def write_reviews(C,run,cfg,progress):
    import numpy as np
    import pandas as pd
    state,ms=clean_panel(STATE_PANEL,C.var_names); boundary,mb=clean_panel(BOUNDARY_PANEL,C.var_names)
    json_write(run/"tables/missing_markers.json",{"state_panel":ms,"boundary_panel":mb})
    allpan={**state,**boundary}; genes=list(dict.fromkeys(g for gs in allpan.values() for g in gs))
    bx=res_tag(cfg["baseline_resolution"]); baseline="lineage_L2_"+bx
    allsummary=[]; baseline_mean=baseline_frac=None
    for r in cfg["leiden_resolutions"]:
        tag=res_tag(r); key="lineage_L2_"+tag; folder=run/"tables"/tag; folder.mkdir(exist_ok=True)
        progress("MARKER_REVIEW_"+tag)
        labels=C.obs[key].astype(str); order=list(C.obs[key].cat.categories)
        mean,frac,raw,n=aggregate_matrix_markers(C.X,C.layers["log1p"],labels,genes,C.var_names,order)
        for fn,df in [("marker_mean_log1p.tsv",mean),("marker_detection_fraction.tsv",frac),("marker_mean_raw_counts.tsv",raw)]:
            df.to_csv(folder/fn,sep="\t")
        ev=pd.concat([d.rename_axis(index="cluster_uid",columns="gene").stack().rename(k)
            for d,k in [(mean,"mean_log1p_CP10k"),(frac,"raw_detection_fraction"),(raw,"mean_raw_counts")]],axis=1).reset_index()
        ev.to_csv(folder/"marker_evidence_long.tsv",sep="\t",index=False)
        summary=cluster_sources(C,key,folder)
        dotplot(state,mean,frac,run/"figures"/f"{tag}_classic_states_relative",True,cfg["figure_dpi"])
        dotplot(boundary,mean,frac,run/"figures"/f"{tag}_lineage_checks_relative",True,cfg["figure_dpi"])
        umap_plot(C,key,run/"figures"/f"UMAP_{tag}",cfg)
        top,topmap=discover_markers(C,key,tag,run,cfg)
        summary["top_marker_candidates"]=summary.cluster_uid.map(topmap).fillna("")
        summary["resolution"]=r
        summary["marker_anchor"]=""
        summary["state_description"]=""
        summary["proposed_level2_annotation"]=""
        summary["identity_review"]="PENDING"
        summary["review_decision"]="PENDING"
        summary["review_note"]=""
        summary["confirmed_other_lineage_action"]="Unknown in level2_annotation; level1 unchanged"
        summary.to_csv(folder/"cluster_review_template.tsv",sep="\t",index=False)
        allsummary.append(summary)
        if not top.empty:
            # 6 cluster anchors per figure, all cluster rows for specificity comparison.
            uids=[u for u in order if u in topmap]
            for page in range(0,len(uids),6):
                panel={u:top.loc[top.group.astype(str).eq(u),"names"].head(5).tolist() for u in uids[page:page+6]}
                ng=list(dict.fromkeys(g for gs in panel.values() for g in gs))
                nm,nf,_,_=aggregate_matrix_markers(C.X,C.layers["log1p"],labels,ng,C.var_names,order)
                dotplot(panel,nm,nf,run/"figures"/f"{tag}_FindMarkers_dotplot_{page//6+1:02d}",True,cfg["figure_dpi"])
                nm.to_csv(folder/f"FindMarkers_panel_{page//6+1:02d}_mean_log1p.tsv",sep="\t")
                nf.to_csv(folder/f"FindMarkers_panel_{page//6+1:02d}_detection_fraction.tsv",sep="\t")
        if key==baseline:
            baseline_mean=mean;baseline_frac=frac
            dotplot(state,mean,frac,run/"figures/baseline_classic_states_mean_log1p",False,cfg["figure_dpi"])
            dotplot(boundary,mean,frac,run/"figures/baseline_lineage_checks_mean_log1p",False,cfg["figure_dpi"])
        del mean,frac,raw,ev,top;gc.collect()
    pd.concat(allsummary,ignore_index=True).to_csv(run/"tables/all_resolution_review_candidates.tsv",sep="\t",index=False)
    for r in cfg["leiden_resolutions"]:
        tag=res_tag(r)
        if tag!=bx:
            pd.crosstab(C.obs[baseline],C.obs["lineage_L2_"+tag]).to_csv(run/"tables"/f"{bx}_vs_{tag}_membership.tsv",sep="\t")
    # Every partition gets its own sample-stratified marker summaries.
    for r in cfg["leiden_resolutions"]:
        tag=res_tag(r); write_sample_marker_stats(C,"lineage_L2_"+tag,genes,run/"tables"/tag)
    boundary_evidence_all_partitions(C,{cfg["core_panel_key"]:state[cfg["core_panel_key"]],**boundary},run,cfg)
    for key in ("study","tissue","condition","disease_context","L2_total_raw_counts","L2_detected_genes","L2_pct_counts_mt"):
        if key in C.obs: umap_plot(C,key,run/"figures"/("UMAP_"+key),cfg)
    mapping=C.obs[["parent_cell_id","parent_row",L1_KEY,BATCH]+["lineage_L2_"+res_tag(r) for r in cfg["leiden_resolutions"]]].copy()
    mapping.insert(0,"cell_id",C.obs_names)
    mapping.to_csv(run/"tables/cell_cluster_mapping.tsv.gz",sep="\t",index=False)
    # No formal level2_annotation exists before human review.
    require("level2_annotation" not in C.obs,"Unexpected formal L2 labels in Round1.")
    return pd.concat(allsummary,ignore_index=True)


def save_h5ad_verified(C,path,raw_hash):
    import numpy as np
    import h5py
    try: from anndata.io import read_elem
    except ImportError: from anndata._io.specs import read_elem
    path=Path(path); require(not path.exists(),f"Refusing to overwrite {path}")
    tmp=path.with_name(path.stem+".partial.h5ad")
    if tmp.exists(): tmp.rename(path.with_name(path.stem+".incomplete_"+stamp()+".h5ad"))
    kw={"compression":"lzf"}
    if "convert_strings_to_categoricals" in inspect.signature(C.write_h5ad).parameters:
        kw["convert_strings_to_categoricals"]=False
    C.write_h5ad(tmp,**kw)
    with h5py.File(tmp,"r") as f:
        require(h5_sparse_hash(f["X"])==raw_hash,"Saved local raw counts differ from extracted source.")
        if "log1p" in C.layers:
            require(h5_sparse_hash(f["layers/log1p"])==sparse_hash(C.layers["log1p"]),"Saved log1p differs.")
        obs=read_elem(f["obs"]);var=read_elem(f["var"])
        require(obs.index.equals(C.obs_names) and var.index.equals(C.var_names),"Saved cell/gene order differs.")
        require(obs[L1_KEY].astype(str).eq(TARGET).all(),"Saved parent L1 changed.")
        for key in C.obs:
            if key.startswith("lineage_L2_"):
                require(np.array_equal(obs[key].astype(str).to_numpy(),C.obs[key].astype(str).to_numpy()),"Saved partition mismatch.")
        for key in C.obsm:
            require(np.array_equal(np.asarray(read_elem(f["obsm"][key])),np.asarray(C.obsm[key]),equal_nan=True),"Saved coordinate mismatch.")
        for key in C.obsp:
            require(h5_sparse_hash(f["obsp"][key])==sparse_hash(C.obsp[key]),"Saved graph mismatch.")
    os.link(tmp,path);tmp.unlink();readonly(path)
    return sha_file(path)



def worker(run):
    # Called in a fresh CPU-only worker process created by the scheduler.
    global TARGET, EXPECTED_N, UID_PREFIX, STATE_PANEL, BOUNDARY_PANEL
    run=Path(run).resolve();task=json.loads((run/"task.json").read_text());cfg=task["config"];threads=int(task["threads"])
    TARGET=task["lineage"]; EXPECTED_N=int(task["expected_n_cells"]); UID_PREFIX=task["uid_prefix"]
    require(TARGET in EXPECTED_L1 and TARGET not in ("Macrophage","Schwann cell","Unknown"), "Excluded lineage requested.")
    require(EXPECTED_N==EXPECTED_L1[TARGET], "Invalid lineage count in task.")
    STATE_PANEL=task["state_panel"];BOUNDARY_PANEL=task["boundary_panel"]
    cfg["core_panel_key"]=task["core_panel_key"]
    require(sha_file(Path(__file__).resolve())==task["core_sha256"], "Worker code snapshot changed.")
    require(sha_file(run/"audit/patched_function.py")==PATCH_SHA256, "HVG patch changed.")
    check_source_snapshot(task)

    require(os.environ.get("CUDA_VISIBLE_DEVICES")=="","Worker CUDA is not hidden.")
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import omicverse as ov
    import torch
    import numba
    import matplotlib
    from threadpoolctl import threadpool_limits,threadpool_info
    matplotlib.use("Agg",force=True)
    require(sc.__version__=="1.11.5","Use the existing Scanpy 1.11.5 environment; do not silently change the HVG implementation.")
    ov.settings.cpu_init();sc.settings.n_jobs=threads;sc.settings.autoshow=False;sc.settings.autosave=False
    torch.set_num_threads(threads)
    try:torch.set_num_interop_threads(1)
    except RuntimeError:pass
    require(not torch.cuda.is_available(),"CUDA still visible; abort CPU-only run.")
    numba.set_num_threads(min(threads,numba.config.NUMBA_NUM_THREADS))
    np.random.seed(cfg["seed"]);torch.manual_seed(cfg["seed"])
    start=time.time()
    def progress(s):
        print(f"[{datetime.now().isoformat(timespec='seconds')}] {s}",flush=True)
        json_write(run/"STATUS.json",{"status":s,"time":datetime.now().isoformat(),"L1_modified":False,"level2_frozen":False})
    versions={}
    for p in ("omicverse","scanpy","anndata","numpy","scipy","pandas","scikit-learn","torch","numba","igraph","pynndescent","umap-learn","h5py"):
        try:versions[p]=metadata.version(p)
        except metadata.PackageNotFoundError:versions[p]="unavailable"
    json_write(run/"audit/environment.json",{"versions":versions,"threads":threads,"CUDA_VISIBLE_DEVICES":"", "torch_cuda_available":False,
        "torch_threads":torch.get_num_threads(),"numba_threads":numba.get_num_threads(),"threadpools":threadpool_info()})
    checkpoint=run/"data/lineage_L2_embedding_checkpoint.h5ad"
    ckmeta=run/"audit/embedding_checkpoint.json"
    if checkpoint.exists() and ckmeta.exists():
        progress("VERIFY_AND_REUSE_LOCAL_EMBEDDING_CHECKPOINT")
        meta=json.loads(ckmeta.read_text());integrity=json.loads((run/"audit/input_integrity.json").read_text())
        require(meta["sha256"]==sha_file(checkpoint),"Embedding checkpoint SHA mismatch.")
        check_source_snapshot(task)
        C=sc.read_h5ad(checkpoint)
        require(C.n_obs==EXPECTED_N and C.obs[L1_KEY].astype(str).eq(TARGET).all(),"Wrong checkpoint cells.")
        cfg=json.loads(C.uns["round1_config_json"])
        require(sparse_hash(C.X)==integrity["selected_raw_sha256"],"Checkpoint raw mismatch.")
        lock_sparse(C.X)
    else:
        progress("READ_ONLY_FROZEN_INPUT_AND_LINEAGE_EXTRACTION")
        rawpath=run/"data/lineage_raw_input.h5ad"
        rawmeta=run/"audit/raw_checkpoint.json"
        if rawpath.exists() and rawmeta.exists():
            meta=json.loads(rawmeta.read_text()); require(sha_file(rawpath)==meta["sha256"],"Raw checkpoint changed.")
            integrity=json.loads((run/"audit/input_integrity.json").read_text())
            C=sc.read_h5ad(rawpath)
            require(C.n_obs==EXPECTED_N and C.obs[L1_KEY].astype(str).eq(TARGET).all(),"Wrong raw checkpoint lineage.")
            require(sparse_hash(C.X)==integrity["selected_raw_sha256"],"Raw checkpoint counts changed.")
            lock_sparse(C.X)
        else:
            C,integrity=read_source(task,run)
            if rawpath.exists():
                verify_existing_raw_checkpoint(rawpath,C,integrity)
            else:
                save_h5ad_verified(C,rawpath,integrity["selected_raw_sha256"])
            json_write(rawmeta,{"sha256":sha_file(rawpath)})
        with threadpool_limits(limits=threads):
            C=make_embedding(C,run,cfg,threads,progress)
        if checkpoint.exists():
            checkpoint.rename(run/"data"/("unverified_embedding_"+stamp()+".h5ad"))
        progress("SAVE_LOCAL_EMBEDDING_CHECKPOINT")
        ckhash=save_h5ad_verified(C,checkpoint,integrity["selected_raw_sha256"])
        json_write(ckmeta,{"sha256":ckhash,"input_source_sha256":integrity["source_sha256"]})
    # Gene counts for testing: canonicalize a temporary log1p copy, never the raw matrix.
    C.var["n_cells_detected_lineage"]=np.asarray(C.layers["log1p"].getnnz(axis=0)).ravel()
    frozen_labels=C.obs[L1_KEY].astype(str).copy();frozen_ids=C.obs_names.copy()
    with threadpool_limits(limits=threads):
        summary=write_reviews(C,run,cfg,progress)
    require(C.obs_names.equals(frozen_ids) and C.obs[L1_KEY].astype(str).equals(frozen_labels),"L1 labels/order changed.")
    require(sparse_hash(C.X)==integrity["selected_raw_sha256"],"Raw counts changed.")
    C.uns["round1_config_json"]=json.dumps(cfg,ensure_ascii=False,default=str)
    C.uns["level1_frozen"]=True
    C.uns["L2_started"]=True
    C.uns["expression_layout"]={"X":"unchanged full-gene raw counts of selected lineage cells",
        "layers_log1p":"CP10k and log1p recomputed from a separate raw-count copy",
        "parent_reference_umap":"display reference only, not used to build local graph"}
    C.uns["level2_frozen"]=False
    C.uns["level2_annotation_applied"]=False
    C.uns["status"]="PAUSED_FOR_LINEAGE_L2_ROUND1_REVIEW"
    C.uns["level2_other_lineage_policy"]="If confirmed by review: level2_annotation=Unknown; keep level1_annotation unchanged. Shared programs or isolated markers are not sufficient. No automatic label changes."
    C.uns["naming_style"]="marker anchor + supported state/program + parent cell type; no fixed state count"
    C.uns["source_raw_unchanged_sha256"]=integrity["selected_raw_sha256"]
    final=run/"data/lineage_L2_round1_review.h5ad"
    if final.exists():
        # A prior partially-finished worker must never overwrite a published review file.
        final.rename(run/"data"/("previous_review_"+stamp()+".h5ad"))
    progress("SAVE_FULL_GENE_LINEAGE_REVIEW_OBJECT")
    finalhash=save_h5ad_verified(C,final,integrity["selected_raw_sha256"])
    progress("VERIFY_FROZEN_MOTHER_STILL_UNCHANGED")
    check_source_snapshot(task)
    offline_report(run,summary,cfg,integrity)
    result={"status":"PAUSED_FOR_LINEAGE_L2_ROUND1_REVIEW","lineage":TARGET,"uid_prefix":UID_PREFIX,"analysis_file":str(final),"analysis_sha256":finalhash,
        "n_cells":C.n_obs,"n_genes":C.n_vars,"level1_key":L1_KEY,"level1_value":TARGET,
        "level1_modified":False,"level2_annotation_applied":False,"level2_frozen":False,
        "partition_counts":{res_tag(r):int(C.obs['lineage_L2_'+res_tag(r)].nunique()) for r in cfg['leiden_resolutions']},
        "source_sha256":integrity["source_sha256"],"source_stat_unchanged":True,"batch_full_source_hash_check":"performed_by_scheduler","elapsed_seconds":time.time()-start,
        "baseline_for_review_only":res_tag(cfg["baseline_resolution"]),"actual_harmony_class":cfg.get("actual_harmony_class"),
        "threads":threads,"unknown_policy":"review-confirmed off-lineage or unresolved identity -> L2 Unknown; never modify L1"}
    json_write(run/"RESULTS.json",result)
    (run/"PROJECT_STATE_LINEAGE_L2_ROUND1.md").write_text(
        f"# New {TARGET} L2 Round1 branch\n\n"+json.dumps(result,ensure_ascii=False,indent=2)+
        "\n\nFrozen mother/old PROJECT_STATE are not edited. No L2 labels applied.\n",encoding="utf-8")
    progress("PACKAGING_REVIEW_ONLY_NO_EXPRESSION_H5AD")
    archive=run/"lineage_L2_round1_review_materials.zip"
    # Results, all tables/figures, code and patch; exclude full data h5ad, stale temp files, active log.
    members=[]
    for p in run.rglob("*"):
        if p.is_file() and p!=archive and "data" not in p.relative_to(run).parts and p.name not in ("run.log","worker.log","STATUS.json","review_checksums.json") and not p.name.endswith(".tmp"):
            members.append(p)
    checksum={str(p.relative_to(run)):sha_file(p) for p in members}
    json_write(run/"review_checksums.json",checksum);members.append(run/"review_checksums.json")
    with zipfile.ZipFile(archive,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=4) as z:
        for p in members:z.write(p,p.relative_to(run))
    with zipfile.ZipFile(archive) as z:require(z.testzip() is None,"Review ZIP CRC error.")
    json_write(run/"STATUS.json",result)
    print("\n首轮完成，未修改母图谱或写入正式L2。",flush=True)
    print("审核包：",archive,flush=True)
    print("分析对象：",final,flush=True)
    print("各分辨率计算簇数（不是已确认细胞亚型数）：",result["partition_counts"],flush=True)



def source_stat(path):
    s=Path(path).stat()
    return {"size":s.st_size,"mtime_ns":s.st_mtime_ns,"inode":s.st_ino,"device":s.st_dev}


def check_source_snapshot(task):
    require(source_stat(task["source"])==task["source_stat"], "Frozen source file changed during batch; stop.")
    require(sha_file(task["source_manifest"])==task["manifest_sha256"], "Frozen source manifest changed.")
    ticket=Path(task["batch_verification"])
    require(sha_file(ticket)==task["batch_verification_sha256"], "Batch verification ticket changed.")
    info=json.loads(ticket.read_text())
    require(info["source_sha256"]==task["verified_source_sha256"] and info["source_stat"]==task["source_stat"],
            "Invalid verified source ticket.")


def verify_existing_raw_checkpoint(path,C,integrity):
    import h5py
    try: from anndata.io import read_elem
    except ImportError: from anndata._io.specs import read_elem
    with h5py.File(path,"r") as f:
        require(h5_sparse_hash(f["X"])==integrity["selected_raw_sha256"],"Old raw checkpoint differs.")
        o=read_elem(f["obs"]);v=read_elem(f["var"])
        require(o.index.equals(C.obs_names) and v.index.equals(C.var_names),"Raw checkpoint identities differ.")
        require(o[L1_KEY].astype(str).eq(TARGET).all(),"Raw checkpoint has wrong L1.")


def write_sample_marker_stats(C,key,genes,folder):
    import numpy as np
    import pandas as pd
    codes,levels=pd.factorize(pd.MultiIndex.from_arrays([C.obs[key].astype(str),C.obs[BATCH].astype(str)]),sort=True)
    ids=np.asarray([f"CS_{i:05d}" for i in codes])
    cm,cf,_,cn=aggregate_matrix_markers(C.X,C.layers["log1p"],ids,genes,C.var_names)
    cm.to_csv(folder/"cluster_sample_marker_mean_log1p.tsv.gz",sep="\t")
    cf.to_csv(folder/"cluster_sample_marker_detection_fraction.tsv.gz",sep="\t")
    rows=[{"group_id":f"CS_{i:05d}","cluster_uid":str(u),"sample_uid":str(s),"n_cells":int(cn[f"CS_{i:05d}"])}
          for i,(u,s) in enumerate(levels)]
    pd.DataFrame(rows).to_csv(folder/"cluster_sample_groups.tsv",sep="\t",index=False)


def boundary_evidence_all_partitions(C,panels,run,cfg):
    import numpy as np
    import pandas as pd
    records={};avail=set(C.var_names)
    for name,genes in panels.items():
        gs=list(dict.fromkeys(g for g in genes if g in avail))
        if not gs:continue
        x=C.X[:,C.var_names.get_indexer(gs)].tocsr(copy=True)
        x.sum_duplicates();x.eliminate_zeros()
        for threshold in (1,2):
            records[name+f"__detected_genes_ge{threshold}"]=np.asarray((x>=threshold).sum(axis=1)).ravel().astype(np.uint8)
    percell=pd.DataFrame(records,index=C.obs_names)
    for r in cfg["leiden_resolutions"]:
        key="lineage_L2_"+res_tag(r);percell.insert(0,key,C.obs[key].astype(str).to_numpy())
    percell.to_csv(run/"tables/boundary_multi_marker_counts_per_cell.tsv.gz",sep="\t",index_label="cell_id")
    rows=[]
    for r in cfg["leiden_resolutions"]:
        key="lineage_L2_"+res_tag(r)
        for uid,ix in C.obs.groupby(key,observed=True,sort=False).indices.items():
            for name in panels:
                for t in (1,2):
                    col=name+f"__detected_genes_ge{t}"
                    if col not in percell:continue
                    vals=percell[col].iloc[ix].to_numpy()
                    rows.append({"resolution":r,"cluster_uid":str(uid),"panel":name,"raw_count_cutoff":t,
                        "n_cells":len(ix),"mean_n_detected_panel_genes":float(vals.mean()),
                        "fraction_at_least_2_panel_genes":float((vals>=2).mean()),
                        "interpretation":"DESCRIPTIVE_ONLY_NOT_AN_IDENTITY_CLASSIFIER"})
    pd.DataFrame(rows).to_csv(run/"tables/boundary_multi_marker_by_cluster.tsv",sep="\t",index=False)


def offline_report(run,summary,cfg,integrity):
    counts=summary.groupby("resolution").size().to_dict()
    notes=f'''# {TARGET} L2 CPU Round1 — 候选审核
输入：冻结 level1_annotation == {TARGET}; {EXPECTED_N:,} 个细胞，全部基因保留。
母图谱和巨噬细胞分支未修改；未应用正式level2_annotation。
各分辨率簇数：{counts}；基线{cfg['baseline_resolution']}仅为审核锚点。

1. 每个分辨率均使用自身FindMarkers与自身sample分层marker表；不得混用编号。
2. 命名风格：特征marker + 有证据支持的状态描述 + 细胞类型。名称/亚群数不预设。
3. 经典marker和FindMarkers共同支持人工审核，不以最高module score自动分配。
4. 经审核确认属于其他谱系或身份不明时，只把新的level2_annotation记为Unknown；冻结L1不动。
5. 单个标记、HLA、IFN、应激、细胞毒或ECM共享程序本身不等于另一细胞类型。
6. 少于2细胞的sample只不参与HVG方差估计，细胞仍保留在后续分析。
7. 未重新运行Scrublet，无细胞删除，无疾病推断、CoVarNet、CellChat或轨迹。

流程：raw工作副本→Pearson residual HVG {cfg['n_hvg']}，逐sample clip=sqrt(n)
→全基因CP10k+log1p→临时HVG矩阵scale（max=10）→PCA {cfg['n_pcs']}
→CPU Harmony(theta=2)→neighbors 30→Leiden 0.4/0.8/1.2→UMAP。
FindMarkers：全基因中至少3细胞检出的基因，Scanpy Wilcoxon，layer=log1p，use_raw=False。
FDR/效应筛选仅用于画特征marker图；不改聚类或L1，不用来自动删小群。
DotPlot：含零细胞的平均log1p；点大小=raw counts>0；相对图是每基因内0–1缩放。
临床字段、donor缺失、NR及PROVISIONAL均保留；覆盖样本数不等于独立患者数。
L1筛选后同一细胞的全基因log1p与母图谱本应可一致，本程序按同一公式在工作副本重算。
单study状态不自动删除，也不能宣称跨研究稳定；采样/疾病不能指导聚类。
所有status均为待审核，不是已验证的正式L2亚型。
冻结输入SHA256: {integrity['source_sha256']}
'''
    (run/"REVIEW_CN.md").write_text(notes,encoding="utf-8")
    cols=[c for c in ["resolution","cluster_uid","n_cells","n_sample_uid","n_study","n_donor_uid_v2",
                       "top_study","top_study_fraction","top_marker_candidates"] if c in summary]
    docs=["<!doctype html><meta charset='utf-8'><style>body{font:15px sans-serif;margin:24px}img{max-width:100%}pre{white-space:pre-wrap}td,th{padding:5px;border:1px solid #ccc}</style>",
          "<pre>"+html.escape(notes)+"</pre>",summary[cols].to_html(index=False,escape=True)]
    for p in sorted((run/"figures").glob("*.png")):
        docs += ["<h2>"+html.escape(p.stem)+"</h2>","<img src='data:image/png;base64,"+base64.b64encode(p.read_bytes()).decode()+"'>"]
    (run/"LINEAGE_L2_ROUND1_OFFLINE.html").write_text("\n".join(docs),encoding="utf-8")


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--worker",required=True)
    args=parser.parse_args()
    try:
        worker(args.worker)
    except BaseException as exc:
        json_write(Path(args.worker)/"STATUS.json",{"status":"FAILED_NOT_FROZEN","error":repr(exc),
            "traceback":traceback.format_exc(),"time":datetime.now().isoformat(),
            "L1_modified":False,"macrophage_modified":False})
        traceback.print_exc();sys.exit(1)
