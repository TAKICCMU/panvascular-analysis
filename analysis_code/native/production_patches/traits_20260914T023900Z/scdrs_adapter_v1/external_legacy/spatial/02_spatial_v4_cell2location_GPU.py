from __future__ import annotations

import gc
import json
import os
import re
import warnings
import platform
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import torch
from scipy.io import mmread

from cell2location.models import Cell2location, RegressionModel
from cell2location.utils.filtering import filter_genes


ROOT = Path(os.environ.get(
    "SPATIAL_V4_ROOT",
    "/home/anzhen_vhd/scRNAseq/codex/spatial_validation_v4_gpu_final",
))
SOURCE_ROOT = Path(os.environ.get(
    "SPATIAL_V3_SOURCE_ROOT",
    "/home/anzhen_vhd/scRNAseq/codex/spatial_validation_v3",
))
ATLAS_PATH = Path(
    "/home/anzhen_vhd/scRNAseq/codex/Level2_celltype_final_v1_metadata_v2.h5ad"
)
MANIFEST_PATH = SOURCE_ROOT / "00_inputs" / "spatial_sample_manifest_curated.tsv"
REF_OUT = ROOT / "01_reference"
MAP_OUT = ROOT / "02_cell2location"
AUDIT_OUT = ROOT / "00_audit"
for directory in (ROOT, REF_OUT, MAP_OUT, AUDIT_OUT):
    directory.mkdir(parents=True, exist_ok=True)

USE_GPU = True
RANDOM_SEED = int(os.environ.get("SPATIAL_RANDOM_SEED", "20260814"))
MAX_CELLS_PER_STATE = int(os.environ.get("SPATIAL_MAX_CELLS_PER_STATE", "1500"))
MIN_CELLS_PER_STATE = int(os.environ.get("SPATIAL_MIN_CELLS_PER_STATE", "40"))
MAX_STUDIES_PER_STATE = int(os.environ.get("SPATIAL_MAX_STUDIES_PER_STATE", "20"))
REF_MAX_EPOCHS = int(os.environ.get("SPATIAL_REF_EPOCHS", "250"))
SPATIAL_MAX_EPOCHS = int(os.environ.get("SPATIAL_MAPPING_EPOCHS", "30000"))
N_CELLS_PER_LOCATION = int(os.environ.get("SPATIAL_N_CELLS_PER_LOCATION", "12"))
DETECTION_ALPHA = float(os.environ.get("SPATIAL_DETECTION_ALPHA", "20"))
POSTERIOR_SAMPLES = int(os.environ.get("SPATIAL_POSTERIOR_SAMPLES", "1000"))

np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(RANDOM_SEED)
try:
    import scvi
    scvi.settings.seed = RANDOM_SEED
except Exception:
    scvi = None


def trainer_device_kwargs() -> dict:
    """Return Lightning 2/scvi-tools device arguments."""
    if USE_GPU:
        if not torch.cuda.is_available():
            raise RuntimeError(
                "USE_GPU=True but torch.cuda.is_available() is False. "
                "Use the GPU-enabled Jupyter kernel or set USE_GPU=False."
            )
        torch.set_float32_matmul_precision("high")
        return {"accelerator": "gpu"}
    return {"accelerator": "cpu"}

FORCED_GENES = {
    "KCNK3", "CKB", "ZNF589", "SCARB1", "LRRC10B",
    "COL8A1", "ELN", "LOX", "LOXL1", "DCN", "BGN", "LUM", "THBS2",
    "TIMP1", "SERPINH1", "RUNX2", "ALPL", "IBSP", "BGLAP", "SPP1",
    "MGP", "BMP2", "ENPP1", "TREM2", "LPL", "APOE", "APOC1", "GPNMB",
    "LGALS3", "CTSB", "CD68", "IL1B", "TNF", "CCL2", "CXCL8", "NFKBIA",
    "S100A8", "S100A9", "CTSS", "VCAM1", "SELE", "ICAM1", "VWF", "KDR",
    "CXCL2",
}


def save_history(model, stem: str) -> None:
    """Save ELBO/training histories when supported by the installed version."""
    history = getattr(model, "history", None)
    if history is None:
        return
    rows = []
    if isinstance(history, dict):
        for metric, values in history.items():
            try:
                array = np.asarray(values).ravel()
            except Exception:
                continue
            for epoch, value in enumerate(array):
                rows.append({"metric": str(metric), "epoch": epoch, "value": value})
    elif hasattr(history, "items"):
        for metric, values in history.items():
            array = np.asarray(values).ravel()
            for epoch, value in enumerate(array):
                rows.append({"metric": str(metric), "epoch": epoch, "value": value})
    if rows:
        pd.DataFrame(rows).to_csv(AUDIT_OUT / f"{stem}__training_history.tsv.gz", sep="\t", index=False, compression="gzip")


def gpu_audit() -> dict:
    audit = {
        "python": sys.version,
        "platform": platform.platform(),
        "torch": torch.__version__,
        "cuda_available": bool(torch.cuda.is_available()),
        "random_seed": RANDOM_SEED,
        "reference_max_epochs": REF_MAX_EPOCHS,
        "spatial_max_epochs": SPATIAL_MAX_EPOCHS,
        "posterior_samples": POSTERIOR_SAMPLES,
        "gene_cap": None,
    }
    if torch.cuda.is_available():
        free, total = torch.cuda.mem_get_info()
        audit.update({
            "gpu_name": torch.cuda.get_device_name(0),
            "gpu_free_GiB_at_start": free / 1024**3,
            "gpu_total_GiB": total / 1024**3,
        })
        minimum_free = float(os.environ.get("SPATIAL_MIN_FREE_GPU_GB", "10"))
        if free / 1024**3 < minimum_free:
            raise RuntimeError(
                f"Only {free / 1024**3:.1f} GiB GPU memory is free; "
                f"at least {minimum_free:.1f} GiB was requested. Stop the other GPU job "
                "or lower SPATIAL_MIN_FREE_GPU_GB only after checking memory use."
            )
    (AUDIT_OUT / "GPU_and_parameter_audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(audit, indent=2, ensure_ascii=False))
    return audit


def clean_label(value: str) -> str:
    value = re.sub(r"[^0-9A-Za-z_]+", "_", str(value).strip())
    return re.sub(r"_+", "_", value).strip("_")


def make_unique(values) -> list[str]:
    counts: dict[str, int] = {}
    output = []
    for value in map(str, values):
        count = counts.get(value, 0)
        output.append(value if count == 0 else f"{value}__dup{count}")
        counts[value] = count + 1
    return output


def normalize_gene_names(adata: ad.AnnData) -> ad.AnnData:
    candidates = ["gene_symbol_std", "gene_symbol", "gene_symbols"]
    names = None
    for column in candidates:
        if column in adata.var.columns:
            values = adata.var[column].astype(str)
            if values.notna().mean() > 0.8:
                names = values
                break
    if names is None:
        names = pd.Series(adata.var_names, index=adata.var_names, dtype=str)
    adata.var_names = pd.Index(make_unique(names.str.upper()))
    return adata


def valid_string(series: pd.Series) -> pd.Series:
    values = series.astype("object").astype(str).str.strip()
    return ~values.str.lower().isin({"", "nan", "none", "na", "unknown", "unassigned"})


def select_state_key(obs: pd.DataFrame) -> str:
    candidates = [
        "Level2_CoVarNet_node_v1",
        "Level2_celltype_for_enrichment_v1",
        "Level2_celltype_final_v1",
    ]
    key = next((x for x in candidates if x in obs.columns), None)
    if key is None:
        raise KeyError(f"No final Level2 field found. Tried: {candidates}")
    return key


def select_level1_key(obs: pd.DataFrame) -> str:
    candidates = ["Level1_celltype_final_v3", "Level1_celltype_final"]
    key = next((x for x in candidates if x in obs.columns), None)
    if key is None:
        raise KeyError(f"No final Level1 field found. Tried: {candidates}")
    return key


def select_study_key(obs: pd.DataFrame) -> str:
    candidates = ["study_canonical", "study"]
    key = next((x for x in candidates if x in obs.columns), None)
    if key is None:
        raise KeyError(f"No study field found. Tried: {candidates}")
    return key


def build_spatial_reference_labels(obs: pd.DataFrame) -> tuple[pd.Series, dict]:
    state_key = select_state_key(obs)
    level1_key = select_level1_key(obs)
    state = obs[state_key].astype("object").astype(str).str.strip()
    level1 = obs[level1_key].astype("object").astype(str).str.strip()
    keep_state = valid_string(obs[state_key])
    labels = state.where(keep_state, "Other | " + level1)
    labels = labels.where(valid_string(labels), "Other | Unknown")
    return labels, {"state_key": state_key, "level1_key": level1_key}


def balanced_indices(
    obs: pd.DataFrame,
    label_key: str,
    study_key: str,
    maximum: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    selected: list[np.ndarray] = []
    for _, state_idx in obs.groupby(label_key, observed=True).indices.items():
        state_idx = np.asarray(state_idx, dtype=int)
        if len(state_idx) <= maximum:
            selected.append(state_idx)
            continue
        state_studies = obs.iloc[state_idx][study_key].astype(str)
        study_names = state_studies.value_counts().index[:MAX_STUDIES_PER_STATE]
        quota = max(1, maximum // max(1, len(study_names)))
        chosen = []
        for study in study_names:
            local = state_idx[state_studies.to_numpy() == study]
            chosen.extend(rng.choice(local, min(quota, len(local)), replace=False).tolist())
        if len(chosen) < maximum:
            remaining = np.setdiff1d(state_idx, np.asarray(chosen), assume_unique=False)
            chosen.extend(
                rng.choice(remaining, min(maximum - len(chosen), len(remaining)), replace=False).tolist()
            )
        selected.append(np.asarray(chosen[:maximum], dtype=int))
    return np.sort(np.concatenate(selected))


def extract_reference() -> tuple[ad.AnnData, str]:
    atlas = ad.read_h5ad(ATLAS_PATH, backed="r")
    labels, keys = build_spatial_reference_labels(atlas.obs)
    study_key = select_study_key(atlas.obs)
    obs = pd.DataFrame(index=atlas.obs_names)
    obs["spatial_reference_state"] = labels.to_numpy()
    obs["reference_study"] = (
        atlas.obs[study_key].astype("object").fillna("missing_study").astype(str).to_numpy()
    )
    counts = obs["spatial_reference_state"].value_counts()
    allowed = counts[counts >= MIN_CELLS_PER_STATE].index
    keep = obs["spatial_reference_state"].isin(allowed).to_numpy()
    eligible = np.flatnonzero(keep)
    local_obs = obs.iloc[eligible].reset_index(drop=True)
    chosen_local = balanced_indices(
        local_obs, "spatial_reference_state", "reference_study",
        MAX_CELLS_PER_STATE, RANDOM_SEED,
    )
    chosen = eligible[chosen_local]

    # Older anndata versions keep backed Raw slicing lazy. Calling
    # Raw.to_adata() then incorrectly carries full-atlas obsm arrays into the
    # cell subset. Read only raw counts in bounded row blocks and construct a
    # fresh AnnData with no inherited obsm/uns payload.
    raw_var = atlas.raw.var.copy()
    raw_var_names = atlas.raw.var_names.copy()
    count_blocks = []
    block_size = 2500
    for start in range(0, len(chosen), block_size):
        stop = min(start + block_size, len(chosen))
        block_rows = np.asarray(chosen[start:stop], dtype=int)
        block = atlas.raw.X[block_rows, :]
        if not sp.issparse(block):
            block = sp.csr_matrix(block)
        count_blocks.append(block.tocsr().astype(np.float32))
        if start == 0 or stop == len(chosen) or (start // block_size) % 10 == 0:
            print(f"Reference raw extraction: {stop:,}/{len(chosen):,}")
    raw_counts = sp.vstack(count_blocks, format="csr")
    del count_blocks
    atlas.file.close()

    reference_obs = obs.iloc[chosen].copy()
    reference_obs.index = pd.Index(obs.index[chosen])
    raw_var.index = raw_var_names
    ref = ad.AnnData(X=raw_counts, obs=reference_obs, var=raw_var)
    ref = normalize_gene_names(ref)
    ref.obs["spatial_reference_state_original"] = ref.obs["spatial_reference_state"].astype(str)
    ref.obs["spatial_reference_state_safe"] = ref.obs["spatial_reference_state_original"].map(clean_label)

    label_map = (
        ref.obs[["spatial_reference_state_original", "spatial_reference_state_safe"]]
        .drop_duplicates()
        .sort_values("spatial_reference_state_original")
    )
    duplicated_safe = label_map["spatial_reference_state_safe"].duplicated(keep=False)
    if duplicated_safe.any():
        bad = label_map.loc[duplicated_safe]
        raise ValueError(f"Safe label collision detected:\n{bad.to_string(index=False)}")
    label_map.to_csv(REF_OUT / "spatial_reference_state_label_map.tsv", sep="\t", index=False)
    ref.obs["reference_study"] = ref.obs["reference_study"].astype(str)

    selected_gene_result = filter_genes(
        ref,
        cell_count_cutoff=5,
        cell_percentage_cutoff2=0.03,
        nonz_mean_cutoff=1.12,
    )
    result = np.asarray(
        selected_gene_result.to_numpy()
        if isinstance(selected_gene_result, (pd.Index, pd.Series))
        else selected_gene_result
    )
    if result.dtype == bool:
        selected_genes = set(ref.var_names[result])
    elif np.issubdtype(result.dtype, np.integer):
        selected_genes = set(ref.var_names[result])
    else:
        selected_genes = {str(x).upper() for x in result if str(x).upper() in ref.var_names}
    selected_genes.update(gene for gene in FORCED_GENES if gene in ref.var_names)
    ref = ref[:, sorted(selected_genes)].copy()
    if not sp.issparse(ref.X):
        ref.X = sp.csr_matrix(ref.X)
    ref.X = ref.X.astype(np.float32)

    audit = {
        **keys,
        "study_key": study_key,
        "n_cells": int(ref.n_obs),
        "n_genes": int(ref.n_vars),
        "n_states": int(ref.obs["spatial_reference_state_safe"].nunique()),
        "forced_genes_present": sorted(FORCED_GENES.intersection(ref.var_names)),
        "gene_selection_rule": "cell2location.filter_genes plus forced target/pathology genes; no fixed maximum-gene cap",
    }
    (REF_OUT / "reference_audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    ref.write_h5ad(REF_OUT / "final_atlas_reference_downsampled_counts.h5ad", compression="gzip")
    return ref, "reference_study"


def fit_or_load_reference() -> pd.DataFrame:
    signature_path = REF_OUT / "final_atlas_reference_signatures.tsv.gz"
    if signature_path.exists() and signature_path.stat().st_size > 1000:
        print("Reusing reference signatures:", signature_path)
        return pd.read_csv(signature_path, sep="\t", index_col=0)

    ref, batch_key = extract_reference()
    RegressionModel.setup_anndata(
        adata=ref,
        batch_key=batch_key,
        labels_key="spatial_reference_state_safe",
    )
    model = RegressionModel(ref)
    model.train(
        max_epochs=REF_MAX_EPOCHS,
        **trainer_device_kwargs(),
    )
    save_history(model, "reference_regression")
    ref = model.export_posterior(
        ref,
        sample_kwargs={"num_samples": POSTERIOR_SAMPLES, "batch_size": 2500},
    )
    model.save(str(REF_OUT / "RegressionModel_final_atlas"), overwrite=True)
    candidates = [key for key in ref.varm.keys() if "means_per_cluster" in key]
    if not candidates:
        raise KeyError(f"No cell2location signature matrix in ref.varm: {list(ref.varm.keys())}")
    signatures = pd.DataFrame(ref.varm[candidates[0]], index=ref.var_names)
    factor_names = ref.uns.get("mod", {}).get("factor_names")
    if factor_names is not None and len(factor_names) == signatures.shape[1]:
        signatures.columns = list(map(str, factor_names))
    signatures.to_csv(signature_path, sep="\t", compression="gzip")
    ref.write_h5ad(REF_OUT / "reference_with_posterior.h5ad", compression="gzip")
    del model, ref
    gc.collect()
    return signatures


def read_features(path: Path) -> pd.Index:
    table = pd.read_csv(path, sep="\t", header=None)
    names = table.iloc[:, 1] if table.shape[1] > 1 else table.iloc[:, 0]
    return pd.Index(make_unique(names.astype(str).str.upper()), name=None)


def read_positions(path: Path) -> pd.DataFrame:
    first = pd.read_csv(path, header=None, nrows=1)
    if str(first.iloc[0, 0]).lower() in {"barcode", "barcodes"}:
        table = pd.read_csv(path)
        table = table.iloc[:, :6]
    else:
        table = pd.read_csv(path, header=None).iloc[:, :6]
    table.columns = [
        "barcode", "in_tissue", "array_row", "array_col",
        "pxl_row_in_fullres", "pxl_col_in_fullres",
    ]
    table["barcode"] = table["barcode"].astype(str)
    return table


def read_one_section(row: pd.Series) -> ad.AnnData:
    matrix = mmread(row["matrix"]).tocsr().T.astype(np.float32)
    barcodes = (
        pd.read_csv(row["barcodes"], sep="\t", header=None)
        .iloc[:, 0]
        .astype(str)
        .to_numpy()
    )
    section = ad.AnnData(X=matrix)
    section.obs_names = pd.Index(barcodes, name=None)
    section.var_names = read_features(Path(row["features"]))
    section.var_names_make_unique()
    section.obs.index.name = None
    section.var.index.name = None

    positions = read_positions(Path(row["positions"])).set_index("barcode")
    common = section.obs_names.intersection(positions.index)
    if len(common) == 0:
        stripped_obs = pd.Index(section.obs_names.str.replace(r"-1$", "", regex=True))
        stripped_pos = positions.copy()
        stripped_pos.index = stripped_pos.index.str.replace(r"-1$", "", regex=True)
        overlap = stripped_obs.intersection(stripped_pos.index)
        old_by_stripped = dict(zip(stripped_obs, section.obs_names))
        common = pd.Index([old_by_stripped[x] for x in overlap])
        positions = stripped_pos.loc[overlap].copy()
        positions.index = common
    section = section[common].copy()
    positions = positions.loc[common]
    if "in_tissue" in positions:
        in_tissue = positions["in_tissue"].astype(int).to_numpy() == 1
        section = section[in_tissue].copy()
        positions = positions.iloc[np.flatnonzero(in_tissue)]

    for column in positions.columns:
        section.obs[column] = positions[column].to_numpy()
    section.obsm["spatial"] = positions[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(float)
    for column in [
        "GSE", "sample_id", "section_id", "biological_unit_id", "metadata_status",
        "condition", "tissue", "dataset_role",
    ]:
        section.obs[column] = str(row[column])
    section.obs["barcode_original"] = section.obs_names.astype(str)
    section.obs_names = pd.Index(
        [f"{row['section_id']}::{x}" for x in section.obs_names],
        name=None,
    )
    section.obs.index.name = None
    section.var.index.name = None
    section.layers["counts"] = section.X.copy()
    section.uns["source_files"] = {key: str(row[key]) for key in ["matrix", "barcodes", "features", "positions"]}
    return section


def load_one_gse(gse: str, manifest: pd.DataFrame) -> ad.AnnData:
    sections = []
    for _, row in manifest[(manifest["GSE"] == gse) & manifest["include_mapping"]].iterrows():
        section = read_one_section(row)
        minimum_counts = 500 if section.n_vars > 10000 else 100
        minimum_genes = 200 if section.n_vars > 10000 else 50
        counts = np.asarray(section.X.sum(axis=1)).ravel()
        genes = np.asarray((section.X > 0).sum(axis=1)).ravel()
        keep = (counts >= minimum_counts) & (genes >= minimum_genes)
        if keep.sum() < 50:
            keep = counts > 0
        section = section[keep].copy()
        sections.append(section)
        print("Loaded", row["section_id"], section.shape)
    if not sections:
        raise ValueError(f"No included sections for {gse}")
    result = sc.concat(sections, join="outer", fill_value=0, merge="same")
    result.var_names_make_unique()
    result.obs.index.name = None
    result.var.index.name = None
    result.layers["counts"] = result.X.copy()
    return result


def map_one_gse(gse: str, manifest: pd.DataFrame, signatures: pd.DataFrame) -> Path:
    output = MAP_OUT / f"{gse}__final_atlas_cell2location.h5ad"
    if output.exists() and output.stat().st_size > 10000:
        print("Reusing mapped output:", output)
        return output
    spatial = load_one_gse(gse, manifest)
    spatial.obs["total_counts_original"] = np.asarray(spatial.layers["counts"].sum(axis=1)).ravel()
    spatial.obs["n_genes_original"] = np.asarray((spatial.layers["counts"] > 0).sum(axis=1)).ravel()
    common = spatial.var_names.intersection(signatures.index)
    if len(common) < 1000:
        raise ValueError(f"Only {len(common)} common genes for {gse}")
    spatial = spatial[:, common].copy()
    spatial.obs.index.name = None
    spatial.var.index.name = None
    signature_common = signatures.loc[common].copy()

    Cell2location.setup_anndata(adata=spatial, batch_key="sample_id")
    model = Cell2location(
        spatial,
        cell_state_df=signature_common,
        N_cells_per_location=N_CELLS_PER_LOCATION,
        detection_alpha=DETECTION_ALPHA,
    )
    started = time.time()
    model.train(
        max_epochs=SPATIAL_MAX_EPOCHS,
        batch_size=None,
        train_size=1,
        **trainer_device_kwargs(),
    )
    save_history(model, f"{gse}__cell2location")
    spatial = model.export_posterior(
        spatial,
        sample_kwargs={
            "num_samples": POSTERIOR_SAMPLES,
            "batch_size": spatial.n_obs,
        },
    )
    model.save(str(MAP_OUT / f"{gse}__model"), overwrite=True)
    spatial.obs.index.name = None
    spatial.var.index.name = None
    spatial.write_h5ad(output, compression="gzip")
    pd.DataFrame([{
        "GSE": gse,
        "n_spots": int(spatial.n_obs),
        "n_genes": int(spatial.n_vars),
        "max_epochs": SPATIAL_MAX_EPOCHS,
        "posterior_samples": POSTERIOR_SAMPLES,
        "elapsed_seconds": time.time() - started,
        "q05_available": "q05_cell_abundance_w_sf" in spatial.obsm,
        "mean_available": "means_cell_abundance_w_sf" in spatial.obsm,
    }]).to_csv(AUDIT_OUT / f"{gse}__mapping_audit.tsv", sep="\t", index=False)
    del model, spatial
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return output


if not MANIFEST_PATH.exists():
    raise FileNotFoundError(
        f"Run 01_spatial_v3_input_audit.py first; manifest not found: {MANIFEST_PATH}"
    )
manifest = pd.read_csv(MANIFEST_PATH, sep="\t")
gpu_audit()
for column in ["include_mapping", "include_within_section", "include_group_contrast"]:
    manifest[column] = manifest[column].astype(str).str.lower().isin({"true", "1", "yes"})

signatures = fit_or_load_reference()
mapped = []
failures = []
requested_gses = sorted(manifest.loc[manifest["include_mapping"], "GSE"].unique())
only_gses = {
    value.strip()
    for value in os.environ.get("SPATIAL_ONLY_GSES", "").split(",")
    if value.strip()
}
if only_gses:
    requested_gses = [gse for gse in requested_gses if gse in only_gses]
    if not requested_gses:
        raise ValueError(f"SPATIAL_ONLY_GSES did not match the manifest: {sorted(only_gses)}")
    print("Restricted mapping run:", requested_gses)
for gse in requested_gses:
    try:
        mapped.append(str(map_one_gse(gse, manifest, signatures)))
    except Exception as exc:
        failures.append({"GSE": gse, "error": repr(exc)})
        print("FAILED", gse, repr(exc))
pd.DataFrame({"mapped_h5ad": mapped}).to_csv(MAP_OUT / "mapped_h5ad_manifest.tsv", sep="\t", index=False)
pd.DataFrame(failures, columns=["GSE", "error"]).to_csv(
    MAP_OUT / "mapping_failures.tsv", sep="\t", index=False
)
if failures:
    raise RuntimeError(f"{len(failures)} GSE mappings failed; inspect mapping_failures.tsv")
print("All requested final-atlas cell2location mappings completed.")
