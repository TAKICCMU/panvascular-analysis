from __future__ import annotations

import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr


ROOT = Path(os.environ.get(
    "SPATIAL_V4_ROOT",
    "/path/to/panvas_project/spatial_validation_v4_gpu_final",
))
OUT = Path(os.environ.get(
    "CELL2LOCATION_QC_OUT",
    str(ROOT / "06_cell2location_model_quality_audit_v1"),
))
OUT.mkdir(parents=True, exist_ok=True)
MAP_ROOT = ROOT / "02_cell2location"
AUDIT_ROOT = ROOT / "00_audit"
SIGNATURE_PATH = ROOT / "01_reference" / "final_atlas_reference_signatures.tsv.gz"

CHUNK_SIZE = int(os.environ.get("CELL2LOCATION_QC_CHUNK_SIZE", "128"))
EPS = 1e-8


def abundance_frame(adata: ad.AnnData, key: str) -> pd.DataFrame:
    value = adata.obsm[key]
    if isinstance(value, pd.DataFrame):
        result = value.copy()
        result.index = adata.obs_names
        result.columns = list(map(str, result.columns))
        return result
    names = adata.uns.get("mod", {}).get("factor_names")
    if names is None or len(names) != value.shape[1]:
        raise ValueError(f"Cannot resolve factor names for {key}")
    return pd.DataFrame(np.asarray(value), index=adata.obs_names, columns=list(map(str, names)))


def poisson_deviance(y: np.ndarray, mu: np.ndarray) -> float:
    mu = np.maximum(mu, EPS)
    positive = y > 0
    term = np.zeros_like(y, dtype=np.float64)
    term[positive] = y[positive] * np.log(y[positive] / mu[positive])
    return float(2.0 * np.sum(term - (y - mu)))


def training_history_summary(path: Path) -> dict:
    table = pd.read_csv(path, sep="\t")
    table["value"] = pd.to_numeric(table["value"], errors="coerce")
    candidates = table[table["metric"].astype(str).str.contains("elbo|loss", case=False, regex=True)]
    if candidates.empty:
        candidates = table
    metric = candidates["metric"].value_counts().index[0]
    x = candidates[candidates["metric"] == metric].sort_values("epoch")
    values = x["value"].to_numpy(float)
    finite = np.isfinite(values)
    if finite.sum() < 5:
        return {
            "history_file": str(path), "metric": metric,
            "n_epochs_recorded": int(len(values)), "finite_fraction": float(finite.mean()),
            "history_usable": False,
        }
    values = values[finite]
    n_tail = max(5, int(math.ceil(len(values) * 0.1)))
    tail = values[-n_tail:]
    axis = np.arange(len(tail), dtype=float)
    slope = float(np.polyfit(axis, tail, 1)[0]) if len(tail) > 1 else np.nan
    scale = max(float(np.ptp(values)), abs(float(np.median(tail))), EPS)
    return {
        "history_file": str(path),
        "metric": metric,
        "n_epochs_recorded": int(len(values)),
        "finite_fraction": float(finite.mean()),
        "first_value": float(values[0]),
        "last_value": float(values[-1]),
        "tail_median": float(np.median(tail)),
        "tail_CV": float(np.std(tail) / max(abs(np.mean(tail)), EPS)),
        "tail_normalized_slope": slope / scale,
        "history_usable": True,
    }


if not SIGNATURE_PATH.exists():
    raise FileNotFoundError(SIGNATURE_PATH)
signatures = pd.read_csv(SIGNATURE_PATH, sep="\t", index_col=0)
signatures.index = signatures.index.astype(str).str.upper()
signatures.columns = signatures.columns.astype(str)

mapped_paths = sorted(MAP_ROOT.glob("*__final_atlas_cell2location.h5ad"))
if not mapped_paths:
    raise FileNotFoundError(f"No mapped h5ad files under {MAP_ROOT}")

history_rows = []
for path in sorted(AUDIT_ROOT.glob("*__training_history.tsv.gz")):
    result = training_history_summary(path)
    stem = path.name.replace("__training_history.tsv.gz", "")
    result["dataset"] = stem.replace("__cell2location", "")
    history_rows.append(result)
history_table = pd.DataFrame(history_rows)
history_table.to_csv(OUT / "01_training_history_summary.tsv", sep="\t", index=False)

dataset_rows = []
factor_rows = []
for path in mapped_paths:
    gse = path.name.split("__", 1)[0]
    print("Auditing", gse, flush=True)
    spatial = ad.read_h5ad(path)
    for key in ("q05_cell_abundance_w_sf", "means_cell_abundance_w_sf"):
        if key not in spatial.obsm:
            raise KeyError(f"{key} absent from {path}")
    q05 = abundance_frame(spatial, "q05_cell_abundance_w_sf")
    mean = abundance_frame(spatial, "means_cell_abundance_w_sf")
    common_factors = q05.columns.intersection(mean.columns).intersection(signatures.columns)
    if len(common_factors) < 20:
        raise ValueError(f"Only {len(common_factors)} common factors for {gse}")
    q05 = q05[common_factors]
    mean = mean[common_factors]

    posterior_values = np.concatenate([q05.to_numpy().ravel(), mean.to_numpy().ravel()])
    finite_fraction = float(np.isfinite(posterior_values).mean())
    nonnegative_fraction = float((posterior_values[np.isfinite(posterior_values)] >= 0).mean())
    mean_ge_q05_fraction = float((mean.to_numpy() + 1e-10 >= q05.to_numpy()).mean())
    total_mean = mean.sum(axis=1).to_numpy(float)
    original_library = pd.to_numeric(
        spatial.obs.get("total_counts_original", pd.Series(index=spatial.obs_names, dtype=float)),
        errors="coerce",
    ).to_numpy(float)
    original_genes = pd.to_numeric(
        spatial.obs.get("n_genes_original", pd.Series(index=spatial.obs_names, dtype=float)),
        errors="coerce",
    ).to_numpy(float)
    rho_total_counts = spearmanr(total_mean, original_library, nan_policy="omit").statistic
    rho_detected_genes = spearmanr(total_mean, original_genes, nan_policy="omit").statistic

    per_factor_rho = []
    for factor in common_factors:
        rho = spearmanr(q05[factor], mean[factor], nan_policy="omit").statistic
        per_factor_rho.append(rho)
        factor_rows.append({
            "dataset": gse,
            "factor": factor,
            "q05_mean_spearman_rho": rho,
            "q05_positive_fraction": float((q05[factor].to_numpy() > 0).mean()),
            "mean_positive_fraction": float((mean[factor].to_numpy() > 0).mean()),
            "median_q05": float(np.median(q05[factor])),
            "median_mean": float(np.median(mean[factor])),
        })

    # Posterior-mean reconstruction audit on mapping genes.
    common_genes = pd.Index(spatial.var_names.astype(str).str.upper()).intersection(signatures.index)
    if len(common_genes) < 1000:
        raise ValueError(f"Only {len(common_genes)} reconstruction genes for {gse}")
    positions = pd.Index(spatial.var_names.astype(str).str.upper()).get_indexer(common_genes)
    X = spatial.X[:, positions]
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X = X.tocsr().astype(np.float64)
    S = signatures.loc[common_genes, common_factors].to_numpy(float)
    A = mean[common_factors].to_numpy(float)
    global_gene_counts = np.asarray(X.sum(axis=0)).ravel()
    global_gene_p = global_gene_counts / max(global_gene_counts.sum(), EPS)
    observed_gene_sum = np.zeros(len(common_genes), dtype=float)
    predicted_gene_sum = np.zeros(len(common_genes), dtype=float)
    cosine_values = []
    dev_model = 0.0
    dev_null = 0.0
    for start in range(0, spatial.n_obs, CHUNK_SIZE):
        stop = min(start + CHUNK_SIZE, spatial.n_obs)
        y = X[start:stop].toarray()
        mu = A[start:stop] @ S.T
        library = y.sum(axis=1)
        mu_sum = mu.sum(axis=1)
        mu *= (library / np.maximum(mu_sum, EPS))[:, None]
        null_mu = library[:, None] * global_gene_p[None, :]
        observed_gene_sum += y.sum(axis=0)
        predicted_gene_sum += mu.sum(axis=0)
        numerator = np.sum(y * mu, axis=1)
        denominator = np.sqrt(np.sum(y * y, axis=1) * np.sum(mu * mu, axis=1))
        cosine_values.extend((numerator / np.maximum(denominator, EPS)).tolist())
        dev_model += poisson_deviance(y, mu)
        dev_null += poisson_deviance(y, null_mu)
    gene_mean_rho = spearmanr(observed_gene_sum, predicted_gene_sum).statistic
    deviance_explained = 1 - dev_model / dev_null if dev_null > 0 else np.nan

    row = {
        "dataset": gse,
        "mapped_h5ad": str(path),
        "n_spots": int(spatial.n_obs),
        "n_mapping_genes": int(spatial.n_vars),
        "n_reconstruction_genes": int(len(common_genes)),
        "n_factors": int(len(common_factors)),
        "posterior_finite_fraction": finite_fraction,
        "posterior_nonnegative_fraction": nonnegative_fraction,
        "mean_ge_q05_fraction": mean_ge_q05_fraction,
        "median_q05_mean_factor_rho": float(np.nanmedian(per_factor_rho)),
        "minimum_q05_mean_factor_rho": float(np.nanmin(per_factor_rho)),
        "total_abundance_vs_original_counts_rho": rho_total_counts,
        "total_abundance_vs_detected_genes_rho": rho_detected_genes,
        "posterior_mean_reconstruction_gene_rho": gene_mean_rho,
        "posterior_mean_reconstruction_median_spot_cosine": float(np.nanmedian(cosine_values)),
        "posterior_mean_reconstruction_deviance_explained": deviance_explained,
    }
    row["gate_finite_nonnegative"] = finite_fraction == 1 and nonnegative_fraction == 1
    row["gate_posterior_order"] = mean_ge_q05_fraction >= 0.999
    row["gate_q05_mean_stability"] = row["median_q05_mean_factor_rho"] >= 0.9
    row["gate_reconstruction"] = (
        np.isfinite(gene_mean_rho) and gene_mean_rho >= 0.5
        and np.isfinite(deviance_explained) and deviance_explained > 0
    )
    row["all_numeric_gates_pass"] = bool(
        row["gate_finite_nonnegative"]
        and row["gate_posterior_order"]
        and row["gate_q05_mean_stability"]
        and row["gate_reconstruction"]
    )
    dataset_rows.append(row)
    del spatial, X

dataset_table = pd.DataFrame(dataset_rows)
factor_table = pd.DataFrame(factor_rows)
dataset_table.to_csv(OUT / "02_dataset_model_quality_summary.tsv", sep="\t", index=False)
factor_table.to_csv(OUT / "03_factor_posterior_stability.tsv.gz", sep="\t", index=False, compression="gzip")

fig, axes = plt.subplots(1, 3, figsize=(8.0, 2.6))
plot_specs = [
    ("median_q05_mean_factor_rho", "q05-mean factor rho", (0, 1)),
    ("posterior_mean_reconstruction_gene_rho", "Gene-mean reconstruction rho", (0, 1)),
    ("posterior_mean_reconstruction_deviance_explained", "Deviance explained", None),
]
for ax, (column, label, ylim) in zip(axes, plot_specs):
    ax.bar(dataset_table["dataset"], dataset_table[column], color="#2C5F8A", width=0.72)
    ax.axhline(0, color="#1D2733", lw=0.6)
    ax.set_ylabel(label, fontsize=7)
    ax.tick_params(axis="x", rotation=60, labelsize=6)
    ax.tick_params(axis="y", labelsize=6)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout(w_pad=1.2)
fig.savefig(OUT / "cell2location_model_quality_summary.pdf", bbox_inches="tight")
fig.savefig(OUT / "cell2location_model_quality_summary.svg", bbox_inches="tight")
fig.savefig(OUT / "cell2location_model_quality_summary.tiff", dpi=600, bbox_inches="tight")
plt.close(fig)

manifest = {
    "status": "completed",
    "generated_at": datetime.now(timezone.utc).isoformat(),
    "root": str(ROOT),
    "n_datasets": int(len(dataset_table)),
    "n_all_numeric_gates_pass": int(dataset_table["all_numeric_gates_pass"].sum()),
    "quality_gates": {
        "finite_nonnegative": "all posterior values finite and nonnegative",
        "posterior_order": "mean abundance >= q05 abundance for >=99.9% of entries",
        "posterior_stability": "median factor-wise Spearman(q05, mean) >=0.90",
        "reconstruction": "gene-mean rho >=0.50 and posterior-mean deviance explained >0",
    },
    "scope_note": (
        "The reconstruction is a posterior-mean model audit using estimated cell-state "
        "abundances and reference signatures. It is not labelled as a formal posterior "
        "predictive p-value."
    ),
}
(OUT / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
print(json.dumps(manifest, indent=2, ensure_ascii=False))
print("Output:", OUT)
