from __future__ import annotations

import json
import math
import os
import platform
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import norm
import statsmodels
import statsmodels.api as sm


# =============================================================================
# Configuration
# =============================================================================

ROOT = Path(
    "/home/anzhen_vhd/scRNAseq/codex/"
    "spatial_validation_v5_locked71_gpu_final"
)
OUT = Path(
    "/home/anzhen_vhd/scRNAseq/codex/"
    "spatial_validation_v5_locked71_spatial_statistics_audit_v1"
)
OUT.mkdir(parents=True, exist_ok=True)

RANDOM_SEED = int(os.environ.get("SPATIAL_AUDIT_RANDOM_SEED", "20260826"))
N_MORAN_PERMUTATIONS = int(os.environ.get("SPATIAL_AUDIT_MORAN_PERMUTATIONS", "999"))
N_MATCHED_GENES = int(os.environ.get("SPATIAL_AUDIT_MATCHED_GENES", "100"))
K_NEIGHBORS = int(os.environ.get("SPATIAL_AUDIT_K_NEIGHBORS", "6"))
MIN_SPOTS = int(os.environ.get("SPATIAL_AUDIT_MIN_SPOTS", "40"))
MIN_POSITIVE_SPOTS = int(os.environ.get("SPATIAL_AUDIT_MIN_POSITIVE_SPOTS", "10"))
POSTERIOR_KEY = "q05_cell_abundance_w_sf"

MAPPED_DIR = ROOT / "02_cell2location"
LABEL_MAP = ROOT / "01_reference" / "spatial_reference_state_label_map.tsv"
ROSTER_CANDIDATES = [
    ROOT / "05_all_public_Tier1_spatial_validation_v1" / "00_audit" / "Tier1_roster_used.tsv",
    ROOT / "05_all_public_Tier1_spatial_validation_v1" / "00_audit" / "Tier1_roster_redaction_audit.tsv",
]


def bh(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce").to_numpy(float)
    out = np.full(len(p), np.nan)
    keep = np.isfinite(p)
    if not keep.any():
        return pd.Series(out, index=values.index)
    x = p[keep]
    order = np.argsort(x)
    ranked = x[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.clip(adjusted, 0, 1)
    out[np.where(keep)[0]] = restored
    return pd.Series(out, index=values.index)


def first_existing(paths: list[Path]) -> Path:
    for path in paths:
        if path.exists():
            return path
    raise FileNotFoundError("None of the required files exists:\n" + "\n".join(map(str, paths)))


def choose_column(frame: pd.DataFrame, candidates: list[str], label: str) -> str:
    for name in candidates:
        if name in frame.columns:
            return name
    raise KeyError(f"Cannot resolve {label}. Tried {candidates}; observed={list(frame.columns)}")


def factor_names(adata: ad.AnnData, key: str) -> list[str]:
    value = adata.obsm[key]
    names = adata.uns.get("mod", {}).get("factor_names")
    if names is not None and len(names) == value.shape[1]:
        return list(map(str, names))
    if isinstance(value, pd.DataFrame):
        columns = list(map(str, value.columns))
        prefix = key.replace("_cell", "cell", 1) + "_"
        if all(column.startswith(prefix) for column in columns):
            return [column[len(prefix):] for column in columns]
        return columns
    raise ValueError(f"Cannot resolve factor names for {key}")


def abundance_frame(adata: ad.AnnData, key: str) -> pd.DataFrame:
    value = adata.obsm[key]
    values = value.to_numpy(float) if isinstance(value, pd.DataFrame) else np.asarray(value, dtype=float)
    return pd.DataFrame(values, index=adata.obs_names, columns=factor_names(adata, key))


def resolve_counts(adata: ad.AnnData):
    candidates = []
    if "counts" in adata.layers:
        candidates.append(("layers[counts]", adata.layers["counts"]))
    candidates.append(("X", adata.X))
    for label, matrix in candidates:
        sample = matrix[: min(100, adata.n_obs), : min(100, adata.n_vars)]
        sample = sample.toarray() if sparse.issparse(sample) else np.asarray(sample)
        if np.isfinite(sample).all() and np.nanmin(sample) >= 0:
            fraction_integer = np.mean(np.isclose(sample, np.round(sample), atol=1e-6))
            if fraction_integer >= 0.98:
                return label, matrix
    raise ValueError("No nonnegative integer-like spatial count matrix was found")


def matrix_column(matrix, row_indices: np.ndarray, column_index: int) -> np.ndarray:
    value = matrix[row_indices, column_index]
    if sparse.issparse(value):
        return value.toarray().ravel().astype(float)
    return np.asarray(value).ravel().astype(float)


def row_sums(matrix) -> np.ndarray:
    value = matrix.sum(axis=1)
    return np.asarray(value).ravel().astype(float)


def block_groups(coords: np.ndarray, bins: int = 5) -> np.ndarray:
    groups = np.zeros(len(coords), dtype=int)
    for axis in range(2):
        x = coords[:, axis]
        edges = np.quantile(x, np.linspace(0, 1, bins + 1))
        edges = np.unique(edges)
        if len(edges) <= 2:
            labels = np.zeros(len(x), dtype=int)
        else:
            labels = np.digitize(x, edges[1:-1], right=False)
        groups += labels * (bins ** axis)
    return groups


def zscore(x: np.ndarray) -> np.ndarray:
    sd = np.nanstd(x, ddof=1)
    if not np.isfinite(sd) or sd <= 0:
        return np.zeros_like(x, dtype=float)
    return (x - np.nanmean(x)) / sd


def fit_count_models(y: np.ndarray, context: np.ndarray, library: np.ndarray, coords: np.ndarray) -> dict:
    x = zscore(np.log1p(np.maximum(context, 0)))
    if np.std(x) <= 0 or np.sum(y > 0) < MIN_POSITIVE_SPOTS:
        raise ValueError("Insufficient target expression or context variation")
    design = sm.add_constant(x, has_constant="add")
    offset = np.log(np.maximum(library, 1))
    groups = block_groups(coords)

    poisson = sm.GLM(y, design, family=sm.families.Poisson(), offset=offset)
    try:
        pfit = poisson.fit(
            maxiter=100,
            disp=0,
            cov_type="cluster",
            cov_kwds={"groups": groups, "use_correction": True},
        )
    except Exception:
        pfit = poisson.fit(maxiter=100, disp=0, cov_type="HC1")

    mean_y = float(np.mean(y))
    var_y = float(np.var(y, ddof=1)) if len(y) > 1 else mean_y
    alpha = min(max((var_y - mean_y) / max(mean_y * mean_y, 1e-8), 1e-4), 10.0)
    nb = sm.GLM(
        y, design, family=sm.families.NegativeBinomial(alpha=alpha), offset=offset
    )
    try:
        nfit = nb.fit(
            maxiter=100,
            disp=0,
            cov_type="cluster",
            cov_kwds={"groups": groups, "use_correction": True},
        )
    except Exception:
        nfit = nb.fit(maxiter=100, disp=0, cov_type="HC1")

    return {
        "poisson_beta": float(pfit.params[1]),
        "poisson_se": float(pfit.bse[1]),
        "poisson_p": float(pfit.pvalues[1]),
        "poisson_dispersion": float(np.sum(np.square(pfit.resid_pearson)) / max(pfit.df_resid, 1)),
        "nb_alpha": alpha,
        "nb_beta": float(nfit.params[1]),
        "nb_se": float(nfit.bse[1]),
        "nb_p": float(nfit.pvalues[1]),
        "pearson_residuals": np.asarray(nfit.resid_pearson, dtype=float),
    }


def neighbor_edges(coords: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    k_use = min(k + 1, len(coords))
    indices = cKDTree(coords).query(coords, k=k_use)[1]
    if indices.ndim == 1:
        indices = indices[:, None]
    neighbors = indices[:, 1:] if indices.shape[1] > 1 else indices[:, :0]
    sources = np.repeat(np.arange(len(coords)), neighbors.shape[1])
    targets = neighbors.reshape(-1)
    keep = sources != targets
    return sources[keep], targets[keep]


def moran_statistic(values: np.ndarray, sources: np.ndarray, targets: np.ndarray) -> float:
    centered = values - np.mean(values)
    denominator = np.sum(np.square(centered))
    if denominator <= 0 or len(sources) == 0:
        return np.nan
    return float(len(values) / len(sources) * np.sum(centered[sources] * centered[targets]) / denominator)


def moran_permutation(values: np.ndarray, coords: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    sources, targets = neighbor_edges(coords, K_NEIGHBORS)
    observed = moran_statistic(values, sources, targets)
    if not np.isfinite(observed):
        return np.nan, np.nan
    extreme = 0
    for _ in range(N_MORAN_PERMUTATIONS):
        permuted = rng.permutation(values)
        value = moran_statistic(permuted, sources, targets)
        extreme += int(abs(value) >= abs(observed))
    return observed, (1 + extreme) / (N_MORAN_PERMUTATIONS + 1)


def matched_gene_roster(matrix, var_names: pd.Index, targets: set[str]) -> pd.DataFrame:
    n_spots = matrix.shape[0]
    totals = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
    if sparse.issparse(matrix):
        detected = np.asarray((matrix > 0).sum(axis=0)).ravel().astype(float)
    else:
        detected = np.sum(np.asarray(matrix) > 0, axis=0).astype(float)
    features = pd.DataFrame(
        {
            "gene": var_names.astype(str),
            "mean_count": totals / max(n_spots, 1),
            "detection_fraction": detected / max(n_spots, 1),
        }
    )
    features["log_mean"] = np.log1p(features["mean_count"])
    p = np.clip(features["detection_fraction"], 1e-5, 1 - 1e-5)
    features["logit_detection"] = np.log(p / (1 - p))
    eligible = (
        (features["detection_fraction"] >= 0.005)
        & (features["gene"].str.match(r"^[A-Za-z0-9.-]+$"))
        & (~features["gene"].str.upper().str.startswith(("MT-", "RPL", "RPS", "HB")))
        & (~features["gene"].isin(targets))
    )
    pool = features[eligible].copy()
    rows = []
    for target in sorted(targets):
        hit = features[features["gene"] == target]
        if hit.empty:
            continue
        target_row = hit.iloc[0]
        scale_mean = max(pool["log_mean"].std(), 1e-8)
        scale_detect = max(pool["logit_detection"].std(), 1e-8)
        distance = np.sqrt(
            np.square((pool["log_mean"] - target_row["log_mean"]) / scale_mean)
            + np.square((pool["logit_detection"] - target_row["logit_detection"]) / scale_detect)
        )
        selected = pool.assign(match_distance=distance).nsmallest(N_MATCHED_GENES, "match_distance")
        for rank, (_, row) in enumerate(selected.iterrows(), start=1):
            rows.append(
                {
                    "target_gene": target,
                    "matched_gene": row["gene"],
                    "match_rank": rank,
                    "match_distance": row["match_distance"],
                    "target_mean_count": target_row["mean_count"],
                    "matched_mean_count": row["mean_count"],
                    "target_detection_fraction": target_row["detection_fraction"],
                    "matched_detection_fraction": row["detection_fraction"],
                }
            )
    return pd.DataFrame(rows)


def linear_association_z(y: np.ndarray, x: np.ndarray) -> float:
    keep = np.isfinite(y) & np.isfinite(x)
    y = y[keep]
    x = x[keep]
    if len(y) < MIN_SPOTS or np.std(x) <= 0 or np.std(y) <= 0:
        return np.nan
    x = zscore(x)
    y = zscore(y)
    rho = float(np.dot(x, y) / max(len(x) - 1, 1))
    rho = float(np.clip(rho, -0.999999, 0.999999))
    return float(np.arctanh(rho) * math.sqrt(max(len(x) - 3, 1)))


def resolve_roster() -> pd.DataFrame:
    path = first_existing(ROSTER_CANDIDATES)
    x = pd.read_csv(path, sep="\t")
    gene_col = choose_column(x, ["gene_symbol", "gene", "gene_symbol_match"], "Tier 1 gene")
    state_col = choose_column(x, ["cell_state", "cell_context"], "Tier 1 cell state")
    level_col = choose_column(x, ["annotation_level", "context_level"], "annotation level")
    result = pd.DataFrame(
        {
            "gene": x[gene_col].astype(str),
            "cell_state": x[state_col].astype(str),
            "annotation_level": x[level_col].astype(str),
        }
    ).drop_duplicates()
    result["pair_id"] = result["gene"] + " | " + result["cell_state"]
    return result


def context_vector(
    abundance: pd.DataFrame,
    roster_row: pd.Series,
    original_to_safe: dict[str, str],
) -> np.ndarray:
    level = (
        str(roster_row["annotation_level"])
        .lower()
        .replace(" ", "")
        .replace("_", "")
        .replace("-", "")
    )
    state = roster_row["cell_state"]
    if "level1" in level:
        original = [name for name in original_to_safe if name.startswith(state + " | ")]
        safe = [original_to_safe[name] for name in original if original_to_safe[name] in abundance.columns]
        if not safe:
            raise KeyError(f"No locked Level 2 factors found below Level 1 context {state}")
        return abundance[safe].sum(axis=1).to_numpy(float)
    safe = original_to_safe.get(state, state)
    if safe not in abundance.columns:
        raise KeyError(f"Level 2 factor absent: original={state}, safe={safe}")
    return abundance[safe].to_numpy(float)


def main() -> None:
    rng = np.random.default_rng(RANDOM_SEED)
    mapped = sorted(MAPPED_DIR.glob("*__final_atlas_cell2location.h5ad"))
    if len(mapped) != 7:
        raise FileNotFoundError(f"Expected 7 mapped h5ad files under {MAPPED_DIR}; found {len(mapped)}")
    labels = pd.read_csv(LABEL_MAP, sep="\t")
    original_to_safe = dict(
        zip(
            labels["spatial_reference_state_original"].astype(str),
            labels["spatial_reference_state_safe"].astype(str),
        )
    )
    roster = resolve_roster()
    target_genes = set(roster["gene"])

    dataset_audit = []
    target_rows = []
    match_rosters = []
    match_section_rows = []
    failures = []

    for path in mapped:
        gse = path.name.split("__", 1)[0]
        print("\nProcessing", gse, path, flush=True)
        adata = ad.read_h5ad(path)
        if POSTERIOR_KEY not in adata.obsm or "spatial" not in adata.obsm:
            raise KeyError(f"Required posterior or spatial coordinates absent from {path}")
        count_source, matrix = resolve_counts(adata)
        var_names = pd.Index(adata.var_names.astype(str))
        if var_names.has_duplicates:
            raise ValueError(f"Duplicate spatial gene identifiers in {gse}")
        gene_index = {gene: i for i, gene in enumerate(var_names)}
        abundance = abundance_frame(adata, POSTERIOR_KEY)
        section_col = choose_column(
            adata.obs,
            ["section_id", "spatial_section_uid", "sample_uid", "sample_id", "library_id", "sample"],
            "spatial section",
        )
        unit_col = next(
            (x for x in ["biological_unit_id", "donor_uid", "donor_id", "patient", "donor"] if x in adata.obs),
            section_col,
        )
        library = row_sums(matrix)
        match_roster = matched_gene_roster(matrix, var_names, target_genes)
        match_roster.insert(0, "GSE", gse)
        match_rosters.append(match_roster)
        matches = match_roster.groupby("target_gene")["matched_gene"].apply(list).to_dict()

        dataset_audit.append(
            {
                "GSE": gse,
                "n_spots": adata.n_obs,
                "n_genes": adata.n_vars,
                "n_sections": adata.obs[section_col].astype(str).nunique(),
                "section_field": section_col,
                "biological_unit_field": unit_col,
                "count_source": count_source,
                "n_factors": abundance.shape[1],
            }
        )

        section_values = adata.obs[section_col].astype(str).to_numpy()
        unit_values = adata.obs[unit_col].astype(str).to_numpy()
        coords_all = np.asarray(adata.obsm["spatial"], dtype=float)
        for section in sorted(pd.unique(section_values)):
            idx = np.where(section_values == section)[0]
            if len(idx) < MIN_SPOTS:
                continue
            unit = pd.Series(unit_values[idx]).mode().iloc[0]
            coords = coords_all[idx]
            lib = library[idx]
            for _, row in roster.iterrows():
                gene = row["gene"]
                pair_id = row["pair_id"]
                if gene not in gene_index:
                    failures.append({"GSE": gse, "section_id": section, "pair_id": pair_id, "error": "gene_absent"})
                    continue
                try:
                    context_all = context_vector(abundance, row, original_to_safe)
                    context = context_all[idx]
                    y = matrix_column(matrix, idx, gene_index[gene])
                    model = fit_count_models(y, context, lib, coords)
                    moran_i, moran_p = moran_permutation(model.pop("pearson_residuals"), coords, rng)
                    target_rows.append(
                        {
                            "GSE": gse,
                            "section_id": section,
                            "biological_unit_id": f"{gse}::{unit}",
                            "pair_id": pair_id,
                            "gene": gene,
                            "annotation_level": row["annotation_level"],
                            "cell_state": row["cell_state"],
                            "n_spots": len(idx),
                            "n_positive_spots": int(np.sum(y > 0)),
                            **model,
                            "nb_residual_moran_I": moran_i,
                            "nb_residual_moran_permutation_p": moran_p,
                        }
                    )

                    candidate_genes = [gene] + matches.get(gene, [])
                    candidate_indices = [gene_index[x] for x in candidate_genes if x in gene_index]
                    candidate_genes = [x for x in candidate_genes if x in gene_index]
                    sub = matrix[idx][:, candidate_indices]
                    if sparse.issparse(sub):
                        sub = sub.toarray()
                    sub = np.asarray(sub, dtype=float)
                    lognorm = np.log1p(sub / np.maximum(lib[:, None], 1) * 1e4)
                    context_scaled = zscore(np.log1p(np.maximum(context, 0)))
                    for column, candidate in enumerate(candidate_genes):
                        match_section_rows.append(
                            {
                                "GSE": gse,
                                "section_id": section,
                                "biological_unit_id": f"{gse}::{unit}",
                                "pair_id": pair_id,
                                "target_gene": gene,
                                "tested_gene": candidate,
                                "is_target_gene": candidate == gene,
                                "association_z": linear_association_z(lognorm[:, column], context_scaled),
                                "n_spots": len(idx),
                            }
                        )
                except Exception as exc:
                    failures.append(
                        {"GSE": gse, "section_id": section, "pair_id": pair_id, "error": repr(exc)}
                    )
        del adata

    target = pd.DataFrame(target_rows)
    target["nb_residual_moran_FDR"] = bh(target["nb_residual_moran_permutation_p"])
    target["poisson_nb_direction_agreement"] = np.sign(target["poisson_beta"]) == np.sign(target["nb_beta"])
    target.to_csv(OUT / "01_target_NB_Poisson_residual_Moran_by_section.tsv.gz", sep="\t", index=False, compression="gzip")

    target_summary = (
        target.groupby("pair_id", as_index=False)
        .agg(
            n_sections=("section_id", "nunique"),
            n_datasets=("GSE", "nunique"),
            n_biological_units=("biological_unit_id", "nunique"),
            poisson_nb_direction_agreement_fraction=("poisson_nb_direction_agreement", "mean"),
            median_poisson_dispersion=("poisson_dispersion", "median"),
            median_nb_residual_moran_I=("nb_residual_moran_I", "median"),
            fraction_residual_moran_FDR_lt_0_05=("nb_residual_moran_FDR", lambda x: np.mean(x < 0.05)),
        )
    )
    target_summary.to_csv(OUT / "02_target_model_sensitivity_summary.tsv", sep="\t", index=False)

    match_roster_all = pd.concat(match_rosters, ignore_index=True)
    match_roster_all.to_csv(OUT / "03_expression_detection_matched_gene_roster.tsv.gz", sep="\t", index=False, compression="gzip")
    match_sections = pd.DataFrame(match_section_rows)
    match_sections.to_csv(OUT / "04_matched_gene_associations_by_section.tsv.gz", sep="\t", index=False, compression="gzip")

    unit_meta = (
        match_sections.dropna(subset=["association_z"])
        .groupby(["pair_id", "target_gene", "tested_gene", "is_target_gene", "biological_unit_id"], as_index=False)
        .agg(unit_z=("association_z", lambda x: np.sum(x) / math.sqrt(len(x))))
    )
    gene_meta = (
        unit_meta.groupby(["pair_id", "target_gene", "tested_gene", "is_target_gene"], as_index=False)
        .agg(
            n_biological_units=("biological_unit_id", "nunique"),
            meta_z=("unit_z", lambda x: np.sum(x) / math.sqrt(len(x))),
        )
    )
    gene_meta.to_csv(OUT / "05_matched_gene_meta_null_distribution.tsv.gz", sep="\t", index=False, compression="gzip")

    empirical_rows = []
    for pair_id, g in gene_meta.groupby("pair_id", sort=False):
        target_hit = g[g["is_target_gene"]]
        null = g[~g["is_target_gene"]]
        if target_hit.empty or null.empty:
            continue
        z = float(target_hit.iloc[0]["meta_z"])
        empirical_rows.append(
            {
                "pair_id": pair_id,
                "target_gene": target_hit.iloc[0]["target_gene"],
                "n_biological_units": int(target_hit.iloc[0]["n_biological_units"]),
                "target_meta_z": z,
                "n_matched_genes": len(null),
                "matched_null_median_meta_z": float(null["meta_z"].median()),
                "matched_null_95pct_abs_meta_z": float(null["meta_z"].abs().quantile(0.95)),
                "empirical_P_two_sided": float((1 + np.sum(null["meta_z"].abs() >= abs(z))) / (len(null) + 1)),
            }
        )
    empirical = pd.DataFrame(empirical_rows)
    empirical["matched_gene_FDR_two_sided"] = bh(empirical["empirical_P_two_sided"])
    empirical.to_csv(OUT / "06_matched_gene_empirical_summary.tsv", sep="\t", index=False)

    pd.DataFrame(dataset_audit).to_csv(OUT / "00_dataset_input_audit.tsv", sep="\t", index=False)
    pd.DataFrame(failures).to_csv(OUT / "07_failures.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "quantity": "section-level random-module empirical P",
                "source": "existing 1,000 lineage-matched module permutations",
                "minimum_resolution": 1 / 1001,
            },
            {
                "quantity": "NB residual Moran permutation P",
                "source": f"{N_MORAN_PERMUTATIONS} coordinate-preserving residual-label permutations",
                "minimum_resolution": 1 / (N_MORAN_PERMUTATIONS + 1),
            },
            {
                "quantity": "matched-gene empirical P",
                "source": f"up to {N_MATCHED_GENES} genes matched for mean expression and detection",
                "minimum_resolution": 1 / (N_MATCHED_GENES + 1),
            },
        ]
    ).to_csv(OUT / "08_P_value_provenance.tsv", sep="\t", index=False)

    manifest = {
        "status": "completed",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "input_root": str(ROOT),
        "output_root": str(OUT),
        "random_seed": RANDOM_SEED,
        "n_moran_permutations": N_MORAN_PERMUTATIONS,
        "n_matched_genes": N_MATCHED_GENES,
        "posterior": POSTERIOR_KEY,
        "primary_target_sensitivity": "negative-binomial versus Poisson direction and effect",
        "spatial_residual_test": "Moran's I of negative-binomial Pearson residuals on six-nearest-neighbour graph",
        "matched_null": "genes matched within each spatial dataset for mean count and detection fraction",
        "python": sys.version,
        "platform": platform.platform(),
        "anndata": ad.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "statsmodels": statsmodels.__version__,
    }
    (OUT / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("\nSPATIAL STATISTICAL AUDIT COMPLETE")
    print("Target models:", len(target))
    print("Matched-gene section associations:", len(match_sections))
    print("Failures:", len(failures))
    print("Output:", OUT)


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("default")
        main()
