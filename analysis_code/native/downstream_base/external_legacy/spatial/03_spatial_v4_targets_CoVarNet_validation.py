from __future__ import annotations

import itertools
import json
import math
import re
import warnings
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import binomtest, norm, spearmanr

try:
    import statsmodels.api as sm
except ImportError as exc:
    raise ImportError("statsmodels is required for the spatial count models") from exc

try:
    from sklearn.neighbors import NearestNeighbors
except ImportError as exc:
    raise ImportError("scikit-learn is required for spatial-neighbour statistics") from exc


ROOT = Path(os.environ.get(
    "SPATIAL_V4_ROOT",
    "/home/anzhen_vhd/scRNAseq/codex/spatial_validation_v4_gpu_final",
))
SOURCE_ROOT = Path(os.environ.get(
    "SPATIAL_V3_SOURCE_ROOT",
    "/home/anzhen_vhd/scRNAseq/codex/spatial_validation_v3",
))
MAP_OUT = ROOT / "02_cell2location"
INPUT = ROOT / "00_inputs"
OUT = ROOT / "03_validation_targets_CoVarNet"
TABLE_OUT = OUT / "tables"
SPOT_OUT = OUT / "spot_level"
for directory in (OUT, TABLE_OUT, SPOT_OUT):
    directory.mkdir(parents=True, exist_ok=True)

MANIFEST_PATH = SOURCE_ROOT / "00_inputs" / "spatial_sample_manifest_curated.tsv"
LABEL_MAP_PATH = ROOT / "01_reference" / "spatial_reference_state_label_map.tsv"
RANDOM_SEED = 20260810
N_NULL = 1000
K_NEIGHBORS = 6

MODULE_WEIGHT_CANDIDATES = [
    INPUT / "module_node_weights_K14.tsv",
    SOURCE_ROOT / "00_inputs" / "module_node_weights_K14.tsv",
    Path("/home/anzhen_vhd/scRNAseq/codex/CoVarNet_final_v3/stable_NMF/module_node_weights_K14.tsv"),
    Path("/home/anzhen_vhd/scRNAseq/codex/covarnet/stable_NMF/module_node_weights_K14.tsv"),
]

PATHOLOGY_SIGNATURES = {
    "fibrous_cap_core": [
        "COL8A1", "ELN", "LOX", "LOXL1", "DCN", "BGN", "LUM", "THBS2",
        "TIMP1", "SERPINH1",
    ],
    "fibrous_cap_nonstructural_sensitivity": [
        "SERPINE1", "THBS1", "TGFBI", "CTGF", "TIMP1", "MMP2", "MMP14", "CXCL12",
    ],
    "calcification": ["RUNX2", "ALPL", "IBSP", "BGLAP", "SPP1", "MGP", "BMP2", "ENPP1"],
    "lipid_foam": ["TREM2", "LPL", "APOE", "APOC1", "GPNMB", "LGALS3", "CTSB", "CD68"],
    "inflammation": ["IL1B", "TNF", "CCL2", "CXCL8", "NFKBIA", "S100A8", "S100A9", "CTSS"],
    "endothelial_activation": ["VCAM1", "SELE", "ICAM1", "VWF", "KDR", "CXCL2"],
}

# These hypotheses are fixed from the cross-disease pseudobulk target screen.
# Spatial transcriptomics tests regional localisation and context dependence;
# it does not by itself establish cell-intrinsic differential expression.
CANDIDATE_HYPOTHESES = [
    {"gene": "CKB", "cell_context": "Macrophage", "context_level": "Level1", "pathology": "lipid_foam", "expected": -1, "tier": "primary_level1"},
    {"gene": "CKB", "cell_context": "Fibroblast", "context_level": "Level1", "pathology": "fibrous_cap_core", "expected": -1, "tier": "primary_level1"},
    {"gene": "CKB", "cell_context": "Macrophage | SPP1 lipid-associated macrophage", "context_level": "Level2", "pathology": "lipid_foam", "expected": -1, "tier": "primary_level2"},
    {"gene": "KCNK3", "cell_context": "Fibroblast", "context_level": "Level1", "pathology": "fibrous_cap_core", "expected": 1, "tier": "secondary_level1"},
    {"gene": "ZNF589", "cell_context": "Fibroblast", "context_level": "Level1", "pathology": "fibrous_cap_core", "expected": -1, "tier": "secondary_level1"},
]

# CM01 and CM02 are the two high-confidence pan-vascular modules from the final
# stability analysis. CM03 is retained as an explicitly exploratory candidate.
PRIMARY_MODULE_PATHOLOGY = [
    {"module": "CM01", "signature": "fibrous_cap_core", "expected": 1},
    {"module": "CM01", "signature": "calcification", "expected": 1},
    {"module": "CM02", "signature": "inflammation", "expected": 1},
    {"module": "CM03", "signature": "endothelial_activation", "expected": 1},
]


def module_priority(module: str) -> str:
    if module in {"CM01", "CM02"}:
        return "high_confidence_pan_vascular"
    if module == "CM03":
        return "candidate_rank_sensitive"
    return "internal_comparator"


def bh_fdr(values) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    output = np.full(len(p), np.nan)
    finite = np.isfinite(p)
    if not finite.any():
        return output
    pf = p[finite]
    order = np.argsort(pf)
    ranked = pf[order] * len(pf) / np.arange(1, len(pf) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    restored = np.empty_like(ranked)
    restored[order] = np.minimum(ranked, 1)
    output[finite] = restored
    return output


def zscore(values) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    result = np.full(values.shape, np.nan)
    finite = np.isfinite(values)
    if finite.sum() < 3:
        return result
    sd = np.nanstd(values[finite])
    if sd <= 0:
        result[finite] = 0
    else:
        result[finite] = (values[finite] - np.nanmean(values[finite])) / sd
    return result


def random_effects(beta, se) -> dict:
    beta = np.asarray(beta, dtype=float)
    se = np.asarray(se, dtype=float)
    keep = np.isfinite(beta) & np.isfinite(se) & (se > 0)
    beta, se = beta[keep], se[keep]
    k = len(beta)
    if k == 0:
        return {key: np.nan for key in ["beta", "se", "z", "p", "tau2", "I2"]} | {"k": 0}
    w = 1 / se**2
    fixed = np.sum(w * beta) / np.sum(w)
    q = np.sum(w * (beta - fixed) ** 2)
    c = np.sum(w) - np.sum(w**2) / np.sum(w)
    tau2 = max(0.0, (q - (k - 1)) / c) if k > 1 and c > 0 else 0.0
    wr = 1 / (se**2 + tau2)
    pooled = np.sum(wr * beta) / np.sum(wr)
    pooled_se = math.sqrt(1 / np.sum(wr))
    z = pooled / pooled_se
    p = 2 * norm.sf(abs(z))
    i2 = max(0.0, (q - (k - 1)) / q) * 100 if k > 1 and q > 0 else 0.0
    return {"k": k, "beta": pooled, "se": pooled_se, "z": z, "p": p, "tau2": tau2, "I2": i2}


def fixed_effect(beta, se) -> tuple[float, float]:
    beta = np.asarray(beta, dtype=float)
    se = np.asarray(se, dtype=float)
    keep = np.isfinite(beta) & np.isfinite(se) & (se > 0)
    if not keep.any():
        return np.nan, np.nan
    weights = 1 / se[keep] ** 2
    return float(np.sum(weights * beta[keep]) / weights.sum()), float(np.sqrt(1 / weights.sum()))


def strip_cell2location_prefix(value: str) -> str:
    prefixes = [
        "q05cell_abundance_w_sf_", "meanscell_abundance_w_sf_",
        "q05_cell_abundance_w_sf_", "means_cell_abundance_w_sf_",
    ]
    value = str(value)
    for prefix in prefixes:
        if value.startswith(prefix):
            return value[len(prefix):]
    return value


def abundance_frame(adata: ad.AnnData, posterior: str, safe_to_original: dict) -> pd.DataFrame:
    key = {"q05": "q05_cell_abundance_w_sf", "mean": "means_cell_abundance_w_sf"}[posterior]
    if key not in adata.obsm:
        raise KeyError(f"{key} absent from {list(adata.obsm.keys())}")
    values = adata.obsm[key]
    if isinstance(values, pd.DataFrame):
        result = values.copy()
    else:
        names = adata.uns.get("mod", {}).get("factor_names")
        if names is None or len(names) != values.shape[1]:
            raise ValueError("Cannot recover cell2location factor names")
        result = pd.DataFrame(values, index=adata.obs_names, columns=names)
    safe = [strip_cell2location_prefix(x) for x in result.columns]
    result.columns = [safe_to_original.get(x, x) for x in safe]
    if result.columns.duplicated().any():
        result = result.T.groupby(level=0).sum().T
    return result.astype(float)


def lineage_from_state(value: str) -> str:
    value = str(value)
    if " | " in value:
        return value.split(" | ", 1)[0]
    if value.startswith("Other_"):
        return value.replace("Other_", "", 1)
    return value


def lineage_abundance(abundance: pd.DataFrame) -> pd.DataFrame:
    lineages: dict[str, list[str]] = {}
    for state in abundance.columns:
        lineages.setdefault(lineage_from_state(state), []).append(state)
    return pd.DataFrame(
        {lineage: abundance[states].sum(axis=1) for lineage, states in lineages.items()},
        index=abundance.index,
    )


def normalized_gene_matrix(adata: ad.AnnData, genes: list[str]) -> tuple[pd.DataFrame, pd.Series]:
    matrix = adata.layers["counts"] if "counts" in adata.layers else adata.X
    if "total_counts_original" in adata.obs:
        library = pd.to_numeric(adata.obs["total_counts_original"], errors="coerce").to_numpy(float)
    else:
        library = np.asarray(matrix.sum(axis=1)).ravel().astype(float)
    positions = {str(g).upper(): i for i, g in enumerate(adata.var_names)}
    valid = [gene for gene in genes if gene in positions]
    if not valid:
        return pd.DataFrame(index=adata.obs_names), pd.Series(library, index=adata.obs_names)
    columns = [positions[gene] for gene in valid]
    block = matrix[:, columns]
    if sp.issparse(block):
        block = block.toarray()
    block = np.asarray(block, dtype=float)
    block = np.log1p(block * (1e4 / np.maximum(library, 1))[:, None])
    return pd.DataFrame(block, index=adata.obs_names, columns=valid), pd.Series(library, index=adata.obs_names)


def score_signatures(expression: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
    output = pd.DataFrame(index=expression.index)
    coverage = []
    for name, genes in PATHOLOGY_SIGNATURES.items():
        valid = [gene for gene in genes if gene in expression.columns]
        coverage.append({"signature": name, "n_requested": len(genes), "n_valid": len(valid), "valid_genes": ",".join(valid)})
        if len(valid) < 3:
            output[name] = np.nan
            continue
        standardized = np.column_stack([zscore(expression[gene]) for gene in valid])
        output[name] = np.nanmean(standardized, axis=1)
    return output, coverage


def spatial_blocks(obs: pd.DataFrame) -> np.ndarray:
    if {"array_row", "array_col"}.issubset(obs.columns):
        row = pd.to_numeric(obs["array_row"], errors="coerce").fillna(0).astype(int)
        col = pd.to_numeric(obs["array_col"], errors="coerce").fillna(0).astype(int)
        return (row // 5).astype(str).str.cat((col // 5).astype(str), sep="_").to_numpy()
    coords = np.asarray(obs.index.to_series().rank()).reshape(-1, 1)
    return (coords.ravel() // 50).astype(int)


def fit_spatial_count_model(
    counts: np.ndarray,
    library: np.ndarray,
    abundance: np.ndarray,
    pathology: np.ndarray,
    blocks: np.ndarray,
) -> list[dict]:
    finite = np.isfinite(counts) & np.isfinite(library) & np.isfinite(abundance) & np.isfinite(pathology)
    if finite.sum() < 100 or np.sum(counts[finite] > 0) < 20:
        return []
    a = zscore(abundance[finite])
    p = zscore(pathology[finite])
    design = pd.DataFrame({"constant": 1.0, "cell_abundance": a, "pathology": p, "interaction": a * p})
    try:
        fit = sm.GLM(
            counts[finite],
            design,
            family=sm.families.Poisson(link=sm.families.links.Log()),
            offset=np.log(np.maximum(library[finite], 1)),
        ).fit(cov_type="cluster", cov_kwds={"groups": blocks[finite]})
    except Exception as exc:
        warnings.warn(f"Cluster-robust Poisson failed: {exc!r}")
        return []
    rows = []
    for term in ["cell_abundance", "pathology", "interaction"]:
        rows.append({
            "term": term,
            "beta": float(fit.params[term]),
            "se": float(fit.bse[term]),
            "z": float(fit.tvalues[term]),
            "p": float(fit.pvalues[term]),
            "n_spots": int(finite.sum()),
            "n_positive_spots": int(np.sum(counts[finite] > 0)),
            "model": "Poisson_loglink_spatial_block_cluster_robust",
        })
    return rows


def coordinates(adata: ad.AnnData, index: np.ndarray) -> np.ndarray:
    if "spatial" in adata.obsm:
        return np.asarray(adata.obsm["spatial"])[index]
    obs = adata.obs.iloc[index]
    for columns in [["pxl_col_in_fullres", "pxl_row_in_fullres"], ["array_col", "array_row"]]:
        if set(columns).issubset(obs.columns):
            return obs[columns].to_numpy(float)
    raise ValueError("No spatial coordinates found")


def neighbour_lag(values: np.ndarray, coords: np.ndarray, k: int = K_NEIGHBORS) -> np.ndarray:
    n = len(coords)
    use_k = min(k + 1, n)
    indices = NearestNeighbors(n_neighbors=use_k).fit(coords).kneighbors(return_distance=False)
    if indices.shape[1] > 1:
        indices = indices[:, 1:]
    return np.nanmean(values[indices], axis=1)


def matrix_correlations(values: np.ndarray) -> np.ndarray:
    standardized = np.column_stack([zscore(values[:, j]) for j in range(values.shape[1])])
    standardized = np.nan_to_num(standardized)
    denominator = max(1, standardized.shape[0] - 1)
    corr = (standardized.T @ standardized) / denominator
    return np.clip(corr, -1, 1)


def upper_mean(matrix: np.ndarray, indices: list[int]) -> float:
    if len(indices) < 2:
        return np.nan
    block = matrix[np.ix_(indices, indices)]
    return float(np.nanmean(block[np.triu_indices(len(indices), 1)]))


def spatial_cross_matrix(values: np.ndarray, coords: np.ndarray) -> np.ndarray:
    standardized = np.column_stack([zscore(values[:, j]) for j in range(values.shape[1])])
    standardized = np.nan_to_num(standardized)
    lagged = neighbour_lag(standardized, coords)
    n = max(1, len(coords) - 1)
    cross = standardized.T @ lagged / n
    return (cross + cross.T) / 2


def moran_i(values: np.ndarray, coords: np.ndarray) -> float:
    z = np.nan_to_num(zscore(values))
    lag = neighbour_lag(z[:, None], coords).ravel()
    denominator = np.sum(z**2)
    return float(np.sum(z * lag) / denominator) if denominator > 0 else np.nan


def module_scores(abundance: pd.DataFrame, weights: pd.DataFrame) -> pd.DataFrame:
    common = abundance.columns.intersection(weights.index)
    if len(common) < 2:
        raise ValueError("Too few CoVarNet states matched cell2location abundance")
    positive = abundance[common].to_numpy()
    positive_values = positive[positive > 0]
    pseudocount = np.median(positive_values) * 0.01 if positive_values.size else 1e-4
    log_abundance = np.log(positive + max(pseudocount, 1e-8))
    clr = log_abundance - log_abundance.mean(axis=1, keepdims=True)
    result = {}
    for module in weights.columns:
        w = weights.loc[common, module].to_numpy(float)
        if w.sum() <= 0:
            result[module] = np.nan
        else:
            result[module] = clr @ (w / w.sum())
    return pd.DataFrame(result, index=abundance.index)


def sample_matched_null(
    module_states: list[str],
    available_states: list[str],
    rng: np.random.Generator,
) -> list[str]:
    forbidden = set(module_states)
    selected = []
    for lineage in [lineage_from_state(x) for x in module_states]:
        pool = [x for x in available_states if lineage_from_state(x) == lineage and x not in forbidden and x not in selected]
        if not pool:
            pool = [x for x in available_states if x not in forbidden and x not in selected]
        if not pool:
            return []
        selected.append(str(rng.choice(pool)))
    return selected


def resolve_module_weights() -> Path:
    hits = [path for path in MODULE_WEIGHT_CANDIDATES if path.exists()]
    if not hits:
        raise FileNotFoundError(
            "Upload module_node_weights_K14.tsv to "
            f"{INPUT}; searched: {[str(x) for x in MODULE_WEIGHT_CANDIDATES]}"
        )
    return hits[0]


manifest = pd.read_csv(MANIFEST_PATH, sep="\t")
label_map = pd.read_csv(LABEL_MAP_PATH, sep="\t")
safe_to_original = dict(zip(label_map["spatial_reference_state_safe"], label_map["spatial_reference_state_original"]))
weight_path = resolve_module_weights()
weights = pd.read_csv(weight_path, sep="\t").set_index("cell_state")
weights = weights.apply(pd.to_numeric, errors="coerce").fillna(0)
weights.to_csv(INPUT / "module_node_weights_K14.used.tsv", sep="\t")

all_required_genes = sorted(
    {item["gene"] for item in CANDIDATE_HYPOTHESES}
    | set(itertools.chain.from_iterable(PATHOLOGY_SIGNATURES.values()))
)

candidate_rows = []
detection_rows = []
module_rows = []
coverage_rows = []
module_match_rows = []
spot_files = []
rng = np.random.default_rng(RANDOM_SEED)

mapped_paths = sorted(MAP_OUT.glob("GSE*__final_atlas_cell2location.h5ad"))
if not mapped_paths:
    raise FileNotFoundError(f"No mapped h5ad files in {MAP_OUT}; run script 01 first")

for path in mapped_paths:
    print("Validating", path.name)
    adata = ad.read_h5ad(path)
    # Curated donor/condition metadata can be corrected after deconvolution
    # without rerunning the expensive cell2location models.
    manifest_by_section = manifest.drop_duplicates("section_id").set_index("section_id")
    section_ids = adata.obs["section_id"].astype(str)
    for metadata_column in [
        "biological_unit_id", "metadata_status", "condition", "tissue", "dataset_role"
    ]:
        if metadata_column in manifest_by_section.columns:
            mapped_metadata = section_ids.map(manifest_by_section[metadata_column])
            if mapped_metadata.notna().any():
                adata.obs[metadata_column] = mapped_metadata.fillna(
                    adata.obs[metadata_column].astype("object").astype(str)
                ).to_numpy()
    expression, library = normalized_gene_matrix(adata, all_required_genes)
    pathology_scores, coverage = score_signatures(expression)
    gse = str(adata.obs["GSE"].iloc[0])
    coverage_rows.extend([{**row, "GSE": gse} for row in coverage])

    for posterior in ["q05", "mean"]:
        abundance = abundance_frame(adata, posterior, safe_to_original)
        lineage = lineage_abundance(abundance)
        scores = module_scores(abundance, weights)
        for module in weights.columns:
            weighted_states = weights.index[weights[module] > 0].tolist()
            mapped_states = [state for state in weighted_states if state in abundance.columns]
            module_match_rows.append({
                "GSE": gse,
                "posterior": posterior,
                "module": module,
                "n_weighted_states": len(weighted_states),
                "n_mapped_states": len(mapped_states),
                "mapped_weight_fraction": float(weights.loc[mapped_states, module].sum() / weights[module].sum())
                if weights[module].sum() > 0 else np.nan,
                "missing_states": ";".join(sorted(set(weighted_states) - set(mapped_states))),
            })

        for section_id, positions_idx in adata.obs.groupby("section_id", observed=True).indices.items():
            idx = np.asarray(positions_idx, dtype=int)
            obs = adata.obs.iloc[idx]
            section_coords = coordinates(adata, idx)
            section_abundance = abundance.iloc[idx]
            section_lineage = lineage.iloc[idx]
            section_scores = scores.iloc[idx]
            section_pathology = pathology_scores.iloc[idx]
            section_expression = expression.iloc[idx]
            section_library = library.iloc[idx]
            blocks = spatial_blocks(obs)
            metadata = {
                "GSE": gse,
                "section_id": str(section_id),
                "sample_id": str(obs["sample_id"].iloc[0]),
                "biological_unit_id": str(obs["biological_unit_id"].iloc[0]),
                "metadata_status": str(obs["metadata_status"].iloc[0]),
                "condition": str(obs["condition"].iloc[0]),
                "tissue": str(obs["tissue"].iloc[0]),
                "posterior": posterior,
            }

            # Candidate target expression, cell-context abundance and pathology interaction.
            raw_matrix = adata.layers["counts"] if "counts" in adata.layers else adata.X
            gene_pos = {str(g).upper(): j for j, g in enumerate(adata.var_names)}
            for hypothesis in CANDIDATE_HYPOTHESES:
                gene = hypothesis["gene"]
                context = hypothesis["cell_context"]
                pathology = hypothesis["pathology"]
                context_frame = section_abundance if hypothesis["context_level"] == "Level2" else section_lineage
                if gene not in gene_pos or context not in context_frame or pathology not in section_pathology:
                    continue
                y = raw_matrix[idx, gene_pos[gene]]
                y = y.toarray().ravel() if sp.issparse(y) else np.asarray(y).ravel()
                detect = float(np.mean(y > 0))
                detection_rows.append({
                    **metadata, **hypothesis,
                    "n_spots": len(idx), "n_positive_spots": int(np.sum(y > 0)),
                    "detection_fraction": detect,
                    "mean_logCP10K": float(section_expression[gene].mean()) if gene in section_expression else np.nan,
                })
                model_rows = fit_spatial_count_model(
                    y, section_library.to_numpy(float),
                    context_frame[context].to_numpy(float),
                    section_pathology[pathology].to_numpy(float), blocks,
                )
                for model_row in model_rows:
                    candidate_rows.append({**metadata, **hypothesis, **model_row})

            # All-module projection and pre-specified spatial coherence/null tests.
            available_states = list(section_abundance.columns)
            log_values = np.log1p(section_abundance.to_numpy(float))
            corr = matrix_correlations(log_values)
            spatial_cross = spatial_cross_matrix(log_values, section_coords)
            state_index = {state: j for j, state in enumerate(available_states)}
            for module in weights.columns:
                ordered = weights[module].sort_values(ascending=False)
                top_states = [state for state in ordered.index if ordered[state] > 0 and state in state_index][:10]
                top_idx = [state_index[state] for state in top_states]
                coherence = upper_mean(corr, top_idx)
                cross_coherence = upper_mean(spatial_cross, top_idx)
                null_coherence = []
                null_cross = []
                for _ in range(N_NULL):
                    null_states = sample_matched_null(top_states, available_states, rng)
                    if len(null_states) != len(top_states):
                        continue
                    null_idx = [state_index[state] for state in null_states]
                    null_coherence.append(upper_mean(corr, null_idx))
                    null_cross.append(upper_mean(spatial_cross, null_idx))
                null_coherence = np.asarray(null_coherence, dtype=float)
                null_cross = np.asarray(null_cross, dtype=float)
                null_coherence = null_coherence[np.isfinite(null_coherence)]
                null_cross = null_cross[np.isfinite(null_cross)]
                coherence_null_sd = float(np.nanstd(null_coherence))
                cross_null_sd = float(np.nanstd(null_cross))
                score = section_scores[module].to_numpy(float)
                row = {
                    **metadata,
                    "module": module,
                    "module_priority": module_priority(module),
                    "n_top_nodes_mapped": len(top_states),
                    "top_nodes_mapped": ";".join(top_states),
                    "within_spot_coherence": coherence,
                    "spatial_lag_coherence": cross_coherence,
                    "within_spot_null_mean": float(np.nanmean(null_coherence)),
                    "within_spot_null_sd": coherence_null_sd,
                    "within_spot_null_z": (coherence - float(np.nanmean(null_coherence))) / coherence_null_sd
                    if coherence_null_sd > 0 else np.nan,
                    "within_spot_empirical_p": float((1 + np.sum(null_coherence >= coherence)) / (1 + len(null_coherence))),
                    "spatial_lag_null_mean": float(np.nanmean(null_cross)),
                    "spatial_lag_null_sd": cross_null_sd,
                    "spatial_lag_null_z": (cross_coherence - float(np.nanmean(null_cross))) / cross_null_sd
                    if cross_null_sd > 0 else np.nan,
                    "spatial_lag_empirical_p": float((1 + np.sum(null_cross >= cross_coherence)) / (1 + len(null_cross))),
                    "moran_I": moran_i(score, section_coords),
                    "median_module_score": float(np.nanmedian(score)),
                    "p90_module_score": float(np.nanquantile(score, 0.9)),
                }
                for signature in PATHOLOGY_SIGNATURES:
                    x = section_pathology[signature].to_numpy(float)
                    finite = np.isfinite(score) & np.isfinite(x)
                    rho, p_value = spearmanr(score[finite], x[finite]) if finite.sum() >= 50 else (np.nan, np.nan)
                    row[f"rho__{signature}"] = rho
                    row[f"p__{signature}"] = p_value
                module_rows.append(row)

            if posterior == "q05":
                spot = pd.DataFrame({
                    "spot_id": obs.index.astype(str),
                    "GSE": gse,
                    "section_id": str(section_id),
                    "sample_id": metadata["sample_id"],
                    "biological_unit_id": metadata["biological_unit_id"],
                    "metadata_status": metadata["metadata_status"],
                    "condition": metadata["condition"],
                    "tissue": metadata["tissue"],
                    "x": section_coords[:, 0], "y": section_coords[:, 1],
                    "library_size": section_library.to_numpy(float),
                })
                for gene in {x["gene"] for x in CANDIDATE_HYPOTHESES}:
                    spot[f"expr__{gene}"] = section_expression[gene].to_numpy(float) if gene in section_expression else np.nan
                for signature in PATHOLOGY_SIGNATURES:
                    spot[f"pathology__{signature}"] = section_pathology[signature].to_numpy(float)
                for lineage_name in ["Fibroblast", "Fibromyocyte", "Smooth muscle cell", "Pericyte", "Macrophage", "T cell", "Natural killer cell"]:
                    spot[f"abundance__{lineage_name}"] = section_lineage[lineage_name].to_numpy(float) if lineage_name in section_lineage else 0
                for exact_state in sorted({
                    item["cell_context"] for item in CANDIDATE_HYPOTHESES
                    if item["context_level"] == "Level2"
                }):
                    spot[f"abundance__{exact_state}"] = section_abundance[exact_state].to_numpy(float) if exact_state in section_abundance else 0
                for module in weights.columns:
                    spot[f"module__{module}"] = section_scores[module].to_numpy(float)
                safe_section = re.sub(r"[^0-9A-Za-z_.-]+", "_", str(section_id))
                spot_path = SPOT_OUT / f"{safe_section}__spot_scores.tsv.gz"
                spot.to_csv(spot_path, sep="\t", index=False, compression="gzip")
                spot_files.append(str(spot_path))
    del adata

candidate = pd.DataFrame(candidate_rows)
detection = pd.DataFrame(detection_rows)
modules = pd.DataFrame(module_rows)
coverage = pd.DataFrame(coverage_rows)
module_match = pd.DataFrame(module_match_rows)
candidate.to_csv(TABLE_OUT / "candidate_target_spatial_models_by_section.tsv", sep="\t", index=False)
detection.to_csv(TABLE_OUT / "candidate_target_detection_by_section.tsv", sep="\t", index=False)
modules.to_csv(TABLE_OUT / "CoVarNet_all14_spatial_validation_by_section.tsv", sep="\t", index=False)
coverage.to_csv(TABLE_OUT / "pathology_signature_gene_coverage.tsv", sep="\t", index=False)
module_match.to_csv(TABLE_OUT / "CoVarNet_module_state_mapping_audit.tsv", sep="\t", index=False)
pd.DataFrame({"spot_score_file": sorted(set(spot_files))}).to_csv(
    SPOT_OUT / "spot_score_manifest.tsv", sep="\t", index=False
)

# Capture audit is unique per target/section and separates detection from
# positive-spot intensity. It prevents low-capture sections from driving a
# biological interpretation.
if not detection.empty:
    detection_unique = detection.drop_duplicates(
        ["gene", "GSE", "section_id", "posterior"]
    ).copy()
    detection_unique["model_eligible"] = (
        (detection_unique["n_spots"] >= 100)
        & (detection_unique["n_positive_spots"] >= 20)
    )
    detection_unique.to_csv(
        TABLE_OUT / "target_capture_audit_by_section.tsv", sep="\t", index=False
    )
    target_capture_summary = (
        detection_unique[detection_unique["posterior"] == "q05"]
        .groupby(["gene", "GSE"], observed=True)
        .agg(
            n_sections=("section_id", "nunique"),
            n_spots=("n_spots", "sum"),
            n_positive_spots=("n_positive_spots", "sum"),
            median_section_detection=("detection_fraction", "median"),
            median_section_logCP10K=("mean_logCP10K", "median"),
            n_model_eligible_sections=("model_eligible", "sum"),
        )
        .reset_index()
    )
    target_capture_summary.to_csv(
        TABLE_OUT / "target_capture_summary_by_dataset.tsv", sep="\t", index=False
    )

# Candidate meta-analysis: primary q05 posterior; biological units with non-ambiguous metadata.
meta_rows = []
if not candidate.empty:
    primary = candidate[candidate["posterior"] == "q05"].copy()
    for scope, scoped in {
        "confirmed_or_unique_biological_units": primary[~primary["metadata_status"].str.contains("needs_", na=False)],
        "all_sections_sensitivity": primary,
    }.items():
        group_columns = ["gene", "cell_context", "context_level", "pathology", "expected", "tier", "term"]
        for keys, block in scoped.groupby(group_columns, observed=True):
            if scope == "confirmed_or_unique_biological_units":
                unit_effects = []
                for _, unit in block.groupby("biological_unit_id", observed=True):
                    beta, se = fixed_effect(unit["beta"], unit["se"])
                    unit_effects.append((beta, se))
                beta = [x[0] for x in unit_effects]
                se = [x[1] for x in unit_effects]
            else:
                beta, se = block["beta"], block["se"]
            pooled = random_effects(beta, se)
            row = dict(zip(group_columns, keys))
            meta_rows.append({"scope": scope, **row, **pooled})
candidate_meta = pd.DataFrame(meta_rows)
if not candidate_meta.empty:
    candidate_meta["FDR_within_scope_term"] = (
        candidate_meta.groupby(["scope", "term"], observed=True)["p"].transform(bh_fdr)
    )
    candidate_meta["direction_matches_pseudobulk"] = np.where(
        candidate_meta["expected"] == 0, np.nan,
        np.sign(candidate_meta["beta"]) == candidate_meta["expected"],
    )
candidate_meta.to_csv(TABLE_OUT / "candidate_target_spatial_meta_analysis.tsv", sep="\t", index=False)

# q05 is primary; posterior mean is a sensitivity analysis. Report whether
# section-level coefficient directions are retained rather than treating the
# two posteriors as independent evidence.
sensitivity_rows = []
if not candidate.empty:
    merge_keys = [
        "GSE", "section_id", "gene", "cell_context", "context_level",
        "pathology", "term",
    ]
    wide = candidate.pivot_table(
        index=merge_keys, columns="posterior", values="beta", aggfunc="first"
    ).reset_index()
    if {"q05", "mean"}.issubset(wide.columns):
        for keys, block in wide.groupby(
            ["gene", "cell_context", "context_level", "pathology", "term"],
            observed=True,
        ):
            finite = np.isfinite(block["q05"]) & np.isfinite(block["mean"])
            rho, p_value = spearmanr(
                block.loc[finite, "q05"], block.loc[finite, "mean"]
            ) if finite.sum() >= 3 else (np.nan, np.nan)
            sensitivity_rows.append({
                **dict(zip(
                    ["gene", "cell_context", "context_level", "pathology", "term"],
                    keys,
                )),
                "n_sections": int(finite.sum()),
                "spearman_q05_vs_mean": rho,
                "spearman_p": p_value,
                "direction_agreement_fraction": float(np.mean(
                    np.sign(block.loc[finite, "q05"]) == np.sign(block.loc[finite, "mean"])
                )) if finite.any() else np.nan,
            })
pd.DataFrame(sensitivity_rows).to_csv(
    TABLE_OUT / "candidate_target_q05_vs_mean_sensitivity.tsv", sep="\t", index=False
)

# Cross-dataset direction consistency treats each GEO study, rather than each spot or section, as one unit.
candidate_consistency_rows = []
if not candidate.empty:
    interaction = candidate[(candidate["posterior"] == "q05") & (candidate["term"] == "interaction")]
    group_columns = ["gene", "cell_context", "context_level", "pathology", "expected", "tier"]
    for keys, block in interaction.groupby(group_columns, observed=True):
        dataset_effect = block.groupby("GSE", observed=True)["beta"].median().dropna()
        expected = keys[4]
        if expected == 0:
            n_match, p_value = np.nan, np.nan
        else:
            n_match = int(np.sum(np.sign(dataset_effect) == expected))
            p_value = float(binomtest(n_match, len(dataset_effect), 0.5, alternative="greater").pvalue) if len(dataset_effect) else np.nan
        row = dict(zip(group_columns, keys))
        candidate_consistency_rows.append({
            **row,
            "n_datasets": len(dataset_effect),
            "n_datasets_matching_expected_direction": n_match,
            "direction_consistency_fraction": n_match / len(dataset_effect) if len(dataset_effect) and np.isfinite(n_match) else np.nan,
            "one_sided_binomial_p": p_value,
            "dataset_effects": ";".join(f"{gse}:{value:.4g}" for gse, value in dataset_effect.items()),
        })
candidate_consistency = pd.DataFrame(candidate_consistency_rows)
if not candidate_consistency.empty:
    candidate_consistency["FDR"] = bh_fdr(candidate_consistency["one_sided_binomial_p"])
candidate_consistency.to_csv(
    TABLE_OUT / "candidate_target_cross_dataset_direction_consistency.tsv", sep="\t", index=False
)

# Module summary across sections. P-values are empirical null tests; pathology associations are descriptive until donor IDs are confirmed.
module_summary_rows = []
if not modules.empty:
    primary_modules = modules[modules["posterior"] == "q05"].copy()
    for module, block in primary_modules.groupby("module", observed=True):
        dataset_z_within = block.groupby("GSE", observed=True)["within_spot_null_z"].median().dropna()
        dataset_z_lag = block.groupby("GSE", observed=True)["spatial_lag_null_z"].median().dropna()
        stouffer_within = dataset_z_within.sum() / np.sqrt(len(dataset_z_within)) if len(dataset_z_within) else np.nan
        stouffer_lag = dataset_z_lag.sum() / np.sqrt(len(dataset_z_lag)) if len(dataset_z_lag) else np.nan
        row = {
            "module": module,
            "module_priority": block["module_priority"].iloc[0],
            "n_sections": block["section_id"].nunique(),
            "n_datasets": block["GSE"].nunique(),
            "median_within_spot_coherence": block["within_spot_coherence"].median(),
            "fraction_within_spot_empirical_p_lt_0.05": np.mean(block["within_spot_empirical_p"] < 0.05),
            "median_spatial_lag_coherence": block["spatial_lag_coherence"].median(),
            "fraction_spatial_lag_empirical_p_lt_0.05": np.mean(block["spatial_lag_empirical_p"] < 0.05),
            "median_moran_I": block["moran_I"].median(),
            "cross_dataset_stouffer_z_within_spot": stouffer_within,
            "cross_dataset_p_within_spot": float(norm.sf(stouffer_within)) if np.isfinite(stouffer_within) else np.nan,
            "cross_dataset_stouffer_z_spatial_lag": stouffer_lag,
            "cross_dataset_p_spatial_lag": float(norm.sf(stouffer_lag)) if np.isfinite(stouffer_lag) else np.nan,
        }
        for signature in PATHOLOGY_SIGNATURES:
            row[f"median_rho__{signature}"] = block[f"rho__{signature}"].median()
            row[f"positive_fraction__{signature}"] = np.mean(block[f"rho__{signature}"] > 0)
        module_summary_rows.append(row)
module_summary = pd.DataFrame(module_summary_rows)
if not module_summary.empty:
    module_summary["cross_dataset_FDR_within_spot"] = bh_fdr(module_summary["cross_dataset_p_within_spot"])
    module_summary["cross_dataset_FDR_spatial_lag"] = bh_fdr(module_summary["cross_dataset_p_spatial_lag"])
module_summary.to_csv(TABLE_OUT / "CoVarNet_all14_spatial_validation_summary.tsv", sep="\t", index=False)

module_pathology_rows = []
if not modules.empty:
    primary_modules = modules[modules["posterior"] == "q05"]
    for hypothesis in PRIMARY_MODULE_PATHOLOGY:
        block = primary_modules[primary_modules["module"] == hypothesis["module"]]
        column = f"rho__{hypothesis['signature']}"
        dataset_effect = block.groupby("GSE", observed=True)[column].median().dropna()
        n_match = int(np.sum(np.sign(dataset_effect) == hypothesis["expected"]))
        p_value = float(binomtest(n_match, len(dataset_effect), 0.5, alternative="greater").pvalue) if len(dataset_effect) else np.nan
        module_pathology_rows.append({
            **hypothesis,
            "n_datasets": len(dataset_effect),
            "n_datasets_matching_expected_direction": n_match,
            "direction_consistency_fraction": n_match / len(dataset_effect) if len(dataset_effect) else np.nan,
            "one_sided_binomial_p": p_value,
            "median_dataset_rho": float(dataset_effect.median()) if len(dataset_effect) else np.nan,
            "dataset_effects": ";".join(f"{gse}:{value:.4g}" for gse, value in dataset_effect.items()),
        })
module_pathology = pd.DataFrame(module_pathology_rows)
if not module_pathology.empty:
    module_pathology["FDR"] = bh_fdr(module_pathology["one_sided_binomial_p"])
module_pathology.to_csv(
    TABLE_OUT / "CoVarNet_prespecified_pathology_cross_dataset_consistency.tsv", sep="\t", index=False
)

run_summary = {
    "mapped_files": [str(path) for path in mapped_paths],
    "module_weights": str(weight_path),
    "n_candidate_section_models": int(len(candidate)),
    "n_module_section_results": int(len(modules)),
    "n_spot_score_files": int(len(set(spot_files))),
    "high_confidence_modules": ["CM01", "CM02"],
    "candidate_module": "CM03",
    "candidate_targets": ["CKB", "KCNK3", "ZNF589"],
}
(OUT / "validation_run_summary.json").write_text(
    json.dumps(run_summary, indent=2, ensure_ascii=False), encoding="utf-8"
)
print("Spatial validation completed:", OUT)
