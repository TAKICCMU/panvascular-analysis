from __future__ import annotations

import json
import math
import os
import platform
import re
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import binomtest, norm, spearmanr

try:
    import statsmodels.api as sm
except ImportError as exc:
    raise ImportError("Install statsmodels in this Jupyter environment") from exc


# =============================================================================
# 0. Configuration
# =============================================================================

ROOT = Path(os.environ.get(
    "SPATIAL_V4_ROOT",
    "/path/to/panvas_project/spatial_validation_v4_gpu_final",
))
SOURCE_ROOT = Path(os.environ.get(
    "SPATIAL_V3_SOURCE_ROOT",
    "/path/to/panvas_project/spatial_validation_v3",
))
MAP_OUT = ROOT / "02_cell2location"
OUT = ROOT / "05_all_public_Tier1_spatial_validation_v1"
AUDIT_OUT = OUT / "00_audit"
TABLE_OUT = OUT / "01_tables"
FIGURE_OUT = OUT / "02_figures"
SOURCE_DATA_OUT = OUT / "03_source_data"
for directory in (OUT, AUDIT_OUT, TABLE_OUT, FIGURE_OUT, SOURCE_DATA_OUT):
    directory.mkdir(parents=True, exist_ok=True)

# Upload the frozen submitted roster here if available. The embedded roster
# below is an exact fallback for the 22 submitted gene-cell-state pairs.
TIER1_FILE_CANDIDATES = [
    ROOT / "00_inputs" / "Fig3d_Tier1_pairs_redacted.tsv",
    Path("/path/to/panvas_project/Fig3d_Tier1_pairs_redacted.tsv"),
    Path(
        "/path/to/panvas_project/pseudobulk_target_discovery_v1/"
        "edgeR_meta_results_v4/prioritised_target_candidates.tsv"
    ),
]

MANIFEST_PATH = SOURCE_ROOT / "00_inputs" / "spatial_sample_manifest_curated.tsv"
LABEL_MAP_PATH = ROOT / "01_reference" / "spatial_reference_state_label_map.tsv"
PRIMARY_POSTERIOR = "q05"
POSTERIORS = ("q05", "mean")
MIN_SPOTS = 100
MIN_POSITIVE_SPOTS = 20
RANDOM_SEED = 20260823
EXPORT_SPOT_LEVEL = False  # Set True only when journal source data require it.

COLORS = {
    "blue": "#2C5F8A",
    "teal": "#2A9D8F",
    "red": "#C44E52",
    "orange": "#E69F00",
    "grey": "#8B949E",
    "light": "#EAF0F4",
    "ink": "#1D2733",
}

PATHOLOGY_SIGNATURES = {
    "fibrous_cap_core": [
        "COL8A1", "ELN", "LOX", "LOXL1", "DCN", "BGN", "LUM", "THBS2",
        "TIMP1", "SERPINH1",
    ],
    "lipid_foam": [
        "TREM2", "LPL", "APOE", "APOC1", "GPNMB", "LGALS3", "CTSB", "CD68",
    ],
    "endothelial_activation": ["VCAM1", "SELE", "ICAM1", "VWF", "KDR", "CXCL2"],
    "inflammation": ["IL1B", "TNF", "CCL2", "CXCL8", "NFKBIA", "S100A8", "S100A9", "CTSS"],
}

# Frozen submitted Tier 1 roster: 22 pairs, 21 genes.
EMBEDDED_TIER1 = [
    ("HP", "Level1", "Fibroblast", -1),
    ("LRRC10B", "Level2", "Smooth muscle cell | TFF1 valve-associated SMC", -1),
    ("C17ORF58", "Level1", "Fibroblast", -1),
    ("KCNK5", "Level1", "Fibroblast", -1),
    ("CORO6", "Level1", "Fibroblast", -1),
    ("SREBF1", "Level1", "Fibroblast", -1),
    ("ZNF589", "Level1", "Fibroblast", -1),
    ("CCNL1", "Level1", "Fibroblast", -1),
    ("HDAC9", "Level2", "Macrophage | resident-like macrophage", -1),
    ("ATP1B1", "Level2", "Macrophage | resident-like macrophage", -1),
    ("SCARB1", "Level2", "Macrophage | resident-like macrophage", -1),
    ("ARHGEF12", "Level2", "Macrophage | resident-like macrophage", -1),
    ("ARHGEF12", "Level1", "Macrophage", -1),
    ("BAG5", "Level1", "Fibroblast", -1),
    ("SBF2", "Level1", "Fibroblast", -1),
    ("BCAP29", "Level1", "Fibroblast", 1),
    ("SERPINH1", "Level1", "Fibroblast", 1),
    ("DCLRE1B", "Level1", "Fibroblast", 1),
    ("GMDS", "Level1", "Fibroblast", 1),
    ("FHL3", "Level1", "Fibroblast", 1),
    ("FN1", "Level1", "Fibroblast", 1),
    ("KCNK3", "Level1", "Fibroblast", 1),
]


# =============================================================================
# 1. General utilities
# =============================================================================

def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def bh_fdr(values) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    output = np.full(len(p), np.nan)
    finite = np.isfinite(p)
    if not finite.any():
        return output
    pf = p[finite]
    order = np.argsort(pf)
    adjusted = pf[order] * len(pf) / np.arange(1, len(pf) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1)
    output[finite] = restored
    return output


def zscore(values) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    output = np.full(values.shape, np.nan)
    finite = np.isfinite(values)
    if finite.sum() < 3:
        return output
    sd = np.nanstd(values[finite])
    output[finite] = 0 if sd <= 0 else (values[finite] - np.nanmean(values[finite])) / sd
    return output


def fixed_effect(beta, se) -> tuple[float, float]:
    beta = np.asarray(beta, dtype=float)
    se = np.asarray(se, dtype=float)
    keep = np.isfinite(beta) & np.isfinite(se) & (se > 0)
    if not keep.any():
        return np.nan, np.nan
    weights = 1 / se[keep] ** 2
    return float(np.sum(weights * beta[keep]) / weights.sum()), float(np.sqrt(1 / weights.sum()))


def random_effects(beta, se) -> dict:
    beta = np.asarray(beta, dtype=float)
    se = np.asarray(se, dtype=float)
    keep = np.isfinite(beta) & np.isfinite(se) & (se > 0)
    beta, se = beta[keep], se[keep]
    k = len(beta)
    empty = {"k": 0, "beta": np.nan, "se": np.nan, "z": np.nan, "p": np.nan,
             "CI_low": np.nan, "CI_high": np.nan, "tau2": np.nan, "I2": np.nan}
    if k == 0:
        return empty
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
    return {
        "k": k, "beta": pooled, "se": pooled_se, "z": z, "p": p,
        "CI_low": pooled - 1.96 * pooled_se, "CI_high": pooled + 1.96 * pooled_se,
        "tau2": tau2, "I2": i2,
    }


def primary_pathology(parent_lineage: str) -> str:
    if parent_lineage in {"Macrophage", "Monocyte"}:
        return "lipid_foam"
    if parent_lineage == "Endothelial cell":
        return "endothelial_activation"
    if parent_lineage in {
        "Fibroblast", "Fibromyocyte", "Smooth muscle cell", "Pericyte"
    }:
        return "fibrous_cap_core"
    return "inflammation"


def parent_lineage(cell_state: str) -> str:
    return str(cell_state).split(" | ", 1)[0]


def load_tier1_roster() -> tuple[pd.DataFrame, str, int]:
    source = "embedded_frozen_public_roster_20260822"
    table = None
    for path in TIER1_FILE_CANDIDATES:
        if not path.exists():
            continue
        candidate = pd.read_csv(path, sep="\t")
        if "Tier" in candidate.columns:
            candidate = candidate[candidate["Tier"].astype(str).str.startswith("Tier1")].copy()
        required = {"gene_symbol_match", "annotation_level", "cell_state"}
        if required.issubset(candidate.columns):
            table = candidate.copy()
            source = str(path)
            break
    if table is None:
        table = pd.DataFrame(
            EMBEDDED_TIER1,
            columns=["gene_symbol_match", "annotation_level", "cell_state", "expected_direction"],
        )
        table["consensus_direction"] = np.where(
            table["expected_direction"] > 0, "up_in_disease", "down_in_disease"
        )
        table["Tier"] = "Tier1_submitted"
    else:
        table["gene_symbol_match"] = table["gene_symbol_match"].astype(str).str.upper()
        if "meta_logFC" in table.columns:
            table["expected_direction"] = np.sign(pd.to_numeric(table["meta_logFC"], errors="coerce"))
        elif "consensus_direction" in table.columns:
            table["expected_direction"] = table["consensus_direction"].map(
                {"up_in_disease": 1, "down_in_disease": -1}
            )
        else:
            raise ValueError("Tier 1 file lacks meta_logFC and consensus_direction")

    table = table.rename(columns={
        "gene_symbol_match": "gene",
        "annotation_level": "context_level",
        "cell_state": "cell_context",
    })
    table["gene"] = table["gene"].astype(str).str.upper()
    table["context_level"] = table["context_level"].astype(str)
    table["cell_context"] = table["cell_context"].astype(str)
    table["expected_direction"] = pd.to_numeric(table["expected_direction"], errors="coerce")
    submitted_pairs = {(gene, level, state) for gene, level, state, _ in EMBEDDED_TIER1}
    in_submitted_roster = [
        (row.gene, row.context_level, row.cell_context) in submitted_pairs
        for row in table.itertuples(index=False)
    ]
    n_nonroster_removed = int(len(table) - sum(in_submitted_roster))
    table = table.loc[in_submitted_roster].copy()
    table["parent_lineage"] = table["cell_context"].map(parent_lineage)
    table["pathology"] = table["parent_lineage"].map(primary_pathology)
    table["pair_id"] = table["gene"] + " | " + table["cell_context"]
    table = table.drop_duplicates(["gene", "context_level", "cell_context"]).reset_index(drop=True)

    if len(table) != 22 or table["gene"].nunique() != 21:
        raise ValueError(
            f"Frozen roster mismatch: expected 22 pairs/21 genes, observed "
            f"{len(table)} pairs/{table['gene'].nunique()} genes"
        )
    return table, source, n_nonroster_removed


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


def lineage_abundance(abundance: pd.DataFrame) -> pd.DataFrame:
    groups: dict[str, list[str]] = {}
    for state in abundance.columns:
        groups.setdefault(parent_lineage(state), []).append(state)
    return pd.DataFrame(
        {lineage: abundance[states].sum(axis=1) for lineage, states in groups.items()},
        index=abundance.index,
    )


def raw_matrix_and_library(adata: ad.AnnData):
    matrix = adata.layers["counts"] if "counts" in adata.layers else adata.X
    if "total_counts_original" in adata.obs:
        library = pd.to_numeric(adata.obs["total_counts_original"], errors="coerce").to_numpy(float)
    else:
        library = np.asarray(matrix.sum(axis=1)).ravel().astype(float)
    return matrix, library


def gene_position_map(adata: ad.AnnData) -> dict[str, int]:
    return {str(gene).upper(): i for i, gene in enumerate(adata.var_names)}


def normalized_gene_matrix(
    matrix, library: np.ndarray, positions: dict[str, int], genes: list[str], index: np.ndarray,
) -> pd.DataFrame:
    valid = [gene for gene in genes if gene in positions]
    if not valid:
        return pd.DataFrame(index=np.arange(len(index)))
    # scipy sparse matrices interpret matrix[row_array, col_array] as paired
    # element-wise indexing. Select the rectangular spot-by-gene block in two
    # steps so different numbers of spots and genes are always valid.
    row_index = np.asarray(index, dtype=int)
    column_index = np.asarray([positions[gene] for gene in valid], dtype=int)
    block = matrix[row_index, :][:, column_index]
    if sp.issparse(block):
        block = block.toarray()
    block = np.asarray(block, dtype=float)
    lib = np.maximum(library[index], 1)
    block = np.log1p(block * (1e4 / lib)[:, None])
    return pd.DataFrame(block, columns=valid)


def raw_gene_vector(matrix, index: np.ndarray, position: int) -> np.ndarray:
    row_index = np.asarray(index, dtype=int)
    values = matrix[row_index, :][:, int(position)]
    if sp.issparse(values):
        values = values.toarray()
    return np.asarray(values).ravel().astype(float)


def score_signatures(expression: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
    scores = pd.DataFrame(index=expression.index)
    coverage = []
    for name, genes in PATHOLOGY_SIGNATURES.items():
        valid = [gene for gene in genes if gene in expression.columns]
        coverage.append({
            "signature": name, "n_requested": len(genes), "n_valid": len(valid),
            "valid_genes": ";".join(valid),
        })
        if len(valid) < 3:
            scores[name] = np.nan
        else:
            scores[name] = np.nanmean(
                np.column_stack([zscore(expression[gene]) for gene in valid]), axis=1
            )
    return scores, coverage


def spatial_blocks(obs: pd.DataFrame) -> np.ndarray:
    if {"array_row", "array_col"}.issubset(obs.columns):
        row = pd.to_numeric(obs["array_row"], errors="coerce").fillna(0).astype(int)
        col = pd.to_numeric(obs["array_col"], errors="coerce").fillna(0).astype(int)
        return (row // 5).astype(str).str.cat((col // 5).astype(str), sep="_").to_numpy()
    return (np.arange(len(obs)) // 50).astype(str)


def fit_poisson(
    counts: np.ndarray,
    library: np.ndarray,
    abundance: np.ndarray,
    blocks: np.ndarray,
    pathology: np.ndarray | None = None,
) -> list[dict]:
    finite = np.isfinite(counts) & np.isfinite(library) & np.isfinite(abundance)
    if pathology is not None:
        finite &= np.isfinite(pathology)
    if finite.sum() < MIN_SPOTS or np.sum(counts[finite] > 0) < MIN_POSITIVE_SPOTS:
        return []
    if np.unique(blocks[finite]).size < 2:
        return []
    a = zscore(abundance[finite])
    if pathology is None:
        design = pd.DataFrame({"constant": 1.0, "cell_abundance": a})
        model_name = "localization"
        terms = ["cell_abundance"]
    else:
        p = zscore(pathology[finite])
        design = pd.DataFrame({
            "constant": 1.0,
            "cell_abundance": a,
            "pathology": p,
            "interaction": a * p,
        })
        model_name = "pathology_context"
        terms = ["cell_abundance", "pathology", "interaction"]
    try:
        fit = sm.GLM(
            counts[finite], design,
            family=sm.families.Poisson(link=sm.families.links.Log()),
            offset=np.log(np.maximum(library[finite], 1)),
        ).fit(cov_type="cluster", cov_kwds={"groups": blocks[finite]})
    except Exception as exc:
        warnings.warn(f"Poisson model failed: {exc!r}")
        return []
    output = []
    for term in terms:
        output.append({
            "model_scope": model_name,
            "term": term,
            "beta": float(fit.params[term]),
            "se": float(fit.bse[term]),
            "z": float(fit.tvalues[term]),
            "p": float(fit.pvalues[term]),
            "n_spots": int(finite.sum()),
            "n_positive_spots": int(np.sum(counts[finite] > 0)),
            "model": "Poisson_loglink_library_offset_spatial_block_cluster_robust",
        })
    return output


# =============================================================================
# 2. Load frozen hypotheses and input contracts
# =============================================================================

tier1, tier1_source, n_nonroster_removed = load_tier1_roster()
tier1.to_csv(AUDIT_OUT / "Tier1_roster_used.tsv", sep="\t", index=False)
pd.DataFrame([{
    "submitted_pairs_expected": 22,
    "submitted_genes_expected": 21,
    "n_nonroster_rows_removed_from_source": n_nonroster_removed,
    "submitted_pairs_observed": len(tier1),
    "submitted_genes_observed": tier1["gene"].nunique(),
    "status": "PASS",
}]).to_csv(AUDIT_OUT / "Tier1_roster_contract_audit.tsv", sep="\t", index=False)

if not LABEL_MAP_PATH.exists():
    raise FileNotFoundError(LABEL_MAP_PATH)
label_map = pd.read_csv(LABEL_MAP_PATH, sep="\t")
safe_to_original = dict(zip(
    label_map["spatial_reference_state_safe"].astype(str),
    label_map["spatial_reference_state_original"].astype(str),
))

manifest = pd.read_csv(MANIFEST_PATH, sep="\t") if MANIFEST_PATH.exists() else pd.DataFrame()
mapped_paths = sorted(MAP_OUT.glob("GSE*__final_atlas_cell2location.h5ad"))
if not mapped_paths:
    raise FileNotFoundError(f"No mapped h5ad files found under {MAP_OUT}")

print(f"Frozen Tier 1 roster: {len(tier1)} pairs, {tier1['gene'].nunique()} genes")
print(f"Roster source: {tier1_source}")
print(f"Mapped spatial datasets: {len(mapped_paths)}")
print("Reusing cell2location posteriors; deconvolution will NOT be rerun.")


# =============================================================================
# 3. Section-level extraction and models
# =============================================================================

model_rows: list[dict] = []
detection_rows: list[dict] = []
capture_rows: list[dict] = []
coverage_rows: list[dict] = []
context_rows: list[dict] = []

target_genes = sorted(tier1["gene"].unique())
signature_genes = sorted({gene for genes in PATHOLOGY_SIGNATURES.values() for gene in genes})
required_genes = sorted(set(target_genes) | set(signature_genes))

for mapped_path in mapped_paths:
    print("\nProcessing", mapped_path.name)
    adata = ad.read_h5ad(mapped_path)
    gse = str(adata.obs["GSE"].iloc[0]) if "GSE" in adata.obs else mapped_path.name.split("__", 1)[0]

    # Refresh curated section metadata without changing the fitted posterior.
    if not manifest.empty and "section_id" in adata.obs and "section_id" in manifest:
        by_section = manifest.drop_duplicates("section_id").set_index("section_id")
        section_ids = adata.obs["section_id"].astype(str)
        for column in ["biological_unit_id", "metadata_status", "condition", "tissue", "dataset_role"]:
            if column in by_section.columns:
                mapped = section_ids.map(by_section[column])
                if mapped.notna().any():
                    fallback = adata.obs[column].astype(str) if column in adata.obs else pd.Series("NA", index=adata.obs_names)
                    adata.obs[column] = mapped.fillna(fallback).to_numpy()

    matrix, library = raw_matrix_and_library(adata)
    positions = gene_position_map(adata)
    gene_presence = {gene: gene in positions for gene in target_genes}

    for posterior in POSTERIORS:
        abundance = abundance_frame(adata, posterior, safe_to_original)
        lineages = lineage_abundance(abundance)

        for _, hypothesis in tier1.iterrows():
            context_frame = abundance if hypothesis["context_level"] == "Level2" else lineages
            context_rows.append({
                "GSE": gse,
                "posterior": posterior,
                "pair_id": hypothesis["pair_id"],
                "gene": hypothesis["gene"],
                "context_level": hypothesis["context_level"],
                "cell_context": hypothesis["cell_context"],
                "gene_present": bool(gene_presence[hypothesis["gene"]]),
                "cell_context_present": bool(hypothesis["cell_context"] in context_frame.columns),
                "pathology": hypothesis["pathology"],
            })

        if "section_id" not in adata.obs:
            raise KeyError(f"section_id absent from {mapped_path}")
        for section_id, indices in adata.obs.groupby("section_id", observed=True).indices.items():
            idx = np.asarray(indices, dtype=int)
            obs = adata.obs.iloc[idx].copy()
            blocks = spatial_blocks(obs)
            section_expression = normalized_gene_matrix(matrix, library, positions, required_genes, idx)
            pathology_scores, coverage = score_signatures(section_expression)
            coverage_rows.extend([{
                "GSE": gse, "section_id": str(section_id), **row
            } for row in coverage])
            section_abundance = abundance.iloc[idx].reset_index(drop=True)
            section_lineages = lineages.iloc[idx].reset_index(drop=True)

            metadata = {
                "GSE": gse,
                "section_id": str(section_id),
                "sample_id": str(obs["sample_id"].iloc[0]) if "sample_id" in obs else str(section_id),
                "biological_unit_id": str(obs["biological_unit_id"].iloc[0]) if "biological_unit_id" in obs else str(section_id),
                "metadata_status": str(obs["metadata_status"].iloc[0]) if "metadata_status" in obs else "not_available",
                "condition": str(obs["condition"].iloc[0]) if "condition" in obs else "not_available",
                "tissue": str(obs["tissue"].iloc[0]) if "tissue" in obs else "not_available",
                "posterior": posterior,
            }

            for _, hypothesis in tier1.iterrows():
                gene = hypothesis["gene"]
                context_frame = section_abundance if hypothesis["context_level"] == "Level2" else section_lineages
                base = {
                    **metadata,
                    "pair_id": hypothesis["pair_id"],
                    "gene": gene,
                    "context_level": hypothesis["context_level"],
                    "cell_context": hypothesis["cell_context"],
                    "parent_lineage": hypothesis["parent_lineage"],
                    "pathology": hypothesis["pathology"],
                    "expected_direction": hypothesis["expected_direction"],
                }
                if gene not in positions or hypothesis["cell_context"] not in context_frame.columns:
                    continue
                counts = raw_gene_vector(matrix, idx, positions[gene])
                abundance_values = context_frame[hypothesis["cell_context"]].to_numpy(float)
                n_positive = int(np.sum(counts > 0))
                detection_rows.append({
                    **base,
                    "n_spots": len(idx),
                    "n_positive_spots": n_positive,
                    "detection_fraction": float(np.mean(counts > 0)),
                    "mean_logCP10K": float(section_expression[gene].mean()),
                    "median_context_abundance": float(np.nanmedian(abundance_values)),
                    "model_eligible": bool(len(idx) >= MIN_SPOTS and n_positive >= MIN_POSITIVE_SPOTS),
                })

                # Primary test: direct spatial co-localisation with nominated context.
                for result in fit_poisson(
                    counts, library[idx], abundance_values, blocks, pathology=None
                ):
                    model_rows.append({**base, **result})

                # Secondary test: context association conditional on the pre-specified
                # pathology signature, including an abundance-by-pathology interaction.
                signature = hypothesis["pathology"]
                if signature in pathology_scores:
                    for result in fit_poisson(
                        counts, library[idx], abundance_values, blocks,
                        pathology=pathology_scores[signature].to_numpy(float),
                    ):
                        model_rows.append({**base, **result})

            if EXPORT_SPOT_LEVEL and posterior == PRIMARY_POSTERIOR:
                export = pd.DataFrame({
                    "spot_id": obs.index.astype(str),
                    "GSE": gse,
                    "section_id": str(section_id),
                    "library_size": library[idx],
                })
                for gene in target_genes:
                    export[f"logCP10K__{gene}"] = (
                        section_expression[gene].to_numpy(float) if gene in section_expression else np.nan
                    )
                spot_dir = OUT / "04_optional_spot_level"
                spot_dir.mkdir(exist_ok=True)
                safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(section_id))
                export.to_csv(spot_dir / f"{safe}__Tier1_expression.tsv.gz", sep="\t", index=False)

    del adata

models = pd.DataFrame(model_rows)
detection = pd.DataFrame(detection_rows)
coverage = pd.DataFrame(coverage_rows).drop_duplicates()
context_audit = pd.DataFrame(context_rows).drop_duplicates()

models.to_csv(TABLE_OUT / "Tier1_spatial_models_by_section.tsv.gz", sep="\t", index=False, compression="gzip")
detection.to_csv(TABLE_OUT / "Tier1_detection_by_section.tsv", sep="\t", index=False)
coverage.to_csv(AUDIT_OUT / "pathology_signature_gene_coverage.tsv", sep="\t", index=False)
context_audit.to_csv(AUDIT_OUT / "gene_and_cell_context_capture_audit.tsv", sep="\t", index=False)

if context_audit.empty or not context_audit["gene_present"].all():
    missing = context_audit.loc[~context_audit["gene_present"], ["GSE", "gene"]].drop_duplicates()
    warnings.warn(f"Some Tier 1 genes were absent: {missing.to_dict('records')[:10]}")
if not context_audit.empty and not context_audit["cell_context_present"].all():
    missing = context_audit.loc[
        ~context_audit["cell_context_present"], ["GSE", "pair_id", "posterior"]
    ].drop_duplicates()
    warnings.warn(f"Some nominated contexts were absent: {missing.to_dict('records')[:10]}")


# =============================================================================
# 4. Capture summaries and cross-biological-unit meta-analysis
# =============================================================================

if detection.empty:
    raise RuntimeError("No Tier 1 target was captured in the mapped spatial datasets")

capture = (
    detection[detection["posterior"] == PRIMARY_POSTERIOR]
    .groupby(["pair_id", "gene", "context_level", "cell_context", "GSE"], observed=True)
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
capture.to_csv(TABLE_OUT / "Tier1_capture_summary_by_pair_dataset.tsv", sep="\t", index=False)

meta_rows = []
for scope, scoped in {
    "confirmed_or_unique_biological_units": models[
        ~models["metadata_status"].astype(str).str.contains("needs_", na=False)
    ],
    "all_sections_sensitivity": models,
}.items():
    group_columns = [
        "pair_id", "gene", "context_level", "cell_context", "parent_lineage",
        "pathology", "expected_direction", "posterior", "model_scope", "term",
    ]
    for keys, block in scoped.groupby(group_columns, observed=True):
        if scope == "confirmed_or_unique_biological_units":
            effects = []
            for _, unit in block.groupby("biological_unit_id", observed=True):
                effects.append(fixed_effect(unit["beta"], unit["se"]))
            beta, se = zip(*effects) if effects else ([], [])
            n_sections = block["section_id"].nunique()
            n_datasets = block["GSE"].nunique()
            n_biological_units = block["biological_unit_id"].nunique()
        else:
            beta, se = block["beta"], block["se"]
            n_sections = block["section_id"].nunique()
            n_datasets = block["GSE"].nunique()
            n_biological_units = block["biological_unit_id"].nunique()
        pooled = random_effects(beta, se)
        meta_rows.append({
            "scope": scope,
            **dict(zip(group_columns, keys)),
            "n_sections": n_sections,
            "n_datasets": n_datasets,
            "n_biological_units": n_biological_units,
            **pooled,
        })

meta = pd.DataFrame(meta_rows)
if meta.empty:
    raise RuntimeError("No eligible section-level models were available for meta-analysis")

# Multiplicity is controlled across the 22 frozen pairs separately for each
# inferential family (scope, posterior, model and term).
meta["FDR"] = meta.groupby(
    ["scope", "posterior", "model_scope", "term"], observed=True
)["p"].transform(bh_fdr)
# A nominated cell-state localisation signal is expected to be positive and is
# not expected to share the sign of a within-cell-state disease log fold change.
# Pseudobulk direction is compared only with pathology-context terms.
meta["positive_cell_state_localization"] = np.where(
    meta["model_scope"] == "localization", meta["beta"] > 0, np.nan
)
meta["direction_matches_pseudobulk"] = np.where(
    meta["model_scope"] == "pathology_context",
    np.sign(meta["beta"]) == np.sign(meta["expected_direction"]),
    np.nan,
)
meta.to_csv(TABLE_OUT / "Tier1_spatial_meta_analysis.tsv", sep="\t", index=False)


# =============================================================================
# 5. Posterior sensitivity and cross-dataset direction consistency
# =============================================================================

sensitivity_rows = []
merge_keys = [
    "GSE", "section_id", "pair_id", "gene", "context_level", "cell_context",
    "model_scope", "term",
]
wide = models.pivot_table(index=merge_keys, columns="posterior", values="beta", aggfunc="first").reset_index()
if {"q05", "mean"}.issubset(wide.columns):
    for keys, block in wide.groupby(
        ["pair_id", "gene", "context_level", "cell_context", "model_scope", "term"],
        observed=True,
    ):
        finite = np.isfinite(block["q05"]) & np.isfinite(block["mean"])
        rho, p_value = spearmanr(
            block.loc[finite, "q05"], block.loc[finite, "mean"]
        ) if finite.sum() >= 3 else (np.nan, np.nan)
        sensitivity_rows.append({
            **dict(zip(
                ["pair_id", "gene", "context_level", "cell_context", "model_scope", "term"], keys
            )),
            "n_sections": int(finite.sum()),
            "spearman_q05_vs_mean": rho,
            "spearman_p": p_value,
            "direction_agreement_fraction": (
                float(np.mean(np.sign(block.loc[finite, "q05"]) == np.sign(block.loc[finite, "mean"])))
                if finite.any() else np.nan
            ),
        })
sensitivity = pd.DataFrame(sensitivity_rows)
sensitivity.to_csv(TABLE_OUT / "Tier1_q05_vs_mean_sensitivity.tsv", sep="\t", index=False)

consistency_rows = []
primary_models = models[
    (models["posterior"] == PRIMARY_POSTERIOR)
    & (
        ((models["model_scope"] == "localization") & (models["term"] == "cell_abundance"))
        | ((models["model_scope"] == "pathology_context") & (models["term"] == "interaction"))
    )
]
for keys, block in primary_models.groupby(
    ["pair_id", "gene", "context_level", "cell_context", "model_scope", "term", "expected_direction"],
    observed=True,
):
    dataset_effect = block.groupby("GSE", observed=True)["beta"].median().dropna()
    pseudobulk_expected = float(keys[-1])
    expected_for_test = 1.0 if keys[4] == "localization" else pseudobulk_expected
    n_match = int(np.sum(np.sign(dataset_effect) == np.sign(expected_for_test)))
    p_value = (
        float(binomtest(n_match, len(dataset_effect), 0.5, alternative="greater").pvalue)
        if len(dataset_effect) else np.nan
    )
    consistency_rows.append({
        **dict(zip(
            ["pair_id", "gene", "context_level", "cell_context", "model_scope", "term", "expected_direction"],
            keys,
        )),
        "n_datasets": len(dataset_effect),
        "expected_direction_for_consistency_test": expected_for_test,
        "n_datasets_matching_expected_direction": n_match,
        "direction_consistency_fraction": n_match / len(dataset_effect) if len(dataset_effect) else np.nan,
        "one_sided_binomial_p": p_value,
        "dataset_effects": ";".join(f"{name}:{value:.5g}" for name, value in dataset_effect.items()),
    })
consistency = pd.DataFrame(consistency_rows)
if not consistency.empty:
    consistency["FDR"] = consistency.groupby(["model_scope", "term"], observed=True)[
        "one_sided_binomial_p"
    ].transform(bh_fdr)
consistency.to_csv(TABLE_OUT / "Tier1_cross_dataset_direction_consistency.tsv", sep="\t", index=False)


# =============================================================================
# 6. One-row-per-pair evidence table
# =============================================================================

primary_meta = meta[
    (meta["scope"] == "confirmed_or_unique_biological_units")
    & (meta["posterior"] == PRIMARY_POSTERIOR)
]
localization = primary_meta[
    (primary_meta["model_scope"] == "localization") & (primary_meta["term"] == "cell_abundance")
].copy()
localization = localization.rename(columns={
    "beta": "localization_beta", "se": "localization_se", "p": "localization_p",
    "FDR": "localization_FDR", "CI_low": "localization_CI_low",
    "CI_high": "localization_CI_high", "I2": "localization_I2",
    "k": "localization_k_units", "n_datasets": "localization_n_datasets",
    "positive_cell_state_localization": "localization_positive",
})
keep = [
    "pair_id", "gene", "context_level", "cell_context", "parent_lineage", "pathology",
    "expected_direction", "localization_beta", "localization_se", "localization_CI_low",
    "localization_CI_high", "localization_p", "localization_FDR", "localization_I2",
    "localization_k_units", "localization_n_datasets",
    "localization_positive",
]
localization = localization[[column for column in keep if column in localization.columns]]

interaction = primary_meta[
    (primary_meta["model_scope"] == "pathology_context") & (primary_meta["term"] == "interaction")
].copy()
interaction = interaction.rename(columns={
    "beta": "interaction_beta", "se": "interaction_se", "p": "interaction_p",
    "FDR": "interaction_FDR", "CI_low": "interaction_CI_low",
    "CI_high": "interaction_CI_high", "I2": "interaction_I2",
    "direction_matches_pseudobulk": "interaction_direction_matches_pseudobulk",
})
interaction_keep = [
    "pair_id", "interaction_beta", "interaction_se", "interaction_CI_low",
    "interaction_CI_high", "interaction_p", "interaction_FDR", "interaction_I2",
    "interaction_direction_matches_pseudobulk",
]
interaction = interaction[[column for column in interaction_keep if column in interaction.columns]]

capture_pair = (
    capture.groupby(["pair_id", "gene", "context_level", "cell_context"], observed=True)
    .agg(
        capture_n_datasets=("GSE", "nunique"),
        capture_n_sections=("n_sections", "sum"),
        capture_n_positive_spots=("n_positive_spots", "sum"),
        median_dataset_detection=("median_section_detection", "median"),
        n_model_eligible_sections=("n_model_eligible_sections", "sum"),
    )
    .reset_index()
)

pair_summary = tier1.merge(capture_pair, on=["pair_id", "gene", "context_level", "cell_context"], how="left")
pair_summary = pair_summary.merge(localization, on=[
    column for column in ["pair_id", "gene", "context_level", "cell_context", "parent_lineage", "pathology", "expected_direction"]
    if column in localization.columns
], how="left")
pair_summary = pair_summary.merge(interaction, on="pair_id", how="left")
pair_summary["spatial_support_class"] = np.select(
    [
        (pair_summary["localization_FDR"] < 0.05)
        & pair_summary["localization_positive"].fillna(False),
        (pair_summary["localization_p"] < 0.05)
        & pair_summary["localization_positive"].fillna(False),
        pair_summary["n_model_eligible_sections"].fillna(0).eq(0),
    ],
    ["FDR-supported_positive_localization", "nominal_positive_localization", "not_testable_low_capture"],
    default="not_positive_localization_supported",
)
pair_summary.to_csv(TABLE_OUT / "Tier1_pair_evidence_summary.tsv", sep="\t", index=False)


# =============================================================================
# 7. Journal-ready overview figures and source data
# =============================================================================

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "axes.linewidth": 0.7,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})


def save_figure(fig, stem: str):
    fig.savefig(FIGURE_OUT / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(FIGURE_OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIGURE_OUT / f"{stem}.tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.show()
    plt.close(fig)


# A. Capture heatmap: median section detection fraction by gene and dataset.
capture_gene = (
    capture.groupby(["gene", "GSE"], observed=True)["median_section_detection"]
    .median().unstack("GSE")
)
capture_gene = capture_gene.reindex(target_genes)
capture_gene.to_csv(SOURCE_DATA_OUT / "Tier1_capture_heatmap_source.tsv", sep="\t")
fig, ax = plt.subplots(figsize=(max(4.8, 0.65 * capture_gene.shape[1]), 6.5))
image = ax.imshow(capture_gene.to_numpy(float), aspect="auto", cmap="Blues", vmin=0, vmax=max(0.05, np.nanmax(capture_gene.to_numpy(float))))
ax.set_xticks(range(capture_gene.shape[1]), capture_gene.columns, rotation=45, ha="right")
ax.set_yticks(range(capture_gene.shape[0]), capture_gene.index)
ax.set_xlabel("Spatial transcriptomic dataset")
ax.set_ylabel("Tier 1 gene")
cbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
cbar.set_label("Median spot detection fraction")
save_figure(fig, "Tier1_capture_heatmap")


# B. Primary localization forest plot for all 22 gene-state pairs.
forest = pair_summary.copy()
forest["display"] = forest["gene"] + "  |  " + forest["cell_context"]
forest = forest.sort_values(["localization_beta", "gene"], na_position="first").reset_index(drop=True)
forest.to_csv(SOURCE_DATA_OUT / "Tier1_localization_forest_source.tsv", sep="\t", index=False)
fig, ax = plt.subplots(figsize=(7.2, 7.4))
y = np.arange(len(forest))
valid = np.isfinite(forest["localization_beta"])
colors = np.where(
    forest["localization_FDR"].fillna(1) < 0.05,
    COLORS["red"],
    np.where(forest["localization_p"].fillna(1) < 0.05, COLORS["orange"], COLORS["blue"]),
)
ax.axvline(0, color=COLORS["grey"], lw=0.8, ls="--", zorder=0)
ax.errorbar(
    forest.loc[valid, "localization_beta"], y[valid],
    xerr=[
        forest.loc[valid, "localization_beta"] - forest.loc[valid, "localization_CI_low"],
        forest.loc[valid, "localization_CI_high"] - forest.loc[valid, "localization_beta"],
    ],
    fmt="none", ecolor=COLORS["grey"], elinewidth=0.8, capsize=1.5, zorder=1,
)
ax.scatter(forest.loc[valid, "localization_beta"], y[valid], c=colors[valid], s=24, edgecolor="white", linewidth=0.4, zorder=2)
ax.set_yticks(y, forest["display"])
ax.set_xlabel("Spatial association with nominated cell-state abundance, beta (95% CI)")
ax.set_ylabel("")
ax.grid(axis="x", color="#D9E0E5", lw=0.5)
save_figure(fig, "Tier1_localization_forest")


# C. Compact evidence matrix for capture, localization and pathology interaction.
matrix_source = pair_summary[[
    "pair_id", "gene", "cell_context", "median_dataset_detection",
    "localization_beta", "localization_FDR", "interaction_beta", "interaction_FDR",
    "spatial_support_class",
]].copy()
matrix_source.to_csv(SOURCE_DATA_OUT / "Tier1_evidence_matrix_source.tsv", sep="\t", index=False)
plot_matrix = matrix_source.set_index("pair_id")[["localization_beta", "interaction_beta"]]
limit = np.nanquantile(np.abs(plot_matrix.to_numpy(float)), 0.95)
limit = max(float(limit) if np.isfinite(limit) else 1, 0.1)
fig, ax = plt.subplots(figsize=(7.2, 7.1))
image = ax.imshow(plot_matrix.to_numpy(float), aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit)
ax.set_xticks([0, 1], ["Cell-state localization", "Pathology interaction"])
ax.set_yticks(range(len(plot_matrix)), plot_matrix.index)
ax.set_ylabel("")
cbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
cbar.set_label("Meta-analytic beta")
save_figure(fig, "Tier1_spatial_evidence_matrix")


# =============================================================================
# 8. Freeze manifest
# =============================================================================

run_manifest = {
    "status": "completed",
    "generated_at": now_iso(),
    "tier1_source": tier1_source,
    "n_pairs": int(len(tier1)),
    "n_unique_genes": int(tier1["gene"].nunique()),
    "CKB_redacted": True,
    "mapped_h5ad_files": [str(path) for path in mapped_paths],
    "primary_posterior": PRIMARY_POSTERIOR,
    "minimum_spots": MIN_SPOTS,
    "minimum_positive_spots": MIN_POSITIVE_SPOTS,
    "primary_model": "gene counts ~ nominated cell-state abundance + offset(log library size)",
    "secondary_model": "gene counts ~ abundance + pathology + abundance:pathology + offset(log library size)",
    "variance_estimator": "spatial-block cluster-robust standard errors",
    "meta_analysis": "fixed effect within biological unit, DerSimonian-Laird random effects across units",
    "multiple_testing": "Benjamini-Hochberg across 22 frozen pairs within each scope/posterior/model/term family",
    "interpretation_guardrail": (
        "Visium associations indicate regional co-localisation and do not establish "
        "cell-intrinsic differential expression or causal regulation."
    ),
    "software": {
        "python": sys.version,
        "platform": platform.platform(),
        "anndata": ad.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "statsmodels": getattr(sm, "__version__", "unknown"),
    },
}
(OUT / "run_manifest.json").write_text(json.dumps(run_manifest, indent=2, ensure_ascii=False), encoding="utf-8")

print("\nALL PUBLIC TIER 1 SPATIAL VALIDATION COMPLETE")
print("Output:", OUT)
print("Pairs:", len(tier1), "Unique genes:", tier1["gene"].nunique())
print("Section-level model rows:", len(models))
print("Meta-analysis rows:", len(meta))
print("Pair-level summary:", TABLE_OUT / "Tier1_pair_evidence_summary.tsv")
