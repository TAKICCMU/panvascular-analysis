from __future__ import annotations

import gc
import gzip
import json
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite


# =============================================================================
# Configuration
# =============================================================================

ATLAS_CANDIDATES = [
    Path("/home/anzhen_vhd/scRNAseq/codex/Level2_celltype_final_v1_metadata_v2.h5ad"),
    Path("/home/anzhen_vhd/scRNAseq/codex/Level2_celltype_final_v1.h5ad"),
]
ATLAS_H5AD = next((p for p in ATLAS_CANDIDATES if p.exists()), ATLAS_CANDIDATES[-1])
OUT_ROOT = Path("/home/anzhen_vhd/scRNAseq/codex/pseudobulk_target_discovery_v1")
DISEASE_MAP = OUT_ROOT / "sample_disease_map.tsv"

LEVEL1_KEY = "Level1_celltype_final_v3"
LEVEL2_KEY = "Level2_celltype_for_enrichment_v1"
LEVEL2_ELIGIBLE = "Level2_primary_eligible_v1"
BIO_SAMPLE_KEY = "biological_sample_uid"
ANALYSIS_UNIT_KEY = "pseudobulk_unit_uid"
SAMPLE_KEY = BIO_SAMPLE_KEY
MIN_CELLS_PER_SAMPLE_STATE = 20

# User-requested exclusions. Unknown is also excluded because it is not a lineage.
LEVEL1_EXCLUDE = {"Schwann cell", "Epithelial cell", "Unknown"}

OUT_ROOT.mkdir(parents=True, exist_ok=True)


# =============================================================================
# Load the immutable final atlas and prepare a disease-context mapping template
# =============================================================================

metadata_v2_available = ATLAS_CANDIDATES[0].exists()
if "adata" in globals() and all(
    k in adata.obs for k in [LEVEL1_KEY, LEVEL2_KEY, LEVEL2_ELIGIBLE]
) and (not metadata_v2_available or "condition_metadata_source_v2" in adata.obs):
    print("Reusing compatible adata already loaded in this kernel")
else:
    print("Reading:", ATLAS_H5AD)
    adata = ad.read_h5ad(ATLAS_H5AD)

if adata.raw is None:
    raise ValueError("Full-gene raw counts in adata.raw.X are required")

required = {
    SAMPLE_KEY, "study_canonical", "condition_v2", "tissue_v2",
    "sex_v2", "age_v2", LEVEL1_KEY, LEVEL2_KEY, LEVEL2_ELIGIBLE,
}
missing = required - set(adata.obs)
if missing:
    raise KeyError(f"Missing metadata: {sorted(missing)}")

meta_cols = [
    SAMPLE_KEY, "sample_uid", "donor_uid_v2", "study_canonical",
    "condition_v2", "tissue_v2", "sex_v2", "age_v2",
]
cell_meta = adata.obs[[c for c in meta_cols if c in adata.obs]].copy()
for c in cell_meta:
    cell_meta[c] = cell_meta[c].astype(str).replace({"nan": np.nan, "None": np.nan})


def one_value_per_sample(frame: pd.DataFrame, sample: str) -> pd.DataFrame:
    # One biological sample can legitimately contain multiple technical
    # sample_uid values (for example, two libraries/runs from the same tissue).
    # These identifiers are retained as an audit list, whereas biological
    # attributes must remain unique within biological_sample_uid.
    multivalue_identifier_columns = {"sample_uid"}
    rows = []
    for sid, sub in frame.groupby(sample, observed=True, sort=False):
        row = {sample: sid}
        for c in frame.columns:
            if c == sample:
                continue
            values = pd.Series(sub[c]).dropna().astype(str).unique()
            if c in multivalue_identifier_columns:
                row["technical_sample_uids"] = ";".join(sorted(values))
                row["n_technical_sample_uids"] = len(values)
                continue
            if len(values) > 1:
                raise ValueError(
                    f"{sample}={sid} has multiple values for {c}: {values[:10].tolist()}"
                )
            row[c] = values[0] if len(values) else np.nan
        row["n_atlas_cells"] = len(sub)
        rows.append(row)
    return pd.DataFrame(rows)


sample_meta = one_value_per_sample(cell_meta, SAMPLE_KEY)
template = sample_meta.copy()
template["condition_for_pseudobulk"] = template["condition_v2"]
template["disease_context"] = ""
template["include_primary"] = True
template["mapping_note"] = ""
template_path = OUT_ROOT / "sample_disease_map_template.tsv"
template.to_csv(template_path, sep="\t", index=False)

if not DISEASE_MAP.exists():
    raise FileNotFoundError(
        "Disease-specific pseudobulk requires a completed sample mapping.\n"
        f"A template has been written to:\n{template_path}\n\n"
        "Complete disease_context for BOTH disease and matched healthy samples "
        "(for example CAD, AAA, CAVS, TAAD), save it as:\n"
        f"{DISEASE_MAP}\n"
        "Then rerun this same cell. Do not infer disease_context from tissue alone."
    )

disease_map = pd.read_csv(DISEASE_MAP, sep="\t", dtype=str)
map_required = {SAMPLE_KEY, "disease_context", "condition_for_pseudobulk"}
if not map_required.issubset(disease_map):
    raise KeyError(f"{DISEASE_MAP} must contain {sorted(map_required)}")
if disease_map[SAMPLE_KEY].duplicated().any():
    raise ValueError("sample_disease_map.tsv has duplicated biological_sample_uid")

disease_map["disease_context"] = disease_map["disease_context"].astype(str).str.strip()
disease_map["condition_for_pseudobulk"] = (
    disease_map["condition_for_pseudobulk"].astype(str).str.strip().str.lower()
)
invalid_condition = ~disease_map["condition_for_pseudobulk"].isin(["healthy", "disease"])
if invalid_condition.any():
    raise ValueError(
        "Invalid condition_for_pseudobulk values: "
        + "; ".join(sorted(disease_map.loc[invalid_condition, "condition_for_pseudobulk"].unique()))
    )
disease_map = disease_map.loc[
    disease_map["disease_context"].notna()
    & ~disease_map["disease_context"].isin(["", "nan", "NA"])
].copy()
if "include_primary" in disease_map:
    disease_map = disease_map.loc[
        disease_map["include_primary"].astype(str).str.lower().isin(["true", "1", "yes", "y"])
    ]

sample_meta = sample_meta.merge(
    disease_map[[SAMPLE_KEY, "disease_context", "condition_for_pseudobulk"]],
    on=SAMPLE_KEY, how="left",
)
sample_meta["condition_v2_original"] = sample_meta["condition_v2"]
sample_meta["condition_v2"] = sample_meta["condition_for_pseudobulk"].where(
    sample_meta["condition_for_pseudobulk"].notna(), sample_meta["condition_v2"]
)

# The inferential unit is donor x tissue x condition. Including condition keeps
# legitimate paired specimens separate (for example, diseased and non-diseased
# carotid tissue from the same GSE159677 donor), while still collapsing multiple
# technical libraries or repeated pieces from the same donor/tissue/condition.
sample_meta[ANALYSIS_UNIT_KEY] = (
    sample_meta["donor_uid_v2"].astype(str)
    + "::TISSUE::"
    + sample_meta["tissue_v2"].astype(str)
    + "::CONDITION::"
    + sample_meta["condition_v2"].astype(str)
)


def collapse_analysis_units(frame: pd.DataFrame) -> pd.DataFrame:
    invariant = [
        "study_canonical", "condition_v2", "disease_context", "tissue_v2",
        "sex_v2", "age_v2", "donor_uid_v2",
    ]
    rows = []
    for uid, sub in frame.groupby(ANALYSIS_UNIT_KEY, observed=True, sort=False):
        row = {ANALYSIS_UNIT_KEY: uid}
        for c in invariant:
            values = sub[c].dropna().astype(str).unique()
            if len(values) > 1:
                raise ValueError(f"{uid} has multiple {c} values: {values[:10].tolist()}")
            row[c] = values[0] if len(values) else np.nan
        row["biological_sample_uids"] = ";".join(sorted(sub[SAMPLE_KEY].astype(str).unique()))
        row["n_biological_samples_merged"] = sub[SAMPLE_KEY].nunique()
        row["technical_sample_uids"] = ";".join(sorted({
            x for value in sub["technical_sample_uids"].dropna().astype(str)
            for x in value.split(";") if x
        }))
        row["n_atlas_cells"] = int(sub["n_atlas_cells"].sum())
        rows.append(row)
    return pd.DataFrame(rows)


analysis_meta = collapse_analysis_units(
    sample_meta.loc[sample_meta["disease_context"].notna()].copy()
)
sample_meta.to_csv(OUT_ROOT / "sample_metadata_with_disease_context.tsv", sep="\t", index=False)
analysis_meta.to_csv(OUT_ROOT / "pseudobulk_analysis_unit_metadata.tsv", sep="\t", index=False)

bio_to_unit = (
    sample_meta.loc[sample_meta["disease_context"].notna()]
    .set_index(SAMPLE_KEY)[ANALYSIS_UNIT_KEY]
    .to_dict()
)
mapped_samples = set(sample_meta.loc[sample_meta["disease_context"].notna(), SAMPLE_KEY])
print("Mapped biological samples:", len(mapped_samples), "/", len(sample_meta))
print("Independent donor x tissue units:", len(analysis_meta))
print(analysis_meta["disease_context"].value_counts(dropna=False))


# =============================================================================
# Gene identifiers and duplicate-feature projection
# =============================================================================


def ensembl_ids(raw_var: pd.DataFrame, raw_names: pd.Index) -> pd.Series:
    for c in ["gene_ids", "gene_id_clean", "gene_id", "ensembl_id"]:
        if c in raw_var:
            s = raw_var[c].astype("string").str.replace(r"\..*$", "", regex=True)
            if s.str.upper().str.startswith("ENSG").mean() > 0.5:
                return s
    s = pd.Series(raw_names.astype(str), dtype="string").str.replace(
        r"\..*$", "", regex=True
    )
    if s.str.upper().str.startswith("ENSG").mean() <= 0.5:
        raise ValueError("No Ensembl IDs found in adata.raw.var")
    return s


def gene_symbols(raw_var: pd.DataFrame, raw_names: pd.Index) -> pd.Series:
    for c in ["gene_symbol_std", "gene_symbol", "gene_symbols", "symbol", "gene_name"]:
        if c in raw_var:
            s = raw_var[c].astype("string").str.strip().str.upper()
            if (s.notna() & (s != "")).mean() > 0.5:
                return s
    names = pd.Series(raw_names.astype(str), dtype="string")
    return names.str.upper() if names.str.startswith("ENSG").mean() < 0.5 else pd.Series(
        [pd.NA] * len(names), dtype="string"
    )


raw_var = adata.raw.var.copy()
raw_names = pd.Index(adata.raw.var_names)
ensg = ensembl_ids(raw_var, raw_names)
symbols = gene_symbols(raw_var, raw_names)
valid_gene = ensg.notna() & ensg.str.upper().str.startswith("ENSG")
raw_idx = np.flatnonzero(valid_gene.to_numpy())
ensg_valid = ensg.iloc[raw_idx].astype(str).to_numpy()
symbol_valid = symbols.iloc[raw_idx].astype("string").to_numpy()
unique_ensg = pd.Index(pd.unique(ensg_valid), name="gene_id")
codes = pd.Categorical(ensg_valid, categories=unique_ensg).codes
projection = sp.csr_matrix(
    (np.ones(len(codes), np.int8), (np.arange(len(codes)), codes)),
    shape=(len(codes), len(unique_ensg)),
)
symbol_map = (
    pd.DataFrame({"gene_id": ensg_valid, "gene_symbol": symbol_valid})
    .dropna(subset=["gene_symbol"])
    .query("gene_symbol != ''")
    .drop_duplicates("gene_id")
    .set_index("gene_id")["gene_symbol"]
)
features = pd.DataFrame({"gene_id": unique_ensg})
features["gene_symbol"] = features["gene_id"].map(symbol_map)


# =============================================================================
# Aggregate raw counts by biological sample x cell state
# =============================================================================


def write_matrix_market_gz(matrix, path: Path):
    tmp = Path(str(path).removesuffix(".gz"))
    mmwrite(tmp, matrix)
    with tmp.open("rb") as src, gzip.open(path, "wb", compresslevel=6) as dst:
        shutil.copyfileobj(src, dst)
    tmp.unlink()


def export_level(level: str, state_series: pd.Series, eligible: np.ndarray):
    out = OUT_ROOT / level
    out.mkdir(exist_ok=True)
    count_file = out / "counts_genes_by_pseudobulk.mtx.gz"
    unit_file = out / "pseudobulk_unit_metadata.tsv"
    feature_file = out / "features.tsv"
    if count_file.exists() and unit_file.exists() and feature_file.exists():
        print("Reusing completed export:", level)
        return

    state = state_series.astype(str).to_numpy()
    biological_samples = adata.obs[SAMPLE_KEY].astype(str).to_numpy()
    samples = pd.Series(biological_samples).map(bio_to_unit).to_numpy(object)
    mapped = pd.notna(samples)
    use = eligible & mapped
    pair_counts = (
        pd.DataFrame({"state": state[use], "sample": samples[use]})
        .groupby(["state", "sample"], observed=True)
        .size().rename("n_cells").reset_index()
    )
    pair_counts = pair_counts.loc[
        pair_counts["n_cells"] >= MIN_CELLS_PER_SAMPLE_STATE
    ].copy()
    valid_pairs = set(zip(pair_counts["state"], pair_counts["sample"]))

    blocks, unit_rows = [], []
    Xraw = adata.raw.X
    for i, state_name in enumerate(sorted(pd.unique(state[use]))):
        state_idx = np.flatnonzero(use & (state == state_name))
        if len(state_idx) == 0:
            continue
        state_samples = samples[state_idx]
        keep_sample_names = sorted({
            sid for sid in pd.unique(state_samples)
            if (state_name, sid) in valid_pairs
        })
        if len(keep_sample_names) == 0:
            continue
        sample_to_code = {x: j for j, x in enumerate(keep_sample_names)}
        keep_cell = np.array([x in sample_to_code for x in state_samples])
        state_idx = state_idx[keep_cell]
        state_samples = state_samples[keep_cell]
        sample_codes = np.array([sample_to_code[x] for x in state_samples], dtype=int)

        design = sp.csr_matrix(
            (np.ones(len(sample_codes), np.int8), (np.arange(len(sample_codes)), sample_codes)),
            shape=(len(sample_codes), len(keep_sample_names)),
        )
        block = Xraw[state_idx]
        block = block.tocsr() if sp.issparse(block) else sp.csr_matrix(block)
        agg = (design.T @ block[:, raw_idx] @ projection).tocsr()
        agg.eliminate_zeros()
        blocks.append(agg)

        pair_lookup = pair_counts.set_index(["state", "sample"])["n_cells"]
        sm = analysis_meta.set_index(ANALYSIS_UNIT_KEY)
        for sid in keep_sample_names:
            row = sm.loc[sid].to_dict()
            row.update({
                "pseudobulk_id": f"{level}::{state_name}::{sid}",
                "annotation_level": level,
                "cell_state": state_name,
                ANALYSIS_UNIT_KEY: sid,
                "n_cells": int(pair_lookup.loc[(state_name, sid)]),
            })
            unit_rows.append(row)
        print(level, f"{i + 1}: {state_name} -> {len(keep_sample_names)} units")

    if not blocks:
        raise RuntimeError(f"No pseudobulk units passed filters for {level}")
    units_by_genes = sp.vstack(blocks, format="csr")
    units = pd.DataFrame(unit_rows)
    if units_by_genes.shape[0] != len(units):
        raise RuntimeError("Pseudobulk matrix/metadata row mismatch")
    if units["pseudobulk_id"].duplicated().any():
        raise RuntimeError("Duplicated pseudobulk IDs")

    features.to_csv(feature_file, sep="\t", index=False)
    units.to_csv(unit_file, sep="\t", index=False)
    write_matrix_market_gz(units_by_genes.T.tocoo(), count_file)
    audit = {
        "annotation_level": level,
        "matrix_orientation": "genes_by_pseudobulk_units",
        "n_genes": int(units_by_genes.shape[1]),
        "n_units": int(units_by_genes.shape[0]),
        "n_cell_states": int(units["cell_state"].nunique()),
        "minimum_cells_per_sample_state": MIN_CELLS_PER_SAMPLE_STATE,
        "raw_count_source": "adata.raw.X",
    }
    (out / "export_manifest.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    del units_by_genes, blocks
    gc.collect()


l1_series = adata.obs[LEVEL1_KEY].astype(str)
l1_eligible = ~l1_series.isin(LEVEL1_EXCLUDE).to_numpy()
export_level("Level1", l1_series, l1_eligible)

l2_series = adata.obs[LEVEL2_KEY].astype(str)
l2_eligible = adata.obs[LEVEL2_ELIGIBLE].astype(bool).to_numpy()
export_level("Level2", l2_series, l2_eligible)

print("\nPseudobulk export completed:", OUT_ROOT)
print("Next run 05b_pseudobulk_edgeR_meta_analysis.R on the server.")
