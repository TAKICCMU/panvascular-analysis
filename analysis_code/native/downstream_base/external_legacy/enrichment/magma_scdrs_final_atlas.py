from __future__ import annotations

import gc
import hashlib
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr, wilcoxon
from pandas.api.types import is_categorical_dtype, is_object_dtype, is_string_dtype


# =============================================================================
# 0. Configuration
# =============================================================================

ATLAS_H5AD = Path("/path/to/panvas_project/Level2_celltype_final_v1.h5ad")
GWAS_ROOT = Path("/path/to/magma_scdrs")
OUT_ROOT = GWAS_ROOT / "final_atlas_MAGMA_scDRS_v1"
COMMON = OUT_ROOT / "00_common_input"

LEVEL1_KEY = "Level1_celltype_final_v3"
LEVEL2_KEY = "Level2_celltype_for_enrichment_v1"
LEVEL2_ELIGIBLE_KEY = "Level2_primary_eligible_v1"
L1_TEST_KEY = "Level1_scDRS_group_v1"
L2_TEST_KEY = "Level2_scDRS_group_v1"

SCDRS_ENV = "scdrs_env"
CONDA = Path("/path/to/miniconda3/bin/conda")
MAGMA = Path("/path/to/scrnaseq/magma")
N_CTRL = 1000
CHUNK_SIZE = 20_000
MIN_GROUP_CELLS = 500
MIN_GROUP_STUDIES = 3
MIN_GROUP_STUDY_CELLS = 20
MIN_SAMPLE_STATE_CELLS = 20
RANDOM_SEED = 20260806

# r0.1 is the primary LD-aware Q_SNP-locus-filtered GWAS; r0.6 changes the
# LD-proxy threshold as a sensitivity analysis. The old
# primary_bifactor_G_LDr20.1Q_20260731 directory is an analysis-output folder,
# not a fourth independent set of FUMA gene-level inputs.
SCENARIOS = {
    "primary_r0.1": GWAS_ROOT / "r0.1",
    "primary_r0.6": GWAS_ROOT / "r0.6",
    "leave_one_CAD_r0.1": GWAS_ROOT / "LOTr0.1",
}

for d in [OUT_ROOT, COMMON]:
    d.mkdir(parents=True, exist_ok=True)

if not CONDA.exists():
    raise FileNotFoundError(CONDA)
if not MAGMA.exists():
    magma_in_path = shutil.which("magma")
    if magma_in_path is None:
        raise FileNotFoundError("MAGMA executable was not found")
    MAGMA = Path(magma_in_path)


def clean_token(x: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", str(x)).strip("_")


def bh_fdr(values) -> np.ndarray:
    p = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(float)
    out = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    pv = p[valid]
    if len(pv) == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    restored = np.empty_like(q)
    restored[order] = q
    out[valid] = restored
    return out


def file_md5(path: Path, block=8 * 1024 * 1024) -> str:
    h = hashlib.md5()
    with path.open("rb") as handle:
        while True:
            part = handle.read(block)
            if not part:
                break
            h.update(part)
    return h.hexdigest()


def run_command(args, log_file: Path, cwd: Path | None = None):
    args = [str(x) for x in args]
    print("$", " ".join(args))
    env = os.environ.copy()
    # Prevent the Jupyter inline backend from leaking into the older scDRS env.
    env["MPLBACKEND"] = "Agg"
    completed = subprocess.run(
        args,
        cwd=str(cwd) if cwd else None,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    log_file.parent.mkdir(parents=True, exist_ok=True)
    log_file.write_text(completed.stdout, encoding="utf-8")
    print(completed.stdout[-4000:])
    if completed.returncode != 0:
        raise subprocess.CalledProcessError(completed.returncode, args)


def scdrs_command(*args):
    return [
        CONDA, "run", "--no-capture-output", "-n", SCDRS_ENV, "scdrs", *args
    ]


def resolve_one(root: Path, filename: str) -> Path:
    if not root.exists():
        raise FileNotFoundError(root)
    hits = [
        p for p in root.rglob(filename)
        if OUT_ROOT not in p.parents and not any(x.startswith("._") for x in p.parts)
    ]
    if len(hits) == 0:
        raise FileNotFoundError(f"No {filename} below {root}")
    if len(hits) > 1:
        preferred = [p for p in hits if "FUMA" in str(p) or p.parent == root]
        if len(preferred) == 1:
            return preferred[0]
        raise RuntimeError(
            f"Ambiguous {filename} below {root}:\n" + "\n".join(map(str, hits))
        )
    return hits[0]


# =============================================================================
# 1. Load and validate the locked final atlas
# =============================================================================

if "adata" in globals() and all(
    k in adata.obs for k in [LEVEL1_KEY, LEVEL2_KEY, LEVEL2_ELIGIBLE_KEY]
):
    print("Reusing compatible adata already loaded in this Jupyter kernel")
else:
    print("Reading:", ATLAS_H5AD)
    adata = ad.read_h5ad(ATLAS_H5AD)

if adata.raw is None:
    raise ValueError("adata.raw is absent; full-gene raw counts are required")
if adata.raw.n_obs != adata.n_obs:
    raise ValueError("adata.raw and adata have different cell counts")

required_obs = {
    LEVEL1_KEY, LEVEL2_KEY, LEVEL2_ELIGIBLE_KEY,
    "study_canonical", "biological_sample_uid", "condition_v2",
}
missing = required_obs - set(adata.obs)
if missing:
    raise KeyError(f"Missing final-atlas metadata: {sorted(missing)}")

obs = adata.obs
l1_counts = obs[LEVEL1_KEY].astype(str).value_counts()
l1_valid = set(l1_counts[l1_counts >= MIN_GROUP_CELLS].index) - {"Unknown"}
l2_eligible = obs[LEVEL2_ELIGIBLE_KEY].astype(bool)
l2_counts = obs.loc[l2_eligible, LEVEL2_KEY].astype(str).value_counts()
l2_valid = set(l2_counts[l2_counts >= MIN_GROUP_CELLS].index)

l1_test = obs[LEVEL1_KEY].astype(str).where(
    obs[LEVEL1_KEY].astype(str).isin(l1_valid), "Excluded_low_confidence_or_rare"
)
l2_test = obs[LEVEL2_KEY].astype(str).where(
    l2_eligible & obs[LEVEL2_KEY].astype(str).isin(l2_valid),
    "Excluded_nonstable_or_Level1_only",
)

print("Atlas:", adata.shape, "raw:", adata.raw.shape)
print("Level1 groups tested:", len(l1_valid), sorted(l1_valid))
print("Stable Level2 groups tested:", len(l2_valid))


# =============================================================================
# 2. Build one reusable full-gene raw-count scDRS input
# =============================================================================

SCDRS_H5AD = COMMON / "final_atlas_full_gene_raw_counts_scDRS.h5ad"
SCDRS_COV = COMMON / "final_atlas_scDRS.cov"
SCDRS_MANIFEST = COMMON / "scDRS_input_manifest.json"


def choose_symbols(raw_var: pd.DataFrame, raw_names: pd.Index):
    names = pd.Series(raw_names.astype(str), dtype="string")
    if names.str.upper().str.startswith("ENSG").mean() < 0.5:
        return names.str.strip().str.upper(), "raw.var_names"
    candidates = [
        "gene_symbol_std", "gene_symbol", "gene_symbols", "symbol",
        "gene_name", "feature_name",
    ]
    scored = []
    for c in candidates:
        if c in raw_var:
            s = raw_var[c].astype("string").str.strip()
            scored.append((float((s.notna() & (s != "")).mean()), c))
    if not scored:
        raise ValueError("Raw var names are Ensembl IDs but no gene-symbol column exists")
    _, source = max(scored)
    return raw_var[source].astype("string").str.strip().str.upper(), source


if not (SCDRS_H5AD.exists() and SCDRS_COV.exists() and SCDRS_MANIFEST.exists()):
    # A failed H5AD write can leave a partial file. It is safe to replace only
    # this versioned, pipeline-generated output when its manifest is absent.
    if SCDRS_H5AD.exists() and not SCDRS_MANIFEST.exists():
        print("Removing incomplete pipeline-generated H5AD:", SCDRS_H5AD)
        SCDRS_H5AD.unlink()
    raw_var = adata.raw.var.copy()
    symbols, symbol_source = choose_symbols(raw_var, pd.Index(adata.raw.var_names))
    valid = symbols.notna() & ~symbols.isin(["", "NA", "NAN", "NONE"])
    valid_idx = np.flatnonzero(valid.to_numpy())
    symbol_array = symbols.iloc[valid_idx].astype(str).to_numpy()
    unique_symbols = pd.Index(pd.unique(symbol_array), name="gene_symbol")

    X = adata.raw.X
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
    X = X[:, valid_idx]
    if len(unique_symbols) != len(symbol_array):
        codes = pd.Categorical(symbol_array, categories=unique_symbols).codes
        projection = sp.csr_matrix(
            (np.ones(len(codes), np.float32), (np.arange(len(codes)), codes)),
            shape=(len(codes), len(unique_symbols)),
        )
        X = (X @ projection).tocsr()

    keep_obs = [
        LEVEL1_KEY, LEVEL2_KEY, LEVEL2_ELIGIBLE_KEY, "study_canonical",
        "sample_uid", "biological_sample_uid", "condition_v2", "tissue_v2",
        "sex_v2", "age_v2", "total_counts", "n_genes_by_counts", "pct_counts_mt",
    ]
    sc_obs = adata.obs[[c for c in keep_obs if c in adata.obs]].copy()
    sc_obs[L1_TEST_KEY] = l1_test.astype(str).to_numpy(dtype=object)
    sc_obs[L2_TEST_KEY] = l2_test.astype(str).to_numpy(dtype=object)
    for c in sc_obs:
        dtype = sc_obs[c].dtype
        if is_categorical_dtype(dtype) or is_string_dtype(dtype) or is_object_dtype(dtype):
            values = sc_obs[c].astype(object).to_numpy()
            values[pd.isna(values)] = "NA"
            # Force a legacy NumPy object array. This avoids pandas nullable
            # StringArray categories, which older anndata readers cannot write.
            sc_obs[c] = pd.Series(values, index=sc_obs.index, dtype=object)

    n_counts = np.asarray(X.sum(axis=1)).ravel()
    n_genes = np.asarray(X.getnnz(axis=1)).ravel()
    sc_obs["scdrs_log1p_counts"] = np.log1p(n_counts)
    sc_obs["scdrs_log1p_genes"] = np.log1p(n_genes)

    sc_adata = ad.AnnData(
        X=X,
        obs=sc_obs,
        var=pd.DataFrame(index=unique_symbols),
    )
    sc_adata.uns["raw_count_source"] = "locked final atlas adata.raw.X"
    sc_adata.uns["annotation_version"] = "Level1_final_v3; Level2_final_v1"
    print("Writing reusable scDRS H5AD:", SCDRS_H5AD, sc_adata.shape)
    sc_adata.write_h5ad(
        SCDRS_H5AD,
        compression="lzf",
        convert_strings_to_categoricals=False,
    )

    cov = pd.DataFrame(index=sc_adata.obs_names)
    cov.index.name = "index"
    cov["const"] = 1.0
    for name, x in {
        "log1p_counts": np.log1p(n_counts),
        "log1p_genes": np.log1p(n_genes),
        "pct_mt": pd.to_numeric(
            sc_obs.get("pct_counts_mt", pd.Series(0, index=sc_obs.index)),
            errors="coerce",
        ).fillna(0).to_numpy(),
    }.items():
        sd = np.std(x)
        cov[name] = (x - np.mean(x)) / sd if sd > 0 else 0.0
    cov.to_csv(SCDRS_COV, sep="\t")
    SCDRS_MANIFEST.write_text(json.dumps({
        "atlas": str(ATLAS_H5AD),
        "atlas_n_obs": int(adata.n_obs),
        "raw_n_vars": int(adata.raw.n_vars),
        "scdrs_n_vars": int(sc_adata.n_vars),
        "gene_symbol_source": symbol_source,
        "level1_groups": sorted(l1_valid),
        "level2_groups": sorted(l2_valid),
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    del sc_adata, X, cov
    gc.collect()
else:
    print("Reusing:", SCDRS_H5AD)


# =============================================================================
# 3. Regenerate final-atlas MAGMA specificity matrices
# =============================================================================

SPEC_ROOT = COMMON / "MAGMA_specificity"
SPEC_ROOT.mkdir(exist_ok=True)


def stable_ensembl_ids(raw_var: pd.DataFrame, raw_names: pd.Index) -> pd.Series:
    candidates = ["gene_ids", "gene_id_clean", "gene_id", "ensembl_id"]
    for c in candidates:
        if c in raw_var:
            s = raw_var[c].astype("string").str.replace(r"\..*$", "", regex=True)
            if s.str.upper().str.startswith("ENSG").mean() > 0.5:
                return s
    names = pd.Series(raw_names.astype(str), dtype="string").str.replace(
        r"\..*$", "", regex=True
    )
    if names.str.upper().str.startswith("ENSG").mean() <= 0.5:
        raise ValueError("No Ensembl gene IDs were found for MAGMA specificity")
    return names


def make_group_mapping(series: pd.Series, eligible: np.ndarray, prefix: str):
    tab = pd.DataFrame({"group": series.astype(str), "eligible": eligible})
    audit = (
        pd.DataFrame({
            "group": series.astype(str),
            "study": adata.obs["study_canonical"].astype(str),
            "eligible": eligible,
        })
        .query("eligible")
        .groupby("group", observed=True)
        .agg(n_cells=("group", "size"), n_studies=("study", "nunique"))
        .reset_index()
    )
    audit = audit.query(
        "n_cells >= @MIN_GROUP_CELLS and n_studies >= @MIN_GROUP_STUDIES"
    ).sort_values("group")
    audit["variable"] = [f"{prefix}{i:03d}" for i in range(1, len(audit) + 1)]
    return audit[["variable", "group", "n_cells", "n_studies"]]


L1_MAP = make_group_mapping(
    adata.obs[LEVEL1_KEY], adata.obs[LEVEL1_KEY].astype(str).isin(l1_valid).to_numpy(), "L1_"
)
L2_MAP = make_group_mapping(
    adata.obs[LEVEL2_KEY],
    (adata.obs[LEVEL2_ELIGIBLE_KEY].astype(bool) & adata.obs[LEVEL2_KEY].astype(str).isin(l2_valid)).to_numpy(),
    "L2_",
)
L1_MAP.to_csv(SPEC_ROOT / "Level1_variable_mapping.tsv", sep="\t", index=False)
L2_MAP.to_csv(SPEC_ROOT / "Level2_variable_mapping.tsv", sep="\t", index=False)


def build_specificity(level_name: str, group_key: str, mapping: pd.DataFrame):
    cell_file = SPEC_ROOT / f"{level_name}_cell_weighted_specificity.txt"
    study_file = SPEC_ROOT / f"{level_name}_study_balanced_specificity.txt"
    if cell_file.exists() and study_file.exists():
        print("Reusing specificity:", level_name)
        return cell_file, study_file

    group_to_code = dict(zip(mapping["group"], range(len(mapping))))
    group_labels = adata.obs[group_key].astype(str).map(group_to_code).fillna(-1).to_numpy(int)
    study_labels, study_names = pd.factorize(
        adata.obs["study_canonical"].astype(str), sort=True
    )
    n_group, n_study, n_gene = len(mapping), len(study_names), adata.raw.n_vars
    sums = np.zeros((n_group, n_gene), dtype=np.float32)
    counts = np.zeros(n_group, dtype=np.int64)
    gs_sums = np.zeros((n_group * n_study, n_gene), dtype=np.float32)
    gs_counts = np.zeros(n_group * n_study, dtype=np.int64)

    Xraw = adata.raw.X
    for start in range(0, adata.n_obs, CHUNK_SIZE):
        stop = min(start + CHUNK_SIZE, adata.n_obs)
        codes = group_labels[start:stop]
        keep = codes >= 0
        if not keep.any():
            continue
        block = Xraw[start:stop]
        block = block.tocsr() if sp.issparse(block) else sp.csr_matrix(block)
        block = block[keep].astype(np.float32)
        row_total = np.asarray(block.sum(axis=1)).ravel()
        row_total[row_total <= 0] = 1
        block = sp.diags(1e6 / row_total) @ block
        block.data = np.log1p(block.data)

        c = codes[keep]
        d = sp.csr_matrix(
            (np.ones(len(c), np.float32), (np.arange(len(c)), c)),
            shape=(len(c), n_group),
        )
        sums += (d.T @ block).toarray().astype(np.float32)
        counts += np.bincount(c, minlength=n_group)

        gs = c * n_study + study_labels[start:stop][keep]
        dgs = sp.csr_matrix(
            (np.ones(len(gs), np.float32), (np.arange(len(gs)), gs)),
            shape=(len(gs), n_group * n_study),
        )
        # Convert only group-study pairs present in this chunk to dense form.
        # This avoids a transient (groups x studies x genes) dense allocation.
        observed_gs = np.unique(gs)
        gs_chunk = (dgs.T @ block).tocsr()
        gs_sums[observed_gs] += gs_chunk[observed_gs].toarray().astype(np.float32)
        gs_counts += np.bincount(gs, minlength=n_group * n_study)
        if start % (CHUNK_SIZE * 10) == 0:
            print(level_name, f"{stop:,}/{adata.n_obs:,}")

    means_cell = sums / np.maximum(counts[:, None], 1)
    gs_sums = gs_sums.reshape(n_group, n_study, n_gene)
    gs_counts = gs_counts.reshape(n_group, n_study)
    valid_gs = gs_counts >= MIN_GROUP_STUDY_CELLS
    means_study = np.zeros((n_group, n_gene), dtype=np.float32)
    for g in range(n_group):
        valid_studies = np.flatnonzero(valid_gs[g])
        if len(valid_studies) < MIN_GROUP_STUDIES:
            raise ValueError(
                f"{mapping.iloc[g]['group']} has <{MIN_GROUP_STUDIES} adequately sampled studies"
            )
        per_study = gs_sums[g, valid_studies] / gs_counts[g, valid_studies, None]
        means_study[g] = per_study.mean(axis=0)

    ensg = stable_ensembl_ids(adata.raw.var, pd.Index(adata.raw.var_names))
    valid_gene = ensg.notna() & ensg.str.upper().str.startswith("ENSG")
    ensg_valid = ensg[valid_gene].astype(str).to_numpy()

    def write_matrix(means, output):
        x = means[:, valid_gene.to_numpy()]
        frame = pd.DataFrame(x.T, columns=mapping["variable"].tolist())
        frame.insert(0, "GENE", ensg_valid)
        # Duplicate Ensembl features are averaged after cell-level count aggregation.
        frame = frame.groupby("GENE", sort=False, as_index=False).sum(numeric_only=True)
        vals = frame[mapping["variable"]].to_numpy(float)
        denominator = vals.sum(axis=1)
        keep = denominator > 0
        frame = frame.loc[keep].copy()
        vals = vals[keep]
        denominator = denominator[keep]
        specificity = vals / denominator[:, None]
        out = pd.DataFrame(specificity, columns=mapping["variable"].tolist())
        out.insert(0, "Average", vals.mean(axis=1))
        out.insert(0, "GENE", frame["GENE"].to_numpy())
        out.to_csv(output, sep="\t", index=False, float_format="%.8g")

    write_matrix(means_cell, cell_file)
    write_matrix(means_study, study_file)
    del sums, gs_sums, means_cell, means_study
    gc.collect()
    return cell_file, study_file


L1_CELL, L1_STUDY = build_specificity(LEVEL1_KEY, LEVEL1_KEY, L1_MAP)
L2_CELL, L2_STUDY = build_specificity(LEVEL2_KEY, LEVEL2_KEY, L2_MAP)
SPECIFICATIONS = {
    "Level1_cell_weighted": (L1_CELL, L1_MAP),
    "Level1_study_balanced": (L1_STUDY, L1_MAP),
    "Level2_cell_weighted": (L2_CELL, L2_MAP),
    "Level2_study_balanced": (L2_STUDY, L2_MAP),
}


# =============================================================================
# 4. Run scDRS and MAGMA for primary and sensitivity GWAS inputs
# =============================================================================

scenario_inputs = {}
for name, root in SCENARIOS.items():
    scenario_inputs[name] = {
        "genes_out": resolve_one(root, "magma.genes.out"),
        "genes_raw": resolve_one(root, "magma.genes.raw"),
    }
pd.DataFrame([
    {"scenario": k, **{x: str(y) for x, y in v.items()}}
    for k, v in scenario_inputs.items()
]).to_csv(OUT_ROOT / "scenario_input_manifest.tsv", sep="\t", index=False)


def valid_table(path: Path, minimum_bytes=100):
    return path.exists() and path.stat().st_size >= minimum_bytes


def parse_magma_gsa(path: Path, mapping: pd.DataFrame) -> pd.DataFrame:
    tab = pd.read_csv(path, sep=r"\s+", comment="#")
    tab = tab.loc[tab["VARIABLE"].astype(str) != "Average"].copy()
    tab = tab.merge(mapping, left_on="VARIABLE", right_on="variable", how="left")
    tab["FDR_BH"] = bh_fdr(tab["P"])
    return tab.sort_values("P")


def sample_balanced_scdrs(score_file: Path, output: Path):
    score = pd.read_csv(score_file, sep="\t", index_col=0)
    score_column = next(
        (c for c in ["zscore", "norm_score", "raw_score"] if c in score), None
    )
    if score_column is None:
        raise KeyError(f"No usable score column in {score_file}: {list(score.columns)}")
    values = pd.to_numeric(score[score_column], errors="coerce").reindex(adata.obs_names)
    base = pd.DataFrame({
        "score": values.to_numpy(),
        "sample": adata.obs["biological_sample_uid"].astype(str).to_numpy(),
        "study": adata.obs["study_canonical"].astype(str).to_numpy(),
        "Level1": l1_test.to_numpy(),
        "Level2": l2_test.to_numpy(),
    })
    rows, sample_tables = [], []
    for level, key, excluded in [
        ("Level1", "Level1", "Excluded_low_confidence_or_rare"),
        ("Level2", "Level2", "Excluded_nonstable_or_Level1_only"),
    ]:
        x = base.loc[base[key] != excluded, ["score", "sample", "study", key]].copy()
        agg = (
            x.groupby([key, "sample", "study"], observed=True)["score"]
            .agg(["median", "size"])
            .reset_index()
            .rename(columns={key: "cell_state", "median": "sample_median_score", "size": "n_cells"})
        )
        agg = agg.loc[agg["n_cells"] >= MIN_SAMPLE_STATE_CELLS]
        agg["level"] = level
        sample_tables.append(agg)
        for state, sub in agg.groupby("cell_state", observed=True):
            val = sub["sample_median_score"].to_numpy(float)
            try:
                stat, p = wilcoxon(val, alternative="greater") if len(val) >= 10 else (np.nan, np.nan)
            except ValueError:
                stat, p = np.nan, np.nan
            study_med = sub.groupby("study")["sample_median_score"].median()
            rows.append({
                "level": level,
                "cell_state": state,
                "n_biological_samples": sub["sample"].nunique(),
                "n_studies": sub["study"].nunique(),
                "median_sample_score": np.median(val) if len(val) else np.nan,
                "positive_sample_fraction": np.mean(val > 0) if len(val) else np.nan,
                "positive_study_fraction": np.mean(study_med > 0) if len(study_med) else np.nan,
                "wilcoxon_greater_p": p,
            })
    result = pd.DataFrame(rows)
    result["FDR_BH"] = result.groupby("level")["wilcoxon_greater_p"].transform(bh_fdr)
    result.to_csv(output, sep="\t", index=False)
    pd.concat(sample_tables, ignore_index=True).to_csv(
        output.with_name(output.stem + "__per_sample.tsv.gz"),
        sep="\t", index=False, compression="gzip",
    )


all_summary = []
gene_priority_tables = []

for scenario, inputs in scenario_inputs.items():
    run = OUT_ROOT / scenario
    input_dir, scdrs_dir, magma_dir, table_dir, log_dir = [
        run / x for x in ["01_input", "02_scDRS", "03_MAGMA", "04_tables", "logs"]
    ]
    for d in [input_dir, scdrs_dir, magma_dir, table_dir, log_dir]:
        d.mkdir(parents=True, exist_ok=True)

    genes = pd.read_csv(inputs["genes_out"], sep=r"\s+")
    required = {"SYMBOL", "ZSTAT", "P"}
    if not required.issubset(genes):
        raise KeyError(f"{inputs['genes_out']} lacks {sorted(required - set(genes))}")
    genes["SYMBOL"] = genes["SYMBOL"].astype("string").str.strip().str.upper()
    genes["P"] = pd.to_numeric(genes["P"], errors="coerce")
    genes["ZSTAT"] = pd.to_numeric(genes["ZSTAT"], errors="coerce")
    genes["FDR_BH"] = bh_fdr(genes["P"])
    genes["Bonferroni_0.05"] = genes["P"] < 0.05 / genes["P"].notna().sum()
    gp = genes[["SYMBOL", "P", "ZSTAT", "FDR_BH", "Bonferroni_0.05"]].copy()
    gp["scenario"] = scenario
    gene_priority_tables.append(gp)

    trait = clean_token("PanVascular_G_" + scenario)
    zfile = input_dir / f"{trait}_gene_zscore.tsv"
    gsfile = input_dir / f"{trait}.gs"
    gene_z = genes.loc[
        genes["SYMBOL"].notna() & np.isfinite(genes["ZSTAT"]), ["SYMBOL", "ZSTAT"]
    ].sort_values("ZSTAT", ascending=False).drop_duplicates("SYMBOL")
    gene_z.columns = ["GENE", trait]
    gene_z.to_csv(zfile, sep="\t", index=False)

    score_dir = scdrs_dir / "score"
    l1_down = scdrs_dir / "downstream_Level1"
    l2_down = scdrs_dir / "downstream_Level2"
    for d in [score_dir, l1_down, l2_down]:
        d.mkdir(exist_ok=True)

    if not valid_table(gsfile):
        run_command(scdrs_command(
            "munge-gs", "--out-file", gsfile, "--zscore-file", zfile,
            "--weight", "zscore", "--n-min", "100", "--n-max", "1000",
        ), log_dir / "01_scdrs_munge.log", run)

    score_file = score_dir / f"{trait}.full_score.gz"
    if not valid_table(score_file, 1000):
        run_command(scdrs_command(
            "compute-score", "--h5ad-file", SCDRS_H5AD,
            "--h5ad-species", "human", "--gs-file", gsfile,
            "--gs-species", "human", "--cov-file", SCDRS_COV,
            "--ctrl-match-opt", "mean_var", "--weight-opt", "vs",
            "--flag-filter-data", "True", "--flag-raw-count", "True",
            "--n-ctrl", str(N_CTRL), "--flag-return-ctrl-raw-score", "False",
            "--flag-return-ctrl-norm-score", "True", "--out-folder", score_dir,
        ), log_dir / "02_scdrs_compute_score.log", run)
    if not valid_table(score_file, 1000):
        candidates = sorted(score_dir.glob("*.full_score.gz"))
        if len(candidates) != 1:
            raise RuntimeError(f"Cannot identify scDRS score in {score_dir}: {candidates}")
        score_file = candidates[0]

    group_outputs = {}
    for level, key, outdir in [
        ("Level1", L1_TEST_KEY, l1_down), ("Level2", L2_TEST_KEY, l2_down)
    ]:
        candidates = sorted(outdir.glob(f"*.scdrs_group.{key}"))
        if len(candidates) != 1 or not valid_table(candidates[0]):
            run_command(scdrs_command(
                "perform-downstream", "--h5ad-file", SCDRS_H5AD,
                "--score-file", score_file, "--out-folder", outdir,
                "--group-analysis", key, "--flag-filter-data", "True",
                "--flag-raw-count", "True",
            ), log_dir / f"03_scdrs_{level}.log", run)
            candidates = sorted(outdir.glob(f"*.scdrs_group.{key}"))
        if len(candidates) != 1:
            raise RuntimeError(f"Expected one {level} scDRS group file: {candidates}")
        tab = pd.read_csv(candidates[0], sep="\t")
        tab["assoc_FDR_BH"] = bh_fdr(tab["assoc_mcp"])
        tab["hetero_FDR_BH"] = bh_fdr(tab["hetero_mcp"])
        tab = tab.loc[~tab["group"].astype(str).str.startswith("Excluded_")]
        tab.to_csv(table_dir / f"scDRS_{level}_with_FDR.tsv", sep="\t", index=False)
        group_outputs[level] = tab

    sample_sens = table_dir / "scDRS_sample_balanced_sensitivity.tsv"
    if not valid_table(sample_sens):
        sample_balanced_scdrs(score_file, sample_sens)

    magma_results = {}
    for spec_name, (spec_file, mapping) in SPECIFICATIONS.items():
        prefix = magma_dir / spec_name / f"MAGMA_{spec_name}"
        prefix.parent.mkdir(parents=True, exist_ok=True)
        gsa = Path(str(prefix) + ".gsa.out")
        if not valid_table(gsa):
            run_command([
                MAGMA, "--gene-results", inputs["genes_raw"],
                "--gene-covar", spec_file, "--model", "direction=pos",
                "condition=Average", "--out", prefix,
            ], log_dir / f"04_MAGMA_{spec_name}.log", run)
        parsed = parse_magma_gsa(gsa, mapping)
        parsed.to_csv(table_dir / f"MAGMA_{spec_name}_with_FDR.tsv", sep="\t", index=False)
        magma_results[spec_name] = parsed

    for level in ["Level1", "Level2"]:
        s = group_outputs[level].copy()
        m = magma_results[f"{level}_cell_weighted"].copy()
        merged = s.merge(m, left_on="group", right_on="group", how="outer", suffixes=("_scDRS", "_MAGMA"))
        merged.to_csv(table_dir / f"{level}_scDRS_MAGMA_primary_integrated.tsv", sep="\t", index=False)
        finite = merged.loc[
            np.isfinite(pd.to_numeric(merged.get("assoc_mcz"), errors="coerce"))
            & np.isfinite(pd.to_numeric(merged.get("BETA_STD"), errors="coerce"))
        ]
        rho, p_rho = spearmanr(finite["assoc_mcz"], finite["BETA_STD"]) if len(finite) >= 3 else (np.nan, np.nan)
        all_summary.append({
            "scenario": scenario, "level": level,
            "n_scDRS_groups": len(s), "n_MAGMA_groups": len(m),
            "n_scDRS_FDR_lt_0.05": int((s["assoc_FDR_BH"] < 0.05).sum()),
            "n_MAGMA_FDR_lt_0.05": int((m["FDR_BH"] < 0.05).sum()),
            "n_both_FDR_lt_0.05": int(((merged["assoc_FDR_BH"] < 0.05) & (merged["FDR_BH"] < 0.05)).sum()),
            "scDRS_MAGMA_spearman_r": rho, "scDRS_MAGMA_spearman_p": p_rho,
        })

    (run / "completed.json").write_text(json.dumps({
        "scenario": scenario,
        "genes_out": str(inputs["genes_out"]),
        "genes_raw": str(inputs["genes_raw"]),
        "score_file": str(score_file),
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    print("COMPLETED:", scenario)


# =============================================================================
# 5. Cross-scenario tables for reporting and pseudobulk integration
# =============================================================================

summary = pd.DataFrame(all_summary)
summary.to_csv(OUT_ROOT / "MAGMA_scDRS_cross_scenario_summary.tsv", sep="\t", index=False)

gene_master = pd.concat(gene_priority_tables, ignore_index=True)
wide = gene_master.pivot_table(
    index="SYMBOL", columns="scenario", values=["P", "FDR_BH", "Bonferroni_0.05"],
    aggfunc="first",
)
wide.columns = [f"{metric}__{scenario}" for metric, scenario in wide.columns]
wide = wide.reset_index()
wide["primary_Bonferroni"] = wide.get("Bonferroni_0.05__primary_r0.1", False).fillna(False).astype(bool)
wide["primary_FDR"] = wide.get("FDR_BH__primary_r0.1", 1.0) < 0.05
wide["retained_leave_one_CAD_nominal"] = wide.get("P__leave_one_CAD_r0.1", 1.0) < 0.05
wide["retained_leave_one_CAD_FDR"] = wide.get("FDR_BH__leave_one_CAD_r0.1", 1.0) < 0.05
wide.to_csv(OUT_ROOT / "GWAS_gene_priority_master.tsv", sep="\t", index=False)

print("\nALL ANALYSES COMPLETED")
print("Output:", OUT_ROOT)
print(summary)
