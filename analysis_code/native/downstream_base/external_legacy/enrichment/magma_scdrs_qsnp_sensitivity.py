from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import tarfile
import time
import zipfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon


# =============================================================================
# 0. Fixed paths and analysis contract
# =============================================================================

GWAS_ROOT = Path("/path/to/magma_scdrs")
OLD_RESULT_ROOT = GWAS_ROOT / "final_atlas_MAGMA_scDRS_v1"
COMMON = OLD_RESULT_ROOT / "00_common_input"
INPUT_ROOT = GWAS_ROOT / "QSNP_sensitivity_inputs_v1"
INPUT_ARCHIVE = GWAS_ROOT / "MAGMA_scDRS_unfiltered_exactQ_FUMA_inputs.zip"
OUT_ROOT = GWAS_ROOT / "final_atlas_MAGMA_scDRS_QSNP_sensitivity_v1"

SCDRS_H5AD = COMMON / "final_atlas_full_gene_raw_counts_scDRS.h5ad"
SCDRS_COV = COMMON / "final_atlas_scDRS.cov"
SPEC_ROOT = COMMON / "MAGMA_specificity"

LEVEL1_GROUP_KEY = "Level1_scDRS_group_v1"
LEVEL2_GROUP_KEY = "Level2_scDRS_group_v1"
EXPECTED_LEVEL1_GROUPS = 16
EXPECTED_LEVEL2_GROUPS = 71

CONDA = Path("/path/to/miniconda3/bin/conda")
SCDRS_ENV = os.environ.get("PANVAS_SCDRS_ENV", "scdrs_env")
MAGMA = Path(os.environ.get("PANVAS_MAGMA", "/path/to/scrnaseq/magma"))

N_CTRL = int(os.environ.get("PANVAS_SCDRS_N_CTRL", "1000"))
MIN_SAMPLE_STATE_CELLS = int(os.environ.get("PANVAS_MIN_SAMPLE_STATE_CELLS", "20"))
RANDOM_SEED = int(os.environ.get("PANVAS_RANDOM_SEED", "20260824"))

SCENARIOS = {
    "unfiltered": {
        "input_dir": INPUT_ROOT / "unfiltered",
        "trait": "PanVascular_G_unfiltered",
    },
    "exact_Q_only": {
        "input_dir": INPUT_ROOT / "exact_Q_only",
        "trait": "PanVascular_G_exact_Q_only",
    },
}

SPECIFICATIONS = {
    "Level1_cell_weighted": (
        SPEC_ROOT / "Level1_celltype_final_v3_cell_weighted_specificity.txt",
        SPEC_ROOT / "Level1_variable_mapping.tsv",
    ),
    "Level1_study_balanced": (
        SPEC_ROOT / "Level1_celltype_final_v3_study_balanced_specificity.txt",
        SPEC_ROOT / "Level1_variable_mapping.tsv",
    ),
    "Level2_cell_weighted": (
        SPEC_ROOT / "Level2_celltype_for_enrichment_v1_cell_weighted_specificity.txt",
        SPEC_ROOT / "Level2_variable_mapping.tsv",
    ),
    "Level2_study_balanced": (
        SPEC_ROOT / "Level2_celltype_for_enrichment_v1_study_balanced_specificity.txt",
        SPEC_ROOT / "Level2_variable_mapping.tsv",
    ),
}

OUT_ROOT.mkdir(parents=True, exist_ok=True)


# =============================================================================
# 1. Utilities
# =============================================================================

def sha256(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def clean_token(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", str(value)).strip("_")


def valid_file(path: Path, minimum_bytes: int = 100) -> bool:
    return path.exists() and path.stat().st_size >= minimum_bytes


def bh_fdr(values) -> np.ndarray:
    p = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(float)
    result = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    if valid.sum() == 0:
        return result
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    restored = np.empty_like(q)
    restored[order] = q
    result[valid] = restored
    return result


def safe_extract_zip(archive: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    root = destination.resolve()
    with zipfile.ZipFile(archive) as handle:
        for member in handle.infolist():
            target = (destination / member.filename).resolve()
            if root not in target.parents and target != root:
                raise RuntimeError(f"Unsafe ZIP member: {member.filename}")
        handle.extractall(destination)


def run_streaming(args, log_file: Path, cwd: Path | None = None) -> None:
    args = [str(x) for x in args]
    log_file.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"
    env["PYTHONHASHSEED"] = str(RANDOM_SEED)
    env["PYTHONUNBUFFERED"] = "1"
    print("$", " ".join(args), flush=True)
    with log_file.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            args,
            cwd=str(cwd) if cwd else None,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="", flush=True)
            log.write(line)
            log.flush()
        return_code = process.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, args)


def scdrs_command(*args) -> list[str | Path]:
    return [
        CONDA,
        "run",
        "--no-capture-output",
        "-n",
        SCDRS_ENV,
        "scdrs",
        *args,
    ]


def parse_magma_gsa(path: Path, mapping_path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep=r"\s+", comment="#")
    if "VARIABLE" not in table or "P" not in table:
        raise ValueError(f"Unexpected MAGMA GSA columns in {path}: {list(table.columns)}")
    table = table.loc[table["VARIABLE"].astype(str) != "Average"].copy()
    mapping = pd.read_csv(mapping_path, sep="\t")
    table = table.merge(mapping, left_on="VARIABLE", right_on="variable", how="left")
    if table["group"].isna().any():
        bad = table.loc[table["group"].isna(), "VARIABLE"].tolist()
        raise ValueError(f"Unmapped MAGMA variables in {path}: {bad[:10]}")
    table["FDR_BH"] = bh_fdr(table["P"])
    table["Bonferroni_0.05"] = (
        pd.to_numeric(table["P"], errors="coerce") < 0.05 / max(1, len(table))
    )
    return table.sort_values("P")


def locate_group_output(folder: Path, group_key: str) -> Path:
    candidates = sorted(folder.glob(f"*.scdrs_group.{group_key}"))
    candidates = [path for path in candidates if valid_file(path)]
    if len(candidates) != 1:
        raise RuntimeError(
            f"Expected exactly one scDRS group output for {group_key} below "
            f"{folder}; found: {candidates}"
        )
    return candidates[0]


def parse_scdrs_group(path: Path, excluded_prefix: str) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    required = {"group", "assoc_mcp", "assoc_mcz", "hetero_mcp", "hetero_mcz"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"Missing scDRS columns in {path}: {sorted(missing)}")
    table = table.loc[~table["group"].astype(str).str.startswith(excluded_prefix)].copy()
    table["assoc_FDR_BH"] = bh_fdr(table["assoc_mcp"])
    table["hetero_FDR_BH"] = bh_fdr(table["hetero_mcp"])
    return table.sort_values("assoc_mcp")


# =============================================================================
# 2. Input extraction and strict preflight audit
# =============================================================================

if not INPUT_ROOT.exists() and INPUT_ARCHIVE.exists():
    print("Extracting FUMA inputs:", INPUT_ARCHIVE)
    safe_extract_zip(INPUT_ARCHIVE, INPUT_ROOT)

required_common = [SCDRS_H5AD, SCDRS_COV, CONDA, MAGMA]
for specificity, mapping in SPECIFICATIONS.values():
    required_common.extend([specificity, mapping])

missing_common = [str(path) for path in required_common if not valid_file(path, 1)]
if missing_common:
    raise FileNotFoundError(
        "Existing reusable MAGMA/scDRS inputs are incomplete:\n" + "\n".join(missing_common)
    )

input_manifest_rows = []
for scenario, config in SCENARIOS.items():
    for filename in ["magma.genes.out", "magma.genes.raw", "params.config"]:
        path = config["input_dir"] / filename
        if not valid_file(path):
            raise FileNotFoundError(
                f"Missing {path}. Upload/extract the supplied FUMA input ZIP first."
            )
        input_manifest_rows.append(
            {
                "scenario": scenario,
                "file": filename,
                "path": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )

print("Auditing reusable scDRS H5AD metadata without loading the expression matrix...")
backed = ad.read_h5ad(SCDRS_H5AD, backed="r")
required_obs = {
    LEVEL1_GROUP_KEY,
    LEVEL2_GROUP_KEY,
    "biological_sample_uid",
    "study_canonical",
}
missing_obs = required_obs - set(backed.obs.columns)
if missing_obs:
    backed.file.close()
    raise KeyError(f"Reusable scDRS H5AD lacks: {sorted(missing_obs)}")

obs = backed.obs[
    [LEVEL1_GROUP_KEY, LEVEL2_GROUP_KEY, "biological_sample_uid", "study_canonical"]
].copy()
obs.index = obs.index.astype(str)
backed.file.close()

level1_groups = sorted(
    x for x in obs[LEVEL1_GROUP_KEY].astype(str).unique()
    if not x.startswith("Excluded_")
)
level2_groups = sorted(
    x for x in obs[LEVEL2_GROUP_KEY].astype(str).unique()
    if not x.startswith("Excluded_")
)
if len(level1_groups) != EXPECTED_LEVEL1_GROUPS:
    raise ValueError(
        f"Expected {EXPECTED_LEVEL1_GROUPS} Level 1 groups, found {len(level1_groups)}"
    )
if len(level2_groups) != EXPECTED_LEVEL2_GROUPS:
    raise ValueError(
        f"Expected {EXPECTED_LEVEL2_GROUPS} locked Level 2 groups, found {len(level2_groups)}"
    )

pd.DataFrame(input_manifest_rows).to_csv(
    OUT_ROOT / "scenario_input_manifest.tsv", sep="\t", index=False
)
(OUT_ROOT / "preflight_audit.json").write_text(
    json.dumps(
        {
            "status": "passed",
            "scdrs_h5ad": str(SCDRS_H5AD),
            "scdrs_cov": str(SCDRS_COV),
            "n_cells": int(len(obs)),
            "n_Level1_groups": len(level1_groups),
            "n_Level2_groups": len(level2_groups),
            "Level1_groups": level1_groups,
            "Level2_groups": level2_groups,
            "N_CTRL": N_CTRL,
            "random_seed_environment": RANDOM_SEED,
        },
        indent=2,
        ensure_ascii=False,
    )
    + "\n",
    encoding="utf-8",
)
print(f"Preflight passed: {len(obs):,} cells; {len(level1_groups)} Level 1; {len(level2_groups)} Level 2")


# =============================================================================
# 3. Sample-balanced scDRS sensitivity
# =============================================================================

def sample_balanced_scdrs(score_file: Path, output: Path) -> None:
    header = pd.read_csv(score_file, sep="\t", nrows=0)
    index_column = header.columns[0]
    score_column = next(
        (column for column in ["zscore", "norm_score", "raw_score"] if column in header),
        None,
    )
    if score_column is None:
        raise KeyError(f"No usable score column in {score_file}: {list(header.columns)}")
    score = pd.read_csv(
        score_file,
        sep="\t",
        usecols=[index_column, score_column],
        index_col=0,
    )
    score.index = score.index.astype(str)
    aligned = pd.to_numeric(score[score_column], errors="coerce").reindex(obs.index)
    if aligned.notna().mean() < 0.999:
        raise ValueError(
            f"Only {aligned.notna().mean():.3%} of scDRS scores align to H5AD cells"
        )

    base = obs.copy()
    base["score"] = aligned.to_numpy(float)
    rows = []
    per_sample_tables = []
    for level, key in [("Level1", LEVEL1_GROUP_KEY), ("Level2", LEVEL2_GROUP_KEY)]:
        data = base.loc[
            ~base[key].astype(str).str.startswith("Excluded_"),
            [key, "biological_sample_uid", "study_canonical", "score"],
        ].copy()
        aggregate = (
            data.groupby(
                [key, "biological_sample_uid", "study_canonical"], observed=True
            )["score"]
            .agg(["median", "size"])
            .reset_index()
            .rename(
                columns={
                    key: "cell_state",
                    "median": "sample_median_score",
                    "size": "n_cells",
                }
            )
        )
        aggregate = aggregate.loc[aggregate["n_cells"] >= MIN_SAMPLE_STATE_CELLS].copy()
        aggregate["level"] = level
        per_sample_tables.append(aggregate)

        for state, subset in aggregate.groupby("cell_state", observed=True):
            values = subset["sample_median_score"].to_numpy(float)
            try:
                _, p_value = (
                    wilcoxon(values, alternative="greater")
                    if len(values) >= 10 else (np.nan, np.nan)
                )
            except ValueError:
                p_value = np.nan
            study_medians = subset.groupby("study_canonical")["sample_median_score"].median()
            rows.append(
                {
                    "level": level,
                    "cell_state": state,
                    "n_biological_samples": subset["biological_sample_uid"].nunique(),
                    "n_studies": subset["study_canonical"].nunique(),
                    "median_sample_score": np.median(values) if len(values) else np.nan,
                    "positive_sample_fraction": np.mean(values > 0) if len(values) else np.nan,
                    "positive_study_fraction": (
                        np.mean(study_medians > 0) if len(study_medians) else np.nan
                    ),
                    "wilcoxon_greater_p": p_value,
                }
            )

    result = pd.DataFrame(rows)
    result["FDR_BH"] = result.groupby("level")["wilcoxon_greater_p"].transform(bh_fdr)
    result.to_csv(output, sep="\t", index=False)
    pd.concat(per_sample_tables, ignore_index=True).to_csv(
        output.with_name(output.stem + "__per_sample.tsv.gz"),
        sep="\t",
        index=False,
        compression="gzip",
    )


# =============================================================================
# 4. Run four MAGMA specificity models and one scDRS score per scenario
# =============================================================================

summary_rows = []
long_rows = []

for scenario, config in SCENARIOS.items():
    print("\n" + "=" * 88)
    print("START SCENARIO:", scenario)
    print("=" * 88)
    started = time.time()
    run_root = OUT_ROOT / scenario
    input_dir = run_root / "01_input"
    scdrs_dir = run_root / "02_scDRS"
    magma_dir = run_root / "03_MAGMA"
    table_dir = run_root / "04_tables"
    log_dir = run_root / "logs"
    for folder in [input_dir, scdrs_dir, magma_dir, table_dir, log_dir]:
        folder.mkdir(parents=True, exist_ok=True)

    genes_out = config["input_dir"] / "magma.genes.out"
    genes_raw = config["input_dir"] / "magma.genes.raw"
    genes = pd.read_csv(genes_out, sep=r"\s+")
    required = {"SYMBOL", "ZSTAT", "P"}
    missing = required - set(genes.columns)
    if missing:
        raise KeyError(f"{genes_out} lacks {sorted(missing)}")
    genes["SYMBOL"] = genes["SYMBOL"].astype("string").str.strip().str.upper()
    genes["ZSTAT"] = pd.to_numeric(genes["ZSTAT"], errors="coerce")
    genes["P"] = pd.to_numeric(genes["P"], errors="coerce")
    genes["FDR_BH"] = bh_fdr(genes["P"])
    genes["Bonferroni_0.05"] = genes["P"] < 0.05 / genes["P"].notna().sum()
    genes.to_csv(table_dir / "MAGMA_gene_results_with_FDR.tsv.gz", sep="\t", index=False)

    trait = clean_token(config["trait"])
    zscore_file = input_dir / f"{trait}_gene_zscore.tsv"
    gene_set_file = input_dir / f"{trait}.gs"
    gene_z = (
        genes.loc[
            genes["SYMBOL"].notna() & np.isfinite(genes["ZSTAT"]),
            ["SYMBOL", "ZSTAT"],
        ]
        .sort_values("ZSTAT", ascending=False)
        .drop_duplicates("SYMBOL")
    )
    gene_z.columns = ["GENE", trait]
    gene_z.to_csv(zscore_file, sep="\t", index=False)

    # MAGMA is quick, so run it before the long scDRS calculation.
    magma_tables = {}
    for specification, (specificity_file, mapping_file) in SPECIFICATIONS.items():
        prefix = magma_dir / specification / f"MAGMA_{specification}"
        prefix.parent.mkdir(parents=True, exist_ok=True)
        gsa_file = Path(str(prefix) + ".gsa.out")
        if not valid_file(gsa_file):
            run_streaming(
                [
                    MAGMA,
                    "--gene-results",
                    genes_raw,
                    "--gene-covar",
                    specificity_file,
                    "--model",
                    "direction=pos",
                    "condition=Average",
                    "--out",
                    prefix,
                ],
                log_dir / f"MAGMA_{specification}.log",
                run_root,
            )
        parsed = parse_magma_gsa(gsa_file, mapping_file)
        output = table_dir / f"MAGMA_{specification}_with_FDR.tsv"
        parsed.to_csv(output, sep="\t", index=False)
        magma_tables[specification] = parsed

    if not valid_file(gene_set_file):
        run_streaming(
            scdrs_command(
                "munge-gs",
                "--out-file",
                gene_set_file,
                "--zscore-file",
                zscore_file,
                "--weight",
                "zscore",
                "--n-min",
                "100",
                "--n-max",
                "1000",
            ),
            log_dir / "scDRS_01_munge_gs.log",
            run_root,
        )

    score_dir = scdrs_dir / "score"
    downstream_dir = scdrs_dir / "downstream_combined"
    score_dir.mkdir(exist_ok=True)
    downstream_dir.mkdir(exist_ok=True)
    score_file = score_dir / f"{trait}.full_score.gz"
    if not valid_file(score_file, 1000):
        run_streaming(
            scdrs_command(
                "compute-score",
                "--h5ad-file",
                SCDRS_H5AD,
                "--h5ad-species",
                "human",
                "--gs-file",
                gene_set_file,
                "--gs-species",
                "human",
                "--cov-file",
                SCDRS_COV,
                "--ctrl-match-opt",
                "mean_var",
                "--weight-opt",
                "vs",
                "--flag-filter-data",
                "True",
                "--flag-raw-count",
                "True",
                "--n-ctrl",
                str(N_CTRL),
                "--flag-return-ctrl-raw-score",
                "False",
                "--flag-return-ctrl-norm-score",
                "True",
                "--out-folder",
                score_dir,
            ),
            log_dir / "scDRS_02_compute_score.log",
            run_root,
        )
    if not valid_file(score_file, 1000):
        candidates = sorted(score_dir.glob("*.full_score.gz"))
        if len(candidates) != 1:
            raise RuntimeError(f"Cannot identify scDRS score in {score_dir}: {candidates}")
        score_file = candidates[0]

    group_outputs_ready = all(
        len([p for p in downstream_dir.glob(f"*.scdrs_group.{key}") if valid_file(p)]) == 1
        for key in [LEVEL1_GROUP_KEY, LEVEL2_GROUP_KEY]
    )
    if not group_outputs_ready:
        run_streaming(
            scdrs_command(
                "perform-downstream",
                "--h5ad-file",
                SCDRS_H5AD,
                "--score-file",
                score_file,
                "--out-folder",
                downstream_dir,
                "--group-analysis",
                f"{LEVEL1_GROUP_KEY},{LEVEL2_GROUP_KEY}",
                "--flag-filter-data",
                "False",
                "--flag-raw-count",
                "False",
            ),
            log_dir / "scDRS_03_perform_downstream_combined.log",
            run_root,
        )

    scdrs_tables = {}
    for level, key in [
        ("Level1", LEVEL1_GROUP_KEY),
        ("Level2", LEVEL2_GROUP_KEY),
    ]:
        group_file = locate_group_output(downstream_dir, key)
        parsed = parse_scdrs_group(group_file, "Excluded_")
        expected = EXPECTED_LEVEL1_GROUPS if level == "Level1" else EXPECTED_LEVEL2_GROUPS
        if len(parsed) != expected:
            raise ValueError(f"{scenario} {level}: expected {expected} groups, found {len(parsed)}")
        parsed.to_csv(table_dir / f"scDRS_{level}_with_FDR.tsv", sep="\t", index=False)
        scdrs_tables[level] = parsed

    sample_sensitivity = table_dir / "scDRS_sample_balanced_sensitivity.tsv"
    if not valid_file(sample_sensitivity):
        sample_balanced_scdrs(score_file, sample_sensitivity)

    for level in ["Level1", "Level2"]:
        scdrs = scdrs_tables[level]
        magma = magma_tables[f"{level}_cell_weighted"]
        integrated = scdrs.merge(magma, on="group", how="outer", suffixes=("_scDRS", "_MAGMA"))
        integrated.to_csv(
            table_dir / f"{level}_scDRS_MAGMA_integrated.tsv",
            sep="\t",
            index=False,
        )
        finite = integrated.loc[
            np.isfinite(pd.to_numeric(integrated["assoc_mcz"], errors="coerce"))
            & np.isfinite(pd.to_numeric(integrated["BETA_STD"], errors="coerce"))
        ]
        rho, rho_p = (
            spearmanr(finite["assoc_mcz"], finite["BETA_STD"])
            if len(finite) >= 3 else (np.nan, np.nan)
        )
        summary_rows.append(
            {
                "scenario": scenario,
                "level": level,
                "n_groups": len(integrated),
                "n_scDRS_FDR_lt_0.05": int((scdrs["assoc_FDR_BH"] < 0.05).sum()),
                "n_MAGMA_FDR_lt_0.05": int((magma["FDR_BH"] < 0.05).sum()),
                "n_both_FDR_lt_0.05": int(
                    (
                        (integrated["assoc_FDR_BH"] < 0.05)
                        & (integrated["FDR_BH"] < 0.05)
                    ).sum()
                ),
                "scDRS_MAGMA_spearman_rho": rho,
                "scDRS_MAGMA_spearman_P": rho_p,
            }
        )
        for _, row in integrated.iterrows():
            long_rows.append(
                {
                    "scenario": scenario,
                    "level": level,
                    "group": row["group"],
                    "MAGMA_BETA_STD": row.get("BETA_STD", np.nan),
                    "MAGMA_P": row.get("P", np.nan),
                    "MAGMA_FDR_BH": row.get("FDR_BH", np.nan),
                    "scDRS_assoc_mcz": row.get("assoc_mcz", np.nan),
                    "scDRS_assoc_P": row.get("assoc_mcp", np.nan),
                    "scDRS_assoc_FDR_BH": row.get("assoc_FDR_BH", np.nan),
                }
            )

    (run_root / "completed.json").write_text(
        json.dumps(
            {
                "status": "completed",
                "scenario": scenario,
                "genes_out": str(genes_out),
                "genes_raw": str(genes_raw),
                "score_file": str(score_file),
                "elapsed_seconds": time.time() - started,
                "N_CTRL": N_CTRL,
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    print("COMPLETED:", scenario)


# =============================================================================
# 5. Cross-scenario sensitivity tables
# =============================================================================

summary = pd.DataFrame(summary_rows)
long_table = pd.DataFrame(long_rows)
summary.to_csv(OUT_ROOT / "MAGMA_scDRS_QSNP_sensitivity_summary.tsv", sep="\t", index=False)
long_table.to_csv(OUT_ROOT / "MAGMA_scDRS_QSNP_sensitivity_long.tsv", sep="\t", index=False)

comparison_rows = []
for level in ["Level1", "Level2"]:
    left = long_table.loc[long_table["scenario"] == "unfiltered"]
    left = left.loc[left["level"] == level].set_index("group")
    right = long_table.loc[long_table["scenario"] == "exact_Q_only"]
    right = right.loc[right["level"] == level].set_index("group")
    common_groups = left.index.intersection(right.index)
    for method, effect, fdr in [
        ("MAGMA", "MAGMA_BETA_STD", "MAGMA_FDR_BH"),
        ("scDRS", "scDRS_assoc_mcz", "scDRS_assoc_FDR_BH"),
    ]:
        x = pd.to_numeric(left.loc[common_groups, effect], errors="coerce")
        y = pd.to_numeric(right.loc[common_groups, effect], errors="coerce")
        valid = np.isfinite(x) & np.isfinite(y)
        rho, rho_p = spearmanr(x[valid], y[valid]) if valid.sum() >= 3 else (np.nan, np.nan)
        left_sig = set(common_groups[pd.to_numeric(left.loc[common_groups, fdr], errors="coerce") < 0.05])
        right_sig = set(common_groups[pd.to_numeric(right.loc[common_groups, fdr], errors="coerce") < 0.05])
        comparison_rows.append(
            {
                "level": level,
                "method": method,
                "n_common_groups": len(common_groups),
                "effect_spearman_rho": rho,
                "effect_spearman_P": rho_p,
                "direction_agreement_fraction": float(np.mean(np.sign(x[valid]) == np.sign(y[valid]))),
                "n_FDR_unfiltered": len(left_sig),
                "n_FDR_exact_Q_only": len(right_sig),
                "n_FDR_shared": len(left_sig & right_sig),
                "FDR_Jaccard": (
                    len(left_sig & right_sig) / len(left_sig | right_sig)
                    if (left_sig | right_sig) else 1.0
                ),
            }
        )

comparison = pd.DataFrame(comparison_rows)
comparison.to_csv(
    OUT_ROOT / "unfiltered_vs_exact_Q_cell_enrichment_concordance.tsv",
    sep="\t",
    index=False,
)


# =============================================================================
# 6. Software provenance and compact review archive
# =============================================================================

def capture(command: list[str | Path]) -> str:
    completed = subprocess.run(
        [str(x) for x in command],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env={**os.environ, "MPLBACKEND": "Agg"},
    )
    return completed.stdout.strip()


versions = {
    "jupyter_python": capture(["python", "--version"]),
    "scdrs_environment": SCDRS_ENV,
    "scdrs_version_probe": capture(
        [CONDA, "run", "-n", SCDRS_ENV, "python", "-c", "import scdrs; print(getattr(scdrs, '__version__', 'unknown'))"]
    ),
    "magma_path": str(MAGMA),
    "magma_version_probe": capture([MAGMA, "--version"]),
}
(OUT_ROOT / "software_and_run_provenance.json").write_text(
    json.dumps(versions, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
)

review_root = OUT_ROOT / "review_bundle"
if review_root.exists():
    shutil.rmtree(review_root)
review_root.mkdir(parents=True)
for filename in [
    "scenario_input_manifest.tsv",
    "preflight_audit.json",
    "MAGMA_scDRS_QSNP_sensitivity_summary.tsv",
    "MAGMA_scDRS_QSNP_sensitivity_long.tsv",
    "unfiltered_vs_exact_Q_cell_enrichment_concordance.tsv",
    "software_and_run_provenance.json",
]:
    shutil.copy2(OUT_ROOT / filename, review_root / filename)

for scenario in SCENARIOS:
    source = OUT_ROOT / scenario
    target = review_root / scenario
    target.mkdir()
    for relative in ["01_input", "03_MAGMA", "04_tables", "logs"]:
        src = source / relative
        if src.exists():
            shutil.copytree(src, target / relative)
    shutil.copy2(source / "completed.json", target / "completed.json")

(review_root / "README.txt").write_text(
    "This compact review bundle excludes the very large per-cell scDRS full_score files.\n"
    "The full score files remain under each scenario/02_scDRS/score on the server.\n",
    encoding="utf-8",
)

archive = OUT_ROOT.with_name(OUT_ROOT.name + "_review_bundle.tar.gz")
with tarfile.open(archive, "w:gz") as handle:
    handle.add(review_root, arcname=review_root.name)

print("\nALL ANALYSES COMPLETED")
print(summary.to_string(index=False))
print("\nUnfiltered versus exact-Q concordance")
print(comparison.to_string(index=False))
print("\nDownload this compact archive:")
print(archive)
print("Archive size MiB:", round(archive.stat().st_size / 1024**2, 2))
