#!/usr/bin/env python3
"""Source-stratified sensitivity analysis for Tier 1 spatial localisation.

This script uses the same two-stage estimator as the primary spatial analysis:
inverse-variance fixed effects across serial sections from one biological sample,
followed by DerSimonian-Laird random effects across biological samples. It also
applies the primary metadata eligibility filter and controls the false discovery
rate separately within each analysis stratum.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


SOURCE_CLASSIFICATION = pd.DataFrame(
    [
        ("GSE241346", "no_identified_atlas_source_overlap", "", "no_identified_source_overlap"),
        ("GSE243179", "no_identified_atlas_source_overlap", "", "no_identified_source_overlap"),
        ("GSE248608", "same_study_documented_donor_overlap", "GSE248610 (single-cell component GSE248609)", "same_project_or_donor_overlap"),
        ("GSE274572", "no_identified_atlas_source_overlap", "", "no_identified_source_overlap"),
        ("GSE293427", "no_identified_atlas_source_overlap", "", "no_identified_source_overlap"),
        ("GSE306288", "same_publication_project_overlap_unresolved", "GSE273980", "same_project_or_donor_overlap"),
        ("GSE314851", "no_identified_atlas_source_overlap", "", "no_identified_source_overlap"),
    ],
    columns=["GSE", "source_tier", "atlas_related_accession", "analysis_group"],
)


def bh(p_values: pd.Series) -> pd.Series:
    values = pd.to_numeric(p_values, errors="coerce").to_numpy(float)
    adjusted = np.full(values.size, np.nan)
    valid = np.isfinite(values)
    if not valid.any():
        return pd.Series(adjusted, index=p_values.index)
    observed = values[valid]
    order = np.argsort(observed)
    ranked = observed[order]
    q = ranked * ranked.size / np.arange(1, ranked.size + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    restored = np.empty_like(q)
    restored[order] = np.clip(q, 0, 1)
    adjusted[np.flatnonzero(valid)] = restored
    return pd.Series(adjusted, index=p_values.index)


def two_sided_normal_p(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0))


def combine_serial_sections(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    group_columns = ["GSE", "biological_unit_id", "analysis_group", "source_tier"]
    for keys, group in frame.groupby(group_columns, dropna=False, sort=False):
        beta = pd.to_numeric(group["beta"], errors="coerce").to_numpy(float)
        se = pd.to_numeric(group["se"], errors="coerce").to_numpy(float)
        keep = np.isfinite(beta) & np.isfinite(se) & (se > 0)
        if not keep.any():
            continue
        weights = 1.0 / np.square(se[keep])
        rows.append(
            {
                **dict(zip(group_columns, keys)),
                "n_sections_in_unit": int(keep.sum()),
                "beta": float(np.sum(weights * beta[keep]) / np.sum(weights)),
                "se": float(math.sqrt(1.0 / np.sum(weights))),
            }
        )
    return pd.DataFrame(rows)


def dersimonian_laird(beta: np.ndarray, se: np.ndarray) -> dict[str, float]:
    keep = np.isfinite(beta) & np.isfinite(se) & (se > 0)
    effects = beta[keep]
    errors = se[keep]
    k = effects.size
    if k == 0:
        return {"k": 0}
    fixed_weights = 1.0 / np.square(errors)
    fixed_effect = float(np.sum(fixed_weights * effects) / np.sum(fixed_weights))
    q = float(np.sum(fixed_weights * np.square(effects - fixed_effect))) if k > 1 else 0.0
    c = float(np.sum(fixed_weights) - np.sum(np.square(fixed_weights)) / np.sum(fixed_weights)) if k > 1 else 0.0
    tau2 = max(0.0, (q - (k - 1)) / c) if k > 1 and c > 0 else 0.0
    weights = 1.0 / (np.square(errors) + tau2)
    effect = float(np.sum(weights * effects) / np.sum(weights))
    meta_se = float(math.sqrt(1.0 / np.sum(weights)))
    z = effect / meta_se
    i2 = max(0.0, (q - (k - 1)) / q) * 100 if q > 0 and k > 1 else 0.0
    return {
        "k": int(k),
        "beta": effect,
        "se": meta_se,
        "z": z,
        "p": two_sided_normal_p(z),
        "CI_low": effect - 1.96 * meta_se,
        "CI_high": effect + 1.96 * meta_se,
        "tau2": tau2,
        "I2": i2,
    }


def run(section_models: Path, primary_meta: Path, output_dir: Path) -> None:
    frame = pd.read_csv(section_models, sep="\t")
    required = {
        "GSE", "section_id", "biological_unit_id", "metadata_status", "pair_id",
        "gene", "context_level", "cell_context", "parent_lineage", "pathology",
        "expected_direction", "posterior", "model_scope", "term", "beta", "se",
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"Section-model input is missing columns: {missing}")

    frame = frame[~frame["metadata_status"].astype(str).str.contains("needs_", na=False)].copy()
    frame = frame.merge(SOURCE_CLASSIFICATION, on="GSE", how="left", validate="many_to_one")
    if frame["analysis_group"].isna().any():
        unknown = sorted(frame.loc[frame["analysis_group"].isna(), "GSE"].unique())
        raise ValueError(f"Unclassified spatial datasets: {unknown}")

    identifiers = [
        "pair_id", "gene", "context_level", "cell_context", "parent_lineage",
        "pathology", "expected_direction", "posterior", "model_scope", "term",
    ]
    strata = ["all_public", "no_identified_source_overlap", "same_project_or_donor_overlap"]
    rows = []
    for keys, group in frame.groupby(identifiers, dropna=False, sort=False):
        for stratum in strata:
            selected = group if stratum == "all_public" else group[group["analysis_group"] == stratum]
            if selected.empty:
                continue
            units = combine_serial_sections(selected)
            result = dersimonian_laird(units["beta"].to_numpy(float), units["se"].to_numpy(float))
            if result["k"] == 0:
                continue
            rows.append(
                {
                    **dict(zip(identifiers, keys)),
                    "source_stratum": stratum,
                    "n_datasets": int(units["GSE"].nunique()),
                    "n_biological_units": int(units["biological_unit_id"].nunique()),
                    "n_sections": int(selected["section_id"].nunique()),
                    **result,
                }
            )

    results = pd.DataFrame(rows)
    results["FDR"] = results.groupby(
        ["source_stratum", "posterior", "model_scope", "term"], dropna=False
    )["p"].transform(bh)
    focused = results[
        (results["posterior"] == "q05")
        & (results["model_scope"] == "localization")
        & (results["term"] == "cell_abundance")
    ].copy()

    primary = pd.read_csv(primary_meta, sep="\t")
    primary = primary[
        (primary["scope"] == "confirmed_or_unique_biological_units")
        & (primary["posterior"] == "q05")
        & (primary["model_scope"] == "localization")
        & (primary["term"] == "cell_abundance")
    ].copy()
    comparison = primary.merge(
        focused[focused["source_stratum"] == "all_public"],
        on="pair_id",
        suffixes=("_primary", "_rebuilt"),
        validate="one_to_one",
    )
    differences = {
        field: float(np.nanmax(np.abs(comparison[f"{field}_primary"] - comparison[f"{field}_rebuilt"])))
        for field in ["beta", "se", "p", "FDR"]
    }
    if len(comparison) != len(primary) or any(value > 1e-10 for value in differences.values()):
        raise RuntimeError(f"All-public stratum did not reproduce the primary analysis: {differences}")

    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_dir / "Tier1_spatial_meta_by_source_independence.tsv.gz", sep="\t", index=False, compression="gzip")
    focused.to_csv(output_dir / "Tier1_q05_localization_by_source_independence.tsv", sep="\t", index=False)
    SOURCE_CLASSIFICATION.to_csv(output_dir / "spatial_source_overlap_audit.tsv", sep="\t", index=False)
    summary = {
        "primary_reproduction_max_absolute_difference": differences,
        "all_public_FDR_lt_0.05": int((focused.query("source_stratum == 'all_public'")["FDR"] < 0.05).sum()),
        "no_overlap_FDR_lt_0.05": int((focused.query("source_stratum == 'no_identified_source_overlap'")["FDR"] < 0.05).sum()),
        "no_overlap_direction_concordance": int((
            np.sign(focused.query("source_stratum == 'all_public'").set_index("pair_id")["beta"])
            == np.sign(focused.query("source_stratum == 'no_identified_source_overlap'").set_index("pair_id")["beta"])
        ).sum()),
    }
    (output_dir / "run_manifest.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--section-models", required=True, type=Path)
    parser.add_argument("--primary-meta", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args.section_models, args.primary_meta, args.output_dir)
