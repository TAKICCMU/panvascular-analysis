#!/usr/bin/env python3
"""Recompute a unified scIB leaderboard from frozen raw metric outputs.

This script does not rerun integration or scIB metrics. It removes the earlier
track-specific scaling, scales every method together within each neighbourhood
size, and distinguishes structural non-applicability from failed metrics.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


K_VALUES = (10, 20, 30)
EXPECTED_RESOLUTIONS = (0.2, 0.5, 1.0, 1.5)

# scIB-conventional grouping. All metrics are already oriented so higher is
# better by the frozen evaluator.
BATCH_METRICS = ("PCR_batch", "iLISI_batch", "ASW_batch", "graph_connectivity")
BIOLOGY_METRICS = ("ASW_label", "NMI_best", "ARI_best", "cLISI_label")

# Metrics available for both low-dimensional embeddings and graph-only methods.
COMMON_BATCH_METRICS = ("iLISI_batch", "graph_connectivity")
COMMON_BIOLOGY_METRICS = ("NMI_best", "ARI_best", "cLISI_label")

# Historical grouping retained only as a sensitivity result.
HISTORICAL_BATCH_METRICS = ("PCR_batch", "iLISI_batch", "ASW_batch")
HISTORICAL_BIOLOGY_METRICS = (
    "graph_connectivity",
    "ASW_label",
    "NMI_best",
    "ARI_best",
    "cLISI_label",
)

LABEL_INFORMED_METHODS = {"scanvi", "scpoli", "cellhint"}
GRAPH_ONLY_METHODS = {"cellhint"}


def sha256(path: Path, block_size: int = 4 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def scale_within_k(
    frame: pd.DataFrame, metrics: tuple[str, ...]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    result = frame.copy()
    bounds: list[dict[str, object]] = []
    for k, positions in result.groupby("K", observed=True, sort=True).groups.items():
        for metric in metrics:
            values = pd.to_numeric(result.loc[positions, metric], errors="coerce")
            finite = values[np.isfinite(values)]
            column = f"{metric}_scaled_unified_within_K"
            result.loc[positions, column] = np.nan
            if finite.empty:
                lo = hi = np.nan
                constant = False
            else:
                lo, hi = float(finite.min()), float(finite.max())
                constant = bool(np.isclose(lo, hi))
                if constant:
                    result.loc[finite.index, column] = 0.5
                else:
                    result.loc[finite.index, column] = (finite - lo) / (hi - lo)
            bounds.append(
                {
                    "K": int(k),
                    "metric": metric,
                    "minimum": lo,
                    "maximum": hi,
                    "constant_metric": constant,
                    "n_finite": int(len(finite)),
                    "n_total": int(len(values)),
                }
            )
    return result, pd.DataFrame(bounds)


def row_mean(frame: pd.DataFrame, metrics: tuple[str, ...]) -> pd.Series:
    columns = [f"{metric}_scaled_unified_within_K" for metric in metrics]
    return frame[columns].mean(axis=1, skipna=True)


def add_metric_contract(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    result["method"] = result["method"].astype(str).str.lower()
    result["uses_input_labels"] = result["method"].isin(LABEL_INFORMED_METHODS)
    result["method_supervision"] = np.where(
        result["uses_input_labels"], "label-informed", "label-free"
    )
    result["output_type"] = np.where(
        result["method"].isin(GRAPH_ONLY_METHODS),
        "graph_only",
        result.get("output_type", "embedding"),
    )

    required_embedding = set(BATCH_METRICS + BIOLOGY_METRICS)
    required_graph = set(COMMON_BATCH_METRICS + COMMON_BIOLOGY_METRICS)
    statuses = []
    missing_required = []
    not_applicable = []
    for _, row in result.iterrows():
        is_graph = str(row["output_type"]) == "graph_only"
        required = required_graph if is_graph else required_embedding
        structural_na = (required_embedding - required_graph) if is_graph else set()
        missing = sorted(
            metric
            for metric in required
            if metric not in result.columns or not np.isfinite(pd.to_numeric(row.get(metric), errors="coerce"))
        )
        statuses.append("complete" if not missing else "incomplete")
        missing_required.append(";".join(missing))
        not_applicable.append(";".join(sorted(structural_na)))
    result["evaluation_status"] = statuses
    result["missing_required_metrics"] = missing_required
    result["structurally_not_applicable_metrics"] = not_applicable
    return result


def summarize_across_k(frame: pd.DataFrame, score_columns: list[str]) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    keys = ["method", "config_id", "method_supervision", "uses_input_labels", "output_type"]
    for key, group in frame.groupby(keys, observed=True, sort=True, dropna=False):
        row = dict(zip(keys, key))
        observed_k = sorted(pd.to_numeric(group["K"], errors="coerce").dropna().astype(int).unique())
        row["K_values_present"] = ";".join(map(str, observed_k))
        row["all_K_present"] = observed_k == list(K_VALUES)
        row["all_rows_complete"] = bool((group["evaluation_status"] == "complete").all())
        row["eligible_for_formal_ranking"] = bool(row["all_K_present"] and row["all_rows_complete"])
        for column in score_columns:
            values = pd.to_numeric(group[column], errors="coerce")
            row[f"mean_{column}"] = float(values.mean()) if values.notna().any() else np.nan
            row[f"min_{column}"] = float(values.min()) if values.notna().any() else np.nan
            row[f"max_{column}"] = float(values.max()) if values.notna().any() else np.nan
        records.append(row)
    summary = pd.DataFrame(records)
    for column in score_columns:
        mean_col = f"mean_{column}"
        rank_col = f"rank_{column}"
        summary[rank_col] = np.nan
        mask = summary["eligible_for_formal_ranking"] & summary[mean_col].notna()
        summary.loc[mask, rank_col] = summary.loc[mask, mean_col].rank(
            ascending=False, method="min"
        )
    return summary.sort_values(
        ["eligible_for_formal_ranking", "mean_unified_overall_score"],
        ascending=[False, False],
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics-tsv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--expected-methods",
        default="pca,harmony,scanorama,liger,scvi,scanvi,scpoli,scalex,cellhint",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(args.metrics_tsv, sep="\t")
    if raw.duplicated(['method', 'config_id', 'K']).any():
        raise RuntimeError('Duplicate method/configuration/K rows would change equal-K weighting')
    required_columns = {
        "method",
        "config_id",
        "output_type",
        "K",
        *BATCH_METRICS,
        *BIOLOGY_METRICS,
    }
    missing_columns = sorted(required_columns.difference(raw.columns))
    if missing_columns:
        raise RuntimeError(f"Raw metric table is missing columns: {missing_columns}")
    observed_k = sorted(pd.to_numeric(raw["K"], errors="coerce").dropna().astype(int).unique())
    if observed_k != list(K_VALUES):
        raise RuntimeError(f"Expected K={K_VALUES}, observed {observed_k}")
    resolution_audit = []
    for _, row in raw.iterrows():
        for metric in ['NMI', 'ARI']:
            cols = [f'{metric}_resolution_{r:g}' for r in EXPECTED_RESOLUTIONS]
            vals = pd.to_numeric(row.reindex(cols), errors='coerce').to_numpy(dtype=float)
            ok = bool(np.isfinite(vals).all() and np.isclose(float(row[f'{metric}_best']), vals.max()))
            resolution_audit.append({'method':row['method'], 'config_id':row['config_id'],
                                     'K':row['K'], 'metric':metric, 'prespecified_resolution_scan_valid':ok})
    pd.DataFrame(resolution_audit).to_csv(args.output_dir/'resolution_scan_audit.tsv', sep='\t',index=False)
    if not all(r['prespecified_resolution_scan_valid'] for r in resolution_audit):
        raise RuntimeError('Missing or inconsistent fixed-resolution NMI/ARI evidence')

    metric_columns = tuple(dict.fromkeys(BATCH_METRICS + BIOLOGY_METRICS))
    scored = add_metric_contract(raw)
    # Exclude an incomplete configuration from scaling all K, so failed outputs
    # cannot alter the bounds for completed methods or receive partial credit.
    for _, group in scored.groupby(['method','config_id']):
        if set(group.K) != set(K_VALUES) or not (group.evaluation_status == 'complete').all():
            scored.loc[group.index,'evaluation_status']='incomplete'
    incomplete_before_scaling = scored.evaluation_status != 'complete'
    scored.loc[incomplete_before_scaling, list(metric_columns)] = np.nan
    scored, bounds = scale_within_k(scored, metric_columns)

    scored["unified_batch_score"] = row_mean(scored, BATCH_METRICS)
    scored["unified_biology_score"] = row_mean(scored, BIOLOGY_METRICS)
    scored["unified_overall_score"] = (
        0.4 * scored["unified_batch_score"] + 0.6 * scored["unified_biology_score"]
    )
    scored["common_graph_batch_score"] = row_mean(scored, COMMON_BATCH_METRICS)
    scored["common_graph_biology_score"] = row_mean(scored, COMMON_BIOLOGY_METRICS)
    scored["common_graph_overall_score"] = (
        0.4 * scored["common_graph_batch_score"]
        + 0.6 * scored["common_graph_biology_score"]
    )
    scored["historical_grouping_batch_score"] = row_mean(scored, HISTORICAL_BATCH_METRICS)
    scored["historical_grouping_biology_score"] = row_mean(scored, HISTORICAL_BIOLOGY_METRICS)
    scored["historical_grouping_overall_score"] = (
        0.4 * scored["historical_grouping_batch_score"]
        + 0.6 * scored["historical_grouping_biology_score"]
    )

    incomplete = scored["evaluation_status"] != "complete"
    formal_scores = [
        "unified_batch_score",
        "unified_biology_score",
        "unified_overall_score",
        "common_graph_batch_score",
        "common_graph_biology_score",
        "common_graph_overall_score",
        "historical_grouping_batch_score",
        "historical_grouping_biology_score",
        "historical_grouping_overall_score",
    ]
    scored.loc[incomplete, formal_scores] = np.nan
    for k, positions in scored.groupby("K", observed=True, sort=True).groups.items():
        for score in ["unified_overall_score", "common_graph_overall_score"]:
            rank = scored.loc[positions, score].rank(ascending=False, method="min")
            scored.loc[positions, f"rank_{score}_within_K"] = rank

    summary = summarize_across_k(scored, formal_scores)
    expected = {x.strip().lower() for x in args.expected_methods.split(",") if x.strip()}
    observed = set(scored["method"])
    missing_methods = sorted(expected.difference(observed))
    unexpected_methods = sorted(observed.difference(expected))

    scored.to_csv(args.output_dir / "unified_scib_metrics_by_K.tsv", sep="\t", index=False)
    summary.to_csv(args.output_dir / "unified_scib_configuration_ranking.tsv", sep="\t", index=False)
    bounds.to_csv(args.output_dir / "unified_metric_scaling_bounds.tsv", sep="\t", index=False)

    formal_complete = bool(not missing_methods and not unexpected_methods and not incomplete.any()
                           and summary['eligible_for_formal_ranking'].all())
    status = {
        "status": "completed" if formal_complete else "incomplete",
        "input": str(args.metrics_tsv.resolve()),
        "input_sha256": sha256(args.metrics_tsv),
        "K": list(K_VALUES),
        "Leiden_resolutions_expected": list(EXPECTED_RESOLUTIONS),
        "scaling_scope": "all included methods and configurations together within K",
        "primary_metric_grouping": {
            "batch": list(BATCH_METRICS),
            "biology": list(BIOLOGY_METRICS),
        },
        "common_graph_metric_grouping": {
            "batch": list(COMMON_BATCH_METRICS),
            "biology": list(COMMON_BIOLOGY_METRICS),
        },
        "composite": "0.4 batch + 0.6 biology",
        "structural_missingness": "excluded only for graph-only metrics declared not applicable",
        "failed_missingness": "marks the row incomplete and suppresses formal scores",
        "expected_methods": sorted(expected),
        "observed_methods": sorted(observed),
        "missing_methods": missing_methods,
        "unexpected_methods": unexpected_methods,
        "n_incomplete_rows": int(incomplete.sum()),
        "formal_ranking_complete": formal_complete,
        "configuration_scaling_eligibility": "Only complete configurations with one row at each K enter scaling; incomplete configurations cannot affect bounds",
        "current_table_interpretation": "formal nine-method comparison" if formal_complete else "provisional completed-method comparison; missing methods can alter scaling and ranking",
    }
    (args.output_dir / "UNIFIED_SCIB_RANKING_STATUS.json").write_text(
        json.dumps(status, indent=2), encoding="utf-8"
    )
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
