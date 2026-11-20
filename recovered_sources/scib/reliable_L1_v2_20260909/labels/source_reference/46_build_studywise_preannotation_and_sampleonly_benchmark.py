#!/usr/bin/env python3
"""Build integration-independent study-wise labels and sample-only scIB inputs.

The biological labels are assigned from prespecified raw-count marker programs
within each study. No integrated embedding, neighbour graph, cluster, UMAP, or
historical annotation participates in label assignment. Benchmark membership is
then sampled solely within sample_uid, independently of every biological label.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


EXPECTED_SHAPE = (1_175_115, 37_336)
SEED = 0
PRIMARY_CAP = 1_500
PILOT_CAP = 500
N_HVGS = 2_000
BATCH_KEY = "sample_uid"
STUDY_KEY = "study_file"
LABEL_KEY = "studywise_marker_level1"
ARGMAX_KEY = "studywise_marker_level1_argmax"
LEGACY_KEY = "legacy_Level1_celltype_final_v3"

# Rare lineages remain available in the complete marker-evidence output but are
# not used to define the primary broad-lineage scIB label. This guards against
# a single ambient transcript spuriously creating a rare class.
PRIMARY_PROGRAMS = [
    "Endothelial cell",
    "Smooth muscle cell",
    "Pericyte",
    "Fibroblast",
    "Fibromyocyte",
    "Macrophage",
    "Monocyte",
    "Dendritic cell",
    "T cell",
    "Natural killer cell",
    "B cell",
    "Plasma cell",
    "Neutrophil",
    "Mast cell",
]


def names_hash(values) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(str(value).encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()


def file_hash(path: Path, block_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def build_studywise_labels(
    source: ad.AnnData,
    marker_dir: Path,
    output_dir: Path,
) -> tuple[pd.DataFrame, dict]:
    manifest_path = marker_dir / "level1_marker_evidence_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    programs = list(manifest["program_order"])
    missing = sorted(set(PRIMARY_PROGRAMS).difference(programs))
    if missing:
        raise RuntimeError(f"Primary marker programs are missing: {missing}")
    expression_path = marker_dir / "level1_marker_log1p_cpm_expression.npy"
    ids_path = marker_dir / "cell_ids.tsv.gz"
    expression = np.load(expression_path, mmap_mode="r")
    marker_ids = pd.read_csv(ids_path, sep="\t")["cell_id"].astype(str)
    if names_hash(marker_ids) != manifest["cell_ids_sha256"]:
        raise RuntimeError("Marker-evidence cell identifiers fail the manifest hash")
    if not np.array_equal(marker_ids.to_numpy(), source.obs_names.astype(str).to_numpy()):
        raise RuntimeError("Marker evidence is not in corrected-atlas cell order")

    marker_genes = list(manifest["marker_gene_order"])
    if expression.shape != (source.n_obs, len(marker_genes)):
        raise RuntimeError("Marker-expression matrix dimensions are inconsistent")
    gene_index = {gene: index for index, gene in enumerate(marker_genes)}
    contract = pd.read_csv(
        marker_dir / "level1_marker_program_contract.tsv", sep="\t"
    ).set_index("program")
    program_gene_positions: list[np.ndarray] = []
    for program in PRIMARY_PROGRAMS:
        genes = str(contract.loc[program, "present_markers"]).split(";")
        positions = np.asarray([gene_index[gene] for gene in genes], dtype=np.int64)
        if len(positions) < 3:
            raise RuntimeError(f"Insufficient markers for {program}")
        program_gene_positions.append(positions)

    # IDF is estimated once from raw-count-derived marker detection, without
    # labels or an integration result. The fixed weights are then applied to
    # every study separately, which keeps the biological label definition
    # comparable across studies while downweighting ubiquitous shared markers.
    detected_cells = np.zeros(len(marker_genes), dtype=np.int64)
    for start in range(0, source.n_obs, 100_000):
        values = np.asarray(expression[start : start + 100_000], dtype=np.float32)
        detected_cells += (values > 0).sum(axis=0)
    idf = np.log((source.n_obs + 1) / (detected_cells + 1)) + 1.0
    pd.DataFrame(
        {
            "marker_gene": marker_genes,
            "detected_cells": detected_cells,
            "detection_fraction": detected_cells / source.n_obs,
            "inverse_detection_frequency_weight": idf,
        }
    ).to_csv(output_dir / "marker_gene_idf_weights.tsv", sep="\t", index=False)

    studies = source.obs[STUDY_KEY].astype(str).to_numpy()
    samples = source.obs[BATCH_KEY].astype(str).to_numpy()
    legacy = (
        source.obs[LEGACY_KEY].astype(str).to_numpy()
        if LEGACY_KEY in source.obs
        else np.repeat("NOT_AVAILABLE", source.n_obs)
    )
    argmax_labels = np.empty(source.n_obs, dtype=object)
    final_labels = np.empty(source.n_obs, dtype=object)
    confidence = np.empty(source.n_obs, dtype=object)
    top_score_all = np.empty(source.n_obs, dtype=np.float32)
    margin_all = np.empty(source.n_obs, dtype=np.float32)
    ratio_all = np.empty(source.n_obs, dtype=np.float32)
    top_detection_all = np.empty(source.n_obs, dtype=np.float32)
    top_n_detected_all = np.empty(source.n_obs, dtype=np.int16)
    study_rows: list[dict] = []

    for study in sorted(pd.unique(studies)):
        positions = np.flatnonzero(studies == study)
        local_expression = np.asarray(expression[positions], dtype=np.float32)
        decision = np.zeros((len(positions), len(PRIMARY_PROGRAMS)), dtype=np.float32)
        detected_counts = np.zeros_like(decision, dtype=np.int16)
        for program_index, gene_positions in enumerate(program_gene_positions):
            values = local_expression[:, gene_positions]
            weights = idf[gene_positions]
            weights = (weights / weights.sum()).astype(np.float32)
            n_detected = (values > 0).sum(axis=1).astype(np.int16)
            # The square-root detection term rewards coherent multi-marker
            # support without forcing low-depth cells to express every marker.
            decision[:, program_index] = (
                (values @ weights)
                * np.sqrt(n_detected / len(gene_positions)).astype(np.float32)
            )
            detected_counts[:, program_index] = n_detected
        observed = detected_counts > 0
        no_marker = ~observed.any(axis=1)
        order = np.argsort(decision, axis=1)
        winner = order[:, -1]
        runner_up = order[:, -2]
        row = np.arange(len(positions))
        top_score = decision[row, winner]
        second_score = decision[row, runner_up]
        margin = top_score - second_score
        ratio = top_score / np.maximum(second_score, 1e-6)
        top_n_detected = detected_counts[row, winner]
        top_program_size = np.asarray(
            [len(program_gene_positions[index]) for index in winner], dtype=np.int16
        )
        top_detection = top_n_detected / top_program_size
        labels = np.asarray(PRIMARY_PROGRAMS, dtype=object)[winner]
        labels[no_marker] = "Unknown"

        # Confidence strata are frozen before any integration result is read.
        # The main scIB label keeps every marker-supported argmax; high/moderate
        # confidence are reported as sensitivity subsets rather than used to
        # determine benchmark membership.
        tiers = np.full(len(positions), "low", dtype=object)
        moderate = (
            (top_n_detected >= 1)
            & (ratio >= 1.20)
            & ~no_marker
        )
        high = (
            (top_n_detected >= 2)
            & (ratio >= 1.50)
            & ~no_marker
        )
        tiers[moderate] = "moderate"
        tiers[high] = "high"
        tiers[no_marker] = "unresolved"

        argmax_labels[positions] = labels
        final_labels[positions] = labels
        confidence[positions] = tiers
        top_score_all[positions] = top_score
        margin_all[positions] = margin
        ratio_all[positions] = ratio
        top_detection_all[positions] = top_detection
        top_n_detected_all[positions] = top_n_detected
        study_rows.append(
            {
                STUDY_KEY: study,
                "n_cells": int(len(positions)),
                "n_samples": int(pd.Series(samples[positions]).nunique()),
                "n_assigned_labels": int(pd.Series(labels).nunique()),
                "unknown_cells": int(no_marker.sum()),
                "high_confidence_fraction": float((tiers == "high").mean()),
                "moderate_or_high_fraction": float(
                    np.isin(tiers, ["moderate", "high"]).mean()
                ),
                "median_top_score": float(np.median(top_score[~no_marker])) if (~no_marker).any() else np.nan,
                "median_top_to_second_ratio": float(np.median(ratio[~no_marker])) if (~no_marker).any() else np.nan,
            }
        )

    table = pd.DataFrame(
        {
            "cell_id": source.obs_names.astype(str),
            "source_row_position": np.arange(source.n_obs, dtype=np.int64),
            STUDY_KEY: studies,
            BATCH_KEY: samples,
            LEGACY_KEY: legacy,
            ARGMAX_KEY: argmax_labels,
            LABEL_KEY: final_labels,
            "studywise_marker_confidence": confidence,
            "marker_top_program_tfidf_score": top_score_all,
            "marker_top_minus_second_tfidf_score": margin_all,
            "marker_top_to_second_score_ratio": ratio_all,
            "marker_top_program_detection_fraction": top_detection_all,
            "marker_top_program_n_markers_detected": top_n_detected_all,
        }
    )
    table_path = output_dir / "studywise_marker_preannotation_all_cells.tsv.gz"
    table.to_csv(table_path, sep="\t", index=False, compression="gzip")
    pd.DataFrame(study_rows).to_csv(
        output_dir / "studywise_preannotation_audit.tsv", sep="\t", index=False
    )
    table.groupby([STUDY_KEY, LABEL_KEY], observed=True).size().rename("n_cells").reset_index().to_csv(
        output_dir / "studywise_preannotation_counts.tsv", sep="\t", index=False
    )
    pd.crosstab(table[LABEL_KEY], table[LEGACY_KEY], margins=True).to_csv(
        output_dir / "studywise_preannotation_vs_legacy.tsv", sep="\t"
    )
    status = {
        "status": "completed",
        "label_key": LABEL_KEY,
        "programs": PRIMARY_PROGRAMS,
        "n_cells": int(len(table)),
        "n_studies": int(table[STUDY_KEY].nunique()),
        "n_samples": int(table[BATCH_KEY].nunique()),
        "unknown_cells": int((table[LABEL_KEY] == "Unknown").sum()),
        "confidence_counts": table["studywise_marker_confidence"].value_counts().to_dict(),
        "assignment_rule": (
            "The same prespecified raw-count marker classifier was applied separately "
            "to each study. Marker expression was weighted by inverse detection "
            "frequency and by the square root of the fraction of programme markers "
            "detected; the highest marker-supported programme was assigned."
        ),
        "independence_contract": (
            "No integrated embedding, neighbour graph, cluster, UMAP, Harmony, "
            "scVI/scANVI, SCALEX, Scanorama, LIGER, scPoli, CellHint, or legacy "
            "annotation was used in assignment. Legacy labels are audit-only."
        ),
        "confidence_contract": {
            "high": "at least 2 top-program markers and top/second score ratio >=1.50",
            "moderate": "at least 1 marker and top/second score ratio >=1.20",
            "low": "marker-supported argmax not satisfying moderate/high rules",
            "unresolved": "no primary-program marker detected",
        },
        "files": {
            "preannotation": str(table_path),
            "marker_manifest": str(manifest_path),
            "marker_manifest_sha256": file_hash(manifest_path),
            "marker_expression_sha256": file_hash(expression_path),
        },
        "generated_at_unix": time.time(),
    }
    (output_dir / "STUDYWISE_PREANNOTATION_COMPLETE.json").write_text(
        json.dumps(status, indent=2), encoding="utf-8"
    )
    return table, status


def sample_only_positions(obs: pd.DataFrame, cap: int, seed: int, eligible: np.ndarray | None = None) -> np.ndarray:
    if eligible is None:
        eligible = np.arange(len(obs), dtype=np.int64)
    frame = pd.DataFrame(
        {
            BATCH_KEY: obs.iloc[eligible][BATCH_KEY].astype(str).to_numpy(),
            "_position": eligible,
        }
    )
    rng = np.random.default_rng(seed)
    selected: list[np.ndarray] = []
    for _, group in frame.groupby(BATCH_KEY, observed=True, sort=True):
        positions = group["_position"].to_numpy(dtype=np.int64)
        if len(positions) > cap:
            positions = rng.choice(positions, size=cap, replace=False)
        selected.append(np.sort(positions))
    return np.sort(np.concatenate(selected))


def write_inputs(
    source: ad.AnnData,
    labels: pd.DataFrame,
    output_dir: Path,
    primary_cap: int,
    pilot_cap: int,
) -> dict:
    obs = source.obs.copy()
    valid_sample = ~obs[BATCH_KEY].astype(str).str.strip().str.lower().isin(
        {"", "unknown", "unassigned", "na", "nan", "none"}
    )
    eligible = np.flatnonzero(valid_sample.to_numpy())
    primary_positions = sample_only_positions(obs, primary_cap, SEED, eligible)
    pilot_positions = sample_only_positions(obs, pilot_cap, SEED, primary_positions)
    if not set(pilot_positions).issubset(set(primary_positions)):
        raise RuntimeError("Pilot sample is not a subset of the primary sample")

    label_values = labels[LABEL_KEY].astype(str).to_numpy()
    confidence_values = labels["studywise_marker_confidence"].astype(str).to_numpy()
    membership = pd.DataFrame(
        {
            "cell_id": obs.index[primary_positions].astype(str),
            "source_row_position": primary_positions,
            "primary_benchmark": True,
            "pilot_benchmark": np.isin(primary_positions, pilot_positions),
            BATCH_KEY: obs.iloc[primary_positions][BATCH_KEY].astype(str).to_numpy(),
            STUDY_KEY: obs.iloc[primary_positions][STUDY_KEY].astype(str).to_numpy(),
            LABEL_KEY: label_values[primary_positions],
            "studywise_marker_confidence": confidence_values[primary_positions],
        }
    )
    membership.to_csv(
        output_dir / "sample_only_benchmark_membership.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    matrix = source.X[primary_positions, :]
    matrix = matrix.tocsr() if sparse.issparse(matrix) else sparse.csr_matrix(matrix)
    matrix.sort_indices()
    primary_obs = obs.iloc[primary_positions].copy()
    primary_obs.index = obs.index[primary_positions].copy()
    primary_obs[LABEL_KEY] = label_values[primary_positions]
    primary_obs["studywise_marker_confidence"] = confidence_values[primary_positions]
    work = ad.AnnData(X=matrix, obs=primary_obs, var=source.var.copy())

    detected = np.asarray(work.X.getnnz(axis=0)).ravel()
    robust = detected / work.n_obs * 100 >= 0.05
    work.var["benchmark_n_cells_detected"] = detected
    work.var["benchmark_percent_cells_detected"] = detected / work.n_obs * 100
    work.var["benchmark_robust"] = robust
    if int(robust.sum()) < N_HVGS:
        raise RuntimeError(f"Only {int(robust.sum())} robust genes remain")
    work = work[:, robust].copy()
    sc.experimental.pp.highly_variable_genes(
        work,
        flavor="pearson_residuals",
        n_top_genes=N_HVGS,
        batch_key=BATCH_KEY,
        check_values=True,
        inplace=True,
    )
    if int(work.var["highly_variable"].sum()) != N_HVGS:
        raise RuntimeError("HVG selection did not return exactly 2,000 genes")
    work.var.to_csv(
        output_dir / "sample_only_hvg_audit.tsv.gz", sep="\t", compression="gzip"
    )
    benchmark = work[:, work.var["highly_variable"].to_numpy()].copy()
    benchmark.layers["counts"] = benchmark.X.copy()
    full_library = pd.to_numeric(benchmark.obs["total_counts"], errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(full_library).all() or np.any(full_library <= 0):
        raise RuntimeError("Invalid full-library size factors")
    benchmark.X = benchmark.layers["counts"].multiply((10_000.0 / full_library)[:, None]).tocsr()
    benchmark.X.data = np.log1p(benchmark.X.data).astype(np.float32, copy=False)

    pilot_ids = set(obs.index[pilot_positions].astype(str))
    pilot_mask = np.fromiter((str(x) in pilot_ids for x in benchmark.obs_names), dtype=bool)
    pilot = benchmark[pilot_mask].copy()
    if pilot.n_obs != len(pilot_positions):
        raise RuntimeError("Pilot membership failed identifier round-trip")

    contract = {
        "contract_version": "2026-09-08.studywise-preannotation.sample-only.v1",
        "source_h5ad": str(source.filename),
        "batch_key": BATCH_KEY,
        "label_key": LABEL_KEY,
        "sampling_keys": [BATCH_KEY],
        "sampling_excludes_biological_label": True,
        "seed": SEED,
        "primary_cap_per_sample_uid": int(primary_cap),
        "pilot_cap_per_sample_uid": int(pilot_cap),
        "n_primary_cells": int(benchmark.n_obs),
        "n_pilot_cells": int(pilot.n_obs),
        "n_samples": int(benchmark.obs[BATCH_KEY].nunique()),
        "n_studies": int(benchmark.obs[STUDY_KEY].nunique()),
        "n_labels": int(benchmark.obs[LABEL_KEY].nunique()),
        "n_hvgs": N_HVGS,
        "hvg_flavor": "pearson_residuals",
        "hvg_batch_key": BATCH_KEY,
        "normalization": "raw HVG counts / audited full-library total_counts * 10,000, then log1p",
        "pca": "not stored in H5AD; authoritative OmicVerse mixed-mode PCA is generated as a hash-locked sidecar",
        "python": platform.python_version(),
        "generated_at_unix": time.time(),
    }
    benchmark.uns["benchmark_contract"] = contract | {"dataset_role": "primary"}
    pilot.uns["benchmark_contract"] = contract | {"dataset_role": "pilot"}
    primary_path = output_dir / "scib_primary_sampleonly_cap1500_hvg2000.h5ad"
    pilot_path = output_dir / "scib_pilot_sampleonly_cap500_hvg2000.h5ad"
    benchmark.write_h5ad(primary_path, compression="lzf")
    pilot.write_h5ad(pilot_path, compression="lzf")

    outputs = []
    for path in [primary_path, pilot_path]:
        check = ad.read_h5ad(path, backed="r")
        outputs.append(
            {
                "path": str(path),
                "shape": [int(check.n_obs), int(check.n_vars)],
                "counts_layer_present": "counts" in check.layers,
                "obs_names_unique": bool(check.obs_names.is_unique),
                "var_names_unique": bool(check.var_names.is_unique),
                "cell_ids_sha256": names_hash(check.obs_names),
                "sha256": file_hash(path),
                "size_bytes": int(path.stat().st_size),
            }
        )
        check.file.close()
    result = {"status": "completed", "configuration": contract, "outputs": outputs}
    (output_dir / "SAMPLE_ONLY_BENCHMARK_INPUTS_COMPLETE.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", required=True, type=Path)
    parser.add_argument("--marker-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--primary-cap", type=int, default=PRIMARY_CAP)
    parser.add_argument("--pilot-cap", type=int, default=PILOT_CAP)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    started = time.time()
    source = ad.read_h5ad(args.counts, backed="r")
    if source.shape != EXPECTED_SHAPE:
        raise RuntimeError(f"Unexpected corrected-atlas shape: {source.shape}")
    for key in [BATCH_KEY, STUDY_KEY, "total_counts"]:
        if key not in source.obs:
            raise KeyError(f"Missing required obs field: {key}")
    labels, label_status = build_studywise_labels(source, args.marker_dir, args.output_dir)
    benchmark_status = write_inputs(
        source, labels, args.output_dir, args.primary_cap, args.pilot_cap
    )
    source.file.close()
    final = {
        "status": "completed",
        "preannotation": label_status,
        "benchmark": benchmark_status,
        "elapsed_seconds": time.time() - started,
    }
    (args.output_dir / "PREANNOTATION_AND_INPUTS_COMPLETE.json").write_text(
        json.dumps(final, indent=2), encoding="utf-8"
    )
    print(json.dumps(final, indent=2), flush=True)


if __name__ == "__main__":
    main()
