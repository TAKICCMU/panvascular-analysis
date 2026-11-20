#!/usr/bin/env python3
"""Build a corrected post-QC full-gene raw-count atlas container.

The locked source object is used only as a cell-identifier-locked count and
metadata container.  All old embeddings, graphs, highly-variable-gene results,
clusters and annotations are discarded.  Previous annotations are retained
only under a ``legacy_`` prefix for benchmarking and later audit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


EXCLUDED_STUDY_FILES = {"GSE198750", "GSE216860"}
EXPECTED_INPUT_CELLS = 1_728_873
EXPECTED_RAW_GENES = 37_336
EXPECTED_OUTPUT_CELLS = 1_175_115
ENSG_RE = re.compile(r"^ENSG\d{11}(?:_PAR_Y)?$")


def sha256(path: Path, block_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit-dir", type=Path, required=True)
    parser.add_argument("--compression", choices=["lzf", "gzip"], default="lzf")
    args = parser.parse_args()
    args.audit_dir.mkdir(parents=True, exist_ok=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.output.with_suffix(args.output.suffix + ".partial")
    if tmp.exists():
        tmp.unlink()

    started = time.time()
    print(f"Opening locked atlas: {args.input}", flush=True)
    source = ad.read_h5ad(args.input, backed="r")
    if source.n_obs != EXPECTED_INPUT_CELLS:
        raise RuntimeError(f"Unexpected input cell count: {source.n_obs}")
    if source.raw is None or source.raw.n_vars != EXPECTED_RAW_GENES:
        raise RuntimeError("The expected full-gene raw-count layer is unavailable")
    if "study_file" not in source.obs:
        raise KeyError("obs['study_file'] is required for exact source exclusion")

    study_file = source.obs["study_file"].astype(str)
    keep_mask = ~study_file.isin(EXCLUDED_STUDY_FILES).to_numpy()
    keep_index = np.flatnonzero(keep_mask)
    excluded = study_file[~keep_mask].value_counts().sort_index()
    if len(keep_index) != EXPECTED_OUTPUT_CELLS:
        raise RuntimeError(
            f"Expected {EXPECTED_OUTPUT_CELLS} retained cells, found {len(keep_index)}"
        )
    if excluded.to_dict() != {"GSE198750": 480405, "GSE216860": 73353}:
        raise RuntimeError(f"Unexpected exclusion counts: {excluded.to_dict()}")

    var = source.raw.var.copy()
    required_var = {"gene_ids", "gene_id_clean", "gene_symbol"}
    missing_var = sorted(required_var.difference(var.columns))
    if missing_var:
        raise KeyError(f"Missing raw.var fields: {missing_var}")
    clean = (
        var["gene_id_clean"]
        .astype("object")
        .where(var["gene_id_clean"].notna(), "")
        .astype(str)
        .str.strip()
        .str.replace(r"\.\d+(?=_PAR_Y$|$)", "", regex=True)
    )
    valid_ensg = clean.map(lambda value: bool(ENSG_RE.fullmatch(value)))
    if clean[valid_ensg].duplicated().any():
        examples = clean[valid_ensg][clean[valid_ensg].duplicated(keep=False)].head(20).tolist()
        raise RuntimeError(f"Non-unique cleaned ENSG identifiers: {examples}")

    # The archived merge deliberately retained features that lacked a recoverable
    # Ensembl identifier by falling back to their standardized gene symbol.  Keep
    # that information rather than silently discarding 735 expression columns.
    display_symbol = pd.Index(var.index.astype(str)).str.strip()
    if display_symbol.has_duplicates or (display_symbol == "").any():
        raise RuntimeError("Raw feature symbols must be non-empty and unique")
    fallback_key = pd.Series(
        [f"SYMBOL_FALLBACK::{symbol}" for symbol in display_symbol],
        index=var.index,
        dtype="object",
    )
    stable_key = clean.where(valid_ensg, fallback_key)
    if stable_key.duplicated().any():
        examples = stable_key[stable_key.duplicated(keep=False)].head(20).tolist()
        raise RuntimeError(f"Non-unique stable feature identifiers: {examples}")

    # Use ordinary object strings here. Nullable StringArray columns require an
    # AnnData opt-in that is not portable to older reader versions.
    var["gene_id_clean"] = clean.where(valid_ensg, "").astype("object")
    var["stable_ensembl_id"] = clean.where(valid_ensg, "").astype("object")
    var["display_gene_symbol"] = display_symbol.to_numpy()
    var["stable_feature_id"] = stable_key.to_numpy()
    var["merge_key"] = stable_key.to_numpy()
    var["gene_mapping_status"] = np.where(valid_ensg, "ensembl", "symbol_fallback")

    mapping_audit = var[
        [
            "display_gene_symbol",
            "stable_ensembl_id",
            "stable_feature_id",
            "gene_mapping_status",
            "gene_ids",
        ]
    ].copy()
    mapping_audit.index.name = "feature_name"
    mapping_audit.to_csv(
        args.audit_dir / "full_gene_identifier_mapping.tsv.gz",
        sep="\t",
        compression="gzip",
    )

    obs = source.obs.iloc[keep_index].copy()
    legacy_columns = [
        col
        for col in obs.columns
        if col.startswith("Level1_")
        or col.startswith("Level2_")
        or col.lower().startswith("leiden")
    ]
    rename = {col: f"legacy_{col}" for col in legacy_columns}
    obs.rename(columns=rename, inplace=True)
    obs["atlas_rebuild_inclusion"] = "retained_after_source_deduplication"
    obs["atlas_rebuild_version"] = "panvascular_atlas_rebuild_v2_20260907"
    if not obs.index.is_unique:
        duplicate_examples = obs.index[obs.index.duplicated()].astype(str)[:20].tolist()
        raise RuntimeError(f"Retained cell identifiers are not unique: {duplicate_examples}")

    print(
        f"Reading {len(keep_index):,} retained rows from raw counts "
        f"({source.raw.n_vars:,} genes)",
        flush=True,
    )
    matrix = source.raw.X[keep_index, :]
    matrix = matrix.tocsr() if sparse.issparse(matrix) else sparse.csr_matrix(matrix)
    matrix.sort_indices()
    if matrix.shape != (EXPECTED_OUTPUT_CELLS, EXPECTED_RAW_GENES):
        raise RuntimeError(f"Unexpected retained count shape: {matrix.shape}")
    if not np.isfinite(matrix.data).all() or (matrix.data < 0).any():
        raise RuntimeError("Count matrix contains negative or non-finite values")
    rounded = np.rint(matrix.data)
    if not np.allclose(matrix.data, rounded, rtol=0, atol=1e-6):
        raise RuntimeError("The purported raw-count layer contains non-integer values")
    matrix.data = rounded.astype(np.float32, copy=False)

    output = ad.AnnData(X=matrix, obs=obs, var=var)
    output.uns["rebuild_provenance"] = {
        "source_h5ad": str(args.input),
        "source_shape": [int(source.n_obs), int(source.n_vars)],
        "source_raw_shape": [int(source.raw.n_obs), int(source.raw.n_vars)],
        "excluded_study_file_values": sorted(EXCLUDED_STUDY_FILES),
        "excluded_cells_by_study_file": {k: int(v) for k, v in excluded.items()},
        "retained_cells": int(output.n_obs),
        "genes": int(output.n_vars),
        "matrix_semantics": "post_QC_raw_UMI_counts",
        "gene_merge_key": (
            "version-stripped Ensembl gene identifier; standardized-symbol "
            "fallback retained when the archived source lacked an Ensembl ID"
        ),
        "n_ensembl_mapped_features": int(valid_ensg.sum()),
        "n_symbol_fallback_features": int((~valid_ensg).sum()),
        "legacy_annotation_columns": sorted(rename.values()),
        "discarded_objects": "source X, obsm, varm, obsp, varp, neighbours, embeddings and uns",
    }
    output.uns["analysis_status"] = "counts_and_metadata_only_before_reintegration"

    print(f"Writing temporary H5AD: {tmp}", flush=True)
    output.write_h5ad(tmp, compression=args.compression)
    del output, matrix
    try:
        source.file.close()
    except Exception:
        pass

    print("Reopening temporary H5AD for contract checks", flush=True)
    check = ad.read_h5ad(tmp, backed="r")
    contract = {
        "shape_ok": tuple(check.shape) == (EXPECTED_OUTPUT_CELLS, EXPECTED_RAW_GENES),
        "obs_names_unique": bool(check.obs_names.is_unique),
        "var_names_unique": bool(check.var_names.is_unique),
        "stable_merge_key_unique": bool(check.var["merge_key"].astype(str).is_unique),
        "gene_mapping_partition_ok": bool(
            check.var["gene_mapping_status"].astype(str).value_counts().to_dict()
            == {"ensembl": 36601, "symbol_fallback": 735}
        ),
        "excluded_study_files_absent": not bool(
            check.obs["study_file"].astype(str).isin(EXCLUDED_STUDY_FILES).any()
        ),
        "old_embeddings_absent": len(check.obsm_keys()) == 0,
        "legacy_annotation_fields_present": all(x in check.obs for x in rename.values()),
    }
    try:
        check.file.close()
    except Exception:
        pass
    if not all(contract.values()):
        raise RuntimeError(f"Temporary H5AD failed contract: {contract}")

    os.replace(tmp, args.output)
    digest = sha256(args.output)
    elapsed = time.time() - started
    result = {
        "status": "completed",
        "output_h5ad": str(args.output),
        "sha256": digest,
        "file_size_bytes": int(args.output.stat().st_size),
        "shape": [EXPECTED_OUTPUT_CELLS, EXPECTED_RAW_GENES],
        "excluded_cells_by_study_file": {k: int(v) for k, v in excluded.items()},
        "contract": contract,
        "elapsed_seconds": elapsed,
    }
    (args.audit_dir / "corrected_counts_build_status.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    pd.DataFrame(
        [{"path": str(args.output), "sha256": digest, "size_bytes": args.output.stat().st_size}]
    ).to_csv(args.audit_dir / "corrected_counts_sha256.tsv", sep="\t", index=False)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
