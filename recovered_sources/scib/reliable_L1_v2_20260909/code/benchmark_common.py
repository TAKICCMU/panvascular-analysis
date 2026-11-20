#!/usr/bin/env python3
"""Shared contracts for the rebuilt pan-vascular scIB benchmark."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import random
import time
from pathlib import Path

import anndata as ad
import numpy as np


BATCH_KEY = os.environ.get("SCIB_BATCH_KEY", "sample_uid")
LABEL_KEY = os.environ.get("SCIB_LABEL_KEY", "studywise_marker_level1")
SEED = 0


def set_seed(seed: int = SEED) -> None:
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch

        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except Exception:
        pass


def obs_names_sha256(names) -> str:
    digest = hashlib.sha256()
    for name in names:
        digest.update(str(name).encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def environment_record(extra_packages: list[str] | None = None) -> dict:
    packages = ["anndata", "numpy", "pandas", "scanpy", "scipy"]
    if extra_packages:
        packages.extend(extra_packages)
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": {name: package_version(name) for name in sorted(set(packages))},
    }


def load_benchmark(path: Path, require_counts: bool = False) -> ad.AnnData:
    adata = ad.read_h5ad(path)
    for key in [BATCH_KEY, LABEL_KEY]:
        if key not in adata.obs:
            raise KeyError(f"Required obs field is missing: {key}")
    if require_counts and "counts" not in adata.layers:
        raise KeyError("Required counts layer is missing")
    if (
        "X_pca" not in adata.obsm
        and os.environ.get("SCIB_ALLOW_PCA_SIDECAR_ONLY", "0") != "1"
    ):
        raise KeyError("Required baseline X_pca is missing")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise RuntimeError("Benchmark identifiers must be unique")
    adata.obs[BATCH_KEY] = adata.obs[BATCH_KEY].astype(str).astype("category")
    adata.obs[LABEL_KEY] = adata.obs[LABEL_KEY].astype(str).astype("category")
    return adata


def save_embedding(
    embedding: np.ndarray,
    adata: ad.AnnData,
    output_dir: Path,
    method: str,
    config_id: str,
    parameters: dict,
    started: float,
    extra_packages: list[str] | None = None,
    extra: dict | None = None,
    *,
    seed: int = SEED,
) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    embedding = np.asarray(embedding, dtype=np.float32)
    if embedding.ndim != 2 or embedding.shape[0] != adata.n_obs:
        raise RuntimeError(
            f"Embedding shape {embedding.shape} is incompatible with {adata.n_obs} cells"
        )
    if not np.isfinite(embedding).all():
        raise RuntimeError("Embedding contains non-finite values")
    stem = f"{method}__{config_id}"
    npy_path = output_dir / f"{stem}.npy"
    tmp_path = output_dir / f"{stem}.npy.partial"
    with tmp_path.open("wb") as handle:
        np.save(handle, embedding, allow_pickle=False)
    tmp_path.replace(npy_path)
    record = {
        "status": "completed",
        "method": method,
        "config_id": config_id,
        "parameters": parameters,
        "input_shape": [int(adata.n_obs), int(adata.n_vars)],
        "embedding_shape": list(embedding.shape),
        "obs_names_sha256": obs_names_sha256(adata.obs_names),
        "batch_key": BATCH_KEY,
        "label_key": LABEL_KEY,
        "n_batches": int(adata.obs[BATCH_KEY].nunique()),
        "n_labels": int(adata.obs[LABEL_KEY].nunique()),
        "seed": int(seed),
        "elapsed_seconds": time.time() - started,
        "environment": environment_record(extra_packages),
    }
    if extra:
        record.update(extra)
    json_path = output_dir / f"{stem}.json"
    json_path.write_text(json.dumps(record, indent=2), encoding="utf-8")
    return npy_path, json_path
