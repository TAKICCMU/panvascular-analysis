from __future__ import annotations

import copy
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tarfile
import time
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


# =============================================================================
# Import the frozen locked-71 contract without running its main analysis
# =============================================================================

PROJECT_ROOT = Path("/path/to/panvas_project")
COMMON_SCRIPT = PROJECT_ROOT / "16_spatial_locked71_full_reanalysis_Jupyter一键运行版.py"
if not COMMON_SCRIPT.exists():
    raise FileNotFoundError(
        f"Upload the companion locked-71 script before running circularity analysis: {COMMON_SCRIPT}"
    )
spec = importlib.util.spec_from_file_location("locked71_spatial_common", COMMON_SCRIPT)
if spec is None or spec.loader is None:
    raise ImportError(f"Cannot import {COMMON_SCRIPT}")
common = importlib.util.module_from_spec(spec)
spec.loader.exec_module(common)


# =============================================================================
# Analysis paths and target-exclusion contract
# =============================================================================

PRIMARY_ROOT = Path(os.environ.get(
    "SPATIAL_LOCKED71_ROOT",
    str(PROJECT_ROOT / "spatial_validation_v5_locked71_gpu_final"),
))
SENSITIVITY_ROOT = Path(os.environ.get(
    "SPATIAL_LOCKED71_CIRCULARITY_ROOT",
    str(PROJECT_ROOT / "spatial_validation_v5_locked71_Tier1_excluded_gpu_v1"),
))
HYBRID_ROOT = Path(os.environ.get(
    "SPATIAL_LOCKED71_CIRCULARITY_HYBRID_ROOT",
    str(PROJECT_ROOT / "spatial_validation_v5_locked71_Tier1_excluded_hybrid_v1"),
))

# Exact submitted Tier 1 genes: 22 gene-cell-state pairs, 21 genes.
TIER1_EXCLUDED_GENES = {
    "HP", "LRRC10B", "C17ORF58", "KCNK5", "CORO6", "SREBF1", "ZNF589",
    "CCNL1", "HDAC9", "ATP1B1", "SCARB1", "ARHGEF12", "BAG5", "SBF2",
    "BCAP29", "SERPINH1", "DCLRE1B", "GMDS", "FHL3", "FN1", "KCNK3",
}
EXPECTED_PAIRS = 22
WAIT_HOURS = float(os.environ.get("SPATIAL_PRIMARY_WAIT_HOURS", "168"))
WAIT_SECONDS = int(os.environ.get("SPATIAL_PRIMARY_POLL_SECONDS", "300"))


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def abundance_frame(adata: ad.AnnData, key: str) -> pd.DataFrame:
    value = adata.obsm[key]
    if isinstance(value, pd.DataFrame):
        result = value.copy()
        result.index = adata.obs_names
        result.columns = list(map(str, result.columns))
        return result
    names = adata.uns.get("mod", {}).get("factor_names")
    if names is None or len(names) != value.shape[1]:
        raise ValueError(f"Cannot resolve factor names for {key}")
    return pd.DataFrame(
        np.asarray(value),
        index=adata.obs_names,
        columns=list(map(str, names)),
    )


def primary_ready() -> tuple[bool, list[str]]:
    missing = []
    for gse in common.expected_gses():
        path = PRIMARY_ROOT / "02_cell2location" / f"{gse}__final_atlas_cell2location.h5ad"
        if not path.exists() or path.stat().st_size < 10000:
            missing.append(str(path))
    tier1 = (
        PRIMARY_ROOT
        / "05_all_public_Tier1_spatial_validation_v1"
        / "01_tables"
        / "Tier1_pair_evidence_summary.tsv"
    )
    if not tier1.exists():
        missing.append(str(tier1))
    return len(missing) == 0, missing


def wait_for_primary() -> None:
    deadline = time.time() + WAIT_HOURS * 3600
    while True:
        ready, missing = primary_ready()
        if ready:
            common.validate_locked71_output(PRIMARY_ROOT)
            print("Primary locked-71 mapping and Tier 1 results are ready.")
            return
        if time.time() >= deadline:
            raise TimeoutError(
                "Primary locked-71 run did not complete before the waiting deadline. "
                "Missing examples:\n" + "\n".join(missing[:10])
            )
        print(
            f"Waiting for primary locked-71 analysis ({len(missing)} checkpoints missing). "
            f"Next check in {WAIT_SECONDS} s.",
            flush=True,
        )
        time.sleep(WAIT_SECONDS)


def hybridise_primary_expression_with_target_excluded_abundance() -> pd.DataFrame:
    for subdir in ("00_audit", "01_reference", "02_cell2location", "scripts", "logs"):
        (HYBRID_ROOT / subdir).mkdir(parents=True, exist_ok=True)

    for name in (
        "spatial_reference_state_label_map.tsv",
        "reference_audit.json",
    ):
        source = SENSITIVITY_ROOT / "01_reference" / name
        if source.exists():
            shutil.copy2(source, HYBRID_ROOT / "01_reference" / name)
    contract = SENSITIVITY_ROOT / "00_audit" / "circularity_exclusion_contract.json"
    if contract.exists():
        shutil.copy2(contract, HYBRID_ROOT / "00_audit" / contract.name)

    rows = []
    for gse in common.expected_gses():
        primary_path = PRIMARY_ROOT / "02_cell2location" / f"{gse}__final_atlas_cell2location.h5ad"
        sensitivity_path = (
            SENSITIVITY_ROOT / "02_cell2location" / f"{gse}__final_atlas_cell2location.h5ad"
        )
        output_path = HYBRID_ROOT / "02_cell2location" / primary_path.name
        if output_path.exists() and output_path.stat().st_size > 10000:
            print("Reusing hybrid object:", output_path)
            rows.append({"GSE": gse, "hybrid_output": str(output_path), "reused": True})
            continue
        print("Hybridising:", gse, flush=True)
        primary = ad.read_h5ad(primary_path)
        sensitivity = ad.read_h5ad(sensitivity_path)
        if set(primary.obs_names) != set(sensitivity.obs_names):
            raise ValueError(f"Spot identity mismatch for {gse}")
        sensitivity = sensitivity[primary.obs_names].copy()
        for key in ("q05_cell_abundance_w_sf", "means_cell_abundance_w_sf"):
            primary_frame = abundance_frame(primary, key)
            sensitivity_frame = abundance_frame(sensitivity, key)
            if set(primary_frame.columns) != set(sensitivity_frame.columns):
                raise ValueError(f"Factor-set mismatch for {gse} {key}")
            if len(sensitivity_frame.columns) != 71:
                raise ValueError(f"Target-excluded abundance has {len(sensitivity_frame.columns)} factors")
            primary.obsm[key] = sensitivity_frame.loc[
                primary.obs_names, primary_frame.columns
            ].copy()
        primary.uns["mod"] = copy.deepcopy(sensitivity.uns.get("mod", {}))
        primary.uns["Tier1_circularity_sensitivity"] = {
            "expression_source": str(primary_path),
            "abundance_source": str(sensitivity_path),
            "excluded_mapping_genes": sorted(TIER1_EXCLUDED_GENES),
            "reference_basis": "71 locked Level 2 states",
        }
        primary.obs.index.name = None
        primary.var.index.name = None
        primary.write_h5ad(output_path, compression="gzip")
        rows.append({
            "GSE": gse,
            "n_spots": int(primary.n_obs),
            "n_expression_genes": int(primary.n_vars),
            "n_target_excluded_factors": 71,
            "expression_source": str(primary_path),
            "abundance_source": str(sensitivity_path),
            "hybrid_output": str(output_path),
            "reused": False,
        })
        del primary, sensitivity
    manifest = pd.DataFrame(rows)
    manifest.to_csv(HYBRID_ROOT / "00_audit" / "hybrid_object_manifest.tsv", sep="\t", index=False)
    return manifest


def compare_primary_and_target_excluded() -> dict:
    primary_path = (
        PRIMARY_ROOT
        / "05_all_public_Tier1_spatial_validation_v1"
        / "01_tables"
        / "Tier1_pair_evidence_summary.tsv"
    )
    sensitivity_path = (
        HYBRID_ROOT
        / "05_all_public_Tier1_spatial_validation_v1"
        / "01_tables"
        / "Tier1_pair_evidence_summary.tsv"
    )
    primary = pd.read_csv(primary_path, sep="\t")
    sensitivity = pd.read_csv(sensitivity_path, sep="\t")
    if len(primary) != EXPECTED_PAIRS or len(sensitivity) != EXPECTED_PAIRS:
        raise ValueError(
            f"Expected {EXPECTED_PAIRS} Tier 1 rows, observed primary={len(primary)}, "
            f"target-excluded={len(sensitivity)}"
        )
    required = [
        "pair_id", "localization_beta", "localization_p", "localization_FDR",
        "interaction_beta", "interaction_p", "interaction_FDR",
    ]
    for name, table in (("primary", primary), ("target-excluded", sensitivity)):
        missing = [column for column in required if column not in table.columns]
        if missing:
            raise KeyError(f"{name} Tier 1 table is missing columns: {missing}")
    optional = ["spatial_support_class"]
    keep = required + [column for column in optional if column in primary.columns and column in sensitivity.columns]
    comparison = primary[keep].merge(
        sensitivity[keep],
        on="pair_id",
        suffixes=("_primary", "_target_excluded"),
        validate="one_to_one",
    )
    for prefix in ("localization", "interaction"):
        comparison[f"{prefix}_direction_agrees"] = (
            np.sign(comparison[f"{prefix}_beta_primary"])
            == np.sign(comparison[f"{prefix}_beta_target_excluded"])
        )
        comparison[f"{prefix}_absolute_beta_change"] = np.abs(
            comparison[f"{prefix}_beta_target_excluded"]
            - comparison[f"{prefix}_beta_primary"]
        )
        denominator = comparison[f"{prefix}_beta_primary"].abs().replace(0, np.nan)
        comparison[f"{prefix}_relative_beta_change"] = (
            comparison[f"{prefix}_absolute_beta_change"] / denominator
        )
    comparison["localization_FDR_support_retained"] = (
        (comparison["localization_FDR_primary"] < 0.05)
        & (comparison["localization_FDR_target_excluded"] < 0.05)
    )

    output = HYBRID_ROOT / "06_circularity_comparison"
    output.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(
        output / "Tier1_locked71_primary_vs_target_excluded.tsv",
        sep="\t",
        index=False,
    )
    summary = {
        "generated_at": now_iso(),
        "n_pairs_compared": int(len(comparison)),
        "reference_basis_primary": "71 locked Level 2 factors",
        "reference_basis_target_excluded": "71 locked Level 2 factors",
        "n_excluded_Tier1_genes": len(TIER1_EXCLUDED_GENES),
        "localization_beta_spearman": float(
            comparison[["localization_beta_primary", "localization_beta_target_excluded"]]
            .corr(method="spearman")
            .iloc[0, 1]
        ),
        "localization_direction_agreement_fraction": float(
            comparison["localization_direction_agrees"].mean()
        ),
        "interaction_direction_agreement_fraction": float(
            comparison["interaction_direction_agrees"].mean()
        ),
        "n_primary_localization_FDR_supported": int(
            (comparison["localization_FDR_primary"] < 0.05).sum()
        ),
        "n_target_excluded_localization_FDR_supported": int(
            (comparison["localization_FDR_target_excluded"] < 0.05).sum()
        ),
        "n_localization_FDR_support_retained": int(
            comparison["localization_FDR_support_retained"].sum()
        ),
    }
    (output / "circularity_sensitivity_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return summary


def write_review_bundle(summary: dict) -> Path:
    completion = {
        "status": "completed",
        "generated_at": now_iso(),
        "primary_root": str(PRIMARY_ROOT),
        "target_excluded_root": str(SENSITIVITY_ROOT),
        "hybrid_root": str(HYBRID_ROOT),
        "n_reference_factors": 71,
        "roster_sha256": common.LOCKED71_SHA256,
        "excluded_Tier1_genes": sorted(TIER1_EXCLUDED_GENES),
        "comparison": summary,
    }
    completion_path = HYBRID_ROOT / "00_audit" / "LOCKED71_CIRCULARITY_COMPLETE.json"
    completion_path.write_text(json.dumps(completion, indent=2, ensure_ascii=False), encoding="utf-8")
    archive = HYBRID_ROOT.parent / f"{HYBRID_ROOT.name}_review_bundle.tar.gz"
    include = [
        SENSITIVITY_ROOT / "00_audit",
        SENSITIVITY_ROOT / "01_reference" / "reference_audit.json",
        SENSITIVITY_ROOT / "01_reference" / "spatial_reference_state_label_map.tsv",
        SENSITIVITY_ROOT / "02_cell2location" / "mapped_h5ad_manifest.tsv",
        SENSITIVITY_ROOT / "02_cell2location" / "mapping_failures.tsv",
        SENSITIVITY_ROOT / "06_cell2location_model_quality_audit_v1",
        HYBRID_ROOT / "00_audit",
        HYBRID_ROOT / "05_all_public_Tier1_spatial_validation_v1" / "00_audit",
        HYBRID_ROOT / "05_all_public_Tier1_spatial_validation_v1" / "01_tables",
        HYBRID_ROOT / "05_all_public_Tier1_spatial_validation_v1" / "02_figures",
        HYBRID_ROOT / "06_circularity_comparison",
    ]
    with tarfile.open(archive, "w:gz") as tar:
        for path in include:
            if path.exists():
                root = SENSITIVITY_ROOT if SENSITIVITY_ROOT in path.parents or path == SENSITIVITY_ROOT else HYBRID_ROOT
                prefix = "target_excluded" if root == SENSITIVITY_ROOT else "hybrid"
                tar.add(path, arcname=str(Path(prefix) / path.relative_to(root)))
    return archive


def main() -> None:
    common.check_required_inputs()
    for root in (SENSITIVITY_ROOT, HYBRID_ROOT):
        if root.resolve() in {common.ORIGINAL_ROOT.resolve(), PRIMARY_ROOT.resolve()}:
            raise RuntimeError("Circularity output roots must not overwrite original or primary results")
        root.mkdir(parents=True, exist_ok=True)

    patched_stage02 = (
        SENSITIVITY_ROOT / "scripts" / "02_spatial_v5_locked71_Tier1_excluded.py"
    )
    common.patch_stage02_locked71(
        common.SOURCE_STAGE02,
        patched_stage02,
        excluded_genes=TIER1_EXCLUDED_GENES,
    )
    common.prepare_root(SENSITIVITY_ROOT, patched_stage02)
    env = common.build_environment(SENSITIVITY_ROOT)
    exclusion_contract = {
        "generated_at": now_iso(),
        "excluded_genes": sorted(TIER1_EXCLUDED_GENES),
        "n_excluded_genes": len(TIER1_EXCLUDED_GENES),
        "n_reference_factors": 71,
        "roster_sha256": common.LOCKED71_SHA256,
        "principle": (
            "All public Tier 1 genes were removed before reference signature estimation "
            "and spatial deconvolution. Gene-expression outcomes are taken from the primary "
            "locked-71 mapped objects; target-excluded abundances are used in the hybrid tests."
        ),
    }
    (SENSITIVITY_ROOT / "00_audit" / "circularity_exclusion_contract.json").write_text(
        json.dumps(exclusion_contract, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print("=" * 100)
    print("LOCKED-71 TIER 1 CIRCULARITY ANALYSIS")
    print("Target-excluded root:", SENSITIVITY_ROOT)
    print("Hybrid root:", HYBRID_ROOT)
    print("CUDA device:", env.get("CUDA_VISIBLE_DEVICES"))
    print("=" * 100)

    common.run_stream(
        [sys.executable, "-u", str(patched_stage02)],
        env,
        SENSITIVITY_ROOT / "logs" / "stage02_locked71_Tier1_excluded.log",
    )
    common.validate_locked71_output(SENSITIVITY_ROOT)

    # The target-excluded fit receives its own quality audit.
    common.run_optional_analysis(
        SENSITIVITY_ROOT / "scripts" / common.QC_SCRIPT.name,
        SENSITIVITY_ROOT,
        env,
        "stage06_locked71_Tier1_excluded_quality.log",
    )
    common.finalize_quality_audit_corrected(SENSITIVITY_ROOT)

    # This permits two-GPU parallel execution: mapping proceeds independently,
    # then waits here until the primary 71-state mapping and Tier 1 table exist.
    wait_for_primary()
    hybridise_primary_expression_with_target_excluded_abundance()

    shutil.copy2(common.TIER1_SCRIPT, HYBRID_ROOT / "scripts" / common.TIER1_SCRIPT.name)
    hybrid_env = common.build_environment(HYBRID_ROOT)
    common.run_optional_analysis(
        HYBRID_ROOT / "scripts" / common.TIER1_SCRIPT.name,
        HYBRID_ROOT,
        hybrid_env,
        "stage05_locked71_Tier1_excluded_hybrid.log",
    )
    summary = compare_primary_and_target_excluded()
    archive = write_review_bundle(summary)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print("Download review bundle:", archive)


if __name__ == "__main__":
    main()
