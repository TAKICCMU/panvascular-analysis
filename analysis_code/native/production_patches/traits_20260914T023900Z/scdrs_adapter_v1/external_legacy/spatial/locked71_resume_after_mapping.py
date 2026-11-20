from __future__ import annotations

import importlib.util
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import pandas as pd


# This script deliberately does not execute Stage 02. It resumes only after the
# seven locked-71 cell2location H5AD files have already been written.
DRIVER = Path(
    "/path/to/panvas_project/"
    "16_spatial_locked71_full_reanalysis_Jupyter一键运行版.py"
)


def load_driver(path: Path):
    if not path.exists():
        raise FileNotFoundError(path)
    spec = importlib.util.spec_from_file_location("locked71_primary_driver", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import driver: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def factor_audit(module, root: Path) -> pd.DataFrame:
    labels_path = root / "01_reference" / "spatial_reference_state_label_map.tsv"
    if not labels_path.exists():
        raise FileNotFoundError(labels_path)
    labels = pd.read_csv(labels_path, sep="\t")
    expected = set(labels["spatial_reference_state_safe"].astype(str))
    rows = []
    mapped_paths = sorted(
        (root / "02_cell2location").glob("*__final_atlas_cell2location.h5ad")
    )
    for path in mapped_paths:
        mapped = ad.read_h5ad(path, backed="r")
        try:
            for key in (
                "q05_cell_abundance_w_sf",
                "means_cell_abundance_w_sf",
            ):
                record = {
                    "GSE": path.name.split("__", 1)[0],
                    "path": str(path),
                    "key": key,
                    "n_spots": int(mapped.n_obs),
                    "key_present": key in mapped.obsm,
                }
                if key in mapped.obsm:
                    names = module.factor_names(mapped, key)
                    observed = set(map(str, names))
                    record.update(
                        {
                            "n_factors": len(names),
                            "factor_set_matches": observed == expected,
                            "missing_factors": " | ".join(sorted(expected - observed)),
                            "unexpected_factors": " | ".join(sorted(observed - expected)),
                        }
                    )
                rows.append(record)
        finally:
            mapped.file.close()
    return pd.DataFrame(rows)


def main() -> None:
    module = load_driver(DRIVER)
    root = module.OUTPUT_ROOT
    audit_dir = root / "00_audit"
    audit_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 100)
    print("LOCKED-71 RESUME AFTER COMPLETED CELL2LOCATION MAPPING")
    print("Driver:", DRIVER)
    print("Output root:", root)
    print("Stage 02 will NOT be run by this script.")
    print("=" * 100)

    mapped_paths = sorted(
        (root / "02_cell2location").glob("*__final_atlas_cell2location.h5ad")
    )
    expected_gses = module.expected_gses()
    observed_gses = [path.name.split("__", 1)[0] for path in mapped_paths]
    print("Expected GSEs:", expected_gses)
    print("Observed GSEs:", observed_gses)
    if observed_gses != expected_gses:
        raise RuntimeError(
            f"Seven-dataset mapping contract failed: expected={expected_gses}, "
            f"observed={observed_gses}"
        )

    detailed = factor_audit(module, root)
    detailed_path = audit_dir / "locked71_factor_contract_preflight.tsv"
    detailed.to_csv(detailed_path, sep="\t", index=False)
    print("Factor preflight:")
    print(
        detailed[
            ["GSE", "key", "key_present", "n_factors", "factor_set_matches"]
        ].to_string(index=False)
    )
    failed = detailed[
        (~detailed["key_present"].fillna(False))
        | (~detailed["factor_set_matches"].fillna(False))
    ]
    if not failed.empty:
        print("\nFAILED FACTOR CONTRACT")
        print(
            failed[
                [
                    "GSE",
                    "key",
                    "n_factors",
                    "missing_factors",
                    "unexpected_factors",
                ]
            ].to_string(index=False)
        )
        raise RuntimeError(
            "Mapped files exist, but at least one posterior does not expose the "
            "exact frozen 71-factor set. Send locked71_factor_contract_preflight.tsv."
        )

    # This writes locked71_mapped_output_contract.tsv on success.
    contract = module.validate_locked71_output(root)
    print("\nLocked-71 output contract PASSED for", len(contract), "datasets")

    env = module.build_environment(root)

    print("\n[Stage 03] Recomputing targets and CoVarNet summaries")
    module.run_optional_analysis(
        root / "scripts" / module.SOURCE_STAGE03.name,
        root,
        env,
        "stage03_locked71_targets_CoVarNet.log",
    )

    print("\n[Stage 05] Running all public Tier 1 spatial validation")
    module.run_optional_analysis(
        root / "scripts" / module.TIER1_SCRIPT.name,
        root,
        env,
        "stage05_locked71_all_public_Tier1.log",
    )

    print("\n[Stage 06] Running cell2location model-quality audit")
    module.run_optional_analysis(
        root / "scripts" / module.QC_SCRIPT.name,
        root,
        env,
        "stage06_locked71_cell2location_quality.log",
    )
    corrected_quality = module.finalize_quality_audit_corrected(root)

    completion = {
        "status": "completed",
        "generated_at": now_iso(),
        "root": str(root),
        "reference_field": module.LOCKED71_FIELD,
        "n_reference_factors": 71,
        "roster_sha256": module.LOCKED71_SHA256,
        "n_mapped_datasets": int(len(contract)),
        "mapped_datasets": contract["GSE"].tolist(),
        "reference_epochs": int(env["SPATIAL_REF_EPOCHS"]),
        "mapping_epochs": int(env["SPATIAL_MAPPING_EPOCHS"]),
        "posterior_samples": int(env["SPATIAL_POSTERIOR_SAMPLES"]),
        "random_seed": int(env["SPATIAL_RANDOM_SEED"]),
        "legacy_figure_stage_run": False,
        "corrected_quality_audit": corrected_quality,
        "resumed_after_completed_stage02": True,
        "interpretation": (
            "All 71 frozen Level 2 states were fitted as separate "
            "cell2location factors. Cells outside the frozen roster were "
            "excluded from reference signature estimation."
        ),
    }
    archive = module.write_completion_bundle(root, completion)
    resume_manifest = audit_dir / "LOCKED71_RESUME_COMPLETE.json"
    resume_manifest.write_text(
        json.dumps(completion, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print("\n", json.dumps(completion, indent=2, ensure_ascii=False), sep="")
    print("Download review bundle:", archive)


if __name__ == "__main__":
    main()
