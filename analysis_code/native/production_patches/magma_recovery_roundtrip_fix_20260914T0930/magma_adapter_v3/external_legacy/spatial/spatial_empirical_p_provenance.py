from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


SUBMISSION = Path(
    "/path/to/desktop/Panvas_Nature Version 1.0/"
    "泛血管最终投稿版与展示/Codex_final_submission_20260826"
)
SPATIAL = SUBMISSION / "Source_Data" / "Spatial_Validation_Final"
OUT = Path(__file__).resolve().parent / "03_spatial_empirical_p_audit"
OUT.mkdir(parents=True, exist_ok=True)


def main() -> None:
    by_section = pd.read_csv(
        SPATIAL / "CoVarNet_all14_spatial_validation_by_section.tsv", sep="\t"
    )
    summary = pd.read_csv(
        SPATIAL / "CoVarNet_all14_spatial_validation_summary.tsv", sep="\t"
    )
    audit_rows = []
    for p_col in ["within_spot_empirical_p", "spatial_lag_empirical_p"]:
        values = pd.to_numeric(by_section[p_col], errors="coerce").dropna()
        scaled = values * 1001.0
        audit_rows.append(
            {
                "quantity": p_col,
                "level": "section",
                "statistical_source": "1,000 lineage-matched random modules within each section",
                "n_values": len(values),
                "minimum": values.min(),
                "maximum": values.max(),
                "expected_minimum_resolution": 1 / 1001,
                "fraction_on_1_over_1001_grid": float(np.mean(np.isclose(scaled, np.round(scaled), atol=1e-8))),
                "may_be_smaller_than_1_over_1001": False,
            }
        )

    for p_col in ["cross_dataset_p_within_spot", "cross_dataset_p_spatial_lag"]:
        values = pd.to_numeric(summary[p_col], errors="coerce").dropna()
        audit_rows.append(
            {
                "quantity": p_col,
                "level": "cross-dataset",
                "statistical_source": "parametric signed Stouffer combination of section-level null z statistics",
                "n_values": len(values),
                "minimum": values.min(),
                "maximum": values.max(),
                "expected_minimum_resolution": np.nan,
                "fraction_on_1_over_1001_grid": np.nan,
                "may_be_smaller_than_1_over_1001": True,
            }
        )
    audit = pd.DataFrame(audit_rows)
    audit.to_csv(OUT / "spatial_P_value_provenance_audit.tsv", sep="\t", index=False)

    priority = summary[summary["module"].isin(["CM01", "CM02"])].copy()
    keep = [
        "module", "n_sections", "n_datasets",
        "median_within_spot_coherence", "fraction_within_spot_empirical_p_lt_0.05",
        "median_spatial_lag_coherence", "fraction_spatial_lag_empirical_p_lt_0.05",
        "median_moran_I", "cross_dataset_stouffer_z_within_spot",
        "cross_dataset_p_within_spot", "cross_dataset_FDR_within_spot",
        "cross_dataset_stouffer_z_spatial_lag", "cross_dataset_p_spatial_lag",
        "cross_dataset_FDR_spatial_lag",
    ]
    priority[keep].to_csv(OUT / "CM01_CM02_parametric_vs_empirical_summary.tsv", sep="\t", index=False)

    model_table = pd.read_csv(SPATIAL / "Tier1_spatial_models_by_section.tsv.gz", sep="\t")
    model_audit = (
        model_table.groupby("model", as_index=False)
        .agg(n_rows=("pair_id", "size"), n_sections=("section_id", "nunique"), n_pairs=("pair_id", "nunique"))
    )
    model_audit["residual_spatial_autocorrelation_available_in_frozen_summary"] = False
    model_audit["negative_binomial_sensitivity_available_in_frozen_summary"] = model_audit["model"].str.contains(
        "negative", case=False, na=False
    )
    model_audit.to_csv(OUT / "Tier1_spatial_model_family_audit.tsv", sep="\t", index=False)

    interpretation = """SPATIAL P-VALUE PROVENANCE AUDIT

1. The per-section empirical P values are generated from 1,000 lineage-matched random
   modules. Their minimum attainable value is therefore 1/1001 = 0.000999.
2. The much smaller cross-dataset P values and FDR values are not empirical random-module
   P values. They are parametric signed Stouffer combinations of section-level null z
   statistics and must be labelled separately.
3. The frozen source tables contain Moran's I for the observed module score, not Moran's I
   for residuals from each target-gene localization model.
4. The target-gene model table contains Poisson log-link models. A negative-binomial or
   quasi-Poisson sensitivity and residual Moran audit require the per-spot expression,
   coordinates and cell2location posterior stored in all seven mapped h5ad files.
5. Expression- and detection-matched random-gene controls likewise cannot be reconstructed
   from section-level coefficient summaries.
"""
    (OUT / "interpretation.txt").write_text(interpretation, encoding="utf-8")
    manifest = {
        "status": "completed_local_component",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "local_components_completed": [
            "empirical P-value resolution audit",
            "parametric versus empirical P-value provenance separation",
            "frozen model-family audit",
        ],
        "server_components_required": [
            "target-model residual Moran's I",
            "negative-binomial or quasi-Poisson sensitivity",
            "expression- and detection-matched random-gene controls",
        ],
    }
    (OUT / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(audit.to_string(index=False))
    print("\n", model_audit.to_string(index=False))
    print("\nOutputs:", OUT)


if __name__ == "__main__":
    main()
