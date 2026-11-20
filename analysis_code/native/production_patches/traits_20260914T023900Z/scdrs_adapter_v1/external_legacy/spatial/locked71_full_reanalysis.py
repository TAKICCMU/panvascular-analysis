from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import time
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import pandas as pd


# =============================================================================
# Frozen analysis contract
# =============================================================================

PROJECT_ROOT = Path("/path/to/panvas_project")
ORIGINAL_ROOT = PROJECT_ROOT / "spatial_validation_v4_gpu_final"
SOURCE_ROOT = PROJECT_ROOT / "spatial_validation_v3"
OUTPUT_ROOT = Path(os.environ.get(
    "SPATIAL_LOCKED71_ROOT",
    str(PROJECT_ROOT / "spatial_validation_v5_locked71_gpu_final"),
))

SOURCE_STAGE02 = ORIGINAL_ROOT / "scripts" / "02_spatial_v4_cell2location_GPU.py"
SOURCE_STAGE03 = ORIGINAL_ROOT / "scripts" / "03_spatial_v4_targets_CoVarNet_validation.py"
SOURCE_STAGE04 = ORIGINAL_ROOT / "scripts" / "04_spatial_v4_publication_figures.py"
TIER1_SCRIPT = Path(os.environ.get(
    "TIER1_SPATIAL_ANALYSIS_SOURCE",
    str(PROJECT_ROOT / "06_all_Tier1_spatial_validation_Jupyter一键运行版.py"),
))
QC_SCRIPT = Path(os.environ.get(
    "CELL2LOCATION_QC_SOURCE",
    str(PROJECT_ROOT / "10_cell2location_model_quality_audit_Jupyter一键运行版.py"),
))

ATLAS_PATH = PROJECT_ROOT / "Level2_celltype_final_v1_metadata_v2.h5ad"
MANIFEST_PATH = SOURCE_ROOT / "00_inputs" / "spatial_sample_manifest_curated.tsv"
LOCKED71_FIELD = "Level2_celltype_for_enrichment_v1"

LOCKED71_STATES = (
    "B cell | IGHD IGHM naive B cell",
    "B cell | immediate-early activated B cell",
    "B cell | memory B cell",
    "Dendritic cell | CD1C cDC2",
    "Dendritic cell | CLEC9A XCR1 cDC1",
    "Dendritic cell | GZMB LILRA4 plasmacytoid dendritic cell",
    "Dendritic cell | OLR1 CCL20 inflammatory cDC2",
    "Endothelial cell | ACKR1 venous endothelial cell",
    "Endothelial cell | CA4 capillary endothelial cell",
    "Endothelial cell | CCL21 lymphatic endothelial cell",
    "Endothelial cell | FCN3 organ-specific capillary endothelial cell",
    "Endothelial cell | GATA4 endocardial-valvular endothelial cell",
    "Endothelial cell | PTGIS stalk-matrix endothelial cell",
    "Endothelial cell | SEMA3G arterial endothelial cell",
    "Endothelial cell | SOD3 mesenchymal-like endothelial cell",
    "Endothelial cell | TSHZ2 quiescent endothelial cell",
    "Fibroblast | C1R ECM-remodelling fibroblast",
    "Fibroblast | CXCL12 perivascular fibroblast",
    "Fibroblast | FBLN1 homeostatic adventitial fibroblast",
    "Fibroblast | ITGA8 SMC-like fibroblast",
    "Fibroblast | LAMA2 tissue-associated fibroblast",
    "Fibroblast | PLA2G2A secretory fibroblast",
    "Fibroblast | POSTN myofibroblast-like fibroblast",
    "Fibroblast | RBMS3 progenitor-like fibroblast",
    "Fibromyocyte | AEBP1 ECM-remodelling fibromyocyte",
    "Fibromyocyte | MYH11 contractile fibromyocyte",
    "Fibromyocyte | S100B neural-crest-like fibromyocyte",
    "Fibromyocyte | SCARA5 fibroblast-like fibromyocyte",
    "Macrophage | CCL3 inflammatory macrophage",
    "Macrophage | FCN1 monocyte-derived macrophage",
    "Macrophage | HLA-DQ antigen-presenting macrophage",
    "Macrophage | LYVE1 FOLR2 resident macrophage",
    "Macrophage | SPP1 lipid-associated macrophage",
    "Macrophage | heat-shock antigen-presenting macrophage",
    "Macrophage | resident-like macrophage",
    "Mast cell | NF-kB activated mast cell",
    "Mast cell | TPSB2 canonical mast cell",
    "Mast cell | quiescent mast cell",
    "Monocyte | C1Q macrophage-like monocyte",
    "Monocyte | FCGR3A non-classical monocyte",
    "Monocyte | FCN1 classical monocyte",
    "Monocyte | S100A8 S100A12 classical monocyte",
    "Monocyte | regulatory monocyte state",
    "Natural killer cell | CCL3 CCL4 activated NK cell",
    "Natural killer cell | FGFBP2 FCGR3A cytotoxic NK cell",
    "Natural killer cell | GZMK NK cell",
    "Natural killer cell | IL7R tissue-resident NK cell",
    "Natural killer cell | proliferating NK cell",
    "Natural killer cell | regulatory NK cell",
    "Neutrophil | S100A8 S100A9 activated neutrophil",
    "Neutrophil | VCAN RETN immature neutrophil",
    "Neutrophil | mature neutrophil",
    "Pericyte | FGF7 contractile pericyte",
    "Pericyte | PDGFRA matrix-associated pericyte",
    "Pericyte | PRKG1 contractile pericyte",
    "Pericyte | RGS5 canonical pericyte",
    "Plasma cell | IgA plasma cell",
    "Plasma cell | IgG plasma cell",
    "Plasma cell | IgG2 FBP1 plasma cell",
    "Smooth muscle cell | COX4I2 pericyte-like SMC",
    "Smooth muscle cell | MYH11 contractile SMC",
    "Smooth muscle cell | PRKG1 mature contractile SMC",
    "Smooth muscle cell | TFF1 valve-associated SMC",
    "Smooth muscle cell | TSHZ2 fibroblast-like SMC",
    "T cell | FOXP3 regulatory T cell",
    "T cell | GZMK cytotoxic-memory CD8 T cell",
    "T cell | IL7R KLRB1 memory CD4 T cell",
    "T cell | NK-like cytotoxic T cell",
    "T cell | ZNF683 tissue-resident cytotoxic T cell",
    "T cell | conventional naive-like T cell",
    "T cell | proliferating T cell",
)

LOCKED71_SHA256 = hashlib.sha256(
    ("\n".join(sorted(LOCKED71_STATES)) + "\n").encode("utf-8")
).hexdigest()
EXPECTED_LOCKED71_SHA256 = "4f58fae84263c600d785a60e4c6231d1b1fd9e2818570f15ef740795dc1f70cf"
if len(LOCKED71_STATES) != 71 or LOCKED71_SHA256 != EXPECTED_LOCKED71_SHA256:
    raise RuntimeError("The embedded locked-71 roster failed its count or checksum contract")


# =============================================================================
# Utilities shared with the circularity runner
# =============================================================================

def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def check_required_inputs() -> None:
    required = [
        SOURCE_STAGE02, SOURCE_STAGE03, SOURCE_STAGE04, TIER1_SCRIPT,
        QC_SCRIPT, ATLAS_PATH, MANIFEST_PATH,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n" + "\n".join(missing))


def expected_gses() -> list[str]:
    manifest = pd.read_csv(MANIFEST_PATH, sep="\t")
    include = manifest["include_mapping"].astype(str).str.lower().isin({"true", "1", "yes"})
    return sorted(manifest.loc[include, "GSE"].astype(str).unique().tolist())


def run_stream(command: list[str], env: dict[str, str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    print("$", " ".join(map(str, command)), flush=True)
    print("Log:", log_path, flush=True)
    with log_path.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            command,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="")
            log.write(line)
            log.flush()
        return_code = process.wait()
    if return_code != 0:
        tail = log_path.read_text(errors="replace").splitlines()[-200:]
        print("\n".join(tail))
        raise RuntimeError(f"Command failed with return code {return_code}: {command}")


def replace_one(source: str, pattern: str, replacement: str, description: str) -> str:
    updated, count = re.subn(pattern, replacement, source, count=1, flags=re.DOTALL)
    if count != 1:
        raise RuntimeError(f"Could not patch exactly one {description}; matched {count}")
    return updated


def patch_stage02_locked71(
    source_path: Path,
    target_path: Path,
    excluded_genes: set[str] | None = None,
) -> None:
    """Create an auditable Stage-02 script with exactly 71 separate factors."""
    source = source_path.read_text(encoding="utf-8")

    source = source.replace(
        "USE_GPU = True",
        'USE_GPU = os.environ.get("SPATIAL_USE_GPU", "1") == "1"',
        1,
    )

    roster_literal = repr(tuple(LOCKED71_STATES))
    state_selector = f'''LOCKED71_STATES = {roster_literal}\nLOCKED71_ROSTER_SHA256 = "{LOCKED71_SHA256}"\nLOCKED71_STATE_FIELD = "{LOCKED71_FIELD}"\n\ndef select_state_key(obs: pd.DataFrame) -> str:\n    if LOCKED71_STATE_FIELD not in obs.columns:\n        raise KeyError(\n            f"Required locked-71 field {{LOCKED71_STATE_FIELD!r}} is absent; "\n            f"available fields include: {{list(obs.columns)}}"\n        )\n    return LOCKED71_STATE_FIELD\n\n\n'''
    source = replace_one(
        source,
        r"def select_state_key\(obs: pd\.DataFrame\) -> str:.*?(?=def select_level1_key)",
        state_selector,
        "state selector",
    )

    label_builder = '''def build_spatial_reference_labels(obs: pd.DataFrame) -> tuple[pd.Series, dict]:
    state_key = select_state_key(obs)
    level1_key = select_level1_key(obs)
    state = obs[state_key].astype("object").astype(str).str.strip()
    locked_set = set(LOCKED71_STATES)
    observed_locked = set(state[state.isin(locked_set)].unique())
    missing = sorted(locked_set.difference(observed_locked))
    if missing:
        raise ValueError(f"Locked-71 states absent from the atlas: {missing}")
    labels = state.where(state.isin(locked_set), pd.NA).astype("object")
    return labels, {
        "state_key": state_key,
        "level1_key": level1_key,
        "locked71_roster_sha256": LOCKED71_ROSTER_SHA256,
        "n_locked71_states_requested": len(LOCKED71_STATES),
        "n_locked71_states_observed": len(observed_locked),
        "non_roster_cells_excluded_from_reference": int(labels.isna().sum()),
        "reference_definition": "71 locked Level 2 states fitted as separate factors; non-roster cells excluded",
    }


'''
    source = replace_one(
        source,
        r"def build_spatial_reference_labels\(obs: pd\.DataFrame\).*?(?=def balanced_indices)",
        label_builder,
        "locked-71 label builder",
    )

    assertion = '''    observed_original = set(label_map["spatial_reference_state_original"].astype(str))
    expected_original = set(LOCKED71_STATES)
    if len(label_map) != 71 or observed_original != expected_original:
        raise RuntimeError(
            "Reference label contract failed: "
            f"n={len(label_map)}, missing={sorted(expected_original - observed_original)}, "
            f"unexpected={sorted(observed_original - expected_original)}"
        )
'''
    needle = '    label_map.to_csv(REF_OUT / "spatial_reference_state_label_map.tsv", sep="\\t", index=False)'
    if needle not in source:
        raise RuntimeError("Could not locate the reference label-map write statement")
    source = source.replace(needle, assertion + needle, 1)

    if excluded_genes:
        excluded = {str(gene).upper() for gene in excluded_genes}
        needle = "    selected_genes.update(gene for gene in FORCED_GENES if gene in ref.var_names)"
        if needle not in source:
            raise RuntimeError("Could not locate the selected-gene update statement")
        injection = needle + "\n" + (
            f"    circularity_excluded_genes = {repr(excluded)}\n"
            "    selected_genes.difference_update(circularity_excluded_genes)\n"
            "    if selected_genes.intersection(circularity_excluded_genes):\n"
            "        raise RuntimeError('Tier 1 genes leaked into the target-excluded reference')"
        )
        source = source.replace(needle, injection, 1)

    target_path.parent.mkdir(parents=True, exist_ok=True)
    target_path.write_text(source, encoding="utf-8")
    compile(source, str(target_path), "exec")


def resolve_module_weights() -> Path:
    candidates = [
        SOURCE_ROOT / "00_inputs" / "module_node_weights_K14.tsv",
        PROJECT_ROOT / "CoVarNet_final_v3" / "stable_NMF" / "module_node_weights_K14.tsv",
        PROJECT_ROOT / "covarnet" / "stable_NMF" / "module_node_weights_K14.tsv",
    ]
    path = next((item for item in candidates if item.exists()), None)
    if path is None:
        raise FileNotFoundError("No CoVarNet module weight file found:\n" + "\n".join(map(str, candidates)))
    return path


def prepare_root(root: Path, stage02_path: Path) -> None:
    if root.resolve() == ORIGINAL_ROOT.resolve():
        raise RuntimeError("The locked-71 run must not overwrite the original 63-factor result")
    for subdir in (
        "00_audit", "00_inputs", "01_reference", "02_cell2location",
        "scripts", "logs",
    ):
        (root / subdir).mkdir(parents=True, exist_ok=True)

    roster = pd.DataFrame({
        "locked71_order": range(1, 72),
        "Level2": list(LOCKED71_STATES),
    })
    roster.to_csv(root / "00_inputs" / "locked71_spatial_reference_roster.tsv", sep="\t", index=False)
    (root / "00_audit" / "locked71_reference_contract.json").write_text(
        json.dumps({
            "generated_at": now_iso(),
            "state_field": LOCKED71_FIELD,
            "n_states": 71,
            "roster_sha256": LOCKED71_SHA256,
            "non_roster_handling": "excluded from reference; no pooled residual factor",
            "source_stage02": str(SOURCE_STAGE02),
            "patched_stage02": str(stage02_path),
        }, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    shutil.copy2(SOURCE_STAGE03, root / "scripts" / SOURCE_STAGE03.name)
    shutil.copy2(SOURCE_STAGE04, root / "scripts" / SOURCE_STAGE04.name)
    shutil.copy2(TIER1_SCRIPT, root / "scripts" / TIER1_SCRIPT.name)
    shutil.copy2(QC_SCRIPT, root / "scripts" / QC_SCRIPT.name)
    shutil.copy2(resolve_module_weights(), root / "00_inputs" / "module_node_weights_K14.tsv")


def build_environment(root: Path) -> dict[str, str]:
    env = os.environ.copy()
    cuda_device = env.get("SPATIAL_CUDA_DEVICE", "0")
    env.update({
        "PYTHONUNBUFFERED": "1",
        "MPLBACKEND": "Agg",
        "CUDA_VISIBLE_DEVICES": cuda_device,
        "PYTORCH_CUDA_ALLOC_CONF": "expandable_segments:True,max_split_size_mb:256",
        "SPATIAL_USE_GPU": env.get("SPATIAL_USE_GPU", "1"),
        "SPATIAL_V4_ROOT": str(root),
        "SPATIAL_V3_SOURCE_ROOT": str(SOURCE_ROOT),
        "SPATIAL_REF_EPOCHS": env.get("SPATIAL_REF_EPOCHS", "250"),
        "SPATIAL_MAPPING_EPOCHS": env.get("SPATIAL_MAPPING_EPOCHS", "30000"),
        "SPATIAL_POSTERIOR_SAMPLES": env.get("SPATIAL_POSTERIOR_SAMPLES", "1000"),
        "SPATIAL_RANDOM_SEED": env.get("SPATIAL_RANDOM_SEED", "20260814"),
        "SPATIAL_MIN_FREE_GPU_GB": env.get("SPATIAL_MIN_FREE_GPU_GB", "10"),
    })
    return env


def factor_names(adata: ad.AnnData, key: str) -> list[str]:
    value = adata.obsm[key]
    # cell2location stores canonical factor names in uns['mod']['factor_names'].
    # DataFrame columns written by export_posterior carry a statistic-specific
    # prefix (for example q05cell_abundance_w_sf_), so they are not the factor
    # identifiers and must not be compared directly with the locked roster.
    names = adata.uns.get("mod", {}).get("factor_names")
    if names is not None and len(names) == value.shape[1]:
        return list(map(str, names))

    if isinstance(value, pd.DataFrame):
        columns = list(map(str, value.columns))
        prefix = key.replace("_cell", "cell", 1) + "_"
        if all(column.startswith(prefix) for column in columns):
            stripped = [column[len(prefix):] for column in columns]
            if len(set(stripped)) == len(stripped):
                return stripped
        return columns

    raise ValueError(f"Cannot recover factor names for {key}")


def validate_locked71_output(root: Path) -> pd.DataFrame:
    audit_path = root / "01_reference" / "reference_audit.json"
    labels_path = root / "01_reference" / "spatial_reference_state_label_map.tsv"
    signatures_path = root / "01_reference" / "final_atlas_reference_signatures.tsv.gz"
    for path in (audit_path, labels_path, signatures_path):
        if not path.exists():
            raise FileNotFoundError(path)

    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    labels = pd.read_csv(labels_path, sep="\t")
    observed = set(labels["spatial_reference_state_original"].astype(str))
    if audit.get("state_key") != LOCKED71_FIELD:
        raise RuntimeError(f"Wrong reference field: {audit.get('state_key')}")
    if int(audit.get("n_states", -1)) != 71:
        raise RuntimeError(f"Reference contains {audit.get('n_states')} factors, expected 71")
    if audit.get("locked71_roster_sha256") != LOCKED71_SHA256:
        raise RuntimeError("Reference roster checksum does not match the frozen 71-state roster")
    if observed != set(LOCKED71_STATES):
        raise RuntimeError("Reference label map does not match the frozen 71-state roster")

    signatures = pd.read_csv(signatures_path, sep="\t", index_col=0, nrows=2)
    if signatures.shape[1] != 71:
        raise RuntimeError(f"Signature matrix has {signatures.shape[1]} factors, expected 71")

    safe_expected = set(labels["spatial_reference_state_safe"].astype(str))
    rows = []
    mapped_paths = sorted((root / "02_cell2location").glob("*__final_atlas_cell2location.h5ad"))
    expected = expected_gses()
    observed_gses = [path.name.split("__", 1)[0] for path in mapped_paths]
    if observed_gses != expected:
        raise RuntimeError(f"Mapped GSE contract failed: expected={expected}, observed={observed_gses}")
    for path in mapped_paths:
        mapped = ad.read_h5ad(path, backed="r")
        record = {"GSE": path.name.split("__", 1)[0], "path": str(path), "n_spots": int(mapped.n_obs)}
        for key in ("q05_cell_abundance_w_sf", "means_cell_abundance_w_sf"):
            if key not in mapped.obsm:
                raise KeyError(f"{key} absent from {path}")
            names = factor_names(mapped, key)
            record[f"{key}__n_factors"] = len(names)
            record[f"{key}__factor_set_matches"] = set(names) == safe_expected
            if len(names) != 71 or set(names) != safe_expected:
                raise RuntimeError(f"{path.name} {key} does not contain the exact 71 factors")
        mapped.file.close()
        rows.append(record)
    result = pd.DataFrame(rows)
    result.to_csv(root / "00_audit" / "locked71_mapped_output_contract.tsv", sep="\t", index=False)
    return result


def run_optional_analysis(script: Path, root: Path, env: dict[str, str], log_name: str) -> None:
    if not script.exists():
        raise FileNotFoundError(script)
    run_stream([sys.executable, "-u", str(script)], env, root / "logs" / log_name)


def finalize_quality_audit_corrected(root: Path) -> dict:
    """Apply the valid gates without treating signature-only deviance as model fit."""
    output = root / "06_cell2location_model_quality_audit_v1"
    history_path = output / "01_training_history_summary.tsv"
    quality_path = output / "02_dataset_model_quality_summary.tsv"
    if not history_path.exists() or not quality_path.exists():
        raise FileNotFoundError(f"Incomplete quality audit under {output}")
    history = pd.read_csv(history_path, sep="\t")
    quality = pd.read_csv(quality_path, sep="\t")
    quality = quality.rename(columns={
        "posterior_mean_reconstruction_deviance_explained":
            "signature_only_deviance_relative_to_global_null",
    })
    quality["gate_signature_gene_rank_concordance"] = (
        pd.to_numeric(quality["posterior_mean_reconstruction_gene_rho"], errors="coerce") >= 0.50
    )
    quality["all_valid_numeric_gates_pass"] = (
        quality["gate_finite_nonnegative"].astype(bool)
        & quality["gate_posterior_order"].astype(bool)
        & quality["gate_q05_mean_stability"].astype(bool)
        & quality["gate_signature_gene_rank_concordance"].astype(bool)
    )
    quality["signature_only_deviance_is_gate"] = False
    quality.to_csv(
        output / "02_dataset_model_quality_summary_corrected.tsv", sep="\t", index=False
    )

    spatial_history = history[history["dataset"] != "reference_regression"].copy()
    spatial_history["tail_stability_gate"] = (
        spatial_history["history_usable"].astype(bool)
        & (pd.to_numeric(spatial_history["finite_fraction"], errors="coerce") == 1)
        & (pd.to_numeric(spatial_history["tail_CV"], errors="coerce") < 0.001)
        & (pd.to_numeric(spatial_history["tail_normalized_slope"], errors="coerce").abs() < 1e-5)
    )
    spatial_history.to_csv(
        output / "01_spatial_training_history_summary_corrected.tsv", sep="\t", index=False
    )
    manifest = {
        "status": "completed",
        "generated_at": now_iso(),
        "n_datasets": int(len(quality)),
        "n_valid_numeric_gates_pass": int(quality["all_valid_numeric_gates_pass"].sum()),
        "n_training_tail_stability_pass": int(spatial_history["tail_stability_gate"].sum()),
        "valid_quality_gates": {
            "finite_nonnegative": "all posterior values finite and nonnegative",
            "posterior_order": "posterior mean >= posterior q05 for >=99.9% of entries",
            "posterior_stability": "median factor-wise Spearman(q05, mean) >=0.90",
            "signature_gene_rank": "observed versus signature-only predicted gene totals rho >=0.50",
            "training_tail": "tail CV <0.001 and absolute normalized slope <1e-5",
        },
        "important_correction": (
            "Signature-only reconstruction omits cell2location background and technical terms. "
            "Its deviance relative to a global gene-frequency null is descriptive and is not a gate."
        ),
    }
    (output / "run_manifest_corrected.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return manifest


def write_completion_bundle(root: Path, completion: dict) -> Path:
    completion_path = root / "00_audit" / "LOCKED71_ANALYSIS_COMPLETE.json"
    completion_path.write_text(json.dumps(completion, indent=2, ensure_ascii=False), encoding="utf-8")
    archive = root.parent / f"{root.name}_review_bundle.tar.gz"
    include = [
        root / "00_audit",
        root / "01_reference" / "reference_audit.json",
        root / "01_reference" / "spatial_reference_state_label_map.tsv",
        root / "02_cell2location" / "mapped_h5ad_manifest.tsv",
        root / "02_cell2location" / "mapping_failures.tsv",
        root / "03_validation_targets_CoVarNet" / "tables",
        root / "05_all_public_Tier1_spatial_validation_v1" / "00_audit",
        root / "05_all_public_Tier1_spatial_validation_v1" / "01_tables",
        root / "05_all_public_Tier1_spatial_validation_v1" / "02_figures",
        root / "06_cell2location_model_quality_audit_v1",
    ]
    with tarfile.open(archive, "w:gz") as tar:
        for path in include:
            if path.exists():
                tar.add(path, arcname=str(path.relative_to(root)))
    return archive


# =============================================================================
# Primary locked-71 analysis
# =============================================================================

def main() -> None:
    check_required_inputs()
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    patched_stage02 = OUTPUT_ROOT / "scripts" / "02_spatial_v5_cell2location_locked71.py"
    patch_stage02_locked71(SOURCE_STAGE02, patched_stage02)
    prepare_root(OUTPUT_ROOT, patched_stage02)
    env = build_environment(OUTPUT_ROOT)

    print("=" * 100)
    print("LOCKED-71 PRIMARY SPATIAL ANALYSIS")
    print("Output root:", OUTPUT_ROOT)
    print("CUDA device:", env.get("CUDA_VISIBLE_DEVICES"))
    print("Roster SHA256:", LOCKED71_SHA256)
    print("=" * 100)

    run_stream(
        [sys.executable, "-u", str(patched_stage02)],
        env,
        OUTPUT_ROOT / "logs" / "stage02_locked71_cell2location.log",
    )
    contract = validate_locked71_output(OUTPUT_ROOT)

    # Recompute module and target summaries from the new 71-factor posteriors.
    run_optional_analysis(
        OUTPUT_ROOT / "scripts" / SOURCE_STAGE03.name,
        OUTPUT_ROOT,
        env,
        "stage03_locked71_targets_CoVarNet.log",
    )

    # The prespecified public 22-pair Tier 1 analysis is the manuscript-facing target screen.
    run_optional_analysis(
        OUTPUT_ROOT / "scripts" / TIER1_SCRIPT.name,
        OUTPUT_ROOT,
        env,
        "stage05_locked71_all_public_Tier1.log",
    )

    # Re-run posterior/training quality checks for the new factor basis.
    run_optional_analysis(
        OUTPUT_ROOT / "scripts" / QC_SCRIPT.name,
        OUTPUT_ROOT,
        env,
        "stage06_locked71_cell2location_quality.log",
    )
    corrected_quality = finalize_quality_audit_corrected(OUTPUT_ROOT)

    # Legacy Figure-5 code contains historical target choices. It is disabled by default;
    # final figures should be rebuilt after reviewing the new 71-state results.
    if os.environ.get("SPATIAL_RUN_LEGACY_FIGURES", "0") == "1":
        run_optional_analysis(
            OUTPUT_ROOT / "scripts" / SOURCE_STAGE04.name,
            OUTPUT_ROOT,
            env,
            "stage04_locked71_legacy_figures.log",
        )

    completion = {
        "status": "completed",
        "generated_at": now_iso(),
        "root": str(OUTPUT_ROOT),
        "reference_field": LOCKED71_FIELD,
        "n_reference_factors": 71,
        "roster_sha256": LOCKED71_SHA256,
        "n_mapped_datasets": int(len(contract)),
        "mapped_datasets": contract["GSE"].tolist(),
        "reference_epochs": int(env["SPATIAL_REF_EPOCHS"]),
        "mapping_epochs": int(env["SPATIAL_MAPPING_EPOCHS"]),
        "posterior_samples": int(env["SPATIAL_POSTERIOR_SAMPLES"]),
        "random_seed": int(env["SPATIAL_RANDOM_SEED"]),
        "legacy_figure_stage_run": os.environ.get("SPATIAL_RUN_LEGACY_FIGURES", "0") == "1",
        "corrected_quality_audit": corrected_quality,
        "interpretation": (
            "All 71 frozen Level 2 states were fitted as separate cell2location factors. "
            "Cells outside the frozen roster were excluded from reference signature estimation."
        ),
    }
    archive = write_completion_bundle(OUTPUT_ROOT, completion)
    print(json.dumps(completion, indent=2, ensure_ascii=False))
    print("Download review bundle:", archive)


if __name__ == "__main__":
    main()
