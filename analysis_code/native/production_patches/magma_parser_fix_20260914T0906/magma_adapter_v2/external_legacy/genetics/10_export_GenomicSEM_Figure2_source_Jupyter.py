from __future__ import annotations

import gzip
import hashlib
import itertools
import json
import math
import re
import shutil
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


# =============================================================================
# Figure 2 source-data export for the pan-vascular Genomic SEM analysis
#
# This script is intentionally read-only with respect to all analysis folders.
# It creates a new export directory and copies/derives only plotting tables.
# Run it in one Jupyter cell with:
#
# exec(open(
#     "/path/to/panvas_analysis/Genomic_SEM_complete_followup/"
#     "10_export_GenomicSEM_Figure2_source_Jupyter.py"
# ).read())
# =============================================================================

FOLLOWUP_ROOT = Path("/path/to/panvas_analysis/Genomic_SEM_complete_followup")
GENOMIC_SEM_ROOT = Path("/path/to/panvas_analysis/Genomic_SEM")
OUT = FOLLOWUP_ROOT / "10_figure2_source_export_v1"
COPIED = OUT / "01_copied_tables"
GWAS_OUT = OUT / "02_genomewide_plot_data"
MODEL_OUT = OUT / "03_model_and_LDSC"
FINEMAP_OUT = OUT / "04_finemapping"

for directory in [OUT, COPIED, GWAS_OUT, MODEL_OUT, FINEMAP_OUT]:
    directory.mkdir(parents=True, exist_ok=True)


def log(message: str) -> None:
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def norm(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(text).lower())


def md5sum(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def safe_relative(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except Exception:
        return str(path)


def read_table(path: Path, nrows: int | None = None, usecols=None) -> pd.DataFrame:
    """Read common Genomic SEM/FUMA text formats without modifying the source."""
    compression = "gzip" if path.suffix == ".gz" else "infer"
    attempts = [
        {"sep": "\t", "compression": compression},
        {"sep": r"\s+", "compression": compression},
        {"sep": ",", "compression": compression},
    ]
    last_error = None
    for kwargs in attempts:
        try:
            frame = pd.read_csv(
                path,
                nrows=nrows,
                usecols=usecols,
                low_memory=False,
                **kwargs,
            )
            if frame.shape[1] > 1:
                return frame
        except Exception as exc:
            last_error = exc
    raise RuntimeError(f"Could not read {path}: {last_error}")


def find_first(candidates: list[Path]) -> Path | None:
    for path in candidates:
        if path.exists() and path.is_file():
            return path
    return None


def recursive_candidates(root: Path, include_terms: list[str], suffixes=None) -> list[Path]:
    suffixes = suffixes or [".tsv", ".txt", ".csv", ".gz", ".rds"]
    hits = []
    if not root.exists():
        return hits
    for path in root.rglob("*"):
        if not path.is_file() or OUT in path.parents:
            continue
        lower = str(path).lower()
        if not all(term.lower() in lower for term in include_terms):
            continue
        if not any(path.name.lower().endswith(s) for s in suffixes):
            continue
        hits.append(path)
    return sorted(hits, key=lambda p: (len(str(p)), str(p)))


def resolve_column(columns, aliases, required=True):
    mapping = {norm(c): c for c in columns}
    for alias in aliases:
        key = norm(alias)
        if key in mapping:
            return mapping[key]
    if required:
        raise KeyError(
            f"Could not resolve any of {aliases}. Available columns: {list(columns)}"
        )
    return None


def copy_small_table(path: Path, destination_dir: Path, label: str, max_mb=80) -> dict:
    size_mb = path.stat().st_size / 1024**2
    record = {
        "label": label,
        "source": str(path),
        "size_mb": round(size_mb, 3),
        "copied": False,
        "destination": "",
        "md5": "",
    }
    if size_mb > max_mb:
        return record
    destination = destination_dir / f"{label}__{path.name}"
    shutil.copy2(path, destination)
    record.update(
        copied=True,
        destination=str(destination),
        md5=md5sum(destination),
    )
    return record


# =============================================================================
# 1. Input discovery and immutable source inventory
# =============================================================================

log("Discovering primary and leave-one-CAD multivariate GWAS files")

primary_gwas = find_first([
    FOLLOWUP_ROOT / "05_bifactor_GWAS/full_bifactor_G/full_bifactor_G_GWAS_compact.tsv.gz",
])
if primary_gwas is None:
    hits = recursive_candidates(FOLLOWUP_ROOT, ["full_bifactor_g", "compact"])
    primary_gwas = hits[0] if hits else None

lot_hits = recursive_candidates(FOLLOWUP_ROOT, ["leave_one_cad", "bifactor", "compact"])
if not lot_hits:
    lot_hits = recursive_candidates(FOLLOWUP_ROOT, ["leave_one_cad", "gwas"])
leave_cad_gwas = lot_hits[0] if lot_hits else None

if primary_gwas is None:
    raise FileNotFoundError(
        "Primary compact Genomic SEM GWAS was not found below " + str(FOLLOWUP_ROOT)
    )

log(f"Primary GWAS: {primary_gwas}")
log(f"Leave-one-CAD GWAS: {leave_cad_gwas or 'NOT FOUND'}")

source_inventory = []
for role, path in [("primary_gwas", primary_gwas), ("leave_one_CAD_gwas", leave_cad_gwas)]:
    if path is None:
        continue
    source_inventory.append({
        "role": role,
        "path": str(path),
        "relative_to_followup": safe_relative(path, FOLLOWUP_ROOT),
        "size_bytes": path.stat().st_size,
        "modified_time": datetime.fromtimestamp(path.stat().st_mtime).isoformat(),
    })

pd.DataFrame(source_inventory).to_csv(
    OUT / "source_file_inventory.tsv", sep="\t", index=False
)


# =============================================================================
# 2. Copy compact LDSC, model, leave-one and fine-mapping tables
# =============================================================================

copy_manifest = []

metadata_candidates = [
    GENOMIC_SEM_ROOT / "gwasfile_info.txt",
    GENOMIC_SEM_ROOT / "gwas_file_info.txt",
]
metadata = find_first(metadata_candidates)
if metadata is not None:
    copy_manifest.append(copy_small_table(metadata, MODEL_OUT, "GWAS_source_metadata"))

search_specs = {
    "LDSC": ["ldsc"],
    "model_comparison": ["model", "comparison"],
    "CFA": ["cfa"],
    "bifactor": ["bifactor"],
    "loading": ["loading"],
    "leave_one_trait": ["leave", "trait"],
    "leave_one_CAD": ["leave", "cad"],
    "Q_SNP": ["q_snp"],
    "QSNP": ["qsnp"],
    "effective_N": ["effective", "n"],
}

seen = set()
for label, terms in search_specs.items():
    for i, path in enumerate(recursive_candidates(FOLLOWUP_ROOT, terms)):
        if path in seen:
            continue
        seen.add(path)
        copy_manifest.append(copy_small_table(path, MODEL_OUT, f"{label}_{i+1:02d}"))

finemap_roots = [
    FOLLOWUP_ROOT / "09_finemapping/susieRSS_mismatch_aware_v3/aggregate",
    FOLLOWUP_ROOT / "09_finemapping",
]
finemap_names = [
    "all_region_summary_with_acceptance.tsv",
    "credible_sets_with_acceptance.tsv",
    "mismatch_corrected_variants_PIP_ge_0.10.tsv",
    "primary_vs_leave_one_CAD_mismatch_corrected_PIP.tsv.gz",
]
for name in finemap_names:
    found = None
    for root in finemap_roots:
        if not root.exists():
            continue
        direct = root / name
        if direct.exists():
            found = direct
            break
        matches = list(root.rglob(name))
        if matches:
            found = matches[0]
            break
    if found is not None:
        copy_manifest.append(copy_small_table(found, FINEMAP_OUT, "finemap"))

pd.DataFrame(copy_manifest).to_csv(
    OUT / "copied_table_manifest.tsv", sep="\t", index=False
)


# =============================================================================
# 3. Genome-wide plotting export
# =============================================================================

CHR_ALIASES = ["CHR", "chromosome", "chrom", "chr"]
BP_ALIASES = ["BP", "POS", "position", "base_pair_location", "basepair"]
SNP_ALIASES = ["SNP", "rsID", "rsid", "ID", "variant_id", "MarkerName"]
P_ALIASES = ["P", "P_G", "p_G", "PVAL", "p_value", "pvalue", "Pvalue"]
BETA_ALIASES = ["BETA", "BETA_G", "beta_G", "Estimate", "est", "b"]
SE_ALIASES = ["SE", "SE_G", "se_G", "StdErr", "stderr"]
Q_P_ALIASES = [
    "Q_SNP_pval", "Q_SNP_P", "Q_SNP_p", "QSNP_P", "Qp", "Q_pval",
    "Q_p", "Q_SNP", "QSNP",
]


def inspect_gwas_columns(path: Path) -> dict:
    head = read_table(path, nrows=5)
    columns = list(head.columns)
    resolved = {
        "chr": resolve_column(columns, CHR_ALIASES),
        "bp": resolve_column(columns, BP_ALIASES),
        "snp": resolve_column(columns, SNP_ALIASES, required=False),
        "p": resolve_column(columns, P_ALIASES),
        "beta": resolve_column(columns, BETA_ALIASES, required=False),
        "se": resolve_column(columns, SE_ALIASES, required=False),
        "q_p": resolve_column(columns, Q_P_ALIASES, required=False),
    }
    return {"columns": columns, "resolved": resolved}


def cleaned_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def export_gwas_plot_data(path: Path, label: str, bin_size=25_000, chunksize=500_000):
    info = inspect_gwas_columns(path)
    resolved = info["resolved"]
    log(f"{label} resolved columns: {resolved}")

    selected_columns = list(dict.fromkeys(v for v in resolved.values() if v is not None))
    delimiter = "\t"
    try:
        iterator = pd.read_csv(
            path,
            sep=delimiter,
            compression="gzip" if path.suffix == ".gz" else "infer",
            usecols=selected_columns,
            chunksize=chunksize,
            low_memory=False,
        )
        first_chunk = next(iterator)
        iterator = itertools.chain([first_chunk], iterator)
    except Exception:
        iterator = pd.read_csv(
            path,
            sep=r"\s+",
            compression="gzip" if path.suffix == ".gz" else "infer",
            usecols=selected_columns,
            chunksize=chunksize,
            low_memory=False,
        )

    plot_parts = []
    gws_parts = []
    q_gws_parts = []
    audit_rows = []

    for chunk_index, chunk in enumerate(iterator, start=1):
        canonical = pd.DataFrame(index=chunk.index)
        canonical["CHR"] = cleaned_numeric(chunk[resolved["chr"]])
        canonical["BP"] = cleaned_numeric(chunk[resolved["bp"]])
        if resolved["snp"] is not None:
            canonical["SNP"] = chunk[resolved["snp"]].astype(str)
        else:
            canonical["SNP"] = (
                canonical["CHR"].astype("Int64").astype(str)
                + ":" + canonical["BP"].astype("Int64").astype(str)
            )
        canonical["P_G"] = cleaned_numeric(chunk[resolved["p"]])
        canonical["BETA_G"] = (
            cleaned_numeric(chunk[resolved["beta"]]) if resolved["beta"] else np.nan
        )
        canonical["SE_G"] = (
            cleaned_numeric(chunk[resolved["se"]]) if resolved["se"] else np.nan
        )
        canonical["P_QSNP"] = (
            cleaned_numeric(chunk[resolved["q_p"]]) if resolved["q_p"] else np.nan
        )

        canonical = canonical[
            canonical["CHR"].between(1, 22)
            & canonical["BP"].gt(0)
            & canonical["P_G"].between(0, 1, inclusive="both")
        ].copy()
        canonical["CHR"] = canonical["CHR"].astype(int)
        canonical["BP"] = canonical["BP"].astype(int)
        canonical.loc[canonical["P_G"] <= 0, "P_G"] = np.nextafter(0, 1)
        canonical.loc[canonical["P_QSNP"] <= 0, "P_QSNP"] = np.nextafter(0, 1)
        canonical["BIN25K"] = canonical["BP"] // bin_size

        keep_signal = canonical["P_G"].le(1e-5)
        if canonical["P_QSNP"].notna().any():
            keep_signal |= canonical["P_QSNP"].le(1e-5)

        idx_g = canonical.groupby(["CHR", "BIN25K"], sort=False)["P_G"].idxmin()
        keep_idx = set(idx_g.tolist())
        if canonical["P_QSNP"].notna().any():
            q_valid = canonical.dropna(subset=["P_QSNP"])
            if not q_valid.empty:
                idx_q = q_valid.groupby(["CHR", "BIN25K"], sort=False)["P_QSNP"].idxmin()
                keep_idx.update(idx_q.tolist())

        keep_background = canonical.index.isin(keep_idx)
        retained = canonical[keep_signal | keep_background].copy()
        plot_parts.append(retained)
        gws_parts.append(canonical[canonical["P_G"].le(5e-8)].copy())
        if canonical["P_QSNP"].notna().any():
            q_gws_parts.append(canonical[canonical["P_QSNP"].le(5e-8)].copy())

        audit_rows.append({
            "chunk": chunk_index,
            "n_valid": len(canonical),
            "n_plot_retained": len(retained),
            "n_G_GWS": int(canonical["P_G"].le(5e-8).sum()),
            "n_Q_GWS": int(canonical["P_QSNP"].le(5e-8).sum()),
        })
        if chunk_index % 5 == 0:
            log(f"{label}: processed {chunk_index:,} chunks")

    plot = pd.concat(plot_parts, ignore_index=True).drop_duplicates(["CHR", "BP", "SNP"])
    plot["minus_log10_P_G"] = -np.log10(plot["P_G"].clip(lower=np.nextafter(0, 1)))
    plot["minus_log10_P_QSNP"] = -np.log10(
        plot["P_QSNP"].clip(lower=np.nextafter(0, 1))
    )
    plot = plot.sort_values(["CHR", "BP", "P_G"])

    gws = pd.concat(gws_parts, ignore_index=True) if gws_parts else pd.DataFrame()
    q_gws = pd.concat(q_gws_parts, ignore_index=True) if q_gws_parts else pd.DataFrame()

    plot_path = GWAS_OUT / f"{label}__Miami_Manhattan_plot_data.tsv.gz"
    gws_path = GWAS_OUT / f"{label}__G_factor_GWS_variants.tsv.gz"
    q_path = GWAS_OUT / f"{label}__QSNP_GWS_variants.tsv.gz"
    audit_path = GWAS_OUT / f"{label}__chunk_audit.tsv"
    columns_path = GWAS_OUT / f"{label}__column_resolution.json"

    plot.to_csv(plot_path, sep="\t", index=False, compression="gzip")
    gws.to_csv(gws_path, sep="\t", index=False, compression="gzip")
    q_gws.to_csv(q_path, sep="\t", index=False, compression="gzip")
    pd.DataFrame(audit_rows).to_csv(audit_path, sep="\t", index=False)
    columns_path.write_text(json.dumps(info, indent=2, ensure_ascii=False))

    summary = {
        "label": label,
        "source": str(path),
        "n_plot_rows": int(len(plot)),
        "n_G_GWS_rows": int(len(gws)),
        "n_Q_GWS_rows": int(len(q_gws)),
        "has_QSNP": bool(resolved["q_p"] is not None),
        "plot_file": str(plot_path),
        "GWS_file": str(gws_path),
        "QSNP_file": str(q_path),
    }
    return summary


gwas_summaries = [export_gwas_plot_data(primary_gwas, "primary_full_bifactor_G")]
if leave_cad_gwas is not None:
    gwas_summaries.append(export_gwas_plot_data(leave_cad_gwas, "leave_one_CAD_bifactor_G"))

pd.DataFrame(gwas_summaries).to_csv(
    GWAS_OUT / "GWAS_plot_export_summary.tsv", sep="\t", index=False
)


# =============================================================================
# 4. Export package manifest and ZIP archive
# =============================================================================

manifest_rows = []
for path in sorted(OUT.rglob("*")):
    if not path.is_file() or path.name.endswith(".zip"):
        continue
    manifest_rows.append({
        "relative_path": safe_relative(path, OUT),
        "size_bytes": path.stat().st_size,
        "md5": md5sum(path),
    })

pd.DataFrame(manifest_rows).to_csv(OUT / "export_manifest.tsv", sep="\t", index=False)

readme = f"""Pan-vascular Genomic SEM Figure 2 source export
Created: {datetime.now().isoformat()}

Contents
01_copied_tables: small source tables discovered without modifying originals
02_genomewide_plot_data: deterministic plotting subset, G-factor GWS hits and Q_SNP GWS hits
03_model_and_LDSC: LDSC, CFA, bifactor and leave-one sensitivity tables
04_finemapping: aggregate SuSiE-RSS credible-set and PIP tables when found

Plot downsampling rule
- retain every variant with G-factor P <= 1e-5 or Q_SNP P <= 1e-5
- additionally retain the minimum G P and minimum Q_SNP P variant in each 25-kb genomic bin
- export every genome-wide-significant G and Q_SNP row separately

This is a plotting/export package. Original analysis files were never moved, deleted or edited.
"""
(OUT / "README.txt").write_text(readme, encoding="utf-8")

archive_base = str(OUT)
archive_path = Path(shutil.make_archive(archive_base, "zip", root_dir=OUT.parent, base_dir=OUT.name))

log("EXPORT COMPLETE")
log(f"Directory: {OUT}")
log(f"ZIP: {archive_path}")
log("Download the ZIP and send it back for the final Python-only Figure 2 build.")
