from __future__ import annotations

import gzip
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


# =============================================================================
# Configuration
# =============================================================================

ROOT = Path("/path/to/panvas_analysis/Genomic_SEM_complete_followup/09_finemapping")
OUT_ROOT = ROOT / "susieRSS_finemapping_v1"
OUT_ROOT.mkdir(parents=True, exist_ok=True)

PRIMARY_GWAS = (
    ROOT
    / "input/primary/full_bifactor_G_FUMA_LDaware_QSNP_locus_filtered_primary_r2_0.10.tsv.gz"
)
LOT_GWAS = (
    ROOT
    / "input/leave_one_CAD/leave_one_CAD_bifactor_G_FUMA_LDaware_QSNP_locus_filtered_primary_r2_0.10.tsv.gz"
)
PRIMARY_LOCI = Path(
    "/path/to/panvas_analysis/Genomic_SEM_complete_followup/08_VEP_functional_annotation/"
    "input/FUMA_LDr20.1Q_bifactor/GenomicRiskLoci.txt"
)
LOT_LOCI = Path(
    "/path/to/panvas_analysis/Genomic_SEM_complete_followup/08_VEP_functional_annotation/"
    "input/FUMA_LOTCAD_LDr20.1Q_bifactor/GenomicRiskLoci.txt"
)
LD_PREFIX = Path(
    "/path/to/panvas_analysis/Genomic_SEM_complete_followup/09_finemapping/"
    "reference/1000G_EUR_GRCh37_MAF01/data_maf0.01_rs_ref"
)
PLINK = Path("/path/to/panvas_analysis/miniconda3/envs/finemapping_tools/bin/plink")
RSCRIPT = Path("/usr/bin/Rscript")
R_LIBRARY = ROOT / "R_library_cran_stable"

SCENARIOS = {
    "primary": {"gwas": PRIMARY_GWAS, "loci": PRIMARY_LOCI},
    "leave_one_CAD": {"gwas": LOT_GWAS, "loci": LOT_LOCI},
}

# FUMA loci are LD-defined. A modest flank protects against edge truncation
# without creating very large, weakly connected chromosome-scale regions.
REGION_PADDING_BP = 100_000
SUMMARY_CHUNK_SIZE = 400_000
MIN_MODEL_VARIANTS = 20
MAX_MODEL_VARIANTS = 20_000
PLINK_THREADS = 8
PLINK_MEMORY_MB = 16_000
RUN_NHAT_SENSITIVITY = True
FORCE_REBUILD_SUMMARY_EXTRACTION = False
FORCE_RERUN_COMPLETED_REGIONS = False

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
PALINDROMIC = {frozenset(("A", "T")), frozenset(("C", "G"))}
REQUIRED_GWAS_COLUMNS = {"SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"}


# =============================================================================
# General helpers
# =============================================================================

def now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def clean_chr(value):
    s = str(value).replace("chr", "").replace("CHR", "").strip()
    try:
        x = int(float(s))
    except Exception:
        return np.nan
    return x if 1 <= x <= 22 else np.nan


def complement_allele(a: str) -> str | None:
    a = str(a).upper()
    if not a or any(x not in COMPLEMENT for x in a):
        return None
    return "".join(COMPLEMENT[x] for x in a)


def is_palindromic(a1: str, a2: str) -> bool:
    return len(a1) == 1 and len(a2) == 1 and frozenset((a1, a2)) in PALINDROMIC


def valid_file(path: Path, min_bytes: int = 20) -> bool:
    return path.exists() and path.stat().st_size >= min_bytes


def run_logged(cmd: list[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as log:
        log.write(f"\n[{now()}] $ {' '.join(map(str, cmd))}\n")
        log.flush()
        result = subprocess.run(
            [str(x) for x in cmd],
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd)


def require_inputs() -> None:
    required = [
        PRIMARY_GWAS,
        LOT_GWAS,
        PRIMARY_LOCI,
        LOT_LOCI,
        Path(str(LD_PREFIX) + ".bed"),
        Path(str(LD_PREFIX) + ".bim"),
        Path(str(LD_PREFIX) + ".fam"),
        PLINK,
        RSCRIPT,
    ]
    missing = [str(x) for x in required if not x.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n" + "\n".join(missing))

    code = f'''
.libPaths(c("{R_LIBRARY}", .libPaths()))
cat(requireNamespace("data.table", quietly=TRUE), "\n")
cat(requireNamespace("susieR", quietly=TRUE), "\n")
if (requireNamespace("susieR", quietly=TRUE)) cat(as.character(packageVersion("susieR")), "\n")
'''
    check = subprocess.run(
        [str(RSCRIPT), "--vanilla", "-e", code], capture_output=True, text=True, check=True
    )
    vals = [x.strip() for x in check.stdout.splitlines() if x.strip()]
    if len(vals) < 2 or vals[0] != "TRUE" or vals[1] != "TRUE":
        raise RuntimeError("data.table or susieR is unavailable in the project R library")
    print("R package audit:", vals, flush=True)


# =============================================================================
# Region construction and one-pass GWAS extraction
# =============================================================================

def build_merged_regions(loci_path: Path, scenario: str) -> pd.DataFrame:
    loci = pd.read_csv(loci_path, sep="\t")
    required = {"GenomicLocus", "rsID", "chr", "pos", "start", "end"}
    if not required.issubset(loci.columns):
        raise ValueError(f"{loci_path} lacks columns: {sorted(required - set(loci.columns))}")

    loci = loci.copy()
    loci["CHR"] = loci["chr"].map(clean_chr)
    for c in ["pos", "start", "end"]:
        loci[c] = pd.to_numeric(loci[c], errors="coerce")
    loci = loci.dropna(subset=["CHR", "pos", "start", "end"]).copy()
    loci["CHR"] = loci["CHR"].astype(int)
    loci["region_start"] = (
        loci[["pos", "start"]].min(axis=1).astype(int) - REGION_PADDING_BP
    ).clip(lower=1)
    loci["region_end"] = (
        loci[["pos", "end"]].max(axis=1).astype(int) + REGION_PADDING_BP
    )
    loci = loci.sort_values(["CHR", "region_start", "region_end"])

    merged = []
    for chromosome, block in loci.groupby("CHR", sort=True):
        current = None
        for row in block.itertuples(index=False):
            member = {
                "GenomicLocus": str(row.GenomicLocus),
                "lead_rsID": str(row.rsID),
                "lead_pos": int(row.pos),
            }
            if current is None or int(row.region_start) > current["end"]:
                if current is not None:
                    merged.append(current)
                current = {
                    "scenario": scenario,
                    "CHR": int(chromosome),
                    "start": int(row.region_start),
                    "end": int(row.region_end),
                    "members": [member],
                }
            else:
                current["end"] = max(current["end"], int(row.region_end))
                current["members"].append(member)
        if current is not None:
            merged.append(current)

    rows = []
    for i, x in enumerate(merged, 1):
        region_id = f"{scenario}__R{i:03d}__chr{x['CHR']}__{x['start']}_{x['end']}"
        rows.append(
            {
                "scenario": scenario,
                "region_id": region_id,
                "CHR": x["CHR"],
                "start": x["start"],
                "end": x["end"],
                "width_bp": x["end"] - x["start"] + 1,
                "n_FUMA_loci": len(x["members"]),
                "FUMA_GenomicLoci": ";".join(m["GenomicLocus"] for m in x["members"]),
                "lead_rsIDs": ";".join(m["lead_rsID"] for m in x["members"]),
                "lead_positions": ";".join(str(m["lead_pos"]) for m in x["members"]),
            }
        )
    return pd.DataFrame(rows)


def extract_region_summaries(gwas_path: Path, regions: pd.DataFrame, scenario_dir: Path) -> None:
    marker = scenario_dir / "SUMMARY_EXTRACTION_COMPLETE.json"
    expected = [scenario_dir / "regions" / rid / "summary_input.tsv.gz" for rid in regions.region_id]
    if (
        not FORCE_REBUILD_SUMMARY_EXTRACTION
        and valid_file(marker)
        and all(valid_file(x) for x in expected)
    ):
        print(f"[{now()}] Reusing extracted GWAS regions for {regions.scenario.iloc[0]}")
        return

    stores: dict[str, list[pd.DataFrame]] = {rid: [] for rid in regions.region_id}
    lookup = {}
    for chromosome, block in regions.groupby("CHR", sort=True):
        b = block.sort_values("start")
        lookup[int(chromosome)] = {
            "starts": b.start.to_numpy(np.int64),
            "ends": b.end.to_numpy(np.int64),
            "ids": b.region_id.to_numpy(str),
        }

    n_rows = 0
    n_selected = 0
    reader = pd.read_csv(
        gwas_path,
        sep="\t",
        compression="gzip",
        chunksize=SUMMARY_CHUNK_SIZE,
        low_memory=False,
    )
    for chunk_no, chunk in enumerate(reader, 1):
        missing = REQUIRED_GWAS_COLUMNS - set(chunk.columns)
        if missing:
            raise ValueError(f"GWAS file lacks required columns: {sorted(missing)}")
        chunk = chunk.copy()
        chunk["CHR"] = chunk["CHR"].map(clean_chr)
        chunk["BP"] = pd.to_numeric(chunk["BP"], errors="coerce")
        chunk = chunk.dropna(subset=["CHR", "BP", "SNP", "A1", "A2"])
        chunk["CHR"] = chunk["CHR"].astype(int)
        chunk["BP"] = chunk["BP"].astype(np.int64)
        for c in ["BETA", "SE", "P", "N"]:
            chunk[c] = pd.to_numeric(chunk[c], errors="coerce")
        if "MAF" in chunk:
            chunk["MAF"] = pd.to_numeric(chunk["MAF"], errors="coerce")
        chunk["SNP"] = chunk["SNP"].astype(str)
        chunk["A1"] = chunk["A1"].astype(str).str.upper()
        chunk["A2"] = chunk["A2"].astype(str).str.upper()
        n_rows += len(chunk)

        for chromosome, idx in chunk.groupby("CHR", sort=False).groups.items():
            if int(chromosome) not in lookup:
                continue
            sub = chunk.loc[idx]
            info = lookup[int(chromosome)]
            bp = sub["BP"].to_numpy(np.int64)
            region_index = np.searchsorted(info["starts"], bp, side="right") - 1
            valid = region_index >= 0
            valid &= np.where(valid, bp <= info["ends"][np.maximum(region_index, 0)], False)
            if not valid.any():
                continue
            selected = sub.iloc[np.flatnonzero(valid)].copy()
            selected["_region_id"] = info["ids"][region_index[valid]]
            for rid, part in selected.groupby("_region_id", sort=False):
                stores[str(rid)].append(part.drop(columns="_region_id"))
                n_selected += len(part)

        if chunk_no % 3 == 0:
            print(
                f"[{now()}] {gwas_path.name}: chunks={chunk_no}, "
                f"rows={n_rows:,}, selected={n_selected:,}",
                flush=True,
            )

    extraction_rows = []
    for row in regions.itertuples(index=False):
        region_dir = scenario_dir / "regions" / row.region_id
        region_dir.mkdir(parents=True, exist_ok=True)
        parts = stores[row.region_id]
        if parts:
            out = pd.concat(parts, ignore_index=True)
            out = out[
                np.isfinite(out["BETA"])
                & np.isfinite(out["SE"])
                & (out["SE"] > 0)
                & np.isfinite(out["P"])
            ].copy()
            out = out.sort_values(["P", "SNP"]).drop_duplicates("SNP", keep="first")
        else:
            out = pd.DataFrame(columns=sorted(REQUIRED_GWAS_COLUMNS))
        out.to_csv(region_dir / "summary_input.tsv.gz", sep="\t", index=False, compression="gzip")
        extraction_rows.append({"region_id": row.region_id, "n_summary_variants": len(out)})

    pd.DataFrame(extraction_rows).to_csv(
        scenario_dir / "summary_extraction_audit.tsv", sep="\t", index=False
    )
    marker.write_text(
        json.dumps(
            {
                "completed_at": now(),
                "gwas": str(gwas_path),
                "n_rows_read": n_rows,
                "n_rows_selected": n_selected,
                "padding_bp": REGION_PADDING_BP,
            },
            indent=2,
        )
    )


# =============================================================================
# Reference alignment and PLINK genotype extraction
# =============================================================================

def read_bim(path: Path) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        names=["CHR_REF", "SNP", "CM", "BP_REF", "REF_A1", "REF_A2"],
        dtype={"SNP": str, "REF_A1": str, "REF_A2": str},
    )


def prepare_reference_for_region(region, region_dir: Path) -> dict:
    summary_path = region_dir / "summary_input.tsv.gz"
    ss = pd.read_csv(summary_path, sep="\t")
    ss["SNP"] = ss["SNP"].astype(str)
    ss["A1"] = ss["A1"].astype(str).str.upper()
    ss["A2"] = ss["A2"].astype(str).str.upper()
    ss = ss.sort_values(["P", "SNP"]).drop_duplicates("SNP", keep="first")

    audit_rows = []
    if len(ss) < MIN_MODEL_VARIANTS:
        return {"status": "INSUFFICIENT_SUMMARY_VARIANTS", "n_summary": len(ss)}

    ids_path = region_dir / "summary_variant_ids.txt"
    ids_path.write_text("\n".join(ss["SNP"].tolist()) + "\n")
    pre = region_dir / "reference_preliminary"
    plink_log = region_dir / "plink_commands.log"

    if not valid_file(Path(str(pre) + ".bim")):
        run_logged(
            [
                str(PLINK),
                "--bfile", str(LD_PREFIX),
                "--chr", str(region.CHR),
                "--from-bp", str(region.start),
                "--to-bp", str(region.end),
                "--extract", str(ids_path),
                "--keep-allele-order",
                "--make-bed",
                "--threads", str(PLINK_THREADS),
                "--memory", str(PLINK_MEMORY_MB),
                "--out", str(pre),
            ],
            plink_log,
        )

    bim = read_bim(Path(str(pre) + ".bim"))
    merged = ss.merge(bim, on="SNP", how="left", validate="one_to_one")
    eligible = []
    for row in merged.itertuples(index=False):
        reason = None
        alignment = None
        ref_effect = None
        ref_other = None
        r1 = str(row.REF_A1).upper() if pd.notna(row.REF_A1) else ""
        r2 = str(row.REF_A2).upper() if pd.notna(row.REF_A2) else ""
        g1 = str(row.A1).upper()
        g2 = str(row.A2).upper()

        if not r1 or not r2:
            reason = "not_in_reference"
        elif is_palindromic(g1, g2):
            reason = "palindromic_AT_or_CG"
        elif {g1, g2} == {r1, r2}:
            alignment = "direct"
            ref_effect, ref_other = g1, g2
        else:
            cg1, cg2 = complement_allele(g1), complement_allele(g2)
            if cg1 is not None and {cg1, cg2} == {r1, r2}:
                alignment = "strand_complement"
                ref_effect, ref_other = cg1, cg2
            else:
                reason = "allele_set_mismatch"

        rec = {
            "SNP": row.SNP,
            "GWAS_A1": g1,
            "GWAS_A2": g2,
            "REF_A1_initial": r1,
            "REF_A2_initial": r2,
            "alignment": alignment,
            "REF_EFFECT_ALLELE": ref_effect,
            "REF_OTHER_ALLELE": ref_other,
            "exclusion_reason": reason,
        }
        audit_rows.append(rec)
        if reason is None:
            eligible.append(rec)

    audit = pd.DataFrame(audit_rows)
    audit.to_csv(region_dir / "allele_alignment_audit.tsv", sep="\t", index=False)
    eligible = pd.DataFrame(eligible)
    if len(eligible) < MIN_MODEL_VARIANTS:
        return {
            "status": "INSUFFICIENT_ALIGNED_VARIANTS",
            "n_summary": len(ss),
            "n_reference": len(bim),
            "n_aligned": len(eligible),
        }

    eligible_ids = region_dir / "eligible_variant_ids.txt"
    eligible_ids.write_text("\n".join(eligible.SNP.tolist()) + "\n")
    force_a1 = region_dir / "force_reference_effect_allele.tsv"
    eligible[["SNP", "REF_EFFECT_ALLELE"]].to_csv(
        force_a1, sep="\t", index=False, header=False
    )

    aligned = region_dir / "reference_aligned"
    if not valid_file(Path(str(aligned) + ".bim")):
        run_logged(
            [
                str(PLINK),
                "--bfile", str(pre),
                "--extract", str(eligible_ids),
                "--a1-allele", str(force_a1), "2", "1",
                "--keep-allele-order",
                "--make-bed",
                "--threads", str(PLINK_THREADS),
                "--memory", str(PLINK_MEMORY_MB),
                "--out", str(aligned),
            ],
            plink_log,
        )

    traw_prefix = region_dir / "reference_aligned_genotypes"
    traw = Path(str(traw_prefix) + ".traw")
    if not valid_file(traw):
        run_logged(
            [
                str(PLINK),
                "--bfile", str(aligned),
                "--keep-allele-order",
                "--recode", "A-transpose",
                "--threads", str(PLINK_THREADS),
                "--memory", str(PLINK_MEMORY_MB),
                "--out", str(traw_prefix),
            ],
            plink_log,
        )

    final_bim = read_bim(Path(str(aligned) + ".bim"))
    model = ss.merge(
        eligible[["SNP", "alignment", "REF_EFFECT_ALLELE", "REF_OTHER_ALLELE"]],
        on="SNP",
        how="inner",
        validate="one_to_one",
    ).merge(
        final_bim[["SNP", "REF_A1", "REF_A2"]], on="SNP", how="inner", validate="one_to_one"
    )
    model.to_csv(region_dir / "model_input.tsv.gz", sep="\t", index=False, compression="gzip")

    return {
        "status": "REFERENCE_READY",
        "n_summary": len(ss),
        "n_reference": len(bim),
        "n_aligned": len(model),
        "traw": str(traw),
    }


# =============================================================================
# R helper for SuSiE-RSS
# =============================================================================

R_HELPER = r'''
args <- commandArgs(trailingOnly=TRUE)
region_dir <- args[[1]]
scenario <- args[[2]]
region_id <- args[[3]]
run_nhat <- identical(args[[4]], "TRUE")
r_lib <- args[[5]]
.libPaths(c(r_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(susieR)
})

model_path <- file.path(region_dir, "model_input.tsv.gz")
traw_path <- file.path(region_dir, "reference_aligned_genotypes.traw")
ss <- fread(model_path)
gt <- fread(traw_path)

if (!all(c("SNP", "COUNTED") %in% names(gt))) {
  stop("Unexpected PLINK .traw columns: ", paste(names(gt)[1:min(10, ncol(gt))], collapse=","))
}

ord <- match(gt$SNP, ss$SNP)
keep_match <- !is.na(ord)
gt <- gt[keep_match]
ord <- ord[keep_match]
ss <- ss[ord]
stopifnot(identical(as.character(gt$SNP), as.character(ss$SNP)))

counted <- toupper(as.character(gt$COUNTED))
expected <- toupper(as.character(ss$REF_EFFECT_ALLELE))
other <- toupper(as.character(ss$REF_OTHER_ALLELE))
sign <- ifelse(counted == expected, 1, ifelse(counted == other, -1, NA_real_))
keep_allele <- is.finite(sign)
gt <- gt[keep_allele]
ss <- ss[keep_allele]
sign <- sign[keep_allele]
counted <- counted[keep_allele]

G <- as.matrix(gt[, 7:ncol(gt), with=FALSE])
storage.mode(G) <- "double"
row_mean <- rowMeans(G, na.rm=TRUE)
for (i in which(rowSums(is.na(G)) > 0)) G[i, is.na(G[i, ])] <- row_mean[i]
row_sd <- apply(G, 1, sd)
keep_poly <- is.finite(row_sd) & row_sd > 0
G <- G[keep_poly, , drop=FALSE]
ss <- ss[keep_poly]
sign <- sign[keep_poly]
counted <- counted[keep_poly]

if (nrow(G) < 20) stop("Fewer than 20 polymorphic aligned variants")

X <- t(G)
X <- scale(X, center=TRUE, scale=TRUE)
R <- crossprod(X) / (nrow(X) - 1)
R <- (R + t(R)) / 2
diag(R) <- 1
z <- (as.numeric(ss$BETA) / as.numeric(ss$SE)) * sign

ok <- is.finite(z) & apply(R, 1, function(x) all(is.finite(x)))
if (!all(ok)) {
  z <- z[ok]
  R <- R[ok, ok, drop=FALSE]
  ss <- ss[ok]
  counted <- counted[ok]
  sign <- sign[ok]
}

p <- length(z)
L <- min(10L, p)
n_ref <- nrow(X)
n_hat <- suppressWarnings(round(median(as.numeric(ss$N), na.rm=TRUE)))
if (!is.finite(n_hat) || n_hat < 100) n_hat <- NA_integer_

effective_rank <- p^2 / sum(R^2)
r_over_B <- effective_rank / n_ref
diagnostic_s_status <- if (p <= 5000L) "attempted" else "not_run_p_gt_5000"
lambda_s <- if (p <= 5000L) {
  tryCatch(
    susieR::estimate_s_rss(z, R, n=ifelse(is.na(n_hat), 100000L, n_hat), method="null-mle"),
    error=function(e) {
      diagnostic_s_status <<- paste0("failed: ", conditionMessage(e))
      NA_real_
    }
  )
} else NA_real_

fit_one <- function(method, n_value=NULL) {
  call_args <- list(
    z=z,
    R=R,
    L=L,
    estimate_residual_variance=FALSE,
    max_iter=200,
    tol=1e-3,
    coverage=0.95,
    min_abs_corr=0.5,
    refine=TRUE
  )
  if (!is.null(n_value)) call_args$n <- n_value
  refined <- TRUE
  fit <- tryCatch(
    do.call(susieR::susie_rss, call_args),
    error=function(e) {
      refined <<- FALSE
      call_args$refine <- FALSE
      do.call(susieR::susie_rss, call_args)
    }
  )

  cs <- fit$sets$cs
  membership <- rep("", p)
  cs_rows <- list()
  if (!is.null(cs) && length(cs) > 0) {
    purity <- fit$sets$purity
    coverage <- fit$sets$coverage
    for (i in seq_along(cs)) {
      idx <- as.integer(cs[[i]])
      cs_name <- if (!is.null(names(cs)) && nzchar(names(cs)[i])) names(cs)[i] else paste0("L", i)
      for (j in idx) membership[j] <- paste(c(membership[j], cs_name)[nzchar(c(membership[j], cs_name))], collapse=";")
      pur <- c(min.abs.corr=NA_real_, mean.abs.corr=NA_real_, median.abs.corr=NA_real_)
      if (!is.null(purity) && nrow(as.data.frame(purity)) >= i) {
        q <- as.numeric(purity[i, ])
        names(q) <- colnames(purity)
        pur[names(pur) %in% names(q)] <- q[names(pur)[names(pur) %in% names(q)]]
      }
      cov_i <- if (!is.null(coverage) && length(coverage) >= i) as.numeric(coverage[i]) else NA_real_
      top_idx <- idx[which.max(fit$pip[idx])]
      cs_rows[[i]] <- data.table(
        scenario=scenario,
        region_id=region_id,
        method=method,
        cs_id=cs_name,
        cs_size=length(idx),
        claimed_coverage=cov_i,
        min_abs_corr=unname(pur["min.abs.corr"]),
        mean_abs_corr=unname(pur["mean.abs.corr"]),
        median_abs_corr=unname(pur["median.abs.corr"]),
        lead_SNP=as.character(ss$SNP[top_idx]),
        lead_PIP=as.numeric(fit$pip[top_idx]),
        variants=paste(as.character(ss$SNP[idx]), collapse=";")
      )
    }
  }

  max_alpha <- apply(fit$alpha, 2, max)
  variants <- copy(ss)
  variants[, `:=`(
    scenario=scenario,
    region_id=region_id,
    method=method,
    Z_aligned=z,
    PLINK_COUNTED_ALLELE=counted,
    alignment_sign=sign,
    PIP=as.numeric(fit$pip),
    max_single_effect_alpha=as.numeric(max_alpha),
    credible_set_ids=membership,
    in_credible_set=nzchar(membership)
  )]

  top <- which.max(fit$pip)
  summary <- data.table(
    scenario=scenario,
    region_id=region_id,
    method=method,
    n_model_variants=p,
    n_reference_samples=n_ref,
    n_hat_median=n_hat,
    L=L,
    converged=isTRUE(fit$converged),
    niter=as.integer(fit$niter),
    refinement_used=refined,
    n_credible_sets=if (is.null(cs)) 0L else length(cs),
    max_PIP=as.numeric(fit$pip[top]),
    top_SNP=as.character(ss$SNP[top]),
    top_P=as.numeric(ss$P[top]),
    estimated_s_null_mle=as.numeric(lambda_s),
    estimated_s_status=diagnostic_s_status,
    LD_effective_rank=effective_rank,
    LD_effective_rank_over_reference_N=r_over_B
  )

  list(
    fit=fit,
    variants=variants,
    credible_sets=if (length(cs_rows)) rbindlist(cs_rows, fill=TRUE) else data.table(),
    summary=summary
  )
}

fits <- list(z_only_primary=fit_one("z_only_primary"))
if (run_nhat && !is.na(n_hat)) {
  fits$median_Nhat_sensitivity <- fit_one("median_Nhat_sensitivity", n_hat)
}

variant_out <- rbindlist(lapply(fits, `[[`, "variants"), fill=TRUE)
cs_out <- rbindlist(lapply(fits, `[[`, "credible_sets"), fill=TRUE)
summary_out <- rbindlist(lapply(fits, `[[`, "summary"), fill=TRUE)

fwrite(variant_out, file.path(region_dir, "variant_PIP.tsv.gz"), sep="\t")
fwrite(cs_out, file.path(region_dir, "credible_sets.tsv"), sep="\t")
fwrite(summary_out, file.path(region_dir, "region_summary.tsv"), sep="\t")
saveRDS(lapply(fits, `[[`, "fit"), file.path(region_dir, "susie_fits.rds"), compress=TRUE)
writeLines(capture.output(sessionInfo()), file.path(region_dir, "R_sessionInfo.txt"))
'''


def write_r_helper() -> Path:
    path = OUT_ROOT / "run_one_region_susieRSS.R"
    path.write_text(R_HELPER)
    return path


# =============================================================================
# Per-scenario execution
# =============================================================================

def write_status(region_dir: Path, payload: dict) -> None:
    payload = dict(payload)
    payload["timestamp"] = now()
    (region_dir / "STATUS.json").write_text(json.dumps(payload, indent=2))


def run_scenario(scenario: str, config: dict, r_helper: Path) -> pd.DataFrame:
    scenario_dir = OUT_ROOT / scenario
    scenario_dir.mkdir(parents=True, exist_ok=True)
    regions = build_merged_regions(config["loci"], scenario)
    regions.to_csv(scenario_dir / "merged_finemapping_regions.tsv", sep="\t", index=False)
    print(
        f"[{now()}] {scenario}: {len(pd.read_csv(config['loci'], sep=chr(9)))} FUMA loci -> "
        f"{len(regions)} merged fine-mapping regions",
        flush=True,
    )
    extract_region_summaries(config["gwas"], regions, scenario_dir)

    for number, region in enumerate(regions.itertuples(index=False), 1):
        region_dir = scenario_dir / "regions" / region.region_id
        done = region_dir / "FINEMAPPING_COMPLETE.json"
        if valid_file(done) and not FORCE_RERUN_COMPLETED_REGIONS:
            print(f"[{now()}] {scenario} {number}/{len(regions)} SKIP: {region.region_id}")
            continue

        print(f"[{now()}] {scenario} {number}/{len(regions)} START: {region.region_id}", flush=True)
        try:
            ref = prepare_reference_for_region(region, region_dir)
            if ref["status"] != "REFERENCE_READY":
                write_status(region_dir, ref)
                print(f"[{now()}] {region.region_id}: {ref['status']}", flush=True)
                continue
            if ref["n_aligned"] > MAX_MODEL_VARIANTS:
                payload = dict(ref)
                payload["status"] = "TOO_MANY_VARIANTS_FOR_DENSE_LD"
                payload["maximum_allowed"] = MAX_MODEL_VARIANTS
                write_status(region_dir, payload)
                print(f"[{now()}] {region.region_id}: TOO MANY VARIANTS ({ref['n_aligned']})")
                continue

            r_log = region_dir / "susieRSS.log"
            run_logged(
                [
                    str(RSCRIPT),
                    "--vanilla",
                    str(r_helper),
                    str(region_dir),
                    scenario,
                    region.region_id,
                    str(RUN_NHAT_SENSITIVITY).upper(),
                    str(R_LIBRARY),
                ],
                r_log,
            )
            if not valid_file(region_dir / "region_summary.tsv"):
                raise RuntimeError("R completed without region_summary.tsv")
            done.write_text(
                json.dumps(
                    {
                        "status": "SUCCESS",
                        "completed_at": now(),
                        **ref,
                    },
                    indent=2,
                )
            )
            write_status(region_dir, {"status": "SUCCESS", **ref})
        except Exception as exc:
            write_status(
                region_dir,
                {"status": "FAILED", "error_type": type(exc).__name__, "error": str(exc)},
            )
            print(f"[{now()}] FAILED {region.region_id}: {type(exc).__name__}: {exc}", flush=True)

    return regions


# =============================================================================
# Aggregation, robustness comparison and figures
# =============================================================================

def aggregate_outputs(all_regions: list[pd.DataFrame]) -> None:
    aggregate = OUT_ROOT / "aggregate"
    figures = OUT_ROOT / "figures"
    aggregate.mkdir(exist_ok=True)
    figures.mkdir(exist_ok=True)

    summaries = []
    credible_sets = []
    variants = []
    status_rows = []

    for scenario in SCENARIOS:
        for region_dir in sorted((OUT_ROOT / scenario / "regions").glob("*")):
            if not region_dir.is_dir():
                continue
            status_path = region_dir / "STATUS.json"
            if status_path.exists():
                x = json.loads(status_path.read_text())
                x["scenario"] = scenario
                x["region_id"] = region_dir.name
                status_rows.append(x)
            if valid_file(region_dir / "region_summary.tsv"):
                summaries.append(pd.read_csv(region_dir / "region_summary.tsv", sep="\t"))
            if valid_file(region_dir / "credible_sets.tsv"):
                try:
                    credible_sets.append(pd.read_csv(region_dir / "credible_sets.tsv", sep="\t"))
                except pd.errors.EmptyDataError:
                    pass
            if valid_file(region_dir / "variant_PIP.tsv.gz"):
                variants.append(pd.read_csv(region_dir / "variant_PIP.tsv.gz", sep="\t"))

    status = pd.DataFrame(status_rows)
    status.to_csv(aggregate / "all_region_status.tsv", sep="\t", index=False)
    summary = pd.concat(summaries, ignore_index=True) if summaries else pd.DataFrame()
    cs = pd.concat(credible_sets, ignore_index=True) if credible_sets else pd.DataFrame()
    var = pd.concat(variants, ignore_index=True) if variants else pd.DataFrame()
    summary.to_csv(aggregate / "all_region_summary.tsv", sep="\t", index=False)
    cs.to_csv(aggregate / "all_credible_sets.tsv", sep="\t", index=False)
    var.to_csv(
        aggregate / "all_variant_PIP.tsv.gz", sep="\t", index=False, compression="gzip"
    )

    if not var.empty:
        high = var[(var["method"] == "z_only_primary") & (var["PIP"] >= 0.1)].copy()
        high.sort_values(["scenario", "PIP"], ascending=[True, False]).to_csv(
            aggregate / "high_priority_variants_PIP_ge_0.10.tsv", sep="\t", index=False
        )
        in_cs = var["in_credible_set"]
        if in_cs.dtype == bool:
            in_cs_mask = in_cs
        else:
            in_cs_mask = in_cs.astype(str).str.upper().isin(["TRUE", "T", "1"])
        cs_variants = var[(var["method"] == "z_only_primary") & in_cs_mask]
        cs_variants.to_csv(
            aggregate / "all_95pct_credible_set_variants.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )

        prim = var[(var.scenario == "primary") & (var.method == "z_only_primary")][
            ["SNP", "CHR", "BP", "PIP", "credible_set_ids"]
        ].rename(columns={"PIP": "PIP_primary", "credible_set_ids": "CS_primary"})
        lot = var[(var.scenario == "leave_one_CAD") & (var.method == "z_only_primary")][
            ["SNP", "PIP", "credible_set_ids"]
        ].rename(columns={"PIP": "PIP_leave_one_CAD", "credible_set_ids": "CS_leave_one_CAD"})
        overlap = prim.merge(lot, on="SNP", how="inner")
        overlap["max_PIP"] = overlap[["PIP_primary", "PIP_leave_one_CAD"]].max(axis=1)
        overlap.sort_values("max_PIP", ascending=False).to_csv(
            aggregate / "primary_vs_leave_one_CAD_variant_PIP.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )

    try:
        import matplotlib.pyplot as plt

        if not summary.empty:
            q = summary[
                (summary["scenario"] == "primary") & (summary["method"] == "z_only_primary")
            ].copy()
            if len(q):
                fig, ax = plt.subplots(figsize=(6.4, 4.8))
                size = 28 + 18 * np.minimum(q["n_credible_sets"].to_numpy(float), 5)
                sc = ax.scatter(
                    q["estimated_s_null_mle"],
                    q["max_PIP"],
                    s=size,
                    c=q["n_credible_sets"],
                    cmap="viridis",
                    alpha=0.8,
                    edgecolor="white",
                    linewidth=0.4,
                )
                ax.set_xlabel("Estimated summary-LD mismatch (s)")
                ax.set_ylabel("Maximum posterior inclusion probability")
                ax.set_title("Fine-mapping reliability across pan-vascular loci", loc="left")
                fig.colorbar(sc, ax=ax, label="No. of 95% credible sets")
                fig.tight_layout()
                fig.savefig(figures / "primary_finemapping_reliability.png", dpi=300)
                fig.savefig(figures / "primary_finemapping_reliability.pdf")
                plt.close(fig)

        if not cs.empty:
            q = cs[(cs.scenario == "primary") & (cs.method == "z_only_primary")].copy()
            if len(q):
                fig, ax = plt.subplots(figsize=(6.2, 4.5))
                bins = np.unique(np.r_[1, np.geomspace(2, max(3, q.cs_size.max()), 14).astype(int)])
                ax.hist(q.cs_size, bins=bins, color="#287271", edgecolor="white")
                ax.set_xscale("log")
                ax.set_xlabel("Variants per 95% credible set")
                ax.set_ylabel("Credible sets")
                ax.set_title("Resolution of pan-vascular fine-mapping", loc="left")
                fig.tight_layout()
                fig.savefig(figures / "primary_credible_set_size_distribution.png", dpi=300)
                fig.savefig(figures / "primary_credible_set_size_distribution.pdf")
                plt.close(fig)
    except Exception as exc:
        print("Figure generation warning:", type(exc).__name__, exc)

    manifest = {
        "completed_at": now(),
        "genome_build": "GRCh37",
        "primary_method": "SuSiE-RSS using aligned Z scores and external 1000G EUR LD",
        "primary_sample_size_mode": "z-only asymptotic model; no N supplied",
        "sensitivity_method": "median per-variant Genomic SEM N_hat supplied to SuSiE-RSS",
        "external_LD_reference_samples": 502,
        "region_padding_bp": REGION_PADDING_BP,
        "palindromic_variants": "excluded",
        "credible_set_coverage": 0.95,
        "credible_set_min_abs_corr": 0.5,
        "max_effects_L": 10,
        "susieR_library": str(R_LIBRARY),
        "n_status_rows": int(len(status)),
        "status_counts": status.get("status", pd.Series(dtype=str)).value_counts().to_dict(),
    }
    (aggregate / "finemapping_run_manifest.json").write_text(json.dumps(manifest, indent=2))

    print("\nFine-mapping aggregation complete")
    print("Output:", OUT_ROOT)
    if not status.empty:
        print(status["status"].value_counts(dropna=False))


# =============================================================================
# Run
# =============================================================================

require_inputs()
r_helper = write_r_helper()
all_regions = []
for scenario, config in SCENARIOS.items():
    all_regions.append(run_scenario(scenario, config, r_helper))
aggregate_outputs(all_regions)
