from __future__ import annotations

import json
import os
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd


# =============================================================================
# Configuration
# =============================================================================

ROOT = Path("/path/to/panvas_analysis/Genomic_SEM_complete_followup/09_finemapping")
SOURCE_ROOT = ROOT / "susieRSS_finemapping_v1"
# Keep the all-failed v2 run untouched as an audit trail. This corrected run
# writes to a new directory so stale FAILED.txt files cannot be mistaken for
# current results.
OUT_ROOT = ROOT / "susieRSS_mismatch_aware_v3"
OUT_ROOT.mkdir(parents=True, exist_ok=True)

CONDA = Path("/path/to/panvas_analysis/miniconda3/bin/conda")
ENV_NAME = "susie016"
R_HELPER = OUT_ROOT / "run_one_region_mismatch_aware.R"

SCENARIOS = ["primary", "leave_one_CAD"]
FORCE_RERUN = False


def now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def valid_file(path: Path, minimum_bytes: int = 20) -> bool:
    return path.exists() and path.stat().st_size >= minimum_bytes


def run_logged(cmd: list[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "2"
    env["OPENBLAS_NUM_THREADS"] = "2"
    env["MKL_NUM_THREADS"] = "2"
    with log_path.open("a", encoding="utf-8") as log:
        log.write(f"\n[{now()}] $ {' '.join(map(str, cmd))}\n")
        log.flush()
        result = subprocess.run(
            [str(x) for x in cmd], stdout=log, stderr=subprocess.STDOUT,
            text=True, env=env
        )
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd)


R_CODE = r'''
args <- commandArgs(trailingOnly=TRUE)
source_dir <- args[[1]]
out_dir <- args[[2]]
scenario <- args[[3]]
region_id <- args[[4]]
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(susieR)
})

if (utils::compareVersion(as.character(packageVersion("susieR")), "0.16.0") < 0) {
  stop("susieR >= 0.16.0 is required")
}
if (!all(c("R_finite", "R_mismatch") %in% names(formals(susieR::susie_rss)))) {
  stop("This susieR build does not expose R_finite and R_mismatch")
}

model_path <- file.path(source_dir, "model_input.tsv.gz")
traw_path <- file.path(source_dir, "reference_aligned_genotypes.traw")
if (!file.exists(model_path) || !file.exists(traw_path)) {
  stop("Missing prepared v1 region inputs")
}

read_input_table <- function(path) {
  if (grepl("\\.gz$", path, ignore.case=TRUE)) {
    # Avoid an unnecessary R.utils dependency in the isolated conda environment.
    return(fread(cmd=paste("gzip -dc", shQuote(path))))
  }
  fread(path)
}

ss <- read_input_table(model_path)
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
alignment_sign <- ifelse(counted == expected, 1, ifelse(counted == other, -1, NA_real_))
keep_allele <- is.finite(alignment_sign)
gt <- gt[keep_allele]
ss <- ss[keep_allele]
alignment_sign <- alignment_sign[keep_allele]
counted <- counted[keep_allele]

G <- as.matrix(gt[, 7:ncol(gt), with=FALSE])
storage.mode(G) <- "double"
row_mean <- rowMeans(G, na.rm=TRUE)
for (i in which(rowSums(is.na(G)) > 0)) G[i, is.na(G[i, ])] <- row_mean[i]
row_sd <- apply(G, 1, sd)
keep_poly <- is.finite(row_sd) & row_sd > 0
G <- G[keep_poly, , drop=FALSE]
ss <- ss[keep_poly]
alignment_sign <- alignment_sign[keep_poly]
counted <- counted[keep_poly]
if (nrow(G) < 20) stop("Fewer than 20 polymorphic aligned variants")

X <- t(G)
X <- scale(X, center=TRUE, scale=TRUE)
R <- crossprod(X) / (nrow(X) - 1)
R <- (R + t(R)) / 2
diag(R) <- 1
z <- (as.numeric(ss$BETA) / as.numeric(ss$SE)) * alignment_sign

ok <- is.finite(z) & apply(R, 1, function(x) all(is.finite(x)))
if (!all(ok)) {
  z <- z[ok]
  R <- R[ok, ok, drop=FALSE]
  ss <- ss[ok]
  counted <- counted[ok]
  alignment_sign <- alignment_sign[ok]
}

p <- length(z)
L <- min(10L, p)
n_ref <- nrow(X)
n_hat <- suppressWarnings(round(median(as.numeric(ss$N), na.rm=TRUE)))
if (!is.finite(n_hat) || n_hat < 100) n_hat <- NA_integer_

safe_scalar <- function(x, default=NA) {
  if (is.null(x) || length(x) != 1L) return(default)
  x
}

flatten_scalars <- function(x, prefix="") {
  out <- list()
  if (is.null(x)) return(out)
  nms <- names(x)
  if (is.null(nms)) nms <- as.character(seq_along(x))
  for (i in seq_along(x)) {
    key <- if (nzchar(prefix)) paste(prefix, nms[[i]], sep=".") else nms[[i]]
    value <- x[[i]]
    if (is.atomic(value) && length(value) == 1L) {
      out[[key]] <- as.character(value)
    } else if (is.list(value) && !inherits(value, "susie")) {
      out <- c(out, flatten_scalars(value, key))
    }
  }
  out
}

extract_cs <- function(fit, method, pip_values=NULL) {
  cs <- fit$sets$cs
  if (is.null(cs) || length(cs) == 0) return(data.table())
  purity <- fit$sets$purity
  coverage <- fit$sets$coverage
  if (is.null(pip_values)) pip_values <- fit$pip
  rows <- vector("list", length(cs))
  for (i in seq_along(cs)) {
    idx <- as.integer(cs[[i]])
    nm <- if (!is.null(names(cs)) && nzchar(names(cs)[i])) names(cs)[i] else paste0("L", i)
    pur <- c(min.abs.corr=NA_real_, mean.abs.corr=NA_real_, median.abs.corr=NA_real_)
    if (!is.null(purity) && nrow(as.data.frame(purity)) >= i) {
      q <- as.numeric(purity[i, ])
      names(q) <- colnames(purity)
      common <- intersect(names(pur), names(q))
      pur[common] <- q[common]
    }
    cov_i <- if (!is.null(coverage) && length(coverage) >= i) as.numeric(coverage[i]) else NA_real_
    top_idx <- idx[which.max(pip_values[idx])]
    rows[[i]] <- data.table(
      scenario=scenario, region_id=region_id, method=method,
      cs_id=nm, cs_size=length(idx), claimed_coverage=cov_i,
      min_abs_corr=unname(pur["min.abs.corr"]),
      mean_abs_corr=unname(pur["mean.abs.corr"]),
      median_abs_corr=unname(pur["median.abs.corr"]),
      lead_SNP=as.character(ss$SNP[top_idx]),
      lead_PIP=as.numeric(pip_values[top_idx]),
      variants=paste(as.character(ss$SNP[idx]), collapse=";")
    )
  }
  rbindlist(rows, fill=TRUE)
}

fit_one <- function(method, n_value=NULL) {
  call_args <- list(
    z=z, R=R, L=L,
    R_finite=n_ref,
    R_mismatch="eb",
    estimate_residual_variance=FALSE,
    max_iter=200,
    tol=1e-3,
    coverage=0.95,
    min_abs_corr=0.5,
    refine=TRUE,
    track_fit=FALSE
  )
  if (!is.null(n_value)) call_args$n <- n_value
  refinement_used <- TRUE
  fit <- tryCatch(
    do.call(susieR::susie_rss, call_args),
    error=function(e) {
      refinement_used <<- FALSE
      call_args$refine <- FALSE
      do.call(susieR::susie_rss, call_args)
    }
  )

  diag_obj <- fit$R_finite_diagnostics
  cs <- fit$sets$cs
  membership <- rep("", p)
  if (!is.null(cs) && length(cs) > 0) {
    for (i in seq_along(cs)) {
      idx <- as.integer(cs[[i]])
      nm <- if (!is.null(names(cs)) && nzchar(names(cs)[i])) names(cs)[i] else paste0("L", i)
      membership[idx] <- vapply(membership[idx], function(old) {
        paste(c(old, nm)[nzchar(c(old, nm))], collapse=";")
      }, character(1))
    }
  }

  variants <- copy(ss)
  variants[, `:=`(
    scenario=scenario, region_id=region_id, method=method,
    Z_aligned=z, PLINK_COUNTED_ALLELE=counted,
    alignment_sign=alignment_sign, PIP=as.numeric(fit$pip),
    credible_set_ids=membership, in_credible_set=nzchar(membership)
  )]

  top <- which.max(fit$pip)
  region_summary <- data.table(
    scenario=scenario, region_id=region_id, method=method,
    n_model_variants=p, n_reference_samples=n_ref, n_hat_median=n_hat,
    L=L, converged=isTRUE(fit$converged), niter=as.integer(fit$niter),
    refinement_used=refinement_used,
    n_credible_sets=if (is.null(cs)) 0L else length(cs),
    max_PIP=as.numeric(fit$pip[top]), top_SNP=as.character(ss$SNP[top]),
    top_P=as.numeric(ss$P[top]),
    Q_art=as.numeric(safe_scalar(diag_obj$Q_art, NA_real_)),
    artifact_flag=as.logical(safe_scalar(diag_obj$artifact_flag, NA)),
    R_sensitivity_flag=as.logical(safe_scalar(diag_obj$R_sensitivity_flag, NA)),
    R_reliability_flag=as.logical(safe_scalar(diag_obj$R_reliability_flag, NA)),
    B=as.numeric(safe_scalar(diag_obj$B, NA_real_)),
    effective_rank=as.numeric(safe_scalar(diag_obj$effective_rank, NA_real_)),
    r_over_B=as.numeric(safe_scalar(diag_obj$r_over_B, NA_real_)),
    lambda_bias=as.numeric(safe_scalar(diag_obj$lambda_bias, NA_real_)),
    B_corrected=as.numeric(safe_scalar(diag_obj$B_corrected, NA_real_)),
    susieR_version=as.character(packageVersion("susieR"))
  )

  diag_flat <- flatten_scalars(diag_obj)
  diag_table <- if (length(diag_flat)) data.table(
    scenario=scenario, region_id=region_id, method=method,
    metric=names(diag_flat), value=unlist(diag_flat, use.names=FALSE)
  ) else data.table()

  ser_variants <- data.table()
  ser_cs <- data.table()
  ser <- diag_obj$ser_model
  if (!is.null(ser)) {
    ser_pip <- if (!is.null(ser$pip)) as.numeric(ser$pip) else as.numeric(ser$alpha[1, ])
    ser_membership <- rep("", p)
    if (!is.null(ser$sets$cs) && length(ser$sets$cs) > 0) {
      for (i in seq_along(ser$sets$cs)) {
        idx <- as.integer(ser$sets$cs[[i]])
        nm <- if (!is.null(names(ser$sets$cs)) && nzchar(names(ser$sets$cs)[i])) names(ser$sets$cs)[i] else paste0("SER", i)
        ser_membership[idx] <- nm
      }
    }
    ser_variants <- copy(ss)
    ser_variants[, `:=`(
      scenario=scenario, region_id=region_id, method=paste0(method, "__SER_diagnostic"),
      PIP=ser_pip, credible_set_ids=ser_membership,
      in_credible_set=nzchar(ser_membership)
    )]
    ser_cs <- extract_cs(ser, paste0(method, "__SER_diagnostic"), ser_pip)
  }

  list(
    fit=fit, variants=variants, credible_sets=extract_cs(fit, method),
    summary=region_summary, diagnostics=diag_table,
    ser_variants=ser_variants, ser_credible_sets=ser_cs
  )
}

methods <- list(finite502_EB_z_only=list(n=NULL))
if (!is.na(n_hat)) methods$finite502_EB_median_Nhat <- list(n=n_hat)

fits <- list()
errors <- list()
for (method in names(methods)) {
  fits[[method]] <- tryCatch(
    fit_one(method, methods[[method]]$n),
    error=function(e) {
      errors[[method]] <<- conditionMessage(e)
      NULL
    }
  )
}

if (!length(fits) || all(vapply(fits, is.null, logical(1)))) {
  stop("All mismatch-aware methods failed: ", paste(unlist(errors), collapse=" | "))
}
fits <- fits[!vapply(fits, is.null, logical(1))]

fwrite(rbindlist(lapply(fits, `[[`, "variants"), fill=TRUE),
       file.path(out_dir, "variant_PIP.tsv.gz"), sep="\t")
fwrite(rbindlist(lapply(fits, `[[`, "credible_sets"), fill=TRUE),
       file.path(out_dir, "credible_sets.tsv"), sep="\t")
fwrite(rbindlist(lapply(fits, `[[`, "summary"), fill=TRUE),
       file.path(out_dir, "region_summary.tsv"), sep="\t")
fwrite(rbindlist(lapply(fits, `[[`, "diagnostics"), fill=TRUE),
       file.path(out_dir, "R_finite_diagnostics_flat.tsv"), sep="\t")
fwrite(rbindlist(lapply(fits, `[[`, "ser_variants"), fill=TRUE),
       file.path(out_dir, "SER_diagnostic_variant_PIP.tsv.gz"), sep="\t")
fwrite(rbindlist(lapply(fits, `[[`, "ser_credible_sets"), fill=TRUE),
       file.path(out_dir, "SER_diagnostic_credible_sets.tsv"), sep="\t")
saveRDS(lapply(fits, `[[`, "fit"), file.path(out_dir, "susie_mismatch_aware_fits.rds"), compress=TRUE)
fwrite(data.table(method=names(errors), error=unlist(errors, use.names=FALSE)),
       file.path(out_dir, "method_errors.tsv"), sep="\t")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo.txt"))
'''

R_HELPER.write_text(R_CODE, encoding="utf-8")


# =============================================================================
# Environment and input checks
# =============================================================================

if not CONDA.exists():
    raise FileNotFoundError(CONDA)
if not SOURCE_ROOT.exists():
    raise FileNotFoundError(SOURCE_ROOT)

check_expr = (
    'library(susieR); '
    'cat(as.character(packageVersion("susieR")), "\\n"); '
    'cat(all(c("R_finite","R_mismatch") %in% names(formals(susieR::susie_rss))), "\\n")'
)
check = subprocess.run(
    [str(CONDA), "run", "-n", ENV_NAME, "Rscript", "--vanilla", "-e", check_expr],
    capture_output=True, text=True, check=True,
)
print("susieR environment audit:\n", check.stdout, flush=True)
if "TRUE" not in check.stdout:
    raise RuntimeError("Run 09d installer first; mismatch-aware susieR is unavailable")


# =============================================================================
# One-region smoke test before launching the full batch
# =============================================================================

smoke_source = sorted(
    p for p in (SOURCE_ROOT / "primary" / "regions").iterdir() if p.is_dir()
)[0]
smoke_out = OUT_ROOT / "primary" / "regions" / smoke_source.name
smoke_out.mkdir(parents=True, exist_ok=True)
smoke_summary = smoke_out / "region_summary.tsv"

if not valid_file(smoke_summary):
    smoke_log = smoke_out / "mismatch_aware_fit.log"
    print(f"[{now()}] Running one-region smoke test: {smoke_source.name}", flush=True)
    try:
        run_logged([
            CONDA, "run", "--no-capture-output", "-n", ENV_NAME,
            "Rscript", "--vanilla", R_HELPER,
            smoke_source, smoke_out, "primary", smoke_source.name,
        ], smoke_log)
    except Exception as exc:
        tail = ""
        if smoke_log.exists():
            tail = "\n".join(
                smoke_log.read_text(encoding="utf-8", errors="replace").splitlines()[-120:]
            )
        raise RuntimeError(
            "The one-region smoke test failed; the full batch was not started.\n"
            f"{exc}\n\nLast log lines:\n{tail}"
        ) from exc

if not valid_file(smoke_summary):
    raise RuntimeError("Smoke test returned without a valid region_summary.tsv")

print(f"[{now()}] Smoke test passed: {smoke_summary}", flush=True)


# =============================================================================
# Refit all prepared regions
# =============================================================================

status_rows = []
for scenario in SCENARIOS:
    source_regions = SOURCE_ROOT / scenario / "regions"
    if not source_regions.exists():
        raise FileNotFoundError(source_regions)
    region_dirs = sorted(p for p in source_regions.iterdir() if p.is_dir())
    print(f"\n[{now()}] {scenario}: {len(region_dirs)} prepared regions", flush=True)

    for i, source_dir in enumerate(region_dirs, 1):
        region_id = source_dir.name
        out_dir = OUT_ROOT / scenario / "regions" / region_id
        out_dir.mkdir(parents=True, exist_ok=True)
        summary_path = out_dir / "region_summary.tsv"

        if not FORCE_RERUN and valid_file(summary_path):
            state = "SKIPPED_COMPLETED"
            error_text = ""
        else:
            try:
                run_logged([
                    CONDA, "run", "--no-capture-output", "-n", ENV_NAME,
                    "Rscript", "--vanilla", R_HELPER,
                    source_dir, out_dir, scenario, region_id,
                ], out_dir / "mismatch_aware_fit.log")
                state = "SUCCESS"
                error_text = ""
            except Exception as exc:
                state = "FAILED"
                error_text = str(exc)
                (out_dir / "FAILED.txt").write_text(error_text, encoding="utf-8")

        status_rows.append({
            "scenario": scenario,
            "region_id": region_id,
            "status": state,
            "region_summary_exists": valid_file(summary_path),
            "error": error_text,
        })
        if i == 1 or i % 10 == 0 or i == len(region_dirs):
            print(f"[{now()}] {scenario}: {i}/{len(region_dirs)} {region_id} {state}", flush=True)


# =============================================================================
# Aggregate outputs and define analysis-ready results
# =============================================================================

aggregate = OUT_ROOT / "aggregate"
aggregate.mkdir(parents=True, exist_ok=True)
status = pd.DataFrame(status_rows)
status.to_csv(aggregate / "all_region_status.tsv", sep="\t", index=False)


def collect(filename: str, minimum_bytes: int = 1) -> pd.DataFrame:
    tables = []
    for scenario in SCENARIOS:
        for path in sorted((OUT_ROOT / scenario / "regions").glob(f"*/{filename}")):
            if path.stat().st_size >= minimum_bytes:
                try:
                    tables.append(pd.read_csv(path, sep="\t"))
                except pd.errors.EmptyDataError:
                    pass
    return pd.concat(tables, ignore_index=True, sort=False) if tables else pd.DataFrame()


region_summary = collect("region_summary.tsv")
credible_sets = collect("credible_sets.tsv")
diagnostics = collect("R_finite_diagnostics_flat.tsv")
variant_pip = collect("variant_PIP.tsv.gz")
ser_variant_pip = collect("SER_diagnostic_variant_PIP.tsv.gz")
ser_cs = collect("SER_diagnostic_credible_sets.tsv")

region_summary.to_csv(aggregate / "all_region_summary.tsv", sep="\t", index=False)
credible_sets.to_csv(aggregate / "all_credible_sets.tsv", sep="\t", index=False)
diagnostics.to_csv(aggregate / "all_R_finite_diagnostics_flat.tsv", sep="\t", index=False)
variant_pip.to_csv(aggregate / "all_variant_PIP.tsv.gz", sep="\t", index=False, compression="gzip")
ser_variant_pip.to_csv(aggregate / "all_SER_diagnostic_variant_PIP.tsv.gz", sep="\t", index=False, compression="gzip")
ser_cs.to_csv(aggregate / "all_SER_diagnostic_credible_sets.tsv", sep="\t", index=False)

if not region_summary.empty:
    region_summary["analysis_ready_region"] = (
        region_summary["converged"].fillna(False).astype(bool)
        & region_summary["R_reliability_flag"].eq(False)
    )
    region_summary.to_csv(aggregate / "all_region_summary_with_acceptance.tsv", sep="\t", index=False)

if not credible_sets.empty and not region_summary.empty:
    accepted = region_summary[["scenario", "region_id", "method", "analysis_ready_region"]]
    cs_review = credible_sets.merge(accepted, on=["scenario", "region_id", "method"], how="left")
    cs_review["credible_set_analysis_ready"] = (
        cs_review["analysis_ready_region"].fillna(False)
        & cs_review["min_abs_corr"].ge(0.5)
    )
    cs_review.to_csv(aggregate / "credible_sets_with_acceptance.tsv", sep="\t", index=False)

if not variant_pip.empty:
    main_method = "finite502_EB_z_only"
    high = variant_pip[(variant_pip["method"] == main_method) & (variant_pip["PIP"] >= 0.10)].copy()
    high.to_csv(aggregate / "mismatch_corrected_variants_PIP_ge_0.10.tsv", sep="\t", index=False)

    p = variant_pip[(variant_pip["scenario"] == "primary") & (variant_pip["method"] == main_method)]
    l = variant_pip[(variant_pip["scenario"] == "leave_one_CAD") & (variant_pip["method"] == main_method)]
    p = p[["SNP", "CHR", "BP", "PIP", "credible_set_ids"]].rename(
        columns={"PIP": "PIP_primary", "credible_set_ids": "CS_primary"}
    )
    l = l[["SNP", "CHR", "BP", "PIP", "credible_set_ids"]].rename(
        columns={"PIP": "PIP_leave_one_CAD", "credible_set_ids": "CS_leave_one_CAD"}
    )
    comparison = p.merge(l, on=["SNP", "CHR", "BP"], how="inner")
    comparison["min_cross_CAD_PIP"] = comparison[["PIP_primary", "PIP_leave_one_CAD"]].min(axis=1)
    comparison.to_csv(
        aggregate / "primary_vs_leave_one_CAD_mismatch_corrected_PIP.tsv.gz",
        sep="\t", index=False, compression="gzip"
    )

manifest = {
    "completed_at": now(),
    "source_prepared_regions": str(SOURCE_ROOT),
    "output": str(OUT_ROOT),
    "susieR_mode": "R_finite=reference sample count; R_mismatch=eb",
    "external_LD_reference_samples_expected": 502,
    "primary_method": "finite502_EB_z_only",
    "sensitivity_method": "finite502_EB_median_Nhat",
    "acceptance_rule": "converged AND R_reliability_flag is FALSE; CS purity >= 0.5",
    "status_counts": status["status"].value_counts(dropna=False).to_dict(),
}
(aggregate / "mismatch_aware_run_manifest.json").write_text(
    json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
)

print("\nCompleted mismatch-aware fine-mapping", flush=True)
print(status.groupby(["scenario", "status"]).size(), flush=True)
if not region_summary.empty:
    print(region_summary.groupby(["scenario", "method"])[
        ["converged", "R_reliability_flag", "n_credible_sets"]
    ].agg({"converged": "sum", "R_reliability_flag": "sum", "n_credible_sets": "sum"}), flush=True)
print("Output:", OUT_ROOT, flush=True)
