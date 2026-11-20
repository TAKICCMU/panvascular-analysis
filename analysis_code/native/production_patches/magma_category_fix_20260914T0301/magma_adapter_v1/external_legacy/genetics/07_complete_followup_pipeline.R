#!/usr/bin/env Rscript

# Complete follow-up pipeline for the pan-vascular Genomic SEM analysis.
#
# Required server inputs:
#   /path/to/panvas_analysis/Genomic_SEM/LDSCoutput.RData
#   /path/to/panvas_analysis/Genomic_SEM/PSYCH_sumstats.RData
#   /path/to/panvas_analysis/Genomic_SEM/gwasfile_info.txt
#   /path/to/panvas_analysis/Genomic_SEM/eur_w_ld_chr/
#   /path/to/panvas_analysis/Genomic_SEM/*.sumstats.gz
#
# Main outputs:
#   /path/to/panvas_analysis/Genomic_SEM_complete_followup/
#
# The original files are read only. This script does not delete, move, or
# overwrite anything under /path/to/panvas_analysis/Genomic_SEM.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(GenomicSEM)
  library(data.table)
  library(Matrix)
  library(psych)
  library(parallel)
})

# ==============================================================================
# 0. Configuration
# ==============================================================================

ROOT <- "/path/to/panvas_analysis/Genomic_SEM"
OUT_ROOT <- "/path/to/panvas_analysis/Genomic_SEM_complete_followup"

# These analyses are required for the current manuscript revision.
RUN_LDSC_AXIS_ANALYSIS <- TRUE
RUN_FULL_CFA_COMPARISON <- TRUE
RUN_ALL_LEAVE_ONE_CFA <- TRUE
RUN_ODD_EVEN_LDSC <- TRUE
RUN_FULL_BIFACTOR_GWAS <- TRUE
RUN_LEAVE_ONE_CAD_GWAS <- TRUE

# Not required before assessing the full and leave-one-CAD results.
# Set TRUE only if all eight SNP-level leave-one-trait analyses are desired.
RUN_ALL_LEAVE_ONE_GWAS <- FALSE

# Every requested GWAS first runs a deterministic pilot and proceeds to the
# complete SNP set only if the pilot passes basic output checks.
RUN_PILOT_FIRST <- TRUE
RUN_FULL_AFTER_PILOT <- TRUE
PILOT_N <- 20000L
PILOT_SEED <- 20260729L

# "chunked" reproduces the fast strategy used by the original userGWAS.Rc:
# split SNPs into batches, run official GenomicSEM::userGWAS on every batch
# with parallel workers, and checkpoint each batch. "direct" submits all SNPs
# in one userGWAS call.
USERGWAS_ENGINE <- "chunked"
USERGWAS_CHUNK_SIZE <- 10000L
USERGWAS_CHUNK_MAX_RETRIES <- 1L

# Large raw userGWAS R objects duplicate the compressed result tables.
SAVE_RAW_USERGWAS_RDS <- FALSE
OVERWRITE_EXISTING_GWAS <- FALSE

N_CORES <- max(1L, floor(parallel::detectCores() / 2))
GWAS_COMPARE_N <- 500000L
GWS_P <- 5e-8
QSNP_GWS_P <- 5e-8

EXPECTED_TRAITS <- c(
  "AA", "Stroke", "AAA", "CAVS",
  "CAD", "TAAD", "PAD", "CA"
)
AORTIC_TRAITS <- c("AA", "AAA", "TAAD")
ATHERO_TRAITS <- c("CAD", "Stroke", "PAD", "CA", "CAVS")

# Confirmed from ASCVD_SEM.R and the N-distribution audit:
# these five traits contain SNP-specific effective sample sizes inferred from
# MAF and logistic SE and bounded at 0.5-1.1 times pooled effective N.
EFFECTIVE_N_TRAITS <- c("Stroke", "CAVS", "CAD", "PAD", "CA")

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(
  OUT_ROOT,
  paste0("pipeline_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)

# ==============================================================================
# 1. General helpers
# ==============================================================================

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

log_message <- function(...) {
  cat(
    "[", timestamp(), "] ",
    paste0(..., collapse = ""),
    "\n",
    sep = ""
  )
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}

write_tsv_gz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(
    as.data.table(x),
    path,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = "gzip"
  )
}

# Compatible with older data.table versions in which
# as.data.table.matrix(..., keep.rownames=) may fail to create the internal
# row-name column before renaming it.
matrix_to_table <- function(x, row_name = "trait") {
  x <- as.matrix(x)
  rn <- rownames(x)
  if (is.null(rn) && nrow(x) == ncol(x) && !is.null(colnames(x))) {
    rn <- colnames(x)
  }
  if (is.null(rn)) rn <- as.character(seq_len(nrow(x)))

  out <- as.data.table(x)
  out[, (row_name) := rn]
  setcolorder(out, c(row_name, setdiff(names(out), row_name)))
  out
}

load_rdata_object <- function(path, object_name) {
  if (!file.exists(path)) stop("Missing input: ", path)
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists(object_name, envir = env, inherits = FALSE)) {
    stop("Object ", object_name, " was not found in ", path)
  }
  env[[object_name]]
}

as_numeric_safe <- function(x) suppressWarnings(as.numeric(as.character(x)))

pick_col <- function(tab, candidates, required = TRUE) {
  hit <- candidates[candidates %in% names(tab)]
  if (length(hit) > 0L) return(hit[1])
  if (required) {
    stop("Missing output column: ", paste(candidates, collapse = "/"))
  }
  NA_character_
}

extract_metric <- function(tab, candidates) {
  hit <- candidates[candidates %in% names(tab)]
  if (length(hit) == 0L || nrow(tab) == 0L) return(NA_real_)
  as_numeric_safe(tab[[hit[1]]][1])
}

pair_index_table <- function(nm) {
  rbindlist(lapply(seq_along(nm), function(i) {
    data.table(
      i = i,
      j = i:length(nm),
      trait1 = nm[i],
      trait2 = nm[i:length(nm)]
    )
  }))
}

subset_ldsc <- function(x, keep_traits) {
  old_traits <- colnames(x$S)
  keep <- match(keep_traits, old_traits)
  if (anyNA(keep)) {
    stop("Traits absent from LDSC object: ",
         paste(keep_traits[is.na(keep)], collapse = ", "))
  }

  pairs <- pair_index_table(old_traits)
  pair_keep <- which(pairs$i %in% keep & pairs$j %in% keep)

  out <- x
  out$S <- x$S[keep, keep, drop = FALSE]
  out$V <- x$V[pair_keep, pair_keep, drop = FALSE]
  if (!is.null(x$I)) out$I <- x$I[keep, keep, drop = FALSE]
  if (!is.null(x$N)) out$N <- x$N[, pair_keep, drop = FALSE]
  if (!is.null(x$S_Stand)) {
    out$S_Stand <- x$S_Stand[keep, keep, drop = FALSE]
  }
  if (!is.null(x$V_Stand)) {
    out$V_Stand <- x$V_Stand[pair_keep, pair_keep, drop = FALSE]
  }
  out
}

factor_line <- function(name, indicators) {
  if (length(indicators) < 2L) {
    stop(name, " has fewer than two indicators.")
  }
  paste0(
    name, " =~ NA*", indicators[1],
    if (length(indicators) > 1L) {
      paste0(" + ", paste(indicators[-1], collapse = " + "))
    } else {
      ""
    }
  )
}

marker_factor_line <- function(name, indicators, preferred_anchors) {
  if (length(indicators) < 2L) {
    stop(name, " has fewer than two indicators.")
  }
  anchor <- preferred_anchors[preferred_anchors %in% indicators][1]
  if (is.na(anchor) || length(anchor) == 0L) anchor <- indicators[1]
  ordered <- c(anchor, setdiff(indicators, anchor))
  list(
    syntax = paste0(
      name, " =~ 1*", ordered[1],
      if (length(ordered) > 1L) {
        paste0(" + ", paste(ordered[-1], collapse = " + "))
      } else {
        ""
      }
    ),
    anchor = anchor
  )
}

make_single_factor_model <- function(keep_traits) {
  factor_line("G", keep_traits)
}

make_correlated_model <- function(keep_traits) {
  aortic <- intersect(AORTIC_TRAITS, keep_traits)
  athero <- intersect(ATHERO_TRAITS, keep_traits)

  lines <- c(
    factor_line("AORTIC", aortic),
    factor_line("ATHERO", athero),
    "AORTIC ~~ ATHERO"
  )
  if ("AAA" %in% keep_traits) {
    lines <- c(lines, "AAA ~~ 1e-6*AAA")
  }
  paste(lines, collapse = "\n")
}

make_bifactor_model <- function(keep_traits, include_snp = FALSE) {
  aortic <- intersect(AORTIC_TRAITS, keep_traits)
  athero <- intersect(ATHERO_TRAITS, keep_traits)

  lines <- c(
    factor_line("G", keep_traits),
    factor_line("AORTIC", aortic),
    factor_line("ATHERO", athero),
    "G ~~ 0*AORTIC",
    "G ~~ 0*ATHERO",
    "AORTIC ~~ 0*ATHERO"
  )
  if ("AAA" %in% keep_traits) {
    lines <- c(lines, "AAA ~~ 1e-6*AAA")
  }
  if (include_snp) lines <- c(lines, "G ~ SNP")
  paste(lines, collapse = "\n")
}

make_bifactor_gwas_spec <- function(keep_traits, include_snp = FALSE) {
  aortic <- intersect(AORTIC_TRAITS, keep_traits)
  athero <- intersect(ATHERO_TRAITS, keep_traits)

  general_line <- marker_factor_line(
    "G", keep_traits,
    preferred_anchors = c("Stroke", "CAD", "PAD", "CA", "CAVS")
  )
  aortic_line <- marker_factor_line(
    "AORTIC", aortic,
    preferred_anchors = c("AAA", "AA", "TAAD")
  )
  athero_line <- marker_factor_line(
    "ATHERO", athero,
    preferred_anchors = c("PAD", "CAD", "Stroke", "CA", "CAVS")
  )

  lines <- c(
    general_line$syntax,
    aortic_line$syntax,
    athero_line$syntax,
    "G ~~ G",
    "AORTIC ~~ AORTIC",
    "ATHERO ~~ ATHERO",
    "G ~~ 0*AORTIC",
    "G ~~ 0*ATHERO",
    "AORTIC ~~ 0*ATHERO"
  )
  if ("AAA" %in% keep_traits) {
    lines <- c(lines, "AAA ~~ 1e-6*AAA")
  }
  if (include_snp) lines <- c(lines, "G ~ SNP")

  list(
    syntax = paste(lines, collapse = "\n"),
    G_anchor = general_line$anchor,
    AORTIC_anchor = aortic_line$anchor,
    ATHERO_anchor = athero_line$anchor
  )
}

# Identification-constrained higher-order negative-control model.
# With only F1 and F2 as second-order indicators, a freely estimated
# higher-order G is not identified. The equal-loading constraint below makes
# the model estimable, but it is not the primary structural model and must not
# be used for the pan-vascular GWAS when its fit is unacceptable.
equal_loading_higher_order_negative_control_model <- paste(
  "F1 =~ 1*AA + l11*AAA + l12*TAAD",
  "F2 =~ 1*CAD + l21*Stroke + l22*PAD + l23*CA + l24*CAVS",
  "G =~ a*F1 + a*F2",
  "G ~~ 1*G",
  "F1 ~~ u1*F1",
  "F2 ~~ u2*F2",
  "F1 ~~ 0*F2",
  "AAA ~~ 1e-6*AAA",
  sep = "\n"
)

run_cfa <- function(covstruc, model, std_lv = TRUE) {
  warnings <- character()
  error_text <- NA_character_

  fit <- tryCatch(
    withCallingHandlers(
      usermodel(
        covstruc = covstruc,
        estimation = "DWLS",
        model = model,
        CFIcalc = TRUE,
        std.lv = std_lv,
        imp_cov = TRUE,
        fix_resid = FALSE,
        toler = 1e-20
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_text <<- conditionMessage(e)
      NULL
    }
  )

  if (is.null(fit)) {
    return(list(
      fit_object = NULL,
      fit = data.table(),
      parameters = data.table(),
      warnings = unique(warnings),
      error = error_text
    ))
  }

  list(
    fit_object = fit,
    fit = as.data.table(fit$modelfit),
    parameters = as.data.table(fit$results),
    warnings = unique(warnings),
    error = error_text
  )
}

diagnose_cfa <- function(result, observed_traits) {
  warning_text <- paste(result$warnings, collapse = " | ")
  info_warning <- grepl(
    "information matrix|not identified|could not compute standard errors",
    warning_text,
    ignore.case = TRUE
  )

  if (is.null(result$fit_object)) {
    return(data.table(
      status = "ERROR",
      error = result$error,
      warnings = warning_text,
      CFI = NA_real_,
      SRMR = NA_real_,
      information_matrix_warning = info_warning,
      n_negative_observed_residuals = NA_integer_,
      negative_residual_traits = NA_character_,
      accepted_for_gwas = FALSE
    ))
  }

  fit <- result$fit
  pars <- copy(result$parameters)
  est_col <- pick_col(
    pars,
    c("Unstand_Est", "est", "Estimate"),
    required = FALSE
  )
  residuals <- data.table()
  if (!is.na(est_col)) {
    pars[, estimate_numeric := as_numeric_safe(get(est_col))]
    residuals <- pars[
      op == "~~" & lhs == rhs & lhs %in% observed_traits &
        is.finite(estimate_numeric) & estimate_numeric < -1e-8
    ]
  }

  cfi <- extract_metric(fit, c("CFI", "cfi"))
  srmr <- extract_metric(fit, c("SRMR", "srmr"))
  fit_ok <- is.finite(cfi) && is.finite(srmr) &&
    cfi >= 0.90 && srmr <= 0.10
  accepted <- fit_ok && !info_warning && nrow(residuals) == 0L

  data.table(
    status = "SUCCESS",
    error = NA_character_,
    warnings = warning_text,
    CFI = cfi,
    SRMR = srmr,
    information_matrix_warning = info_warning,
    n_negative_observed_residuals = nrow(residuals),
    negative_residual_traits = paste(residuals$lhs, collapse = ";"),
    accepted_for_gwas = accepted
  )
}

save_cfa_result <- function(result, model, out_dir, label, observed_traits) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(model, file.path(out_dir, "model_syntax.txt"))
  writeLines(result$warnings, file.path(out_dir, "warnings.txt"))
  if (!is.null(result$fit_object)) {
    write_tsv(result$fit, file.path(out_dir, "model_fit.tsv"))
    write_tsv(result$parameters, file.path(out_dir, "model_parameters.tsv"))
    saveRDS(result$fit_object, file.path(out_dir, "CFA_fit.rds"))
  }
  diagnostic <- diagnose_cfa(result, observed_traits)
  diagnostic[, analysis := label]
  write_tsv(diagnostic, file.path(out_dir, "model_diagnostics.tsv"))
  diagnostic
}

drop_trait_sumstats <- function(sumstats, omitted_trait) {
  if (is.na(omitted_trait) || omitted_trait == "NONE") return(sumstats)

  candidates <- c(
    paste0("beta.", omitted_trait),
    paste0("se.", omitted_trait),
    paste0("beta_", omitted_trait),
    paste0("se_", omitted_trait)
  )
  drop_cols <- intersect(candidates, names(sumstats))
  if (length(drop_cols) != 2L) {
    stop(
      "Expected two beta/se columns for ", omitted_trait,
      "; found: ", paste(drop_cols, collapse = ", ")
    )
  }
  sumstats[, setdiff(names(sumstats), drop_cols), with = FALSE]
}

result_to_table <- function(result) {
  if (is.data.frame(result)) return(as.data.table(result))
  if (is.list(result)) return(rbindlist(result, fill = TRUE))
  stop("Unsupported userGWAS output class: ", paste(class(result), collapse = "/"))
}

standardize_gwas_output <- function(tab) {
  tab <- as.data.table(tab)
  if (all(c("lhs", "op", "rhs") %in% names(tab))) {
    factor_rows <- tab[lhs == "G" & op == "~" & rhs == "SNP"]
    if (nrow(factor_rows) > 0L) tab <- factor_rows
  }

  snp_col <- pick_col(tab, c("SNP", "rsid", "RSID"))
  chr_col <- pick_col(tab, c("CHR", "chr"))
  bp_col <- pick_col(tab, c("BP", "POS", "position"))
  a1_col <- pick_col(tab, c("A1", "EA", "effect_allele"))
  a2_col <- pick_col(tab, c("A2", "NEA", "other_allele"))
  maf_col <- pick_col(tab, c("MAF", "maf"))
  beta_col <- pick_col(tab, c("est", "Unstand_Est", "BETA", "beta"))
  se_col <- pick_col(tab, c("se_c", "SE", "Unstand_SE", "se"))
  p_col <- pick_col(tab, c("Pval_Estimate", "p_value", "P", "p"))
  q_col <- pick_col(tab, c("Q_SNP", "Q"), required = FALSE)
  q_p_col <- pick_col(
    tab,
    c("Q_SNP_pval", "Q_SNP_P", "Q_p", "Q_SNP_p"),
    required = FALSE
  )

  out <- data.table(
    SNP = as.character(tab[[snp_col]]),
    CHR = as_numeric_safe(tab[[chr_col]]),
    BP = as_numeric_safe(tab[[bp_col]]),
    A1 = as.character(tab[[a1_col]]),
    A2 = as.character(tab[[a2_col]]),
    MAF = as_numeric_safe(tab[[maf_col]]),
    BETA = as_numeric_safe(tab[[beta_col]]),
    SE = as_numeric_safe(tab[[se_col]]),
    P = as_numeric_safe(tab[[p_col]])
  )
  out[, Q_SNP := if (!is.na(q_col)) as_numeric_safe(tab[[q_col]]) else NA_real_]
  out[, Q_SNP_P := if (!is.na(q_p_col)) {
    as_numeric_safe(tab[[q_p_col]])
  } else {
    NA_real_
  }]
  out
}

calculate_factor_n <- function(gwas) {
  valid <- is.finite(gwas$MAF) & gwas$MAF >= 0.10 & gwas$MAF <= 0.50 &
    is.finite(gwas$SE) & gwas$SE > 0
  implied_n <- 1 / (
    2 * gwas$MAF[valid] * (1 - gwas$MAF[valid]) * gwas$SE[valid]^2
  )
  implied_n <- implied_n[is.finite(implied_n) & implied_n > 0]
  if (length(implied_n) == 0L) stop("No valid SNPs for factor N_hat.")

  data.table(
    N_hat_mean = mean(implied_n),
    N_hat_median = median(implied_n),
    N_hat_trimmed_mean_01 = mean(implied_n, trim = 0.01),
    N_hat_q025 = quantile(implied_n, 0.025),
    N_hat_q975 = quantile(implied_n, 0.975),
    n_snps_for_N_hat = length(implied_n)
  )
}

write_fuma_outputs <- function(gwas, out_dir, label) {
  n_stats <- calculate_factor_n(gwas)
  n_hat <- n_stats$N_hat_mean[1]
  write_tsv(n_stats, file.path(out_dir, paste0(label, "_N_hat.tsv")))

  valid <- is.finite(gwas$P) & gwas$P > 0 & gwas$P <= 1 &
    is.finite(gwas$BETA) & is.finite(gwas$SE) & gwas$SE > 0 &
    is.finite(gwas$CHR) & is.finite(gwas$BP) &
    !is.na(gwas$SNP) & nzchar(gwas$SNP)

  fuma <- gwas[valid, .(
    SNP, CHR = as.integer(CHR), BP = as.integer(BP),
    A1, A2, MAF, BETA, SE, P, N = n_hat
  )]

  write_tsv_gz(
    fuma,
    file.path(out_dir, paste0(label, "_FUMA_all.tsv.gz"))
  )
  write_tsv(
    fuma[P < GWS_P],
    file.path(out_dir, paste0(label, "_genomewide_significant_SNPs.tsv"))
  )

  q_available <- any(is.finite(gwas$Q_SNP_P))
  if (q_available) {
    q_sig_snps <- gwas[
      is.finite(Q_SNP_P) & Q_SNP_P < QSNP_GWS_P,
      unique(SNP)
    ]
    write_tsv(
      gwas[SNP %in% q_sig_snps],
      file.path(out_dir, paste0(label, "_QSNP_genomewide_significant.tsv"))
    )
    write_tsv_gz(
      fuma[!SNP %in% q_sig_snps],
      file.path(
        out_dir,
        paste0(label, "_FUMA_excluding_QSNP_GWS.tsv.gz")
      )
    )
  }

  data.table(
    analysis = label,
    n_total_rows = nrow(gwas),
    n_valid_fuma_rows = nrow(fuma),
    n_GWS = sum(fuma$P < GWS_P, na.rm = TRUE),
    n_QSNP_GWS = if (q_available) {
      sum(gwas$Q_SNP_P < QSNP_GWS_P, na.rm = TRUE)
    } else {
      NA_integer_
    },
    N_hat = n_hat
  )
}

# ==============================================================================
# 2. Input validation and confirmed prevalence settings
# ==============================================================================

INDEX_PATH <- file.path(ROOT, "gwasfile_info.txt")
LDSC_PATH <- file.path(ROOT, "LDSCoutput.RData")
SUMSTATS_PATH <- file.path(ROOT, "PSYCH_sumstats.RData")

index <- fread(INDEX_PATH)
if ("case" %in% names(index) && !"cases" %in% names(index)) {
  setnames(index, "case", "cases")
}
required_index_cols <- c(
  "phenotype", "cases", "controls", "population.prev"
)
missing_index <- setdiff(required_index_cols, names(index))
if (length(missing_index) > 0L) {
  stop("Missing gwasfile_info columns: ", paste(missing_index, collapse = ", "))
}

traits <- index$phenotype
if (!identical(traits, EXPECTED_TRAITS)) {
  stop(
    "Trait order differs from the validated order.\nObserved: ",
    paste(traits, collapse = ", "),
    "\nExpected: ", paste(EXPECTED_TRAITS, collapse = ", ")
  )
}

sample_prev_counts <- index$cases / (index$cases + index$controls)
sample_prev <- ifelse(
  traits %in% EFFECTIVE_N_TRAITS,
  0.5,
  sample_prev_counts
)
population_prev <- index$population.prev

prevalence_audit <- data.table(
  phenotype = traits,
  cases = index$cases,
  controls = index$controls,
  N_type = ifelse(
    traits %in% EFFECTIVE_N_TRAITS,
    "estimated_snp_specific_effective",
    "total"
  ),
  sample_prev = sample_prev,
  population_prev = population_prev
)
write_tsv(
  prevalence_audit,
  file.path(OUT_ROOT, "00_confirmed_sample_size_and_prevalence.tsv")
)

LDSCoutput <- load_rdata_object(LDSC_PATH, "LDSCoutput")
if (!identical(colnames(LDSCoutput$S), traits)) {
  stop("LDSC trait order does not match gwasfile_info.txt.")
}

write_tsv(
  matrix_to_table(LDSCoutput$S, row_name = "trait"),
  file.path(OUT_ROOT, "00_LDSC_genetic_covariance.tsv")
)
write_tsv(
  matrix_to_table(LDSCoutput$S_Stand, row_name = "trait"),
  file.path(OUT_ROOT, "00_LDSC_genetic_correlation.tsv")
)

# ==============================================================================
# 3. LDSC two-axis structure analysis
# ==============================================================================

run_ldsc_axis_analysis <- function() {
  out <- file.path(OUT_ROOT, "01_ldsc_two_axis")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  R <- LDSCoutput$S_Stand
  VR <- LDSCoutput$V_Stand
  pairs <- pair_index_table(traits)
  pairs[, pair_index := .I]
  pairs[, rg := R[cbind(i, j)]]
  pairs[, rg_se := sqrt(diag(VR))]
  pairs[, rg_z := rg / rg_se]
  pairs[, rg_p := 2 * pnorm(abs(rg_z), lower.tail = FALSE)]
  pairs[, block := fifelse(
    trait1 %in% AORTIC_TRAITS & trait2 %in% AORTIC_TRAITS,
    "within_aortic",
    fifelse(
      trait1 %in% ATHERO_TRAITS & trait2 %in% ATHERO_TRAITS,
      "within_atherosclerosis_related",
      "between_axes"
    )
  )]
  pairs[i == j, block := "diagonal"]
  write_tsv(pairs, file.path(out, "pairwise_rg_with_axis_labels.tsv"))

  off <- pairs[i != j]
  off_idx <- off$pair_index
  y <- off$rg
  Voff <- VR[off_idx, off_idx, drop = FALSE]

  contrast_test <- function(group_a, group_b, label) {
    w <- rep(0, nrow(off))
    w[off$block == group_a] <- 1 / sum(off$block == group_a)
    w[off$block == group_b] <- -1 / sum(off$block == group_b)
    estimate <- sum(w * y)
    se <- sqrt(drop(t(w) %*% Voff %*% w))
    data.table(
      contrast = label,
      estimate = estimate,
      se = se,
      z = estimate / se,
      p = 2 * pnorm(abs(estimate / se), lower.tail = FALSE)
    )
  }

  block_summary <- off[, .(
    n_pairs = .N,
    mean_rg = mean(rg),
    median_rg = median(rg),
    min_rg = min(rg),
    max_rg = max(rg)
  ), by = block]

  contrasts <- rbindlist(list(
    contrast_test(
      "within_aortic", "between_axes",
      "within_aortic_minus_between"
    ),
    contrast_test(
      "within_atherosclerosis_related", "between_axes",
      "within_atherosclerosis_related_minus_between"
    )
  ))

  write_tsv(block_summary, file.path(out, "axis_block_rg_summary.tsv"))
  write_tsv(contrasts, file.path(out, "axis_block_rg_contrasts.tsv"))

  R_smooth <- as.matrix(nearPD(R, corr = TRUE)$mat)
  kmo <- psych::KMO(R_smooth)
  efa1 <- factanal(covmat = R_smooth, factors = 1, rotation = "none")
  efa2 <- factanal(covmat = R_smooth, factors = 2, rotation = "promax")
  write_tsv(
    matrix_to_table(unclass(efa1$loadings), row_name = "trait"),
    file.path(out, "EFA_one_factor_loadings.tsv")
  )
  write_tsv(
    matrix_to_table(unclass(efa2$loadings), row_name = "trait"),
    file.path(out, "EFA_two_factor_promax_loadings.tsv")
  )
  write_tsv(
    data.table(overall_KMO = kmo$MSA),
    file.path(out, "KMO_overall.tsv")
  )
}

# ==============================================================================
# 4. Full candidate-model CFA comparison
# ==============================================================================

run_full_cfa_comparison <- function() {
  out <- file.path(OUT_ROOT, "02_full_CFA_comparison")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  models <- list(
    single_common_factor = list(
      syntax = make_single_factor_model(traits),
      std_lv = TRUE
    ),
    correlated_two_factor = list(
      syntax = make_correlated_model(traits),
      std_lv = TRUE
    ),
    orthogonal_bifactor = list(
      syntax = make_bifactor_model(traits),
      std_lv = TRUE
    ),
    equal_loading_higher_order_negative_control = list(
      syntax = equal_loading_higher_order_negative_control_model,
      std_lv = FALSE
    )
  )

  write_tsv(
    data.table(
      model = c(
        "single_common_factor",
        "correlated_two_factor",
        "orthogonal_bifactor",
        "equal_loading_higher_order_negative_control"
      ),
      role = c(
        "homogeneous one-factor negative control",
        "primary genetic-structure model",
        "pan-vascular general-factor decomposition and GWAS basis",
        "identification-constrained higher-order negative control"
      ),
      eligible_for_primary_GWAS = c(FALSE, FALSE, TRUE, FALSE)
    ),
    file.path(out, "MODEL_ROLES.tsv")
  )

  fit_rows <- list()
  parameter_rows <- list()
  diagnostic_rows <- list()

  for (label in names(models)) {
    model_dir <- file.path(out, label)
    result <- run_cfa(
      LDSCoutput,
      models[[label]]$syntax,
      std_lv = models[[label]]$std_lv
    )
    diagnostic <- save_cfa_result(
      result,
      models[[label]]$syntax,
      model_dir,
      label,
      traits
    )
    if (nrow(result$fit) > 0L) {
      fit_rows[[label]] <- copy(result$fit)[, model := label]
    }
    if (nrow(result$parameters) > 0L) {
      parameter_rows[[label]] <- copy(result$parameters)[, model := label]
    }
    diagnostic_rows[[label]] <- diagnostic
  }

  all_fit <- rbindlist(fit_rows, fill = TRUE)
  all_parameters <- rbindlist(parameter_rows, fill = TRUE)
  all_diagnostics <- rbindlist(diagnostic_rows, fill = TRUE)
  write_tsv(all_fit, file.path(out, "candidate_model_fit.tsv"))
  write_tsv(
    all_parameters,
    file.path(out, "candidate_model_all_parameters.tsv")
  )
  write_tsv(
    all_diagnostics,
    file.path(out, "candidate_model_diagnostics.tsv")
  )

  bif <- all_parameters[model == "orthogonal_bifactor" & op == "=~"]
  loading_col <- pick_col(
    bif,
    c("STD_Genotype", "STD_All", "Unstand_Est")
  )
  bif[, loading := as_numeric_safe(get(loading_col))]
  ss_general <- bif[lhs == "G", sum(loading^2, na.rm = TRUE)]
  ss_specific <- bif[lhs != "G", sum(loading^2, na.rm = TRUE)]
  ecv <- ss_general / (ss_general + ss_specific)
  write_tsv(
    data.table(
      sumsq_general = ss_general,
      sumsq_specific = ss_specific,
      ECV_general = ecv
    ),
    file.path(out, "bifactor_ECV.tsv")
  )
}

# ==============================================================================
# 5. All leave-one-trait CFA analyses
# ==============================================================================

run_all_leave_one_cfa <- function() {
  out <- file.path(OUT_ROOT, "03_all_leave_one_CFA")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  omissions <- c("NONE", traits)
  fit_rows <- list()
  parameter_rows <- list()
  diagnostic_rows <- list()
  axis_rg_rows <- list()

  for (omit in omissions) {
    keep <- if (omit == "NONE") traits else setdiff(traits, omit)
    label <- if (omit == "NONE") "full" else paste0("leave_one_", omit)
    cov_sub <- subset_ldsc(LDSCoutput, keep)

    for (model_type in c("correlated_two_factor", "orthogonal_bifactor")) {
      syntax <- if (model_type == "correlated_two_factor") {
        make_correlated_model(keep)
      } else {
        make_bifactor_model(keep)
      }
      analysis_label <- paste(label, model_type, sep = "__")
      model_dir <- file.path(out, label, model_type)
      result <- run_cfa(cov_sub, syntax, std_lv = TRUE)
      diagnostic <- save_cfa_result(
        result, syntax, model_dir, analysis_label, keep
      )
      diagnostic[, `:=`(
        omitted_trait = ifelse(omit == "NONE", NA_character_, omit),
        model_type = model_type,
        n_aortic_indicators = sum(keep %in% AORTIC_TRAITS),
        n_athero_indicators = sum(keep %in% ATHERO_TRAITS)
      )]
      diagnostic_rows[[analysis_label]] <- diagnostic

      if (nrow(result$fit) > 0L) {
        fit_rows[[analysis_label]] <- copy(result$fit)[, `:=`(
          analysis = label,
          omitted_trait = ifelse(omit == "NONE", NA_character_, omit),
          model_type = model_type
        )]
      }
      if (nrow(result$parameters) > 0L) {
        pars <- copy(result$parameters)
        pars[, `:=`(
          analysis = label,
          omitted_trait = ifelse(omit == "NONE", NA_character_, omit),
          model_type = model_type
        )]
        parameter_rows[[analysis_label]] <- pars

        if (model_type == "correlated_two_factor") {
          axis_rg_rows[[analysis_label]] <- pars[
            lhs == "AORTIC" & op == "~~" & rhs == "ATHERO"
          ]
        }
      }
    }
  }

  all_fit <- rbindlist(fit_rows, fill = TRUE)
  all_parameters <- rbindlist(parameter_rows, fill = TRUE)
  all_diagnostics <- rbindlist(diagnostic_rows, fill = TRUE)
  all_axis_rg <- rbindlist(axis_rg_rows, fill = TRUE)

  write_tsv(all_fit, file.path(out, "all_leave_one_model_fit.tsv"))
  write_tsv(
    all_parameters,
    file.path(out, "all_leave_one_parameters.tsv")
  )
  write_tsv(
    all_diagnostics,
    file.path(out, "all_leave_one_diagnostics.tsv")
  )
  write_tsv(
    all_axis_rg,
    file.path(out, "all_leave_one_axis_correlations.tsv")
  )
}

# ==============================================================================
# 6. Odd/even chromosome cross-validation
# ==============================================================================

run_odd_even_cross_validation <- function() {
  out <- file.path(OUT_ROOT, "04_odd_even_cross_validation")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  gwas_files <- file.path(ROOT, paste0(traits, ".sumstats.gz"))
  if (!all(file.exists(gwas_files))) {
    stop(
      "Missing munged sumstats: ",
      paste(gwas_files[!file.exists(gwas_files)], collapse = ", ")
    )
  }
  ld_dir <- file.path(ROOT, "eur_w_ld_chr")
  if (!dir.exists(ld_dir)) stop("Missing LD directory: ", ld_dir)

  run_half <- function(select_chr, label) {
    cache <- file.path(out, paste0("LDSC_", label, ".rds"))
    if (file.exists(cache)) return(readRDS(cache))

    x <- GenomicSEM::ldsc(
      traits = gwas_files,
      sample.prev = sample_prev,
      population.prev = population_prev,
      ld = ld_dir,
      wld = ld_dir,
      trait.names = traits,
      stand = TRUE,
      select = select_chr,
      ldsc.log = file.path(out, paste0("LDSC_", label))
    )
    saveRDS(x, cache)
    x
  }

  status <- tryCatch({
    odd <- run_half("ODD", "odd_chromosomes")
    even <- run_half("EVEN", "even_chromosomes")

    efa_rows <- list()
    cfa_fit_rows <- list()
    cfa_parameter_rows <- list()

    for (label in c("odd", "even")) {
      x <- if (label == "odd") odd else even
      R_smooth <- as.matrix(nearPD(x$S_Stand, corr = TRUE)$mat)
      efa <- factanal(covmat = R_smooth, factors = 2, rotation = "promax")
      efa_rows[[label]] <- matrix_to_table(
        unclass(efa$loadings),
        row_name = "trait"
      )[, chromosome_half := label]

      for (model_type in c("correlated_two_factor", "orthogonal_bifactor")) {
        syntax <- if (model_type == "correlated_two_factor") {
          make_correlated_model(traits)
        } else {
          make_bifactor_model(traits)
        }
        result <- run_cfa(x, syntax, std_lv = TRUE)
        if (nrow(result$fit) > 0L) {
          cfa_fit_rows[[paste(label, model_type)]] <- copy(result$fit)[, `:=`(
            chromosome_half = label,
            model_type = model_type
          )]
        }
        if (nrow(result$parameters) > 0L) {
          cfa_parameter_rows[[paste(label, model_type)]] <-
            copy(result$parameters)[, `:=`(
              chromosome_half = label,
              model_type = model_type
            )]
        }
      }
    }

    write_tsv(
      rbindlist(efa_rows, fill = TRUE),
      file.path(out, "odd_even_EFA_loadings.tsv")
    )
    write_tsv(
      rbindlist(cfa_fit_rows, fill = TRUE),
      file.path(out, "odd_even_CFA_fit.tsv")
    )
    write_tsv(
      rbindlist(cfa_parameter_rows, fill = TRUE),
      file.path(out, "odd_even_CFA_parameters.tsv")
    )
    data.table(status = "SUCCESS", error = NA_character_)
  }, error = function(e) {
    data.table(status = "ERROR", error = conditionMessage(e))
  })

  write_tsv(status, file.path(out, "odd_even_status.tsv"))
}

# ==============================================================================
# 7. Bifactor G GWAS
# ==============================================================================

run_usergwas_call <- function(
    covstruc,
    snps,
    model,
    label,
    print_warnings = TRUE,
    parallel_run = TRUE,
    cores_run = N_CORES) {
  warnings <- character()
  result <- withCallingHandlers(
    userGWAS(
      covstruc = covstruc,
      SNPs = snps,
      estimation = "DWLS",
      model = model,
      sub = "G~SNP",
      cores = cores_run,
      toler = 1e-20,
      SNPSE = FALSE,
      parallel = parallel_run,
      GC = "standard",
      MPI = FALSE,
      smooth_check = TRUE,
      # Unit-loading identification is required when N_hat is calculated
      # from the resulting factor-GWAS standard errors.
      std.lv = FALSE,
      fix_measurement = TRUE,
      Q_SNP = TRUE,
      printwarn = print_warnings
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  attr(result, "captured_warnings") <- unique(warnings)
  result
}

run_usergwas_chunked <- function(
    covstruc,
    snps,
    model,
    label,
    out_dir) {

  chunk_dir <- file.path(out_dir, "userGWAS_chunks")
  dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

  n_snps <- nrow(snps)
  if (n_snps < 1L) stop("No SNPs supplied to chunked userGWAS.")
  chunk_id <- ceiling(seq_len(n_snps) / USERGWAS_CHUNK_SIZE)
  chunk_indices <- split(seq_len(n_snps), chunk_id)

  model_file <- file.path(out_dir, "GWAS_model.txt")
  if (!file.exists(model_file)) writeLines(model, model_file)
  model_md5 <- unname(tools::md5sum(model_file))

  manifest_path <- file.path(chunk_dir, "chunk_manifest.tsv")
  manifest <- data.table(
    chunk = seq_along(chunk_indices),
    first_row = vapply(chunk_indices, min, integer(1)),
    last_row = vapply(chunk_indices, max, integer(1)),
    n_input = lengths(chunk_indices),
    first_SNP = vapply(
      chunk_indices,
      function(i) as.character(snps$SNP[i[1]]),
      character(1)
    ),
    last_SNP = vapply(
      chunk_indices,
      function(i) as.character(snps$SNP[i[length(i)]]),
      character(1)
    ),
    model_md5 = model_md5,
    status = "PENDING",
    n_output = NA_integer_,
    attempts = 0L,
    elapsed_seconds = NA_real_,
    result_file = file.path(
      chunk_dir,
      sprintf("chunk_%05d.rds", seq_along(chunk_indices))
    ),
    metadata_file = file.path(
      chunk_dir,
      sprintf("chunk_%05d_metadata.tsv", seq_along(chunk_indices))
    ),
    error = NA_character_
  )

  cache_is_valid <- function(row) {
    result_file <- row$result_file
    metadata_file <- row$metadata_file
    if (!file.exists(result_file) || !file.exists(metadata_file)) {
      return(FALSE)
    }
    meta <- tryCatch(fread(metadata_file), error = function(e) NULL)
    if (is.null(meta) || nrow(meta) != 1L) return(FALSE)
    identical(as.character(meta$model_md5[1]), model_md5) &&
      identical(as.integer(meta$n_input[1]), as.integer(row$n_input)) &&
      identical(as.character(meta$first_SNP[1]), row$first_SNP) &&
      identical(as.character(meta$last_SNP[1]), row$last_SNP) &&
      identical(as.character(meta$status[1]), "SUCCESS")
  }

  all_warnings <- character()

  for (k in seq_along(chunk_indices)) {
    idx <- chunk_indices[[k]]
    row <- manifest[k]

    if (!OVERWRITE_EXISTING_GWAS && cache_is_valid(row)) {
      cached <- tryCatch(readRDS(row$result_file), error = function(e) NULL)
      if (!is.null(cached)) {
        cached_n <- tryCatch(
          nrow(standardize_gwas_output(result_to_table(cached))),
          error = function(e) NA_integer_
        )
        if (is.finite(cached_n) && cached_n > 0L) {
          manifest[k, `:=`(
            status = "CACHED_SUCCESS",
            n_output = as.integer(cached_n)
          )]
          write_tsv(manifest, manifest_path)
          log_message(
            label, " | chunk ", k, "/", length(chunk_indices),
            " reused from validated cache."
          )
          next
        }
      }
    }

    chunk_result <- NULL
    chunk_error <- NA_character_
    chunk_warnings <- character()
    started <- Sys.time()

    for (attempt in seq_len(USERGWAS_CHUNK_MAX_RETRIES + 1L)) {
      manifest[k, attempts := attempt]
      # The installed server version can produce an output-name length
      # mismatch when printwarn=FALSE is combined with Q_SNP/smooth_check.
      # Keep printwarn=TRUE, matching the original compiled wrapper. If a
      # parallel task still fails, retry that chunk serially as a local
      # fallback without slowing successful chunks.
      use_parallel <- attempt == 1L
      candidate <- tryCatch(
        run_usergwas_call(
          covstruc = covstruc,
          snps = snps[idx],
          model = model,
          label = paste0(label, "_chunk_", k),
          print_warnings = TRUE,
          parallel_run = use_parallel,
          cores_run = if (use_parallel) N_CORES else 1L
        ),
        error = function(e) e
      )

      if (inherits(candidate, "error")) {
        chunk_error <- conditionMessage(candidate)
        log_message(
          label, " | chunk ", k, "/", length(chunk_indices),
          " failed on attempt ", attempt, ": ", chunk_error
        )
      } else {
        standardized_n <- tryCatch(
          nrow(standardize_gwas_output(result_to_table(candidate))),
          error = function(e) {
            chunk_error <<- conditionMessage(e)
            NA_integer_
          }
        )
        if (is.finite(standardized_n) &&
            standardized_n >= 0.80 * length(idx)) {
          chunk_result <- candidate
          chunk_warnings <- attr(candidate, "captured_warnings")
          break
        }
        if (is.na(chunk_error)) {
          chunk_error <- paste0(
            "Only ", standardized_n, " of ", length(idx),
            " input SNPs produced usable rows."
          )
        }
      }
    }

    elapsed <- as.numeric(
      difftime(Sys.time(), started, units = "secs")
    )

    if (is.null(chunk_result)) {
      manifest[k, `:=`(
        status = "FAILED",
        elapsed_seconds = elapsed,
        error = chunk_error
      )]
      write_tsv(manifest, manifest_path)
      next
    }

    saveRDS(chunk_result, row$result_file, compress = TRUE)
    standardized_n <- nrow(
      standardize_gwas_output(result_to_table(chunk_result))
    )
    meta <- data.table(
      chunk = k,
      model_md5 = model_md5,
      n_input = length(idx),
      n_output = standardized_n,
      first_SNP = row$first_SNP,
      last_SNP = row$last_SNP,
      status = "SUCCESS",
      attempts = manifest$attempts[k],
      execution_mode = if (manifest$attempts[k] == 1L) {
        "parallel"
      } else {
        "serial_retry"
      },
      elapsed_seconds = elapsed,
      R_version = R.version.string,
      GenomicSEM_version = as.character(
        utils::packageVersion("GenomicSEM")
      )
    )
    write_tsv(meta, row$metadata_file)
    manifest[k, `:=`(
      status = "SUCCESS",
      n_output = standardized_n,
      elapsed_seconds = elapsed,
      error = NA_character_
    )]
    write_tsv(manifest, manifest_path)
    all_warnings <- unique(c(all_warnings, chunk_warnings))
    rm(chunk_result)
    gc()
    log_message(
      label, " | chunk ", k, "/", length(chunk_indices),
      " completed in ", round(elapsed, 1), " s."
    )
  }

  failed <- manifest[status == "FAILED"]
  if (nrow(failed) > 0L) {
    stop(
      nrow(failed), " userGWAS chunks failed. ",
      "No combined GWAS or FUMA file was produced. Review ",
      manifest_path
    )
  }

  combined <- vector("list", nrow(manifest))
  for (k in seq_len(nrow(manifest))) {
    chunk_result <- readRDS(manifest$result_file[k])
    combined[[k]] <- standardize_gwas_output(
      result_to_table(chunk_result)
    )
    rm(chunk_result)
    if (k %% 25L == 0L) gc()
  }
  combined <- rbindlist(combined, use.names = TRUE, fill = TRUE)
  if (nrow(combined) < 0.80 * n_snps) {
    stop(
      "Combined chunked GWAS contains only ", nrow(combined),
      " rows for ", n_snps, " input SNPs."
    )
  }

  list(
    gwas = combined,
    warnings = all_warnings,
    manifest = manifest
  )
}

run_one_bifactor_gwas <- function(
    LDSC_full,
    SUMSTATS_full,
    omitted_trait = "NONE") {

  keep <- if (omitted_trait == "NONE") {
    traits
  } else {
    setdiff(traits, omitted_trait)
  }
  label <- if (omitted_trait == "NONE") {
    "full_bifactor_G"
  } else {
    paste0("leave_one_", omitted_trait, "_bifactor_G")
  }
  out <- file.path(OUT_ROOT, "05_bifactor_GWAS", label)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  compact_path <- file.path(out, paste0(label, "_GWAS_compact.tsv.gz"))
  if (file.exists(compact_path) && !OVERWRITE_EXISTING_GWAS) {
    log_message("Skipping existing GWAS: ", label)
    return(data.table(
      analysis = label,
      omitted_trait = omitted_trait,
      status = "SKIPPED_EXISTING",
      output = compact_path
    ))
  }

  cov_sub <- subset_ldsc(LDSC_full, keep)
  snps_sub <- drop_trait_sumstats(SUMSTATS_full, omitted_trait)
  measurement_spec <- make_bifactor_gwas_spec(
    keep,
    include_snp = FALSE
  )
  gwas_spec <- make_bifactor_gwas_spec(
    keep,
    include_snp = TRUE
  )
  measurement_model <- measurement_spec$syntax
  gwas_model <- gwas_spec$syntax
  writeLines(measurement_model, file.path(out, "measurement_model.txt"))
  writeLines(gwas_model, file.path(out, "GWAS_model.txt"))
  write_tsv(
    data.table(
      analysis = label,
      omitted_trait = omitted_trait,
      G_anchor = gwas_spec$G_anchor,
      AORTIC_anchor = gwas_spec$AORTIC_anchor,
      ATHERO_anchor = gwas_spec$ATHERO_anchor,
      std_lv = FALSE
    ),
    file.path(out, "GWAS_factor_scaling.tsv")
  )

  baseline <- run_cfa(cov_sub, measurement_model, std_lv = FALSE)
  baseline_diagnostic <- save_cfa_result(
    baseline,
    measurement_model,
    file.path(out, "baseline_CFA"),
    label,
    keep
  )
  baseline_diagnostic[, `:=`(
    omitted_trait = omitted_trait,
    n_aortic_indicators = sum(keep %in% AORTIC_TRAITS),
    n_athero_indicators = sum(keep %in% ATHERO_TRAITS)
  )]
  write_tsv(
    baseline_diagnostic,
    file.path(out, "baseline_GWAS_gate.tsv")
  )

  if (!isTRUE(baseline_diagnostic$accepted_for_gwas[1])) {
    log_message("Skipping ", label, ": baseline bifactor failed QC.")
    return(data.table(
      analysis = label,
      omitted_trait = omitted_trait,
      status = "SKIPPED_BASELINE_CFA_QC",
      output = NA_character_
    ))
  }

  if (RUN_PILOT_FIRST) {
    set.seed(PILOT_SEED + match(
      omitted_trait,
      c("NONE", traits)
    ))
    n_pilot <- min(PILOT_N, nrow(snps_sub))
    pilot_index <- sample.int(nrow(snps_sub), n_pilot)
    pilot_result <- tryCatch(
      run_usergwas_call(
        cov_sub,
        snps_sub[pilot_index],
        gwas_model,
        paste0(label, "_pilot"),
        print_warnings = TRUE
      ),
      error = function(e) e
    )

    if (inherits(pilot_result, "error")) {
      pilot_qc <- data.table(
        analysis = label,
        status = "PILOT_ERROR",
        error = conditionMessage(pilot_result)
      )
      write_tsv(pilot_qc, file.path(out, "pilot_QC.tsv"))
      return(data.table(
        analysis = label,
        omitted_trait = omitted_trait,
        status = "PILOT_ERROR",
        output = NA_character_
      ))
    }

    pilot_warnings <- attr(pilot_result, "captured_warnings")
    pilot_tab <- standardize_gwas_output(result_to_table(pilot_result))
    pilot_finite <- is.finite(pilot_tab$P) &
      is.finite(pilot_tab$BETA) &
      is.finite(pilot_tab$SE) &
      pilot_tab$SE > 0
    pilot_qc <- data.table(
      analysis = label,
      status = "PILOT_SUCCESS",
      n_requested = n_pilot,
      n_returned = nrow(pilot_tab),
      finite_fraction = mean(pilot_finite),
      n_warnings = length(pilot_warnings),
      warnings = paste(pilot_warnings, collapse = " | ")
    )
    write_tsv(pilot_qc, file.path(out, "pilot_QC.tsv"))
    write_tsv(
      pilot_tab,
      file.path(out, paste0(label, "_pilot_results.tsv"))
    )
    rm(pilot_result, pilot_tab)
    gc()

    if (
      pilot_qc$n_returned[1] < 0.80 * n_pilot ||
      !is.finite(pilot_qc$finite_fraction[1]) ||
      pilot_qc$finite_fraction[1] < 0.90
    ) {
      return(data.table(
        analysis = label,
        omitted_trait = omitted_trait,
        status = "PILOT_FAILED_QC",
        output = NA_character_
      ))
    }
  }

  if (!RUN_FULL_AFTER_PILOT) {
    return(data.table(
      analysis = label,
      omitted_trait = omitted_trait,
      status = "PILOT_ONLY",
      output = NA_character_
    ))
  }

  log_message("Starting full userGWAS: ", label)

  full_run <- tryCatch({
    if (identical(USERGWAS_ENGINE, "chunked")) {
      run_usergwas_chunked(
        covstruc = cov_sub,
        snps = snps_sub,
        model = gwas_model,
        label = label,
        out_dir = out
      )
    } else if (identical(USERGWAS_ENGINE, "direct")) {
      direct_result <- run_usergwas_call(
        cov_sub,
        snps_sub,
        gwas_model,
        label,
        print_warnings = TRUE,
        parallel_run = TRUE,
        cores_run = N_CORES
      )
      list(
        direct_result = direct_result,
        warnings = attr(direct_result, "captured_warnings")
      )
    } else {
      stop("USERGWAS_ENGINE must be 'chunked' or 'direct'.")
    }
  }, error = function(e) e)

  if (inherits(full_run, "error")) {
    error_tab <- data.table(
      analysis = label,
      status = "FULL_GWAS_ERROR",
      error = conditionMessage(full_run)
    )
    write_tsv(error_tab, file.path(out, "full_GWAS_error.tsv"))
    return(data.table(
      analysis = label,
      omitted_trait = omitted_trait,
      status = "FULL_GWAS_ERROR",
      output = NA_character_
    ))
  }

  captured_warnings <- full_run$warnings
  writeLines(
    captured_warnings,
    file.path(out, "full_GWAS_warnings.txt")
  )
  if (identical(USERGWAS_ENGINE, "chunked")) {
    gwas <- full_run$gwas
  } else {
    if (SAVE_RAW_USERGWAS_RDS) {
      saveRDS(
        full_run$direct_result,
        file.path(out, paste0(label, "_raw_userGWAS.rds"))
      )
    }
    gwas <- standardize_gwas_output(
      result_to_table(full_run$direct_result)
    )
  }
  rm(full_run)
  gc()

  write_tsv_gz(gwas, compact_path)
  fuma_qc <- write_fuma_outputs(gwas, out, label)
  fuma_qc[, `:=`(
    omitted_trait = omitted_trait,
    status = "SUCCESS",
    n_captured_warnings = length(captured_warnings)
  )]
  write_tsv(fuma_qc, file.path(out, "full_GWAS_QC.tsv"))

  rm(gwas, snps_sub)
  gc()

  data.table(
    analysis = label,
    omitted_trait = omitted_trait,
    status = "SUCCESS",
    output = compact_path
  )
}

compare_leave_one_with_full <- function(full_path, leave_path, label, out_dir) {
  if (!file.exists(full_path) || !file.exists(leave_path)) return(NULL)

  full <- fread(full_path, select = c("SNP", "BETA", "SE", "P"))
  leave <- fread(leave_path, select = c("SNP", "BETA", "SE", "P"))
  setnames(full, c("BETA", "SE", "P"), c(
    "BETA_full", "SE_full", "P_full"
  ))
  setnames(leave, c("BETA", "SE", "P"), c(
    "BETA_leave", "SE_leave", "P_leave"
  ))

  set.seed(PILOT_SEED)
  n_sample <- min(GWAS_COMPARE_N, nrow(full))
  sampled_snps <- full[sample.int(nrow(full), n_sample), SNP]
  cmp <- merge(
    full[SNP %in% sampled_snps],
    leave[SNP %in% sampled_snps],
    by = "SNP"
  )
  cmp[, `:=`(
    Z_full = BETA_full / SE_full,
    Z_leave = BETA_leave / SE_leave
  )]

  full_gws <- full[P_full < GWS_P]
  leave_gws <- leave[P_leave < GWS_P]
  gws_union <- merge(
    full[SNP %in% union(full_gws$SNP, leave_gws$SNP)],
    leave[SNP %in% union(full_gws$SNP, leave_gws$SNP)],
    by = "SNP",
    all = TRUE
  )
  write_tsv(
    gws_union,
    file.path(out_dir, paste0(label, "_GWS_union_comparison.tsv"))
  )

  summary <- data.table(
    analysis = label,
    n_common_sampled_SNPs = nrow(cmp),
    beta_Pearson_r = cor(
      cmp$BETA_full, cmp$BETA_leave, use = "complete.obs"
    ),
    z_Pearson_r = cor(
      cmp$Z_full, cmp$Z_leave, use = "complete.obs"
    ),
    z_Spearman_r = cor(
      cmp$Z_full, cmp$Z_leave,
      use = "complete.obs", method = "spearman"
    ),
    sign_concordance = mean(
      sign(cmp$BETA_full) == sign(cmp$BETA_leave),
      na.rm = TRUE
    ),
    n_GWS_full = nrow(full_gws),
    n_GWS_leave = nrow(leave_gws),
    n_GWS_shared = length(intersect(full_gws$SNP, leave_gws$SNP))
  )
  write_tsv(
    summary,
    file.path(out_dir, paste0(label, "_comparison_summary.tsv"))
  )
  rm(full, leave, cmp, full_gws, leave_gws, gws_union)
  gc()
  summary
}

run_requested_bifactor_gwas <- function() {
  out <- file.path(OUT_ROOT, "05_bifactor_GWAS")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  PSYCH_sumstats <- as.data.table(
    load_rdata_object(SUMSTATS_PATH, "PSYCH_sumstats")
  )
  expected_sumstat_cols <- unlist(lapply(traits, function(x) {
    c(paste0("beta.", x), paste0("se.", x))
  }))
  missing_sumstat_cols <- setdiff(expected_sumstat_cols, names(PSYCH_sumstats))
  if (length(missing_sumstat_cols) > 0L) {
    stop(
      "PSYCH_sumstats is missing: ",
      paste(missing_sumstat_cols, collapse = ", ")
    )
  }

  requested_omissions <- character()
  if (RUN_FULL_BIFACTOR_GWAS) requested_omissions <- c(
    requested_omissions, "NONE"
  )
  if (RUN_LEAVE_ONE_CAD_GWAS) requested_omissions <- c(
    requested_omissions, "CAD"
  )
  if (RUN_ALL_LEAVE_ONE_GWAS) requested_omissions <- c(
    requested_omissions, traits
  )
  requested_omissions <- unique(requested_omissions)

  status_rows <- list()
  for (omit in requested_omissions) {
    status_rows[[omit]] <- run_one_bifactor_gwas(
      LDSCoutput,
      PSYCH_sumstats,
      omitted_trait = omit
    )
    gc()
  }
  statuses <- rbindlist(status_rows, fill = TRUE)
  write_tsv(statuses, file.path(out, "all_requested_GWAS_status.tsv"))

  full_path <- file.path(
    out,
    "full_bifactor_G",
    "full_bifactor_G_GWAS_compact.tsv.gz"
  )
  comparison_rows <- list()
  for (omit in setdiff(requested_omissions, "NONE")) {
    label <- paste0("leave_one_", omit, "_bifactor_G")
    leave_path <- file.path(
      out,
      label,
      paste0(label, "_GWAS_compact.tsv.gz")
    )
    comparison_rows[[omit]] <- compare_leave_one_with_full(
      full_path,
      leave_path,
      label,
      file.path(out, label)
    )
  }
  comparison_rows <- Filter(Negate(is.null), comparison_rows)
  if (length(comparison_rows) > 0L) {
    write_tsv(
      rbindlist(comparison_rows, fill = TRUE),
      file.path(out, "all_leave_one_vs_full_comparison_summary.tsv")
    )
  }

  fuma_files <- list.files(
    out,
    pattern = "_FUMA_.*\\.tsv\\.gz$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(fuma_files) > 0L) {
    manifest <- data.table(
      file = fuma_files,
      analysis = basename(dirname(fuma_files)),
      recommended_use = fifelse(
        grepl("excluding_QSNP", fuma_files),
        "Q_SNP sensitivity analysis",
        fifelse(
          grepl("full_bifactor_G", fuma_files),
          "primary FUMA analysis",
          "leave-one-trait sensitivity analysis"
        )
      )
    )
    write_tsv(
      manifest,
      file.path(out, "FUMA_UPLOAD_MANIFEST.tsv")
    )
  }
}

# ==============================================================================
# 8. Pipeline execution
# ==============================================================================

main <- function() {
  sink(LOG_FILE, split = TRUE)
  on.exit({
    writeLines(
      capture.output(sessionInfo()),
      file.path(OUT_ROOT, "R_sessionInfo.txt")
    )
    sink()
  }, add = TRUE)

  log_message("Starting complete Genomic SEM follow-up pipeline.")
  log_message("Output root: ", OUT_ROOT)
  log_message("Cores: ", N_CORES)

  step_status <- list()

  run_step <- function(label, enabled, fun) {
    if (!enabled) {
      log_message("Skipping disabled step: ", label)
      return(data.table(
        step = label, status = "DISABLED", error = NA_character_
      ))
    }
    log_message("Starting step: ", label)
    result <- tryCatch({
      fun()
      data.table(step = label, status = "SUCCESS", error = NA_character_)
    }, error = function(e) {
      log_message("Step failed: ", label, " | ", conditionMessage(e))
      data.table(
        step = label,
        status = "ERROR",
        error = conditionMessage(e)
      )
    })
    result
  }

  step_status[["01_ldsc_two_axis"]] <- run_step(
    "01_ldsc_two_axis",
    RUN_LDSC_AXIS_ANALYSIS,
    run_ldsc_axis_analysis
  )
  step_status[["02_full_CFA_comparison"]] <- run_step(
    "02_full_CFA_comparison",
    RUN_FULL_CFA_COMPARISON,
    run_full_cfa_comparison
  )
  step_status[["03_all_leave_one_CFA"]] <- run_step(
    "03_all_leave_one_CFA",
    RUN_ALL_LEAVE_ONE_CFA,
    run_all_leave_one_cfa
  )
  step_status[["04_odd_even"]] <- run_step(
    "04_odd_even",
    RUN_ODD_EVEN_LDSC,
    run_odd_even_cross_validation
  )
  step_status[["05_bifactor_GWAS"]] <- run_step(
    "05_bifactor_GWAS",
    RUN_FULL_BIFACTOR_GWAS ||
      RUN_LEAVE_ONE_CAD_GWAS ||
      RUN_ALL_LEAVE_ONE_GWAS,
    run_requested_bifactor_gwas
  )

  all_status <- rbindlist(step_status, fill = TRUE)
  write_tsv(all_status, file.path(OUT_ROOT, "PIPELINE_STATUS.tsv"))
  log_message("Pipeline finished. Review PIPELINE_STATUS.tsv first.")
}

tryCatch(
  main(),
  error = function(e) {
    message("Fatal pipeline error: ", conditionMessage(e))
    quit(status = 1L)
  }
)
