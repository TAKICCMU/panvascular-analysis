options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(GenomicSEM)
  library(data.table)
  library(parallel)
})

ROOT <- "/path/to/panvas_analysis/Genomic_SEM"
FULL_OUT <- "/path/to/panvas_analysis/Genomic_SEM_bifactor_GWAS/01_full"
OUT <- "/path/to/panvas_analysis/Genomic_SEM_bifactor_GWAS/02_leave_one_CAD"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

AUDITED <- "/path/to/panvas_analysis/Genomic_SEM_validation/00_input_audit/LDSCoutput_audited.RData"
load(if (file.exists(AUDITED)) AUDITED else file.path(ROOT, "LDSCoutput.RData"))
load(file.path(ROOT, "PSYCH_sumstats.RData"))
stopifnot(exists("LDSCoutput"), exists("PSYCH_sumstats"))

pair_index_table <- function(p) {
  rbindlist(lapply(seq_len(p), function(i) {
    data.table(i = i, j = i:p)
  }))
}

subset_ldsc <- function(x, keep_traits) {
  old_traits <- colnames(x$S)
  keep <- match(keep_traits, old_traits)
  stopifnot(!anyNA(keep))
  pairs <- pair_index_table(length(old_traits))
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

all_traits <- colnames(LDSCoutput$S)
keep_traits <- setdiff(all_traits, "CAD")
LDSC_no_CAD <- subset_ldsc(LDSCoutput, keep_traits)

# Keep the exact SNP set used in the full analysis for a direct comparison.
drop_cols <- intersect(c("beta.CAD", "se.CAD"), names(PSYCH_sumstats))
if (length(drop_cols) != 2) {
  stop("Cannot find beta.CAD and se.CAD in PSYCH_sumstats")
}
SUMSTATS_no_CAD <- PSYCH_sumstats[, setdiff(names(PSYCH_sumstats), drop_cols)]

measurement_model <- "
  G =~ 1*Stroke + AA + AAA + TAAD + PAD + CA + CAVS
  AORTIC =~ 1*AAA + AA + TAAD
  ATHERO =~ 1*PAD + Stroke + CA + CAVS
  G ~~ G
  AORTIC ~~ AORTIC
  ATHERO ~~ ATHERO
  G ~~ 0*AORTIC
  G ~~ 0*ATHERO
  AORTIC ~~ 0*ATHERO
  AAA ~~ 1e-6*AAA
"

baseline <- usermodel(
  covstruc = LDSC_no_CAD,
  estimation = "DWLS",
  model = measurement_model,
  CFIcalc = TRUE,
  std.lv = FALSE,
  imp_cov = TRUE,
  fix_resid = FALSE,
  toler = 1e-20
)
fwrite(as.data.table(baseline$modelfit),
       file.path(OUT, "leave_CAD_baseline_fit.tsv"), sep = "\t")
fwrite(as.data.table(baseline$results),
       file.path(OUT, "leave_CAD_baseline_parameters.tsv"), sep = "\t")

fit_check <- as.data.table(baseline$modelfit)
if (fit_check$CFI[1] < 0.95 || fit_check$SRMR[1] > 0.10) {
  stop("Leave-one-CAD bifactor fit is not acceptable. Do not run GWAS.")
}

gwas_model <- paste(measurement_model, "G ~ SNP", sep = "\n")
n_cores <- max(1L, floor(detectCores() / 2))

result <- userGWAS(
  covstruc = LDSC_no_CAD,
  SNPs = SUMSTATS_no_CAD,
  estimation = "DWLS",
  model = gwas_model,
  sub = "G~SNP",
  cores = n_cores,
  toler = 1e-20,
  SNPSE = FALSE,
  parallel = TRUE,
  GC = "standard",
  MPI = FALSE,
  smooth_check = TRUE,
  std.lv = FALSE,
  fix_measurement = TRUE,
  Q_SNP = TRUE,
  printwarn = TRUE
)

saveRDS(result, file.path(OUT, "leave_CAD_bifactor_G_raw.rds"))
gwas <- if (is.data.frame(result)) {
  as.data.table(result)
} else {
  rbindlist(result, fill = TRUE)
}
fwrite(gwas, file.path(OUT, "leave_CAD_bifactor_G_all.tsv"), sep = "\t")

pick_col <- function(tab, candidates, required = TRUE) {
  hit <- candidates[candidates %in% names(tab)]
  if (length(hit)) return(hit[1])
  if (required) stop("Missing output column: ", paste(candidates, collapse = "/"))
  NA_character_
}

beta_col <- pick_col(gwas, c("est", "Unstand_Est", "BETA", "beta"))
se_col <- pick_col(gwas, c("se_c", "SE", "Unstand_SE", "se"))
p_col <- pick_col(gwas, c("Pval_Estimate", "p_value", "P", "p"))
maf_col <- pick_col(gwas, "MAF")

valid_for_n <- is.finite(gwas[[se_col]]) & gwas[[se_col]] > 0 &
  is.finite(gwas[[maf_col]]) & gwas[[maf_col]] >= 0.10 &
  gwas[[maf_col]] < 0.5
N_hat <- mean(
  1 / (2 * gwas[[maf_col]][valid_for_n] *
         (1 - gwas[[maf_col]][valid_for_n]) *
         gwas[[se_col]][valid_for_n]^2),
  na.rm = TRUE
)
fwrite(data.table(
  factor = "G", anchor_trait = "Stroke", omitted_trait = "CAD",
  N_hat = N_hat, n_snps_for_N_hat = sum(valid_for_n)
), file.path(OUT, "leave_CAD_bifactor_G_N_hat.tsv"), sep = "\t")

fuma <- data.table(
  SNP = gwas[[pick_col(gwas, "SNP")]],
  CHR = gwas[[pick_col(gwas, "CHR")]],
  BP = gwas[[pick_col(gwas, "BP")]],
  A1 = gwas[[pick_col(gwas, "A1")]],
  A2 = gwas[[pick_col(gwas, "A2")]],
  MAF = gwas[[maf_col]],
  BETA = gwas[[beta_col]],
  SE = gwas[[se_col]],
  P = gwas[[p_col]],
  N = N_hat
)
fwrite(fuma, file.path(OUT, "leave_CAD_bifactor_G_FUMA.tsv"), sep = "\t")

# Compare with the full bifactor GWAS after that run has been completed.
full_path <- file.path(FULL_OUT, "bifactor_G_full_all.tsv")
if (file.exists(full_path)) {
  full <- fread(full_path)
  snp_full <- pick_col(full, "SNP")
  beta_full <- pick_col(full, c("est", "Unstand_Est", "BETA", "beta"))
  se_full <- pick_col(full, c("se_c", "SE", "Unstand_SE", "se"))
  p_full <- pick_col(full, c("Pval_Estimate", "p_value", "P", "p"))

  cmp <- merge(
    full[, .(
      SNP = get(snp_full),
      beta_full = get(beta_full),
      se_full = get(se_full),
      p_full = get(p_full)
    )],
    gwas[, .(
      SNP = get(pick_col(gwas, "SNP")),
      beta_leave_CAD = get(beta_col),
      se_leave_CAD = get(se_col),
      p_leave_CAD = get(p_col)
    )],
    by = "SNP"
  )
  cmp[, `:=`(
    z_full = beta_full / se_full,
    z_leave_CAD = beta_leave_CAD / se_leave_CAD
  )]

  summary <- data.table(
    n_common_SNPs = nrow(cmp),
    beta_Pearson_r = cor(cmp$beta_full, cmp$beta_leave_CAD,
                         use = "complete.obs"),
    z_Pearson_r = cor(cmp$z_full, cmp$z_leave_CAD,
                      use = "complete.obs"),
    sign_concordance = mean(sign(cmp$beta_full) ==
                              sign(cmp$beta_leave_CAD), na.rm = TRUE),
    n_GWS_full = sum(cmp$p_full < 5e-8, na.rm = TRUE),
    n_GWS_leave_CAD = sum(cmp$p_leave_CAD < 5e-8, na.rm = TRUE)
  )
  fwrite(cmp, file.path(OUT, "full_vs_leave_CAD_SNP_comparison.tsv"),
         sep = "\t")
  fwrite(summary, file.path(OUT, "full_vs_leave_CAD_summary.tsv"),
         sep = "\t")
  print(summary)
}

writeLines(capture.output(sessionInfo()), file.path(OUT, "R_sessionInfo.txt"))
