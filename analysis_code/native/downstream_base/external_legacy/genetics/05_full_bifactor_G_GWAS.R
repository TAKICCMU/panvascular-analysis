options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(GenomicSEM)
  library(data.table)
  library(parallel)
})

ROOT <- "/path/to/panvas_analysis/Genomic_SEM"
OUT <- "/path/to/panvas_analysis/Genomic_SEM_bifactor_GWAS/01_full"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

AUDITED <- "/path/to/panvas_analysis/Genomic_SEM_validation/00_input_audit/LDSCoutput_audited.RData"
load(if (file.exists(AUDITED)) AUDITED else file.path(ROOT, "LDSCoutput.RData"))
load(file.path(ROOT, "PSYCH_sumstats.RData"))
stopifnot(exists("LDSCoutput"), exists("PSYCH_sumstats"))

# Unit-loading identification is used so that implied N is interpretable.
# Stroke is the general-factor anchor and is retained in leave-one-CAD.
# AAA and PAD anchor the two specific factors.
measurement_model <- "
  G =~ 1*Stroke + AA + AAA + TAAD + CAD + PAD + CA + CAVS
  AORTIC =~ 1*AAA + AA + TAAD
  ATHERO =~ 1*PAD + CAD + Stroke + CA + CAVS
  G ~~ G
  AORTIC ~~ AORTIC
  ATHERO ~~ ATHERO
  G ~~ 0*AORTIC
  G ~~ 0*ATHERO
  AORTIC ~~ 0*ATHERO
  AAA ~~ 1e-6*AAA
"

# Confirm that changing only the factor scale preserves the good bifactor fit.
baseline <- usermodel(
  covstruc = LDSCoutput,
  estimation = "DWLS",
  model = measurement_model,
  CFIcalc = TRUE,
  std.lv = FALSE,
  imp_cov = TRUE,
  fix_resid = FALSE,
  toler = 1e-20
)
fwrite(as.data.table(baseline$modelfit),
       file.path(OUT, "baseline_bifactor_fit.tsv"), sep = "\t")
fwrite(as.data.table(baseline$results),
       file.path(OUT, "baseline_bifactor_parameters.tsv"), sep = "\t")

fit_check <- as.data.table(baseline$modelfit)
if (fit_check$CFI[1] < 0.95 || fit_check$SRMR[1] > 0.10) {
  stop("Unit-loading bifactor baseline fit is not acceptable. Do not run GWAS.")
}

gwas_model <- paste(measurement_model, "G ~ SNP", sep = "\n")

# Use "pilot" first. Change to "full" only after checking pilot diagnostics.
RUN_MODE <- "pilot"
set.seed(20260729)
if (RUN_MODE == "pilot") {
  n_test <- min(20000L, nrow(PSYCH_sumstats))
  SNP_INPUT <- PSYCH_sumstats[sample.int(nrow(PSYCH_sumstats), n_test), ]
} else if (RUN_MODE == "full") {
  SNP_INPUT <- PSYCH_sumstats
} else {
  stop("RUN_MODE must be pilot or full")
}

n_cores <- max(1L, floor(detectCores() / 2))
result <- userGWAS(
  covstruc = LDSCoutput,
  SNPs = SNP_INPUT,
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

saveRDS(result, file.path(OUT, paste0("bifactor_G_", RUN_MODE, "_raw.rds")))
gwas <- if (is.data.frame(result)) {
  as.data.table(result)
} else {
  rbindlist(result, fill = TRUE)
}
fwrite(gwas, file.path(OUT, paste0("bifactor_G_", RUN_MODE, "_all.tsv")),
       sep = "\t")

pick_col <- function(candidates, required = TRUE) {
  hit <- candidates[candidates %in% names(gwas)]
  if (length(hit)) return(hit[1])
  if (required) stop("Missing output column: ", paste(candidates, collapse = "/"))
  NA_character_
}

beta_col <- pick_col(c("est", "Unstand_Est", "BETA", "beta"))
se_col <- pick_col(c("se_c", "SE", "Unstand_SE", "se"))
p_col <- pick_col(c("Pval_Estimate", "p_value", "P", "p"))
q_p_col <- pick_col(c("Q_SNP_pval", "Q_SNP_P", "Q_p"), required = FALSE)

maf_col <- pick_col(c("MAF"))
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
  factor = "G", anchor_trait = "Stroke", run_mode = RUN_MODE,
  N_hat = N_hat, n_snps_for_N_hat = sum(valid_for_n)
), file.path(OUT, paste0("bifactor_G_", RUN_MODE, "_N_hat.tsv")),
sep = "\t")

if (RUN_MODE == "full") {
  fuma <- data.table(
    SNP = gwas[[pick_col(c("SNP"))]],
    CHR = gwas[[pick_col(c("CHR"))]],
    BP = gwas[[pick_col(c("BP"))]],
    A1 = gwas[[pick_col(c("A1"))]],
    A2 = gwas[[pick_col(c("A2"))]],
    MAF = gwas[[maf_col]],
    BETA = gwas[[beta_col]],
    SE = gwas[[se_col]],
    P = gwas[[p_col]],
    N = N_hat
  )
  fwrite(fuma, file.path(OUT, "bifactor_G_FUMA_all.tsv"), sep = "\t")

  if (!is.na(q_p_col)) {
    q_sig <- is.finite(gwas[[q_p_col]]) & gwas[[q_p_col]] < 5e-8
    fwrite(gwas[q_sig],
           file.path(OUT, "bifactor_G_QSNP_genomewide_significant.tsv"),
           sep = "\t")
    fwrite(fuma[!q_sig],
           file.path(OUT, "bifactor_G_FUMA_excluding_direct_QSNP.tsv"),
           sep = "\t")
  }
}

writeLines(capture.output(sessionInfo()), file.path(OUT, "R_sessionInfo.txt"))
