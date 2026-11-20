options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(GenomicSEM)
  library(data.table)
})

ROOT <- "/path/to/panvas_analysis/Genomic_SEM"
OUT <- "/path/to/panvas_analysis/Genomic_SEM_validation/03_correlated_leave_one_CFA"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

AUDITED <- "/path/to/panvas_analysis/Genomic_SEM_validation/00_input_audit/LDSCoutput_audited.RData"
load(if (file.exists(AUDITED)) AUDITED else file.path(ROOT, "LDSCoutput.RData"))
stopifnot(exists("LDSCoutput"))

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

make_correlated_model <- function(keep_traits) {
  aortic <- intersect(c("AA", "AAA", "TAAD"), keep_traits)
  athero <- intersect(c("CAD", "Stroke", "PAD", "CA", "CAVS"), keep_traits)
  stopifnot(length(aortic) >= 2, length(athero) >= 2)

  lines <- c(
    paste("AORTIC =~ NA*", aortic[1], " + ",
          paste(aortic[-1], collapse = " + "), sep = ""),
    paste("ATHERO =~ NA*", athero[1], " + ",
          paste(athero[-1], collapse = " + "), sep = ""),
    "AORTIC ~~ r12*ATHERO",
    "AORTIC ~~ 1*AORTIC",
    "ATHERO ~~ 1*ATHERO"
  )
  if ("AAA" %in% keep_traits) lines <- c(lines, "AAA ~~ 1e-6*AAA")
  paste(lines, collapse = "\n")
}

all_traits <- colnames(LDSCoutput$S)
omissions <- c("NONE", all_traits)
fit_rows <- list()
rg_rows <- list()
parameter_rows <- list()

for (omit in omissions) {
  keep <- if (omit == "NONE") all_traits else setdiff(all_traits, omit)
  analysis <- if (omit == "NONE") "full" else paste0("leave_one_", omit)
  model <- make_correlated_model(keep)
  data_sub <- subset_ldsc(LDSCoutput, keep)

  warnings <- character()
  fit <- withCallingHandlers(
    usermodel(
      covstruc = data_sub,
      estimation = "DWLS",
      model = model,
      CFIcalc = TRUE,
      std.lv = TRUE,
      imp_cov = TRUE,
      fix_resid = FALSE,
      toler = 1e-20
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  fit_tab <- as.data.table(fit$modelfit)
  fit_tab[, `:=`(
    analysis = analysis,
    omitted_trait = ifelse(omit == "NONE", NA_character_, omit),
    warnings = paste(unique(warnings), collapse = " | ")
  )]

  par_tab <- as.data.table(fit$results)
  par_tab[, `:=`(
    analysis = analysis,
    omitted_trait = ifelse(omit == "NONE", NA_character_, omit)
  )]
  rg_tab <- par_tab[lhs == "AORTIC" & op == "~~" & rhs == "ATHERO"]

  fit_rows[[analysis]] <- fit_tab
  parameter_rows[[analysis]] <- par_tab
  rg_rows[[analysis]] <- rg_tab
}

fit_all <- rbindlist(fit_rows, fill = TRUE)
par_all <- rbindlist(parameter_rows, fill = TRUE)
rg_all <- rbindlist(rg_rows, fill = TRUE)

est_col <- if ("STD_Genotype" %in% names(rg_all)) "STD_Genotype" else "STD_All"
se_col <- if ("STD_Genotype_SE" %in% names(rg_all)) {
  "STD_Genotype_SE"
} else {
  "Unstand_SE"
}
rg_all[, rg := as.numeric(get(est_col))]
rg_all[, rg_se := as.numeric(get(se_col))]
rg_all[, `:=`(
  rg_CI_low = rg - 1.96 * rg_se,
  rg_CI_high = rg + 1.96 * rg_se
)]

full_rg <- rg_all[analysis == "full", rg][1]
rg_all[, `:=`(
  absolute_change_from_full = rg - full_rg,
  relative_change_from_full = (rg - full_rg) / full_rg
)]

fwrite(fit_all, file.path(OUT, "all_leave_one_model_fit.tsv"), sep = "\t")
fwrite(par_all, file.path(OUT, "all_leave_one_parameters.tsv"), sep = "\t")
fwrite(rg_all, file.path(OUT, "all_leave_one_axis_correlations.tsv"), sep = "\t")

print(fit_all)
print(rg_all[, .(
  analysis, omitted_trait, rg, rg_se, rg_CI_low, rg_CI_high,
  relative_change_from_full
)])
writeLines(capture.output(sessionInfo()), file.path(OUT, "R_sessionInfo.txt"))
