options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(GenomicSEM)
  library(data.table)
})

ROOT <- "/path/to/panvas_analysis/Genomic_SEM"
OUT <- "/path/to/panvas_analysis/Genomic_SEM_validation/02_three_model_CFA"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

AUDITED <- "/path/to/panvas_analysis/Genomic_SEM_validation/00_input_audit/LDSCoutput_audited.RData"
load(if (file.exists(AUDITED)) AUDITED else file.path(ROOT, "LDSCoutput.RData"))
stopifnot(exists("LDSCoutput"))

models <- list(
  single_common_factor = "
    G =~ NA*AA + AAA + TAAD + CAD + Stroke + PAD + CA + CAVS
    G ~~ 1*G
  ",
  correlated_two_factor = "
    AORTIC =~ NA*AA + AAA + TAAD
    ATHERO =~ NA*CAD + Stroke + PAD + CA + CAVS
    AORTIC ~~ ATHERO
    AORTIC ~~ 1*AORTIC
    ATHERO ~~ 1*ATHERO
    AAA ~~ 1e-6*AAA
  ",
  orthogonal_bifactor = "
    G =~ NA*AA + AAA + TAAD + CAD + Stroke + PAD + CA + CAVS
    AORTIC =~ NA*AA + AAA + TAAD
    ATHERO =~ NA*CAD + Stroke + PAD + CA + CAVS
    G ~~ 1*G
    AORTIC ~~ 1*AORTIC
    ATHERO ~~ 1*ATHERO
    G ~~ 0*AORTIC
    G ~~ 0*ATHERO
    AORTIC ~~ 0*ATHERO
    AAA ~~ 1e-6*AAA
  "
)

fit_one <- function(name, syntax) {
  warnings <- character()
  fit <- withCallingHandlers(
    usermodel(
      covstruc = LDSCoutput,
      estimation = "DWLS",
      model = syntax,
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
  fit_tab[, model := name]
  par_tab <- as.data.table(fit$results)
  par_tab[, model := name]

  model_dir <- file.path(OUT, name)
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(fit_tab, file.path(model_dir, "model_fit.tsv"), sep = "\t")
  fwrite(par_tab, file.path(model_dir, "model_parameters.tsv"), sep = "\t")
  writeLines(unique(warnings), file.path(model_dir, "warnings.txt"))

  list(fit = fit_tab, parameters = par_tab, warnings = warnings)
}

fits <- Map(fit_one, names(models), models)
fit_summary <- rbindlist(lapply(fits, `[[`, "fit"), fill = TRUE)
parameter_summary <- rbindlist(lapply(fits, `[[`, "parameters"), fill = TRUE)

fwrite(fit_summary, file.path(OUT, "three_model_fit_comparison.tsv"), sep = "\t")
fwrite(parameter_summary, file.path(OUT, "three_model_all_parameters.tsv"), sep = "\t")

# Bifactor general-factor strength using standardized genetic loadings.
bif <- fits$orthogonal_bifactor$parameters
loading_col <- if ("STD_Genotype" %in% names(bif)) "STD_Genotype" else "STD_All"
loadings <- bif[op == "=~", .(
  factor = lhs, trait = rhs, loading = as.numeric(get(loading_col))
)]

ss_general <- loadings[factor == "G", sum(loading^2, na.rm = TRUE)]
ss_specific <- loadings[factor != "G", sum(loading^2, na.rm = TRUE)]
ecv <- ss_general / (ss_general + ss_specific)

trait_share <- merge(
  loadings[factor == "G", .(trait, general_loading = loading)],
  loadings[factor != "G", .(trait, specific_loading = loading)],
  by = "trait", all.x = TRUE
)
trait_share[, general_common_variance_share :=
              general_loading^2 /
              (general_loading^2 + specific_loading^2)]

fwrite(data.table(
  sumsq_general = ss_general,
  sumsq_specific = ss_specific,
  ECV_general = ecv
), file.path(OUT, "bifactor_ECV.tsv"), sep = "\t")
fwrite(trait_share,
       file.path(OUT, "bifactor_trait_general_variance_share.tsv"), sep = "\t")

print(fit_summary)
print(data.table(ECV_general = ecv))
writeLines(capture.output(sessionInfo()), file.path(OUT, "R_sessionInfo.txt"))
