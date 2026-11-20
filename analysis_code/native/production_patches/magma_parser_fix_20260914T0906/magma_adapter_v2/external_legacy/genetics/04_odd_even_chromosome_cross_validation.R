options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(GenomicSEM)
  library(data.table)
  library(Matrix)
})

ROOT <- "/path/to/panvas_analysis/Genomic_SEM"
OUT <- "/path/to/panvas_analysis/Genomic_SEM_validation/04_odd_even_cross_validation"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

index <- fread(file.path(ROOT, "gwasfile_info.txt"))
traits <- index$phenotype
gwas_files <- file.path(ROOT, paste0(traits, ".sumstats.gz"))
stopifnot(all(file.exists(gwas_files)))

# Fill these values only after checking each original GWAS README.
# TRUE: the N in the munged file is effective case-control N.
# FALSE: the N is total participants or per-SNP total N.
N_IS_EFFECTIVE <- c(
  AA = FALSE,
  Stroke = NA,
  AAA = FALSE,
  CAVS = NA,
  CAD = NA,
  TAAD = FALSE,
  PAD = NA,
  CA = NA
)

if (anyNA(N_IS_EFFECTIVE[traits])) {
  stop(
    "Complete N_IS_EFFECTIVE after the sample-size audit. ",
    "Do not infer effective N only from ismeta=TRUE."
  )
}

sample_prev_counts <- index$cases / (index$cases + index$controls)
sample_prev <- ifelse(N_IS_EFFECTIVE[traits], 0.5, sample_prev_counts)
population_prev <- index$population.prev

ld_dir <- file.path(ROOT, "eur_w_ld_chr")

run_half_ldsc <- function(select_chr, label) {
  out <- ldsc(
    traits = gwas_files,
    sample.prev = sample_prev,
    population.prev = population_prev,
    ld = ld_dir,
    wld = ld_dir,
    trait.names = traits,
    stand = TRUE,
    select = select_chr,
    ldsc.log = file.path(OUT, paste0("LDSC_", label))
  )
  saveRDS(out, file.path(OUT, paste0("LDSC_", label, ".rds")))
  out
}

odd <- run_half_ldsc("ODD", "odd_chromosomes")
even <- run_half_ldsc("EVEN", "even_chromosomes")

efa_two <- function(x, label) {
  R <- as.matrix(nearPD(x$S_Stand, corr = TRUE)$mat)
  efa <- factanal(covmat = R, factors = 2, rotation = "promax")
  tab <- as.data.table(unclass(efa$loadings), keep.rownames = "trait")
  tab[, dataset := label]
  fwrite(tab, file.path(OUT, paste0("EFA_two_factor_", label, ".tsv")), sep = "\t")
  tab
}

model_correlated <- "
  AORTIC =~ NA*AA + AAA + TAAD
  ATHERO =~ NA*CAD + Stroke + PAD + CA + CAVS
  AORTIC ~~ ATHERO
  AORTIC ~~ 1*AORTIC
  ATHERO ~~ 1*ATHERO
  AAA ~~ 1e-6*AAA
"

model_bifactor <- "
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

fit_cfa <- function(x, model, data_label, model_label) {
  fit <- usermodel(
    covstruc = x,
    estimation = "DWLS",
    model = model,
    CFIcalc = TRUE,
    std.lv = TRUE,
    imp_cov = TRUE,
    fix_resid = FALSE,
    toler = 1e-20
  )
  fit_tab <- as.data.table(fit$modelfit)
  fit_tab[, `:=`(dataset = data_label, model = model_label)]
  par_tab <- as.data.table(fit$results)
  par_tab[, `:=`(dataset = data_label, model = model_label)]
  list(fit = fit_tab, parameters = par_tab)
}

# Primary split: discover the two axes on odd chromosomes, confirm on even.
efa_odd <- efa_two(odd, "odd")
confirm_even_cor <- fit_cfa(even, model_correlated, "even", "correlated_two_factor")
confirm_even_bif <- fit_cfa(even, model_bifactor, "even", "orthogonal_bifactor")

# Reverse split is a sensitivity analysis.
efa_even <- efa_two(even, "even")
confirm_odd_cor <- fit_cfa(odd, model_correlated, "odd", "correlated_two_factor")
confirm_odd_bif <- fit_cfa(odd, model_bifactor, "odd", "orthogonal_bifactor")

all_fit <- rbindlist(list(
  confirm_even_cor$fit, confirm_even_bif$fit,
  confirm_odd_cor$fit, confirm_odd_bif$fit
), fill = TRUE)
all_parameters <- rbindlist(list(
  confirm_even_cor$parameters, confirm_even_bif$parameters,
  confirm_odd_cor$parameters, confirm_odd_bif$parameters
), fill = TRUE)

fwrite(rbindlist(list(efa_odd, efa_even), fill = TRUE),
       file.path(OUT, "odd_even_EFA_loadings.tsv"), sep = "\t")
fwrite(all_fit, file.path(OUT, "odd_even_CFA_fit.tsv"), sep = "\t")
fwrite(all_parameters, file.path(OUT, "odd_even_CFA_parameters.tsv"), sep = "\t")

print(all_fit)
writeLines(capture.output(sessionInfo()), file.path(OUT, "R_sessionInfo.txt"))

