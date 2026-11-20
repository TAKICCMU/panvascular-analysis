options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(GenomicSEM)
  library(data.table)
})

ROOT <- "/path/to/panvas_analysis/Genomic_SEM"
OUT <- "/path/to/panvas_analysis/Genomic_SEM_validation/00_input_audit"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

index <- fread(file.path(ROOT, "gwasfile_info.txt"))
traits <- index$phenotype
gwas_files <- file.path(ROOT, paste0(traits, ".sumstats.gz"))
stopifnot(all(file.exists(gwas_files)))

# Complete this vector from the original GWAS documentation.
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
  stop("Complete N_IS_EFFECTIVE before rerunning LDSC.")
}

sample_prev_counts <- index$cases / (index$cases + index$controls)
sample_prev <- ifelse(N_IS_EFFECTIVE[traits], 0.5, sample_prev_counts)

LDSCoutput <- ldsc(
  traits = gwas_files,
  sample.prev = sample_prev,
  population.prev = index$population.prev,
  ld = file.path(ROOT, "eur_w_ld_chr"),
  wld = file.path(ROOT, "eur_w_ld_chr"),
  trait.names = traits,
  stand = TRUE,
  ldsc.log = file.path(OUT, "audited_full_8_trait")
)

save(LDSCoutput, file = file.path(OUT, "LDSCoutput_audited.RData"))
fwrite(data.table(
  phenotype = traits,
  N_is_effective = N_IS_EFFECTIVE[traits],
  sample_prev_used = sample_prev,
  population_prev_used = index$population.prev
), file.path(OUT, "audited_LDSC_prevalence_inputs.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "audited_R_sessionInfo.txt"))

