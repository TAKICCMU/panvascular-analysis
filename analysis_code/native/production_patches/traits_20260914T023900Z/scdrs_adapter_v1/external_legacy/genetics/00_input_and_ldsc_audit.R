options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

ROOT <- "/path/to/panvas_analysis/Genomic_SEM"
OUT <- "/path/to/panvas_analysis/Genomic_SEM_validation/00_input_audit"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

index <- fread(file.path(ROOT, "gwasfile_info.txt"))
required <- c("filename", "phenotype", "cases", "controls",
              "samplesize", "ismeta", "population.prev")
stopifnot(all(required %in% names(index)))

index[, sample_prev_from_counts := cases / (cases + controls)]
index[, effective_N_from_counts := 4 / (1 / cases + 1 / controls)]

# Inspect the N actually stored in each munged LDSC file.
inspect_munged_N <- function(trait) {
  path <- file.path(ROOT, paste0(trait, ".sumstats.gz"))
  if (!file.exists(path)) {
    return(data.table(
      phenotype = trait, munged_file = path, file_exists = FALSE,
      n_rows_read = NA_integer_, N_min = NA_real_, N_q25 = NA_real_,
      N_median = NA_real_, N_mean = NA_real_, N_q75 = NA_real_,
      N_max = NA_real_
    ))
  }

  x <- tryCatch(
    fread(path, select = "N", nrows = 200000, showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(x) || !"N" %in% names(x)) {
    return(data.table(
      phenotype = trait, munged_file = path, file_exists = TRUE,
      n_rows_read = NA_integer_, N_min = NA_real_, N_q25 = NA_real_,
      N_median = NA_real_, N_mean = NA_real_, N_q75 = NA_real_,
      N_max = NA_real_
    ))
  }

  n <- as.numeric(x$N)
  q <- quantile(n, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  data.table(
    phenotype = trait, munged_file = path, file_exists = TRUE,
    n_rows_read = length(n), N_min = q[1], N_q25 = q[2],
    N_median = q[3], N_mean = mean(n, na.rm = TRUE),
    N_q75 = q[4], N_max = q[5]
  )
}

n_audit <- rbindlist(lapply(index$phenotype, inspect_munged_N), fill = TRUE)
audit <- merge(index, n_audit, by = "phenotype", all.x = TRUE)

# This column must be completed from each GWAS README.
# "effective": N is the sum of effective case-control sample sizes.
# "total": N is total participants or per-SNP total participants.
audit[, N_type_confirmed := fifelse(
  ismeta,
  "VERIFY_FROM_GWAS_README",
  "total"
)]
audit[, proposed_sample_prev := fifelse(
  N_type_confirmed == "effective",
  0.5,
  sample_prev_from_counts
)]

fwrite(audit, file.path(OUT, "sample_size_and_prevalence_audit.tsv"), sep = "\t")

cat("\nImportant:\n")
cat("Do not set sample.prev=0.5 solely because ismeta=TRUE.\n")
cat("Use 0.5 only when the N supplied to munge/ldsc is an effective N.\n")
cat("Complete N_type_confirmed using the original GWAS documentation before rerunning LDSC.\n\n")
print(audit[, .(
  phenotype, cases, controls, samplesize, ismeta,
  sample_prev_from_counts, effective_N_from_counts,
  N_median, N_type_confirmed, population.prev
)])

writeLines(capture.output(sessionInfo()), file.path(OUT, "R_sessionInfo.txt"))

