#!/usr/bin/env Rscript

# =============================================================================
# LD-aware Q_SNP locus filtering for the leave-one-CAD bifactor G GWAS
#
# Primary definition:
#   Q_SNP P < 5e-8
#   1000G Phase 3 EUR, GRCh37, MAF >= 0.01
#   LD proxy r2 >= 0.10 within +/-1 Mb
#   Overlapping Q_SNP intervals, or intervals separated by <=250 kb, are merged
#
# Sensitivity definition:
#   Identical workflow with LD proxy r2 >= 0.60
#
# Original inputs and reference files are never modified.
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table", repos = "https://cloud.r-project.org")
  }
  library(data.table)
})

options(scipen = 999)
setDTthreads(max(1L, min(4L, parallel::detectCores(logical = FALSE))))

# =============================================================================
# 1. Fixed local paths
# =============================================================================

FUMA_ROOT <- "/path/to/FUMA"
ROOT <- file.path(FUMA_ROOT, "LOT_CAD")

FUMA_ALL <- file.path(ROOT, "leave_one_CAD_bifactor_G_FUMA_all.tsv.gz")
FUMA_EXACT_Q_FILTERED <- file.path(
  ROOT,
  "leave_one_CAD_bifactor_G_FUMA_excluding_QSNP_GWS.tsv.gz"
)
QSNP_GWS <- file.path(
  ROOT,
  "leave_one_CAD_bifactor_G_QSNP_genomewide_significant.tsv"
)

PLINK_BIN <- file.path(FUMA_ROOT, "software", "plink1.9", "plink")
LD_REF_PREFIX <- file.path(
  FUMA_ROOT,
  "reference",
  "1000G_EUR_GRCh37_MAF01",
  "data_maf0.01_rs_ref"
)

OUT_DIR <- file.path(ROOT, "LD_aware_QSNP_filter")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

QSNP_P_THRESHOLD <- 5e-8
PRIMARY_LD_R2 <- 0.10
SENSITIVITY_LD_R2 <- 0.60
LD_WINDOW_KB <- 1000L
LOCUS_MERGE_KB <- 250L
PLINK_MEMORY_MB <- 8000L
PLINK_THREADS <- 4L
MIN_REFERENCE_COVERAGE <- 0.95
MIN_COORDINATE_MATCH <- 0.95
FORCE_RECOMPUTE_LD <- FALSE

# =============================================================================
# 2. Helpers
# =============================================================================

stop_missing <- function(paths, label) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(label, " missing:\n", paste(missing, collapse = "\n"))
  }
}

write_tsv <- function(x, path) {
  fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
}

write_tsv_gz <- function(x, path) {
  fwrite(
    x,
    path,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = "gzip",
    compressLevel = 6L
  )
}

merge_loci <- function(seed_intervals, merge_gap_bp) {
  x <- copy(seed_intervals)
  if (!nrow(x)) return(x)

  setorder(x, CHR, START, END)
  x[, RUNNING_END := cummax(END), by = CHR]
  x[, NEW_LOCUS := {
    previous_end <- shift(RUNNING_END, fill = NA_integer_)
    seq_len(.N) == 1L |
      (!is.na(previous_end) & START > previous_end + merge_gap_bp)
  }, by = CHR]
  x[, LOCUS_INDEX := cumsum(NEW_LOCUS), by = CHR]

  ans <- x[, .(
    START = min(START),
    END = max(END),
    LOCUS_SIZE_BP = max(END) - min(START) + 1L,
    N_QSNP_SEEDS = uniqueN(SEED),
    N_LD_PROXIES = uniqueN(PROXY),
    QSNP_SEEDS = paste(sort(unique(SEED)), collapse = ",")
  ), by = .(CHR, LOCUS_INDEX)]

  ans[, LOCUS_ID := sprintf("QLOC_chr%s_%04d", CHR, LOCUS_INDEX)]
  setcolorder(
    ans,
    c(
      "LOCUS_ID", "CHR", "START", "END", "LOCUS_SIZE_BP",
      "N_QSNP_SEEDS", "N_LD_PROXIES", "QSNP_SEEDS"
    )
  )
  ans[]
}

mark_variants_in_loci <- function(gwas, loci) {
  hit <- rep(FALSE, nrow(gwas))
  if (!nrow(loci)) return(hit)

  for (chr in sort(unique(loci$CHR))) {
    row_index <- which(gwas$CHR == chr)
    intervals <- loci[CHR == chr][order(START)]
    if (!length(row_index) || !nrow(intervals)) next

    interval_index <- findInterval(gwas$BP[row_index], intervals$START)
    valid_index <- interval_index > 0L
    local_hit <- rep(FALSE, length(row_index))
    local_hit[valid_index] <- (
      gwas$BP[row_index][valid_index] <=
        intervals$END[interval_index[valid_index]]
    )
    hit[row_index] <- local_hit
  }
  hit
}

build_filter_definition <- function(all_links, q_seeds, r2_threshold) {
  links <- all_links[is.finite(R2) & R2 >= r2_threshold]

  self_links <- q_seeds[, .(
    SEED = SNP,
    CHR,
    SEED_BP = BP,
    PROXY = SNP,
    PROXY_BP = BP,
    R2 = 1
  )]
  links <- unique(
    rbindlist(list(links, self_links), use.names = TRUE, fill = TRUE)
  )

  seed_intervals <- links[, .(
    START = min(c(SEED_BP, PROXY_BP), na.rm = TRUE),
    END = max(c(SEED_BP, PROXY_BP), na.rm = TRUE),
    N_PROXIES = uniqueN(PROXY),
    MIN_INCLUDED_R2 = min(R2, na.rm = TRUE),
    MAX_R2 = max(R2, na.rm = TRUE)
  ), by = .(SEED, CHR)]

  loci <- merge_loci(
    seed_intervals = links[, .(
      SEED, CHR,
      START = pmin(SEED_BP, PROXY_BP),
      END = pmax(SEED_BP, PROXY_BP),
      PROXY
    )],
    merge_gap_bp = LOCUS_MERGE_KB * 1000L
  )

  list(
    links = links,
    seed_intervals = seed_intervals,
    loci = loci
  )
}

write_filter_result <- function(
    gwas,
    q_seeds,
    definition,
    r2_threshold,
    label
) {
  result_dir <- file.path(OUT_DIR, label)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  write_tsv_gz(
    definition$links,
    file.path(result_dir, "Q_SNP_seed_proxy_pairs.tsv.gz")
  )
  write_tsv(
    definition$seed_intervals,
    file.path(result_dir, "Q_SNP_seed_LD_intervals.tsv")
  )
  write_tsv(
    definition$loci,
    file.path(result_dir, "Q_SNP_merged_LD_aware_loci.tsv")
  )

  exclude_locus <- mark_variants_in_loci(gwas, definition$loci)
  exclude_exact_q <- gwas$SNP %chin% q_seeds$SNP
  exclude_any <- exclude_locus | exclude_exact_q

  filtered_path <- file.path(
    result_dir,
    paste0(
      "leave_one_CAD_bifactor_G_FUMA_LDaware_QSNP_locus_filtered_",
      label,
      ".tsv.gz"
    )
  )
  excluded_path <- file.path(
    result_dir,
    paste0(
      "leave_one_CAD_bifactor_G_variants_excluded_by_QSNP_loci_",
      label,
      ".tsv.gz"
    )
  )

  write_tsv_gz(gwas[!exclude_any], filtered_path)
  write_tsv_gz(
    gwas[
      exclude_any,
      .(
        SNP, CHR, BP, A1, A2, MAF, BETA, SE, P, N,
        excluded_by_exact_QSNP = exclude_exact_q[exclude_any],
        excluded_by_QSNP_locus = exclude_locus[exclude_any]
      )
    ],
    excluded_path
  )

  qc <- data.table(
    analysis = label,
    LD_r2_threshold = r2_threshold,
    LD_window_kb = LD_WINDOW_KB,
    locus_merge_kb = LOCUS_MERGE_KB,
    n_QSNP_seeds = nrow(q_seeds),
    n_LD_pairs_including_self = nrow(definition$links),
    n_merged_QSNP_loci = nrow(definition$loci),
    n_input_GWAS_variants = nrow(gwas),
    n_exact_QSNP_variants_in_FUMA = sum(exclude_exact_q),
    n_variants_in_QSNP_loci = sum(exclude_locus),
    n_total_excluded_variants = sum(exclude_any),
    n_output_GWAS_variants = sum(!exclude_any),
    n_input_GWS_variants = sum(gwas$P < 5e-8, na.rm = TRUE),
    n_excluded_GWS_variants = sum(exclude_any & gwas$P < 5e-8, na.rm = TRUE),
    n_output_GWS_variants = sum(!exclude_any & gwas$P < 5e-8, na.rm = TRUE),
    output_FUMA_file = filtered_path
  )
  write_tsv(qc, file.path(result_dir, "filter_QC.tsv"))
  qc
}

# =============================================================================
# 3. Preflight checks
# =============================================================================

stop_missing(
  c(FUMA_ALL, FUMA_EXACT_Q_FILTERED, QSNP_GWS),
  "Required Genomic SEM/FUMA input"
)
stop_missing(PLINK_BIN, "PLINK 1.9 executable")
stop_missing(
  paste0(LD_REF_PREFIX, c(".bed", ".bim", ".fam")),
  "1000G Phase 3 EUR GRCh37 reference"
)

plink_version <- system2(PLINK_BIN, "--version", stdout = TRUE)
cat("PLINK: ", paste(plink_version, collapse = " "), "\n", sep = "")
cat("Reference: ", LD_REF_PREFIX, "\n", sep = "")
cat("Output directory: ", OUT_DIR, "\n\n", sep = "")

# =============================================================================
# 4. Read Q_SNP seeds and verify the reference panel
# =============================================================================

q <- fread(QSNP_GWS)
required_q <- c("SNP", "CHR", "BP", "Q_SNP_P")
if (length(setdiff(required_q, names(q)))) {
  stop(
    "Q_SNP file lacks columns: ",
    paste(setdiff(required_q, names(q)), collapse = ", ")
  )
}

q <- unique(
  q[
    is.finite(Q_SNP_P) & Q_SNP_P < QSNP_P_THRESHOLD,
    .(
      SNP = as.character(SNP),
      CHR = as.integer(CHR),
      BP = as.integer(BP),
      Q_SNP_P
    )
  ],
  by = "SNP"
)
if (!nrow(q)) stop("No genome-wide significant Q_SNPs were found.")

bim <- fread(paste0(LD_REF_PREFIX, ".bim"), header = FALSE)
if (ncol(bim) != 6L) stop("The reference BIM file does not have six columns.")
setnames(bim, c("REF_CHR", "SNP", "CM", "REF_BP", "A1", "A2"))
bim <- unique(
  bim[, .(
    SNP = as.character(SNP),
    REF_CHR = as.integer(REF_CHR),
    REF_BP = as.integer(REF_BP)
  )],
  by = "SNP"
)

q_ref <- merge(q, bim, by = "SNP")
reference_coverage <- nrow(q_ref) / nrow(q)
coordinate_match <- if (nrow(q_ref)) {
  mean(q_ref$CHR == q_ref$REF_CHR & q_ref$BP == q_ref$REF_BP)
} else {
  0
}

coverage_qc <- data.table(
  n_QSNP_seeds = nrow(q),
  n_QSNP_seeds_in_reference = nrow(q_ref),
  reference_coverage = reference_coverage,
  coordinate_match_fraction = coordinate_match,
  reference_sample_count = nrow(
    fread(paste0(LD_REF_PREFIX, ".fam"), header = FALSE)
  ),
  reference_variant_count = nrow(bim)
)
write_tsv(
  coverage_qc,
  file.path(OUT_DIR, "01_QSNP_reference_coverage.tsv")
)
write_tsv(
  q[!SNP %chin% q_ref$SNP],
  file.path(OUT_DIR, "01_QSNP_seeds_absent_from_reference.tsv")
)
print(coverage_qc)

if (reference_coverage < MIN_REFERENCE_COVERAGE) {
  stop(
    sprintf(
      "Q_SNP reference coverage is %.2f%%, below the %.2f%% requirement.",
      100 * reference_coverage,
      100 * MIN_REFERENCE_COVERAGE
    )
  )
}
if (coordinate_match < MIN_COORDINATE_MATCH) {
  stop(
    sprintf(
      "Coordinate match is %.2f%%, below the %.2f%% requirement. ",
      100 * coordinate_match,
      100 * MIN_COORDINATE_MATCH
    ),
    "The reference genome build is likely incompatible."
  )
}

q_seed_list <- file.path(OUT_DIR, "02_QSNP_seeds_in_reference.txt")
write_tsv(q_ref[, .(SNP)], q_seed_list)

# PLINK requires a one-column SNP list without a header.
fwrite(
  q_ref[, .(SNP)],
  q_seed_list,
  col.names = FALSE,
  quote = FALSE
)

# =============================================================================
# 5. Calculate all Q_SNP proxy pairs once at r2 >= 0.10
# =============================================================================

plink_out <- file.path(OUT_DIR, "03_QSNP_LD_r2ge0.10")
ld_path <- paste0(plink_out, ".ld.gz")

if (FORCE_RECOMPUTE_LD || !file.exists(ld_path)) {
  plink_args <- c(
    "--bfile", LD_REF_PREFIX,
    "--ld-snp-list", q_seed_list,
    "--r2", "gz",
    "--ld-window", "999999",
    "--ld-window-kb", as.character(LD_WINDOW_KB),
    "--ld-window-r2", as.character(PRIMARY_LD_R2),
    "--memory", as.character(PLINK_MEMORY_MB),
    "--threads", as.character(PLINK_THREADS),
    "--allow-no-sex",
    "--out", plink_out
  )

  cat("\nRunning PLINK LD calculation...\n")
  status <- system2(PLINK_BIN, plink_args)
  if (status != 0L) {
    stop("PLINK failed. Inspect: ", paste0(plink_out, ".log"))
  }
} else {
  cat("\nReusing existing PLINK LD output: ", ld_path, "\n", sep = "")
}

if (!file.exists(ld_path)) {
  stop("Expected PLINK LD output was not created: ", ld_path)
}

ld <- fread(ld_path)
required_ld <- c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2")
if (length(setdiff(required_ld, names(ld)))) {
  stop(
    "Unexpected PLINK LD columns: ",
    paste(names(ld), collapse = ", ")
  )
}
ld <- ld[is.finite(R2) & R2 >= PRIMARY_LD_R2]

q_seed_set <- q_ref$SNP
left_links <- ld[
  SNP_A %chin% q_seed_set,
  .(
    SEED = as.character(SNP_A),
    CHR = as.integer(CHR_A),
    SEED_BP = as.integer(BP_A),
    PROXY = as.character(SNP_B),
    PROXY_BP = as.integer(BP_B),
    R2
  )
]
right_links <- ld[
  SNP_B %chin% q_seed_set,
  .(
    SEED = as.character(SNP_B),
    CHR = as.integer(CHR_B),
    SEED_BP = as.integer(BP_B),
    PROXY = as.character(SNP_A),
    PROXY_BP = as.integer(BP_A),
    R2
  )
]
all_links <- unique(
  rbindlist(list(left_links, right_links), use.names = TRUE)
)

write_tsv_gz(
  all_links,
  file.path(OUT_DIR, "03_all_QSNP_proxy_pairs_r2ge0.10.tsv.gz")
)

# =============================================================================
# 6. Read the complete factor GWAS and create both filtered analyses
# =============================================================================

cat("\nReading leave-one-CAD factor GWAS...\n")
gwas <- fread(FUMA_ALL)
required_gwas <- c(
  "SNP", "CHR", "BP", "A1", "A2", "MAF", "BETA", "SE", "P", "N"
)
if (length(setdiff(required_gwas, names(gwas)))) {
  stop(
    "FUMA input lacks columns: ",
    paste(setdiff(required_gwas, names(gwas)), collapse = ", ")
  )
}

gwas[, `:=`(
  SNP = as.character(SNP),
  CHR = as.integer(CHR),
  BP = as.integer(BP)
)]

if (anyDuplicated(gwas$SNP)) stop("Duplicate SNP IDs detected in FUMA_ALL.")
if (any(!is.finite(gwas$P) | gwas$P <= 0 | gwas$P > 1)) {
  stop("Invalid P values detected in FUMA_ALL.")
}

primary_definition <- build_filter_definition(
  all_links = all_links,
  q_seeds = q,
  r2_threshold = PRIMARY_LD_R2
)
primary_qc <- write_filter_result(
  gwas = gwas,
  q_seeds = q,
  definition = primary_definition,
  r2_threshold = PRIMARY_LD_R2,
  label = "primary_r2_0.10"
)

sensitivity_definition <- build_filter_definition(
  all_links = all_links,
  q_seeds = q,
  r2_threshold = SENSITIVITY_LD_R2
)
sensitivity_qc <- write_filter_result(
  gwas = gwas,
  q_seeds = q,
  definition = sensitivity_definition,
  r2_threshold = SENSITIVITY_LD_R2,
  label = "sensitivity_r2_0.60"
)

# =============================================================================
# 7. Cross-check against the previous exact-Q_SNP-only file
# =============================================================================

exact_filtered <- fread(
  FUMA_EXACT_Q_FILTERED,
  select = c("SNP", "CHR", "BP", "P")
)
exact_filtered[, SNP := as.character(SNP)]

comparison_qc <- data.table(
  version = c(
    "unfiltered",
    "exact_QSNP_only",
    "LDaware_primary_r2_0.10",
    "LDaware_sensitivity_r2_0.60"
  ),
  n_variants = c(
    nrow(gwas),
    nrow(exact_filtered),
    primary_qc$n_output_GWAS_variants,
    sensitivity_qc$n_output_GWAS_variants
  ),
  n_GWS_variants = c(
    sum(gwas$P < 5e-8, na.rm = TRUE),
    sum(exact_filtered$P < 5e-8, na.rm = TRUE),
    primary_qc$n_output_GWS_variants,
    sensitivity_qc$n_output_GWS_variants
  )
)
write_tsv(
  comparison_qc,
  file.path(OUT_DIR, "04_all_filter_versions_comparison.tsv")
)

all_qc <- rbindlist(list(primary_qc, sensitivity_qc), fill = TRUE)
write_tsv(all_qc, file.path(OUT_DIR, "05_all_LD_aware_filter_QC.tsv"))

methods_parameters <- data.table(
  parameter = c(
    "genome_build",
    "ancestry",
    "LD_reference",
    "reference_samples",
    "reference_variants",
    "reference_MAF_minimum",
    "Q_SNP_P_threshold",
    "primary_LD_r2_threshold",
    "sensitivity_LD_r2_threshold",
    "LD_window_kb",
    "locus_merge_kb",
    "PLINK_version",
    "script"
  ),
  value = c(
    "GRCh37",
    "European",
    LD_REF_PREFIX,
    coverage_qc$reference_sample_count,
    coverage_qc$reference_variant_count,
    "0.01",
    QSNP_P_THRESHOLD,
    PRIMARY_LD_R2,
    SENSITIVITY_LD_R2,
    LD_WINDOW_KB,
    LOCUS_MERGE_KB,
    paste(plink_version, collapse = " "),
    file.path(ROOT, "09_leave_one_CAD_LD_aware_QSNP_filter_local.R")
  )
)

write_tsv(
  methods_parameters,
  file.path(OUT_DIR, "06_methods_parameters.tsv")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_DIR, "07_R_sessionInfo.txt")
)

cat("\nLD-aware Q_SNP locus filtering completed successfully.\n\n")
cat("Primary FUMA/MAGMA input:\n")
cat(primary_qc$output_FUMA_file, "\n\n")
cat("Sensitivity FUMA/MAGMA input:\n")
cat(sensitivity_qc$output_FUMA_file, "\n\n")
print(comparison_qc)
print(all_qc)
