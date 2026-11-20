#!/usr/bin/env Rscript

# LD-aware Q_SNP locus filtering for the full bifactor G GWAS.
# The original GWAS files are read-only; all outputs are written to a new folder.

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table", repos = "https://cloud.r-project.org")
  }
  library(data.table)
})

options(scipen = 999)
setDTthreads(max(1L, min(4L, parallel::detectCores(logical = FALSE))))

# =============================================================================
# 1. Configuration: only PLINK_BIN and LD_REF_PREFIX may need to be changed
# =============================================================================

ROOT <- paste0(
  "/path/to/panvas_project/fsdownload/Genomic_SEM_complete_followup/",
  "05_bifactor_GWAS/full_bifactor_G"
)

FUMA_ALL <- file.path(ROOT, "full_bifactor_G_FUMA_all.tsv.gz")
FUMA_EXACT_Q_FILTERED <- file.path(
  ROOT,
  "full_bifactor_G_FUMA_excluding_QSNP_GWS.tsv.gz"
)
QSNP_GWS <- file.path(
  ROOT,
  "full_bifactor_G_QSNP_genomewide_significant.tsv"
)

OUT_DIR <- file.path(ROOT, "LD_aware_QSNP_filter")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# PLINK 1.9 executable. If Sys.which("plink") is empty, replace this with the
# absolute path to the downloaded PLINK 1.9 binary.
PLINK_BIN <- Sys.which("plink")
# Example:
# PLINK_BIN <- "/path/to/panvas_project/software/plink"

# Prefix of a FULL 1000 Genomes Phase 3 EUR, GRCh37 PLINK 1 binary reference.
# The following three files must exist:
#   <LD_REF_PREFIX>.bed
#   <LD_REF_PREFIX>.bim
#   <LD_REF_PREFIX>.fam
LD_REF_PREFIX <- "/path/to/panvas_project/reference/1000G_Phase3_EUR_GRCh37/EUR"

QSNP_P_THRESHOLD <- 5e-8
LD_R2_THRESHOLD <- 0.10
LD_WINDOW_KB <- 1000L
LOCUS_MERGE_KB <- 250L
PLINK_MEMORY_MB <- 8000L
PLINK_THREADS <- 4L

# Stop when fewer than this fraction of Q_SNP seeds are represented in the
# reference panel. A low value means that LD-aware filtering is incomplete.
MIN_REFERENCE_COVERAGE <- 0.80

# =============================================================================
# 2. Helpers and preflight checks
# =============================================================================

stop_missing <- function(paths, label) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(label, " missing:\n", paste(missing, collapse = "\n"))
  }
}

write_tsv_gz <- function(x, path) {
  tmp <- tempfile(pattern = "ldaware_", tmpdir = dirname(path), fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA")
  status <- system(
    sprintf("gzip -c %s > %s", shQuote(tmp), shQuote(path))
  )
  if (status != 0L || !file.exists(path)) {
    stop("Failed to create gzip output: ", path)
  }
}

merge_loci <- function(x, gap_bp) {
  if (!nrow(x)) return(x)
  setorder(x, CHR, START, END)
  x[, RUNNING_END := cummax(END), by = CHR]
  x[, NEW_LOCUS := {
    previous_end <- shift(RUNNING_END, fill = -Inf)
    seq_len(.N) == 1L | START > previous_end + gap_bp
  }, by = CHR]
  x[, LOCUS_INDEX := cumsum(NEW_LOCUS), by = CHR]

  ans <- x[, .(
    START = min(START),
    END = max(END),
    N_QSNP_SEEDS = uniqueN(SEED),
    QSNP_SEEDS = paste(sort(unique(SEED)), collapse = ",")
  ), by = .(CHR, LOCUS_INDEX)]
  ans[, LOCUS_ID := sprintf("QLOC_chr%s_%04d", CHR, LOCUS_INDEX)]
  setcolorder(
    ans,
    c("LOCUS_ID", "CHR", "START", "END", "N_QSNP_SEEDS", "QSNP_SEEDS")
  )
  ans[]
}

stop_missing(
  c(FUMA_ALL, FUMA_EXACT_Q_FILTERED, QSNP_GWS),
  "Required Genomic SEM output"
)

if (!nzchar(PLINK_BIN) || !file.exists(PLINK_BIN)) {
  stop(
    "PLINK 1.9 was not found. Download the macOS PLINK 1.9 binary, then set ",
    "PLINK_BIN to its absolute path."
  )
}

ref_files <- paste0(LD_REF_PREFIX, c(".bed", ".bim", ".fam"))
stop_missing(ref_files, "1000G Phase 3 EUR GRCh37 reference")

cat("PLINK: ", PLINK_BIN, "\n", sep = "")
cat("Reference: ", LD_REF_PREFIX, "\n", sep = "")
cat("Output: ", OUT_DIR, "\n\n", sep = "")

# =============================================================================
# 3. Read Q_SNP seeds and audit genome build/reference coverage
# =============================================================================

q <- fread(QSNP_GWS)
required_q <- c("SNP", "CHR", "BP", "Q_SNP_P")
if (length(setdiff(required_q, names(q)))) {
  stop("Q_SNP file lacks columns: ", paste(setdiff(required_q, names(q)), collapse = ", "))
}

q <- unique(
  q[
    is.finite(Q_SNP_P) & Q_SNP_P < QSNP_P_THRESHOLD,
    .(SNP = as.character(SNP), CHR = as.integer(CHR),
      BP = as.integer(BP), Q_SNP_P)
  ],
  by = "SNP"
)
if (!nrow(q)) stop("No genome-wide significant Q_SNPs were found.")

bim <- fread(
  paste0(LD_REF_PREFIX, ".bim"),
  col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"),
  select = c(1L, 2L, 4L)
)
bim[, `:=`(
  SNP = as.character(SNP),
  CHR = as.integer(CHR),
  BP = as.integer(BP)
)]
bim <- unique(bim, by = "SNP")

q_ref <- merge(
  q,
  bim[, .(SNP, REF_CHR = CHR, REF_BP = BP)],
  by = "SNP"
)
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
  coordinate_match_fraction = coordinate_match
)
fwrite(
  coverage_qc,
  file.path(OUT_DIR, "01_QSNP_reference_coverage.tsv"),
  sep = "\t"
)
fwrite(
  q[!SNP %chin% q_ref$SNP],
  file.path(OUT_DIR, "01_QSNP_seeds_absent_from_reference.tsv"),
  sep = "\t"
)

print(coverage_qc)
if (coordinate_match < 0.95) {
  stop(
    "Fewer than 95% of shared Q_SNPs have identical CHR/BP coordinates. ",
    "The LD panel is probably not GRCh37 or has incompatible variant IDs."
  )
}
if (reference_coverage < MIN_REFERENCE_COVERAGE) {
  stop(
    sprintf(
      "Only %.1f%% of Q_SNP seeds occur in the LD reference; require >= %.1f%%. ",
      100 * reference_coverage,
      100 * MIN_REFERENCE_COVERAGE
    ),
    "Use a fuller 1000G Phase 3 EUR GRCh37 panel."
  )
}

q_seed_list <- file.path(OUT_DIR, "02_QSNP_seeds_in_reference.txt")
fwrite(
  q_ref[, .(SNP)],
  q_seed_list,
  col.names = FALSE,
  quote = FALSE
)

# =============================================================================
# 4. Calculate LD proxies for every Q_SNP seed with PLINK 1.9
# =============================================================================

plink_out <- file.path(OUT_DIR, "03_QSNP_LD")
plink_args <- c(
  "--bfile", LD_REF_PREFIX,
  "--ld-snp-list", q_seed_list,
  "--r2", "gz",
  "--ld-window", "999999",
  "--ld-window-kb", as.character(LD_WINDOW_KB),
  "--ld-window-r2", as.character(LD_R2_THRESHOLD),
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

ld_path <- paste0(plink_out, ".ld.gz")
if (!file.exists(ld_path)) {
  stop("Expected PLINK LD output not found: ", ld_path)
}

ld <- fread(ld_path)
required_ld <- c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2")
if (length(setdiff(required_ld, names(ld)))) {
  stop(
    "Unexpected PLINK LD columns: ",
    paste(names(ld), collapse = ", ")
  )
}
ld <- ld[is.finite(R2) & R2 >= LD_R2_THRESHOLD]

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
links <- unique(rbindlist(list(left_links, right_links), use.names = TRUE))

# Include each seed itself. Seeds missing from the panel are still removed at
# their exact coordinate, but cannot be expanded using LD.
self_links <- q[, .(
  SEED = SNP,
  CHR,
  SEED_BP = BP,
  PROXY = SNP,
  PROXY_BP = BP,
  R2 = 1
)]
links <- unique(rbindlist(list(links, self_links), use.names = TRUE))
write_tsv_gz(
  links,
  file.path(
    OUT_DIR,
    sprintf(
      "03_QSNP_seed_proxy_pairs_r2ge%s.tsv.gz",
      format(LD_R2_THRESHOLD, scientific = FALSE)
    )
  )
)

# =============================================================================
# 5. Build LD-aware Q_SNP loci and merge loci within 250 kb
# =============================================================================

seed_intervals <- links[, .(
  START = min(c(SEED_BP, PROXY_BP), na.rm = TRUE),
  END = max(c(SEED_BP, PROXY_BP), na.rm = TRUE),
  N_PROXIES = uniqueN(PROXY),
  MAX_R2 = max(R2, na.rm = TRUE)
), by = .(SEED, CHR)]

merged_loci <- merge_loci(
  seed_intervals[, .(SEED, CHR, START, END)],
  gap_bp = LOCUS_MERGE_KB * 1000
)
fwrite(
  seed_intervals,
  file.path(OUT_DIR, "04_QSNP_seed_LD_intervals.tsv"),
  sep = "\t"
)
fwrite(
  merged_loci,
  file.path(OUT_DIR, "04_QSNP_merged_LD_aware_loci.tsv"),
  sep = "\t"
)

cat(
  "\nConstructed ", nrow(seed_intervals), " seed intervals and ",
  nrow(merged_loci), " merged Q_SNP loci.\n",
  sep = ""
)

# =============================================================================
# 6. Filter the full factor GWAS without constructing an LD matrix in R
# =============================================================================

cat("\nReading full FUMA GWAS...\n")
gwas <- fread(cmd = paste("gzip -dc", shQuote(FUMA_ALL)))
required_gwas <- c("SNP", "CHR", "BP", "A1", "A2", "MAF", "BETA", "SE", "P", "N")
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

exclude_locus <- rep(FALSE, nrow(gwas))
for (chr in sort(unique(merged_loci$CHR))) {
  row_index <- which(gwas$CHR == chr)
  intervals <- merged_loci[CHR == chr][order(START)]
  if (!length(row_index) || !nrow(intervals)) next

  interval_index <- findInterval(gwas$BP[row_index], intervals$START)
  valid_index <- interval_index > 0L
  local_hit <- rep(FALSE, length(row_index))
  local_hit[valid_index] <- (
    gwas$BP[row_index][valid_index] <= intervals$END[interval_index[valid_index]]
  )
  exclude_locus[row_index] <- local_hit
}

# Always remove exact significant Q_SNP IDs, including seeds absent from the
# reference panel.
exclude_exact_q <- gwas$SNP %chin% q$SNP
exclude_any <- exclude_locus | exclude_exact_q

filtered <- gwas[!exclude_any]
excluded <- gwas[exclude_any]

filtered_path <- file.path(
  OUT_DIR,
  "full_bifactor_G_FUMA_LDaware_QSNP_locus_filtered.tsv.gz"
)
excluded_path <- file.path(
  OUT_DIR,
  "full_bifactor_G_variants_excluded_by_QSNP_loci.tsv.gz"
)

cat("Writing filtered FUMA input...\n")
write_tsv_gz(filtered, filtered_path)
write_tsv_gz(excluded, excluded_path)

# =============================================================================
# 7. Final QC and audit trail
# =============================================================================

exact_filtered_n <- nrow(
  fread(
    cmd = paste("gzip -dc", shQuote(FUMA_EXACT_Q_FILTERED)),
    select = "SNP"
  )
)

qc <- data.table(
  parameter = c(
    "genome_build",
    "LD_reference",
    "Q_SNP_P_threshold",
    "LD_r2_threshold",
    "LD_window_kb",
    "locus_merge_kb",
    "n_QSNP_seeds",
    "n_QSNP_seeds_in_reference",
    "Q_SNP_reference_coverage",
    "coordinate_match_fraction",
    "n_merged_QSNP_loci",
    "n_input_GWAS_variants",
    "n_exact_Q_filtered_variants",
    "n_LDaware_locus_excluded_variants",
    "n_output_GWAS_variants",
    "n_input_GWS_variants",
    "n_output_GWS_variants"
  ),
  value = as.character(c(
    "GRCh37",
    LD_REF_PREFIX,
    QSNP_P_THRESHOLD,
    LD_R2_THRESHOLD,
    LD_WINDOW_KB,
    LOCUS_MERGE_KB,
    nrow(q),
    nrow(q_ref),
    reference_coverage,
    coordinate_match,
    nrow(merged_loci),
    nrow(gwas),
    nrow(gwas) - exact_filtered_n,
    sum(exclude_any),
    nrow(filtered),
    sum(gwas$P < 5e-8, na.rm = TRUE),
    sum(filtered$P < 5e-8, na.rm = TRUE)
  ))
)
fwrite(qc, file.path(OUT_DIR, "05_LD_aware_filter_QC.tsv"), sep = "\t")

cat("\nLD-aware Q_SNP locus filtering completed.\n")
cat("Primary filtered FUMA file:\n", filtered_path, "\n\n", sep = "")
print(qc)
