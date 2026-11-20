#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

INPUT_ROOT <- "/path/to/downloads/spatial_validation_v5_locked71_spatial_statistics_audit_v1"
OUTPUT_ROOT <- "/path/to/panvas_work/reviewer_priority_analyses_20260826/04_spatial_matched_gene_null_corrected"
dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

bh <- function(p) p.adjust(p, method = "BH")

section <- read.delim(
  gzfile(file.path(INPUT_ROOT, "04_matched_gene_associations_by_section.tsv.gz")),
  check.names = FALSE
)
roster <- read.delim(
  gzfile(file.path(INPUT_ROOT, "03_expression_detection_matched_gene_roster.tsv.gz")),
  check.names = FALSE
)

section$is_target_gene <- as.logical(section$is_target_gene)
section$association_z <- suppressWarnings(as.numeric(section$association_z))
section <- section[is.finite(section$association_z), , drop = FALSE]

rank_map <- unique(roster[c("GSE", "target_gene", "matched_gene", "match_rank")])
names(rank_map)[names(rank_map) == "matched_gene"] <- "tested_gene"

null_rows <- merge(
  section[!section$is_target_gene, , drop = FALSE],
  rank_map,
  by = c("GSE", "target_gene", "tested_gene"),
  all.x = TRUE,
  sort = FALSE
)
target_rows <- section[section$is_target_gene, , drop = FALSE]
target_rows$match_rank <- 0L

analysis <- rbind(
  target_rows[names(null_rows)],
  null_rows
)
analysis$match_rank <- suppressWarnings(as.integer(analysis$match_rank))
analysis <- analysis[!is.na(analysis$match_rank), , drop = FALSE]

# Technical sections from one biological unit are averaged before any
# cross-unit statistic, so multiple sections cannot increase its weight.
unit <- aggregate(
  association_z ~ pair_id + target_gene + GSE + biological_unit_id +
    match_rank + is_target_gene,
  data = analysis,
  FUN = mean,
  na.rm = TRUE
)

pair_ids <- unique(unit$pair_id)
null_details <- list()
summary_rows <- list()

for (pair_id in pair_ids) {
  x <- unit[unit$pair_id == pair_id, , drop = FALSE]
  target <- x[x$is_target_gene, c("GSE", "biological_unit_id", "association_z"), drop = FALSE]
  names(target)[3] <- "target_unit_z"
  target <- target[!duplicated(target[c("GSE", "biological_unit_id")]), , drop = FALSE]

  if (nrow(target) < 3L) next

  full_target_z <- sum(target$target_unit_z) / sqrt(nrow(target))
  ranks <- sort(unique(x$match_rank[!x$is_target_gene & x$match_rank > 0L]))
  rows <- list()

  for (rank in ranks) {
    null <- x[!x$is_target_gene & x$match_rank == rank,
              c("GSE", "biological_unit_id", "association_z"), drop = FALSE]
    names(null)[3] <- "null_unit_z"
    null <- null[!duplicated(null[c("GSE", "biological_unit_id")]), , drop = FALSE]
    paired <- merge(target, null, by = c("GSE", "biological_unit_id"), sort = FALSE)
    if (nrow(paired) < 3L) next

    target_z <- sum(paired$target_unit_z) / sqrt(nrow(paired))
    null_z <- sum(paired$null_unit_z) / sqrt(nrow(paired))
    rows[[length(rows) + 1L]] <- data.frame(
      pair_id = pair_id,
      target_gene = x$target_gene[1],
      match_rank = rank,
      n_common_biological_units = nrow(paired),
      n_common_datasets = length(unique(paired$GSE)),
      target_meta_z_same_coverage = target_z,
      null_meta_z = null_z,
      null_at_least_as_extreme = abs(null_z) >= abs(target_z),
      stringsAsFactors = FALSE
    )
  }

  detail <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (!nrow(detail)) next
  null_details[[length(null_details) + 1L]] <- detail

  empirical_p <- (1 + sum(detail$null_at_least_as_extreme)) / (1 + nrow(detail))
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    pair_id = pair_id,
    target_gene = x$target_gene[1],
    n_target_biological_units = nrow(target),
    n_target_datasets = length(unique(target$GSE)),
    target_meta_z_all_units = full_target_z,
    n_balanced_null_replicates = nrow(detail),
    median_common_biological_units = median(detail$n_common_biological_units),
    minimum_common_biological_units = min(detail$n_common_biological_units),
    matched_null_median_meta_z = median(detail$null_meta_z),
    matched_null_95pct_abs_meta_z = as.numeric(quantile(abs(detail$null_meta_z), 0.95)),
    empirical_P_two_sided_balanced = empirical_p,
    stringsAsFactors = FALSE
  )
}

details <- do.call(rbind, null_details)
summary <- do.call(rbind, summary_rows)
summary$matched_gene_FDR_two_sided_balanced <- bh(summary$empirical_P_two_sided_balanced)
summary <- summary[order(summary$matched_gene_FDR_two_sided_balanced,
                         summary$empirical_P_two_sided_balanced), , drop = FALSE]

unit_connection <- gzfile(
  file.path(OUTPUT_ROOT, "01_section_to_biological_unit_associations.tsv.gz"),
  open = "wt"
)
write.table(unit, unit_connection, sep = "\t", quote = FALSE, row.names = FALSE)
close(unit_connection)

null_connection <- gzfile(
  file.path(OUTPUT_ROOT, "02_balanced_rank_matched_null_replicates.tsv.gz"),
  open = "wt"
)
write.table(details, null_connection, sep = "\t", quote = FALSE, row.names = FALSE)
close(null_connection)
write.table(
  summary,
  file.path(OUTPUT_ROOT, "03_balanced_matched_gene_empirical_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

manifest <- c(
  "This balanced analysis supersedes the original pooled 06_matched_gene_empirical_summary.tsv.",
  "Matched genes were selected independently within each spatial dataset.",
  "For null replicate k, match-rank k contributes one matched gene per dataset.",
  "Multiple sections from one biological unit were averaged before meta-analysis.",
  "Each null replicate was compared with the target using exactly the same available biological units.",
  "The empirical two-sided P value is (1 + number of paired null replicates at least as extreme) / (1 + number of usable ranks).",
  paste0("Input: ", INPUT_ROOT),
  paste0("Output: ", OUTPUT_ROOT),
  paste0("Generated: ", format(Sys.time(), tz = "UTC", usetz = TRUE))
)
writeLines(manifest, file.path(OUTPUT_ROOT, "README.txt"))

cat("BALANCED MATCHED-GENE NULL COMPLETE\n")
cat("Pairs:", nrow(summary), "\n")
cat("FDR < 0.05:", sum(summary$matched_gene_FDR_two_sided_balanced < 0.05), "\n")
cat("Output:", OUTPUT_ROOT, "\n")
print(summary, row.names = FALSE)
