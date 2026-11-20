options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: script.R verified_donor_bridge new_output_dir all_atlas|stable_nodes")
ROOT <- normalizePath(args[1], mustWork = TRUE)
OUT <- args[2]
DENOMINATOR <- match.arg(args[3], c("all_atlas", "stable_nodes"))
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
source(file.path(dirname(SCRIPT), "null_statistics.R"))
if (dir.exists(OUT) && length(list.files(OUT, all.files=TRUE, no..=TRUE))) stop("Output exists; do not overwrite")
if (!requireNamespace("jsonlite", quietly=TRUE) || !requireNamespace("digest", quietly=TRUE)) stop("jsonlite/digest required")
bridge_file <- file.path(ROOT, "BRIDGE_COMPLETE.json")
bridge <- jsonlite::fromJSON(bridge_file, simplifyVector=FALSE)
stopifnot(identical(bridge$status, "completed"))
for (field in c("dependency_sha256", "output_sha256")) {
  for (path in names(bridge[[field]])) {
    if (!file.exists(path) || digest::digest(file=path, algo="sha256") != bridge[[field]][[path]]) stop("Bridge SHA mismatch: ", path)
  }
}
gate_file <- file.path(dirname(ROOT), "ANNOTATION_VALIDATION.json")
gate <- jsonlite::fromJSON(gate_file)
stopifnot(isTRUE(gate$publication_ready), !isTRUE(gate$review_only), isTRUE(gate$L1_frozen), isTRUE(gate$L2_frozen))
for (pair in list(c("annotation_file", "annotation_sha256"), c("freeze_file", "freeze_sha256"), c("L1_freeze_file", "L1_freeze_sha256"))) {
  stopifnot(digest::digest(file=gate[[pair[1]]], algo="sha256") == gate[[pair[2]]])
}
N_PERM <- 2000L
SEED <- 20260826L
TOP_N <- 5L
K <- as.integer(bridge$selected_K)
stopifnot(length(K)==1L, is.finite(K), K>=2L)
counts_file <- file.path(ROOT, "CoVarNet_node_counts.tsv")
weights_file <- file.path(ROOT, sprintf("module_node_weights_K%02d.tsv", K))
metadata_file <- file.path(ROOT, "CoVarNet_sample_metadata.tsv")
parent_file <- file.path(ROOT, "CoVarNet_node_parent_map.tsv")
stopifnot(file.exists(counts_file), file.exists(weights_file), file.exists(metadata_file), file.exists(parent_file))
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

counts_df <- read.delim(counts_file, check.names = FALSE)
state_col <- names(counts_df)[1]
states <- as.character(counts_df[[state_col]])
stopifnot(!anyDuplicated(states), !anyDuplicated(names(counts_df)[-1]), length(states)>=TOP_N)
counts <- as.matrix(counts_df[, -1, drop = FALSE])
counts <- validate_count_matrix(counts)
rownames(counts) <- states
if (any(counts < 0L) || any(!is.finite(counts))) stop("Counts must be finite nonnegative integers")

weights_df <- read.delim(weights_file, check.names = FALSE)
stopifnot(!anyDuplicated(weights_df$cell_state), setequal(weights_df$cell_state, states))
rownames(weights_df) <- weights_df$cell_state
weights <- as.matrix(weights_df[states, setdiff(names(weights_df), "cell_state"), drop = FALSE])
storage.mode(weights) <- "double"
stopifnot(ncol(weights)==K, all(is.finite(weights)), all(weights>=0), all(colSums(weights)>0))
TARGET_MODULES <- colnames(weights) # All actual modules, not inherited CM01/CM02.

metadata <- read.delim(metadata_file, check.names = FALSE)
stopifnot(!anyDuplicated(metadata$covarnet_sample_uid), setequal(metadata$covarnet_sample_uid, colnames(counts)))
rownames(metadata) <- metadata$covarnet_sample_uid
metadata <- metadata[colnames(counts), , drop = FALSE]
if (any(is.na(metadata$study_canonical))) stop("Sample metadata did not align to counts")

parent_map <- read.delim(parent_file, check.names=FALSE)
stopifnot(!anyDuplicated(parent_map$cell_state), setequal(parent_map$cell_state, states))
parent <- parent_map$parent_lineage[match(states, parent_map$cell_state)]
stopifnot(!anyNA(parent), all(nzchar(parent)))
names(parent) <- states
lineages <- unique(parent)

sample_totals <- if (DENOMINATOR=="all_atlas") metadata$n_total_atlas_cells else colSums(counts)
stopifnot(length(sample_totals)==ncol(counts), all(is.finite(sample_totals)), all(sample_totals>0), all(colSums(counts)<=sample_totals))
proportions <- sweep(counts, 2, sample_totals, "/")
if (DENOMINATOR=="all_atlas") {
  p <- read.delim(file.path(ROOT,"CoVarNet_node_proportion_of_all_atlas_cells.tsv"), check.names=FALSE, row.names=1)
  stopifnot(max(abs(as.matrix(p[states,colnames(counts),drop=FALSE])-proportions))<1e-12)
}

rank_residual <- function(values, parent_total, study) {
  y <- rank(values, ties.method = "average", na.last = "keep")
  x <- rank(parent_total, ties.method = "average", na.last = "keep")
  fit <- lm(y ~ x + factor(study), na.action = na.exclude)
  as.numeric(residuals(fit))
}

parent_totals <- matrix(0, nrow = length(lineages), ncol = ncol(counts),
                        dimnames = list(lineages, colnames(counts)))
for (lin in lineages) {
  parent_totals[lin, ] <- colSums(counts[parent == lin, , drop = FALSE])
}

residual_qr <- list()
for (lin in lineages) {
  design <- model.matrix(
    ~ rank_parent_total + study,
    data = data.frame(
      rank_parent_total = rank(parent_totals[lin, ] / sample_totals,
                               ties.method = "average"),
      study = factor(metadata$study_canonical)
    )
  )
  residual_qr[[lin]] <- qr(design)
}

partial_residuals <- matrix(NA_real_, nrow = nrow(counts), ncol = ncol(counts),
                            dimnames = dimnames(counts))
within_parent_logit <- matrix(NA_real_, nrow = nrow(counts), ncol = ncol(counts),
                              dimnames = dimnames(counts))
for (i in seq_len(nrow(counts))) {
  lin <- parent[i]
  total <- parent_totals[lin, ]
  partial_residuals[i, ] <- rank_residual(proportions[i, ], total / sample_totals,
                                          metadata$study_canonical)
  n_lin <- sum(parent == lin)
  share <- (counts[i, ] + 0.5) / (total + 0.5 * n_lin)
  share <- pmin(pmax(share, 1e-6), 1 - 1e-6)
  within_parent_logit[i, ] <- qlogis(share)
}

pair_metrics <- function(matrix_values, indices) {
  rows <- list()
  counter <- 1L
  for (a in seq_len(length(indices) - 1L)) {
    for (b in (a + 1L):length(indices)) {
      i <- indices[a]
      j <- indices[b]
      rho <- suppressWarnings(cor(matrix_values[i, ], matrix_values[j, ],
                                  method = "spearman", use = "pairwise.complete.obs"))
      rows[[counter]] <- data.frame(
        state_1 = states[i],
        parent_1 = parent[i],
        state_2 = states[j],
        parent_2 = parent[j],
        cross_parent_pair = parent[i] != parent[j],
        rho = rho
      )
      counter <- counter + 1L
    }
  }
  do.call(rbind, rows)
}

summarise_pairs <- function(x) {
  cross <- x$rho[x$cross_parent_pair & is.finite(x$rho)]
  all <- x$rho[is.finite(x$rho)]
  c(
    median_all_pair_rho = if (length(all)) median(all) else NA_real_,
    positive_all_pair_fraction = if (length(all)) mean(all > 0) else NA_real_,
    median_cross_parent_pair_rho = if (length(cross)) median(cross) else NA_real_,
    positive_cross_parent_pair_fraction = if (length(cross)) mean(cross > 0) else NA_real_,
    n_all_pairs = length(all),
    n_cross_parent_pairs = length(cross)
  )
}

top_indices <- list()
observed_rows <- list()
pair_detail <- list()
for (module in colnames(weights)) {
  idx <- order(weights[, module], decreasing = TRUE)[seq_len(TOP_N)]
  top_indices[[module]] <- idx
  raw_pairs <- pair_metrics(proportions, idx)
  residual_pairs <- pair_metrics(partial_residuals, idx)
  within_pairs <- pair_metrics(within_parent_logit, idx)
  raw_sum <- summarise_pairs(raw_pairs)
  residual_sum <- summarise_pairs(residual_pairs)
  within_sum <- summarise_pairs(within_pairs)
  observed_rows[[module]] <- data.frame(
    module = module,
    top_states = paste(states[idx], collapse = ";"),
    top_parent_lineages = paste(parent[idx], collapse = ";"),
    observed_raw_median_all_pair_rho = raw_sum["median_all_pair_rho"],
    observed_raw_median_cross_parent_pair_rho = raw_sum["median_cross_parent_pair_rho"],
    observed_parent_adjusted_median_all_pair_rho = residual_sum["median_all_pair_rho"],
    observed_parent_adjusted_median_cross_parent_pair_rho = residual_sum["median_cross_parent_pair_rho"],
    observed_within_parent_share_median_all_pair_rho = within_sum["median_all_pair_rho"],
    observed_within_parent_share_median_cross_parent_pair_rho = within_sum["median_cross_parent_pair_rho"]
  )
  if (module %in% TARGET_MODULES) {
    raw_pairs$metric <- paste0("raw_", DENOMINATOR, "_proportion")
    residual_pairs$metric <- "rank_residual_adjusted_for_parent_total_and_study"
    within_pairs$metric <- "within_parent_logit_share"
    detail <- rbind(raw_pairs, residual_pairs, within_pairs)
    detail$module <- module
    pair_detail[[module]] <- detail
  }
}
observed <- do.call(rbind, observed_rows)
write.table(observed, file.path(OUT, "01_observed_parent_adjusted_module_coherence.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(do.call(rbind, pair_detail),
            file.path(OUT, "02_all_modules_parent_adjusted_pair_details.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

set.seed(SEED)
null_rows <- vector("list", N_PERM * length(colnames(weights)))
position <- 1L
target_node_indices <- sort(unique(unlist(top_indices[TARGET_MODULES])))
for (permutation in seq_len(N_PERM)) {
  null_counts <- checked_parent_null(counts, parent)
  null_props <- sweep(null_counts, 2, sample_totals, "/")
  null_parent_adjusted <- matrix(NA_real_, nrow = nrow(counts), ncol = ncol(counts),
                                 dimnames = dimnames(counts))
  for (i in target_node_indices) {
    ranked <- rank(null_props[i, ], ties.method = "average")
    null_parent_adjusted[i, ] <- qr.resid(residual_qr[[parent[i]]], ranked)
  }
  for (module in colnames(weights)) {
    idx <- top_indices[[module]]
    pairs <- pair_metrics(null_props, idx)
    metrics <- summarise_pairs(pairs)
    adjusted_metrics <- if (module %in% TARGET_MODULES) {
      summarise_pairs(pair_metrics(null_parent_adjusted, idx))
    } else {
      c(median_all_pair_rho = NA_real_, median_cross_parent_pair_rho = NA_real_)
    }
    null_rows[[position]] <- data.frame(
      permutation = permutation,
      module = module,
      null_median_all_pair_rho = metrics["median_all_pair_rho"],
      null_positive_all_pair_fraction = metrics["positive_all_pair_fraction"],
      null_median_cross_parent_pair_rho = metrics["median_cross_parent_pair_rho"],
      null_positive_cross_parent_pair_fraction = metrics["positive_cross_parent_pair_fraction"],
      null_parent_adjusted_median_all_pair_rho = adjusted_metrics["median_all_pair_rho"],
      null_parent_adjusted_median_cross_parent_pair_rho = adjusted_metrics["median_cross_parent_pair_rho"]
    )
    position <- position + 1L
  }
  if (permutation == 1L || permutation %% 100L == 0L || permutation == N_PERM) {
    message("Parent-lineage-preserving null: ", permutation, "/", N_PERM)
  }
}
null <- do.call(rbind, null_rows)
gz <- gzfile(file.path(OUT, "03_parent_lineage_preserving_null_distribution.tsv.gz"), "wt")
write.table(null, gz, sep = "\t", row.names = FALSE, quote = FALSE)
close(gz)

summary_rows <- list()
for (module in colnames(weights)) {
  obs <- observed[observed$module == module, ]
  nul <- null[null$module == module, ]
  p_all <- empirical_upper(obs$observed_raw_median_all_pair_rho, nul$null_median_all_pair_rho)
  p_cross <- empirical_upper(obs$observed_raw_median_cross_parent_pair_rho, nul$null_median_cross_parent_pair_rho)
  p_adjusted_all <- empirical_upper(obs$observed_parent_adjusted_median_all_pair_rho, nul$null_parent_adjusted_median_all_pair_rho)
  p_adjusted_cross <- empirical_upper(obs$observed_parent_adjusted_median_cross_parent_pair_rho, nul$null_parent_adjusted_median_cross_parent_pair_rho)
  summary_rows[[module]] <- data.frame(
    module = module,
    n_permutations = N_PERM,
    n_valid_null_all_pairs = sum(is.finite(nul$null_median_all_pair_rho)),
    n_valid_null_cross_parent_pairs = sum(is.finite(nul$null_median_cross_parent_pair_rho)),
    n_valid_null_parent_adjusted_all_pairs = sum(is.finite(nul$null_parent_adjusted_median_all_pair_rho)),
    n_valid_null_parent_adjusted_cross_pairs = sum(is.finite(nul$null_parent_adjusted_median_cross_parent_pair_rho)),
    observed_raw_median_all_pair_rho = obs$observed_raw_median_all_pair_rho,
    null_median_all_pair_rho = median(nul$null_median_all_pair_rho, na.rm = TRUE),
    null_95pct_all_pair_rho = quantile(nul$null_median_all_pair_rho, 0.95, na.rm = TRUE),
    empirical_P_all_pairs = p_all,
    observed_raw_median_cross_parent_pair_rho = obs$observed_raw_median_cross_parent_pair_rho,
    null_median_cross_parent_pair_rho = median(nul$null_median_cross_parent_pair_rho, na.rm = TRUE),
    null_95pct_cross_parent_pair_rho = quantile(nul$null_median_cross_parent_pair_rho, 0.95, na.rm = TRUE),
    empirical_P_cross_parent_pairs = p_cross,
    observed_parent_adjusted_median_all_pair_rho = obs$observed_parent_adjusted_median_all_pair_rho,
    observed_parent_adjusted_median_cross_parent_pair_rho = obs$observed_parent_adjusted_median_cross_parent_pair_rho,
    null_median_parent_adjusted_all_pair_rho = median(nul$null_parent_adjusted_median_all_pair_rho, na.rm = TRUE),
    null_95pct_parent_adjusted_all_pair_rho = quantile(nul$null_parent_adjusted_median_all_pair_rho, 0.95, na.rm = TRUE),
    empirical_P_parent_adjusted_all_pairs = p_adjusted_all,
    null_median_parent_adjusted_cross_parent_pair_rho = median(nul$null_parent_adjusted_median_cross_parent_pair_rho, na.rm = TRUE),
    null_95pct_parent_adjusted_cross_parent_pair_rho = quantile(nul$null_parent_adjusted_median_cross_parent_pair_rho, 0.95, na.rm = TRUE),
    empirical_P_parent_adjusted_cross_parent_pairs = p_adjusted_cross,
    observed_within_parent_share_median_all_pair_rho = obs$observed_within_parent_share_median_all_pair_rho,
    observed_within_parent_share_median_cross_parent_pair_rho = obs$observed_within_parent_share_median_cross_parent_pair_rho
  )
}
summary <- do.call(rbind, summary_rows)
summary$FDR_all_pairs <- p.adjust(summary$empirical_P_all_pairs, method = "BH")
summary$FDR_cross_parent_pairs <- p.adjust(summary$empirical_P_cross_parent_pairs, method = "BH")
summary$FDR_parent_adjusted_all_pairs <- p.adjust(summary$empirical_P_parent_adjusted_all_pairs, method = "BH")
summary$FDR_parent_adjusted_cross_parent_pairs <- p.adjust(summary$empirical_P_parent_adjusted_cross_parent_pairs, method = "BH")
write.table(summary, file.path(OUT, "04_parent_lineage_preserving_null_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(summary[summary$module %in% TARGET_MODULES, ],
            file.path(OUT, "05_all_modules_parent_lineage_null_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

manifest <- list(
  status = "completed",
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  n_permutations = N_PERM,
  seed = SEED,
  n_samples = ncol(counts),
  n_nodes = nrow(counts),
  n_parent_lineages = length(lineages),
  selected_K = K, target_modules = TARGET_MODULES, denominator = DENOMINATOR,
  annotation_validation_sha256 = digest::digest(file=gate_file, algo="sha256"),
  bridge_sha256 = digest::digest(file=bridge_file, algo="sha256"),
  code_sha256 = digest::digest(file=SCRIPT, algo="sha256"),
  helper_sha256 = digest::digest(file=file.path(dirname(SCRIPT),"null_statistics.R"), algo="sha256"),
  nonestimable_observed_statistic = "NA; never convert to small empirical P",
  finite_null_definition = "Report effective null count; p uses finite null distribution only",
  primary_vs_legacy = "all_atlas aligns to main nsNMF denominator; stable_nodes is explicit historical closure sensitivity",
  exact_null_constraints = c(
    "sample-specific parent-lineage totals",
    "global Level 2 node totals within each parent lineage"
  ),
  primary_null_test = "median Spearman correlation among top-five module nodes",
  additional_parent_adjustment = "rank residuals adjusted for parent-lineage total and study"
)
jsonlite::write_json(manifest, file.path(OUT,"PARENT_NULL_COMPLETE.json"), auto_unbox=TRUE, pretty=TRUE, na="null")
json_text <- paste(capture.output(dput(manifest)), collapse = "\n")
writeLines(json_text, file.path(OUT, "run_manifest_dput.txt"))

print(summary[summary$module %in% TARGET_MODULES, ])
cat("Outputs:", OUT, "\n")
