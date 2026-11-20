rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(NMF)
  library(ggplot2)
})


# =============================================================================
# 0. Configuration
# =============================================================================

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop("Usage: stable_reanalysis.R prepared_profile_root")
ROOT <- normalizePath(args[[1L]])
script_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir,"covarnet_review_helpers.R"))
OUT <- file.path(ROOT, "stable_NMF")
if (dir.exists(OUT)) stop("Fresh output required; refuse exists-only historical cache")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

RANK_MIN <- as.integer(Sys.getenv("COVARNET_RANK_MIN", "2"))
RANK_MAX_REQUESTED <- as.integer(Sys.getenv("COVARNET_RANK_MAX", "25"))
RANK_NRUN <- as.integer(Sys.getenv("COVARNET_RANK_NRUN", "50"))
FINAL_NRUN <- as.integer(Sys.getenv("COVARNET_FINAL_NRUN", "200"))
STABILITY_NRUN <- as.integer(Sys.getenv("COVARNET_STABILITY_NRUN", "30"))
N_BOOTSTRAP <- as.integer(Sys.getenv("COVARNET_N_BOOTSTRAP", "30"))
TOP_NODES <- 10L
SEED <- as.integer(Sys.getenv("COVARNET_SEED", "20260806"))

# Optional manual override. Leave as NA to use the pre-specified automatic
# selection below. The script records every decision metric and does not stop
# after the rank survey. The historical K=14 is refitted as sensitivity only.
K_FINAL_OVERRIDE <- NA_integer_
K14_SENSITIVITY_NRUN <- 100L
COPHENETIC_TOLERANCE <- 0.05
RSS_ELBOW_WINDOW <- 2L

workers <- as.integer(Sys.getenv("COVARNET_WORKERS", "2"))
if(!is.finite(workers) || workers<1L) stop("COVARNET_WORKERS must be positive")
nmf_options <- paste0("p", workers)
cat("NMF workers:", workers, "\n")


# =============================================================================
# 1. Read sample-level stable-node frequencies
# =============================================================================

read_matrix <- function(path) {
  x <- fread(path, check.names = FALSE)
  rn <- x[[1L]]
  x[[1L]] <- NULL
  m <- as.matrix(x)
  storage.mode(m) <- "double"
  rownames(m) <- rn
  m
}

prop <- read_matrix(file.path(ROOT, "CoVarNet_node_proportion_of_all_atlas_cells.tsv"))
prop_closed <- read_matrix(file.path(ROOT, "CoVarNet_node_proportion_among_stable_nodes_sensitivity.tsv"))
counts <- read_matrix(file.path(ROOT, "CoVarNet_node_counts.tsv"))
meta <- fread(file.path(ROOT, "CoVarNet_sample_metadata.tsv"))
stopifnot(all(c("covarnet_sample_uid", "covarnet_donor_key") %in% names(meta)))

stopifnot(!anyDuplicated(rownames(prop)), !anyDuplicated(colnames(prop)), !anyDuplicated(meta$covarnet_sample_uid))
stopifnot(identical(dimnames(prop),dimnames(prop_closed)),identical(dimnames(prop),dimnames(counts)))
stopifnot(setequal(colnames(prop),meta$covarnet_sample_uid))
stopifnot(all(is.finite(prop)),all(prop>=0 & prop<=1),all(is.finite(counts)),all(counts>=0),all(counts==round(counts)))
common <- colnames(prop)
if (length(common) < 20L) stop("Too few donor-by-tissue-by-condition samples for stable CoVarNet")
prop <- prop[, common, drop = FALSE]
prop_closed <- prop_closed[rownames(prop), common, drop = FALSE]
counts <- counts[rownames(prop), common, drop = FALSE]
meta <- meta[match(common, covarnet_sample_uid)]
stopifnot(all(meta$covarnet_sample_uid == colnames(prop)))

# Row-wise min-max scaling reproduces the non-negative input used by the original
# CoVarNet/NMF workflow. Raw proportions and logit proportions are retained for
# descriptive transform checks, not independent validation.
row_min <- apply(prop, 1L, min, na.rm = TRUE)
row_max <- apply(prop, 1L, max, na.rm = TRUE)
variable <- is.finite(row_min) & is.finite(row_max) & row_max > row_min
fwrite(data.table(cell_state=rownames(prop),variable_for_NMF=variable,row_min=row_min,row_max=row_max),file.path(OUT,"variable_node_audit.tsv"),sep="\t")
prop <- prop[variable, , drop = FALSE]
prop_closed <- prop_closed[rownames(prop), , drop = FALSE]
counts <- counts[rownames(prop), , drop = FALSE]
mat_nmf <- (prop - row_min[variable]) / (row_max[variable] - row_min[variable])
if(any(!is.finite(mat_nmf))) stop("Nonfinite NMF input")

cat("Stable nodes:", nrow(mat_nmf), " CoVarNet sample units:", ncol(mat_nmf), "\n")
write.table(mat_nmf, file.path(OUT, "NMF_input_minmax.tsv"), sep = "\t", quote = FALSE)


# =============================================================================
# 2. Rank survey
# =============================================================================

rank_max <- min(RANK_MAX_REQUESTED, nrow(mat_nmf) - 1L, ncol(mat_nmf) - 1L)
if (rank_max < RANK_MIN+2L) stop("Insufficient dynamic dimensions for at least three surveyed ranks")
ranks <- seq.int(RANK_MIN, rank_max)
if (length(ranks) < 3L) stop("Insufficient matrix dimensions for rank survey")
rank_file <- file.path(OUT, "rank_survey.rds")

if (!file.exists(rank_file)) {
  set.seed(SEED)
  rank_fit <- nmf(
    mat_nmf,
    rank = ranks,
    method = "nsNMF",
    nrun = RANK_NRUN,
    seed = SEED,
    .options = nmf_options
  )
  saveRDS(rank_fit, rank_file)
} else {
  rank_fit <- readRDS(rank_file)
}

pdf(file.path(OUT, "rank_survey.pdf"), width = 9, height = 7)
plot(rank_fit)
dev.off()

rank_summary <- tryCatch(as.data.frame(summary(rank_fit)), error = function(e) NULL)
if (is.null(rank_summary)) stop("Could not extract NMF rank-survey metrics")

rank_metrics <- as.data.table(rank_summary, keep.rownames = "rank_label")
rank_metrics[, rank := suppressWarnings(as.integer(rank_label))]
if (anyNA(rank_metrics$rank)) {
  rank_metrics[, rank := ranks[seq_len(.N)]]
}
setorder(rank_metrics, rank)

canonical_name <- function(x) {
  gsub("_+", "_", gsub("[^a-z0-9]+", "_", tolower(x)))
}

find_metric <- function(patterns) {
  original <- names(rank_metrics)
  canonical <- canonical_name(original)
  hit <- integer()
  for (pattern in patterns) {
    hit <- which(grepl(pattern, canonical))
    if (length(hit)) break
  }
  hit <- hit[!original[hit] %in% c("rank", "rank_label")]
  if (length(hit)) original[hit[1L]] else NA_character_
}

rescale01 <- function(x, higher_is_better = TRUE) {
  x <- as.numeric(x)
  good <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (!any(good)) return(out)
  lo <- min(x[good])
  hi <- max(x[good])
  if (hi <= lo) {
    out[good] <- 1
  } else {
    out[good] <- (x[good] - lo) / (hi - lo)
  }
  if (!higher_is_better) out[good] <- 1 - out[good]
  out
}

cophenetic_col <- find_metric(c("cophenetic"))
dispersion_col <- find_metric(c("dispersion"))
silhouette_col <- find_metric(c("silhouette_consensus", "silhouette"))
rss_col <- find_metric(c("^rss$", "residual", "residuals"))

metric_map <- data.table(
  role = c("cophenetic", "dispersion", "silhouette_consensus", "rss"),
  detected_column = c(cophenetic_col, dispersion_col, silhouette_col, rss_col)
)
fwrite(metric_map, file.path(OUT, "rank_metric_column_map.tsv"), sep = "\t")

quality_columns <- character()
for (metric_col in c(cophenetic_col, dispersion_col, silhouette_col)) {
  if (!is.na(metric_col)) {
    new_col <- paste0("scaled_", canonical_name(metric_col))
    rank_metrics[, (new_col) := rescale01(get(metric_col), TRUE)]
    quality_columns <- c(quality_columns, new_col)
  }
}

if (length(quality_columns)) {
  quality_matrix <- as.matrix(rank_metrics[, ..quality_columns])
  rank_metrics[, composite_stability_quality := rowMeans(quality_matrix, na.rm = TRUE)]
  rank_metrics[!is.finite(composite_stability_quality), composite_stability_quality := NA_real_]
} else {
  rank_metrics[, composite_stability_quality := NA_real_]
}

# RSS elbow: maximum perpendicular distance from the line joining the first and
# last points of the normalized rank-versus-log(RSS) curve.
RSS_ELBOW_K <- NA_integer_
if (!is.na(rss_col)) {
  rss_values <- as.numeric(rank_metrics[[rss_col]])
  good <- is.finite(rss_values) & rss_values > 0
  if (sum(good) >= 3L) {
    x <- rescale01(rank_metrics$rank[good], TRUE)
    y <- rescale01(log(rss_values[good]), TRUE)
    x1 <- x[1L]; y1 <- y[1L]
    x2 <- x[length(x)]; y2 <- y[length(y)]
    denominator_line <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
    distance <- abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) /
      max(denominator_line, .Machine$double.eps)
    RSS_ELBOW_K <- as.integer(rank_metrics$rank[good][which.max(distance)])
  }
}

if (length(RSS_ELBOW_K) != 1L) stop("RSS elbow selection did not return one rank")
rank_metrics[, rss_elbow_rank := RSS_ELBOW_K]
rank_metrics[, within_rss_elbow_window := if (is.na(RSS_ELBOW_K)) {
  rep(TRUE, .N)
} else {
  abs(rank - RSS_ELBOW_K) <= RSS_ELBOW_WINDOW
}]

if (!is.na(cophenetic_col)) {
  cophenetic_values <- as.numeric(rank_metrics[[cophenetic_col]])
  if (any(is.finite(cophenetic_values))) {
    best_cophenetic <- max(cophenetic_values, na.rm = TRUE)
    rank_metrics[, stable_by_cophenetic :=
                   as.numeric(get(cophenetic_col)) >= best_cophenetic - COPHENETIC_TOLERANCE]
    rank_metrics[is.na(stable_by_cophenetic), stable_by_cophenetic := FALSE]
  } else {
    rank_metrics[, stable_by_cophenetic := TRUE]
  }
} else {
  rank_metrics[, stable_by_cophenetic := TRUE]
}

rank_metrics[, automatic_candidate :=
               stable_by_cophenetic & within_rss_elbow_window &
               is.finite(composite_stability_quality)]

choose_best_rank <- function(dt) {
  candidates <- dt[automatic_candidate == TRUE]
  if (!nrow(candidates)) {
    candidates <- dt[stable_by_cophenetic == TRUE & is.finite(composite_stability_quality)]
  }
  if (nrow(candidates)) {
    setorder(candidates, -composite_stability_quality, rank)
    return(candidates$rank[1L])
  }
  if (!is.na(RSS_ELBOW_K)) return(RSS_ELBOW_K)
  as.integer(round(stats::median(dt$rank)))
}

K_AUTOMATIC <- choose_best_rank(copy(rank_metrics))
K_FINAL <- if (length(K_FINAL_OVERRIDE) == 1L && !is.na(K_FINAL_OVERRIDE)) {
  as.integer(K_FINAL_OVERRIDE)
} else {
  as.integer(K_AUTOMATIC)
}

if (K_FINAL < RANK_MIN || K_FINAL > rank_max) {
  stop("Selected K_FINAL is outside the surveyed range")
}

rank_metrics[, selected_automatically := rank == K_AUTOMATIC]
rank_metrics[, selected_for_primary_analysis := rank == K_FINAL]
rank_metrics[, selection_mode := if (is.na(K_FINAL_OVERRIDE)) "automatic_predefined_rule" else "manual_override"]
fwrite(rank_metrics, file.path(OUT, "rank_metrics_with_decision.tsv"), sep = "\t")

writeLines(c(
  paste0("K_FINAL=", K_FINAL),
  paste0("K_AUTOMATIC=", K_AUTOMATIC),
  paste0("K_FINAL_OVERRIDE=", ifelse(is.na(K_FINAL_OVERRIDE), "NA", K_FINAL_OVERRIDE)),
  paste0("RSS_ELBOW_RANK=", ifelse(is.na(RSS_ELBOW_K), "NA", RSS_ELBOW_K)),
  paste0("SELECTION_MODE=", ifelse(is.na(K_FINAL_OVERRIDE), "automatic_predefined_rule", "manual_override")),
  "RULE=highest composite cophenetic/dispersion/consensus-silhouette quality among ranks within +/-2 of the RSS elbow and within 0.05 of maximum cophenetic stability; documented fallbacks are used if a metric is unavailable",
  "K14=not automatically run in this final-atlas adapter"
), file.path(OUT, "selected_K.txt"))

plot_columns <- quality_columns
if (!is.na(rss_col)) {
  rank_metrics[, scaled_inverse_rss := rescale01(get(rss_col), FALSE)]
  plot_columns <- c(plot_columns, "scaled_inverse_rss")
}
if (length(plot_columns)) {
  rank_plot_data <- melt(
    rank_metrics,
    id.vars = "rank",
    measure.vars = plot_columns,
    variable.name = "metric",
    value.name = "scaled_value"
  )
  p_rank <- ggplot(rank_plot_data, aes(rank, scaled_value, colour = metric)) +
    geom_line(linewidth = 0.65) +
    geom_point(size = 1.5) +
    geom_vline(xintercept = K_FINAL, linetype = 2, linewidth = 0.7) +
    annotate("text", x = K_FINAL, y = 0.03, label = paste0("Primary K=", K_FINAL),
             angle = 90, vjust = -0.4, hjust = 0, size = 3) +
    scale_x_continuous(breaks = ranks) +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "NMF rank", y = "Scaled metric", colour = NULL,
         title = "CoVarNet rank survey and pre-specified rank selection")
  ggsave(file.path(OUT, "rank_selection_diagnostic.pdf"), p_rank, width = 9, height = 5.5)
  ggsave(file.path(OUT, "rank_selection_diagnostic.png"), p_rank, width = 9, height = 5.5,
         dpi = 300)
}

cat("Rank survey completed without interruption.\n")
cat("Automatic K:", K_AUTOMATIC, " Primary K:", K_FINAL, "\n")


# =============================================================================
# 3. Final NMF and module composition
# =============================================================================

fit_file <- file.path(OUT, sprintf("final_nmf_K%02d.rds", K_FINAL))
if (!file.exists(fit_file)) {
  set.seed(SEED + K_FINAL)
  fit <- nmf(
    mat_nmf,
    rank = K_FINAL,
    method = "nsNMF",
    nrun = FINAL_NRUN,
    seed = SEED + K_FINAL,
    .options = nmf_options
  )
  saveRDS(fit, fit_file)
} else {
  fit <- readRDS(fit_file)
}

W <- basis(fit)
H <- coef(fit)
colnames(W) <- sprintf("CM%02d", seq_len(ncol(W)))
rownames(H) <- colnames(W)

# Normalize node weights within module and module coefficients within profile.
# H_norm is relative module coefficient/activity; it is not a cell-abundance estimate.
W_norm <- sweep(W, 2L, colSums(W), "/")
H_norm <- sweep(H, 2L, colSums(H), "/")
W_norm[!is.finite(W_norm)] <- 0
H_norm[!is.finite(H_norm)] <- 0

fwrite(
  as.data.table(W_norm, keep.rownames = "cell_state"),
  file.path(OUT, sprintf("module_node_weights_K%02d.tsv", K_FINAL)), sep = "\t"
)
fwrite(
  as.data.table(t(H_norm), keep.rownames = "covarnet_sample_uid"),
  file.path(OUT, sprintf("sample_module_relative_coefficient_K%02d.tsv", K_FINAL)), sep = "\t"
)

top_nodes <- rbindlist(lapply(seq_len(ncol(W_norm)), function(k) {
  ord <- order(W_norm[, k], decreasing = TRUE)
  n <- min(TOP_NODES, length(ord))
  data.table(
    module = colnames(W_norm)[k], rank = seq_len(n),
    cell_state = rownames(W_norm)[ord[seq_len(n)]],
    weight = W_norm[ord[seq_len(n)], k]
  )
}))
fwrite(top_nodes, file.path(OUT, sprintf("top_nodes_K%02d.tsv", K_FINAL)), sep = "\t")

p_top <- ggplot(top_nodes, aes(x = reorder(cell_state, weight), y = weight, fill = module)) +
  geom_col(width = 0.75, show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~ module, scales = "free_y", ncol = 4) +
  theme_classic(base_size = 9) +
  labs(x = NULL, y = "NMF node weight", title = "Stable cellular modules")
ggsave(file.path(OUT, sprintf("top_nodes_K%02d.pdf", K_FINAL)), p_top,
       width = 12, height = max(7, ceiling(K_FINAL / 4) * 2.5))

if (requireNamespace("pheatmap", quietly = TRUE)) {
  pdf(file.path(OUT, sprintf("module_node_weight_heatmap_K%02d.pdf", K_FINAL)),
      width = 8, height = max(8, nrow(W_norm) * 0.16))
  pheatmap::pheatmap(W_norm, scale = "row", cluster_cols = TRUE,
                     fontsize_row = 6, fontsize_col = 8,
                     main = "NMF node weights (row-scaled)")
  dev.off()
}

# =============================================================================
# 4. Sample phenotype summaries; these explain modules but do not define them
# =============================================================================

sample_modules <- cbind(meta, as.data.table(t(H_norm)))
sample_modules[, module_value_definition := "relative_module_coefficient_H_divided_by_profile_column_sum"]
sample_modules[, dominant_module := colnames(W_norm)[max.col(as.matrix(.SD), ties.method = "first")],
               .SDcols = colnames(W_norm)]
fwrite(sample_modules, file.path(OUT, sprintf("sample_module_metadata_K%02d.tsv", K_FINAL)), sep = "\t")

# Phenotype inference is deliberately delegated to the explicit contrast/model
# module. Original unpaired Kruskal tests pooled repeated donor observations.
writeLines("Descriptive module metadata only. Inference requires approved donor-aware contrasts.",
           file.path(OUT, "PHENOTYPE_INFERENCE_NOT_RUN.txt"))

# =============================================================================
# 5. Pairwise node-network sensitivity on raw and logit proportions
# =============================================================================

# Only descriptive rho is emitted. Repeated donor/tissue profiles are not
# independent observations. A donor bootstrap for NMF stability cannot correct
# ordinary cor.test P values. Any inferential network test needs separate approval.
pairwise_network <- descriptive_profile_spearman

eps <- 0.5 / median(meta$n_total_atlas_cells)
logit_prop <- log((prop + eps) / (1 - prop + eps))
network <- rbind(
  pairwise_network(prop, "raw_proportion", meta$covarnet_donor_key),
  pairwise_network(logit_prop, "logit_proportion", meta$covarnet_donor_key)
)
fwrite(network, file.path(OUT, "stable_node_pairwise_spearman_network.tsv.gz"), sep = "\t")
writeLines(c(
  "Spearman rho is descriptive; P and FDR are intentionally not calculated.",
  "Raw and strictly monotone logit proportions share ranks; these are not independent validations.",
  "Donor-block NMF bootstrap assesses basis stability, not the null distribution of profile rho.",
  "H_norm is a relative module coefficient/activity, not cell abundance."
), file.path(OUT,"COVARNET_INTERPRETATION.txt"))


# =============================================================================
# 6. Leave-one-study and sample-bootstrap module stability
# =============================================================================

# Preserves inherited many-to-one matching as the primary stability description.
# Adds exact one-to-one maximum-total-cosine matching as an explicitly named sensitivity.
compare_basis <- function(reference,candidate,run_type,omitted=NA_character_) {
  compare_basis_review(reference,candidate,run_type,omitted,TOP_NODES)
}

stability_rows <- list()

# Denominator sensitivity: close the composition only within retained stable
# nodes, refit at the selected K, and compare module bases with the primary fit.
closed_min <- apply(prop_closed, 1L, min)
closed_max <- apply(prop_closed, 1L, max)
closed_nmf <- (prop_closed - closed_min) / pmax(closed_max - closed_min, .Machine$double.eps)
closed_nmf[!is.finite(closed_nmf)] <- 0
set.seed(SEED + 500L)
closed_fit <- tryCatch(
  nmf(closed_nmf, K_FINAL, method = "nsNMF", nrun = STABILITY_NRUN,
      seed = SEED + 500L, .options = nmf_options),
  error = function(e) stop("Required denominator sensitivity failed: ", conditionMessage(e))
)
if (!is.null(closed_fit)) {
  stability_rows[["DENOMINATOR_SENSITIVITY"]] <- compare_basis(
    W_norm, basis(closed_fit), "stable_node_denominator_sensitivity", "stable_nodes_only"
  )
}

for (study in unique(meta$study_canonical)) {
  keep <- meta$study_canonical != study
  if (sum(keep) < max(20L, 3L * K_FINAL)) next
  set.seed(SEED + match(study, unique(meta$study_canonical)))
  los_fit <- tryCatch(
    nmf(mat_nmf[, keep, drop = FALSE], K_FINAL, method = "nsNMF",
        nrun = STABILITY_NRUN, seed = SEED, .options = nmf_options),
    error = function(e) stop("Required leave-study fit failed for ", study, ": ", conditionMessage(e))
  )
  if (!is.null(los_fit)) {
    stability_rows[[paste0("LOS_", study)]] <- compare_basis(
      W_norm, basis(los_fit), "leave_one_study", as.character(study)
    )
  }
}

set.seed(SEED + 1000L)
for (b in seq_len(N_BOOTSTRAP)) {
  donor_keys <- unique(meta$covarnet_donor_key)
  draws <- sample(donor_keys, length(donor_keys), replace=TRUE)
  idx <- unlist(lapply(draws, function(d) which(meta$covarnet_donor_key==d)))
  boot <- mat_nmf[, idx, drop = FALSE]
  colnames(boot) <- make.unique(colnames(boot))
  boot_fit <- tryCatch(
    nmf(boot, K_FINAL, method = "nsNMF", nrun = STABILITY_NRUN,
        seed = SEED + b, .options = nmf_options),
    error = function(e) stop("Required bootstrap fit failed at draw ", b, ": ", conditionMessage(e))
  )
  if (!is.null(boot_fit)) {
    stability_rows[[paste0("BOOT_", b)]] <- compare_basis(
      W_norm, basis(boot_fit), "donor_block_bootstrap", as.character(b)
    )
  }
}

stability <- rbindlist(stability_rows, fill = TRUE)
fwrite(stability, file.path(OUT, sprintf("module_stability_K%02d.tsv", K_FINAL)), sep = "\t")
stability_summary <- stability[, .(
  n_runs = .N,
  median_cosine = median(cosine_similarity, na.rm = TRUE),
  q05_cosine = quantile(cosine_similarity, 0.05, na.rm = TRUE),
  median_top_node_jaccard = median(top_node_jaccard, na.rm = TRUE)
), by = .(run_type, matching_mode, reference_module)]
fwrite(stability_summary, file.path(OUT, sprintf("module_stability_summary_K%02d.tsv", K_FINAL)), sep = "\t")
collapse_audit <- stability[, .(
  n_reference_modules=.N,
  n_distinct_matched_modules=uniqueN(matched_module),
  n_unused_candidate_modules=K_FINAL-uniqueN(matched_module),
  n_reference_rows_in_duplicate_matches=sum(duplicate_match),
  module_collapse_detected=any(duplicate_match)
), by=.(run_type,omitted_or_bootstrap,matching_mode)]
fwrite(collapse_audit,file.path(OUT,sprintf("module_matching_collapse_audit_K%02d.tsv",K_FINAL)),sep="\t")

cat("\nStable CoVarNet/NMF analysis completed.\n")
cat("Output:", OUT, "\n")
writeLines(capture.output(sessionInfo()),file.path(OUT,"sessionInfo.txt"))
