options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript g08_covarnet_donor_block_sensitivity.R frozen_input_dir output_dir")
input_dir <- normalizePath(args[[1L]])
out_dir <- args[[2L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(NMF)
})

MODE <- Sys.getenv("G08_COVARNET_MODE", "preflight")
SEED <- 20260906L
selected <- readLines(file.path(input_dir,"selected_K.txt"))
kline <- grep("^K_FINAL=",selected,value=TRUE)
if(length(kline)!=1L) stop("Missing unique selected K")
K <- as.integer(sub("^K_FINAL=","",kline))
if(!is.finite(K) || K<2L) stop("Invalid selected K")
workers <- as.integer(Sys.getenv("COVARNET_WORKERS","2"))
if(!is.finite(workers) || workers<1L) stop("Invalid COVARNET_WORKERS")
rank_file <- function(prefix,extension) sprintf("%s_K%02d.%s",prefix,K,extension)
N_BOOTSTRAP <- 999L
N_PREFLIGHT <- 3L
N_RUN <- 200L
TOP_VALUES <- c(5L, 10L)
RUN_STRUCTURAL <- Sys.getenv("G08_RUN_STRUCTURAL", "0") == "1"
STRUCTURAL_DRAWS <- as.integer(Sys.getenv("G08_STRUCTURAL_DRAWS", "999"))
if (!MODE %in% c("preflight", "full")) stop("G08_COVARNET_MODE must be preflight or full")

read_matrix <- function(path) {
  frame <- read.delim(path, check.names = FALSE)
  default_rows <- identical(rownames(frame), as.character(seq_len(nrow(frame))))
  if (default_rows) {
    row_ids <- as.character(frame[[1L]])
    frame[[1L]] <- NULL
  } else {
    # Files written with col.names=NA are read back with the first field as
    # row.names automatically; no data column should be discarded.
    row_ids <- rownames(frame)
  }
  matrix <- as.matrix(frame)
  storage.mode(matrix) <- "double"
  rownames(matrix) <- row_ids
  matrix
}

sha256 <- function(path) {
  output <- system2("sha256sum", shQuote(path), stdout = TRUE)
  strsplit(output, "[[:space:]]+")[[1L]][1L]
}

normalize_columns <- function(matrix) {
  denominator <- sqrt(colSums(matrix^2))
  out <- sweep(matrix, 2L, pmax(denominator, .Machine$double.eps), "/")
  out[!is.finite(out)] <- 0
  out
}

hungarian_min <- function(cost) {
  n <- nrow(cost)
  m <- ncol(cost)
  if (n > m) stop("hungarian_min requires nrow <= ncol")
  u <- numeric(n + 1L)
  v <- numeric(m + 1L)
  p <- integer(m + 1L)
  way <- integer(m + 1L)
  for (i in seq_len(n)) {
    p[1L] <- i
    j0 <- 1L
    minv <- rep(Inf, m + 1L)
    used <- rep(FALSE, m + 1L)
    repeat {
      used[j0] <- TRUE
      i0 <- p[j0]
      delta <- Inf
      j1 <- NA_integer_
      for (jj in 2:(m + 1L)) {
        if (used[jj]) next
        current <- cost[i0, jj - 1L] - u[i0 + 1L] - v[jj]
        if (current < minv[jj]) {
          minv[jj] <- current
          way[jj] <- j0
        }
        if (minv[jj] < delta) {
          delta <- minv[jj]
          j1 <- jj
        }
      }
      for (jj in seq_len(m + 1L)) {
        if (used[jj]) {
          u[p[jj] + 1L] <- u[p[jj] + 1L] + delta
          v[jj] <- v[jj] - delta
        } else {
          minv[jj] <- minv[jj] - delta
        }
      }
      j0 <- j1
      if (p[j0] == 0L) break
    }
    repeat {
      j1 <- way[j0]
      p[j0] <- p[j1]
      j0 <- j1
      if (j0 == 1L) break
    }
  }
  assignment <- integer(n)
  for (jj in 2:(m + 1L)) if (p[jj] > 0L) assignment[p[jj]] <- jj - 1L
  assignment
}

max_weight_assignment <- function(score) {
  n <- max(nrow(score), ncol(score))
  padded <- matrix(0, n, n)
  padded[seq_len(nrow(score)), seq_len(ncol(score))] <- score
  hungarian_min(max(padded) - padded)[seq_len(nrow(score))]
}

compare_basis <- function(reference, candidate, draw) {
  reference <- normalize_columns(reference)
  candidate <- normalize_columns(candidate[rownames(reference), , drop = FALSE])
  similarity <- t(reference) %*% candidate
  assignment <- max_weight_assignment(similarity)
  rows <- list()
  for (module_index in seq_len(ncol(reference))) {
    matched <- assignment[module_index]
    reference_top <- rownames(reference)[order(reference[, module_index], decreasing = TRUE)[seq_len(10L)]]
    candidate_top <- rownames(candidate)[order(candidate[, matched], decreasing = TRUE)[seq_len(10L)]]
    rows[[module_index]] <- data.frame(
      draw = draw, reference_module = colnames(reference)[module_index],
      matched_component = matched,
      cosine_similarity = similarity[module_index, matched],
      top10_jaccard = length(intersect(reference_top, candidate_top)) /
        length(union(reference_top, candidate_top)),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

median_pairwise_spearman <- function(matrix, state_indices) {
  if (length(state_indices) < 2L) return(NA_real_)
  correlations <- suppressWarnings(cor(t(matrix[state_indices, , drop = FALSE]),
                                       method = "spearman", use = "pairwise.complete.obs"))
  values <- correlations[upper.tri(correlations)]
  values <- values[is.finite(values)]
  if (length(values)) median(values) else NA_real_
}

rank_parent_residuals <- function(proportions, parent, studies) {
  lineages <- unique(parent)
  parent_total <- matrix(0, nrow = length(lineages), ncol = ncol(proportions),
                         dimnames = list(lineages, colnames(proportions)))
  for (lineage in lineages) {
    parent_total[lineage, ] <- colSums(proportions[parent == lineage, , drop = FALSE])
  }
  out <- matrix(NA_real_, nrow = nrow(proportions), ncol = ncol(proportions),
                dimnames = dimnames(proportions))
  for (i in seq_len(nrow(proportions))) {
    frame <- data.frame(
      y = rank(proportions[i, ], ties.method = "average"),
      parent = rank(parent_total[parent[i], ], ties.method = "average"),
      study = factor(studies)
    )
    fit <- lm(y ~ parent + study, data = frame, na.action = na.exclude)
    out[i, ] <- stats::residuals(fit)
  }
  out
}

required <- c(
  "NMF_input_minmax.tsv", rank_file("final_nmf","rds"), rank_file("module_node_weights","tsv"),
  rank_file("sample_module_abundance","tsv"), "CoVarNet_node_proportion_of_all_atlas_cells.tsv",
  "CoVarNet_node_counts.tsv", "CoVarNet_sample_metadata.tsv"
)
missing <- required[!file.exists(file.path(input_dir, required))]
if (length(missing)) stop("Missing CoVarNet frozen inputs: ", paste(missing, collapse = ", "))

nmf_input <- read_matrix(file.path(input_dir, "NMF_input_minmax.tsv"))
proportions <- read_matrix(file.path(input_dir, "CoVarNet_node_proportion_of_all_atlas_cells.tsv"))
counts <- read_matrix(file.path(input_dir, "CoVarNet_node_counts.tsv"))
weights <- read_matrix(file.path(input_dir, rank_file("module_node_weights","tsv")))
module_abundance_frame <- read.delim(file.path(input_dir, rank_file("sample_module_abundance","tsv")), check.names = FALSE)
module_ids <- module_abundance_frame[[1L]]
module_abundance_frame[[1L]] <- NULL
module_abundance <- t(as.matrix(module_abundance_frame))
storage.mode(module_abundance) <- "double"
rownames(module_abundance) <- colnames(module_abundance_frame)
colnames(module_abundance) <- module_ids
metadata <- read.delim(file.path(input_dir, "CoVarNet_sample_metadata.tsv"), check.names = FALSE)
rownames(metadata) <- metadata$covarnet_sample_uid

stopifnot(nrow(nmf_input)>=5L,ncol(nmf_input)>=20L,K<nrow(nmf_input),K<ncol(nmf_input))
stopifnot(!anyDuplicated(rownames(nmf_input)),!anyDuplicated(colnames(nmf_input)))
stopifnot(all(is.finite(nmf_input)),all(nmf_input>=0 & nmf_input<=1),all(is.finite(proportions)),all(proportions>=0 & proportions<=1),all(is.finite(counts)),all(counts>=0),all(counts==round(counts)))
stopifnot(identical(rownames(nmf_input), rownames(proportions)))
stopifnot(identical(colnames(nmf_input), colnames(proportions)))
stopifnot(identical(rownames(nmf_input), rownames(counts)))
stopifnot(identical(colnames(nmf_input), colnames(counts)))
stopifnot(all(colnames(nmf_input) %in% rownames(metadata)))
metadata <- metadata[colnames(nmf_input), , drop = FALSE]
stopifnot(identical(metadata$covarnet_sample_uid, colnames(nmf_input)))

row_min <- apply(proportions, 1L, min, na.rm = TRUE)
row_max <- apply(proportions, 1L, max, na.rm = TRUE)
reconstructed <- sweep(proportions, 1L, row_min, "-")
reconstructed <- sweep(reconstructed, 1L, row_max - row_min, "/")
reconstructed[!is.finite(reconstructed)] <- 0
transformation_error <- max(abs(reconstructed - nmf_input))
if (!is.finite(transformation_error) || transformation_error > 1e-12) {
  stop("Frozen NMF matrix does not match row-wise min-max transformation: ", transformation_error)
}

reference_fit <- readRDS(file.path(input_dir, rank_file("final_nmf","rds")))
reference_basis <- basis(reference_fit)
colnames(reference_basis) <- sprintf("CM%02d", seq_len(ncol(reference_basis)))
reference_export_basis <- sweep(reference_basis,2L,colSums(reference_basis),"/")
stopifnot(all(is.finite(reference_export_basis)),ncol(reference_basis)==K)
weight_error <- max(abs(reference_export_basis - weights[rownames(reference_basis), colnames(reference_basis)]))
if (!is.finite(weight_error) || weight_error > 1e-12) stop("Frozen module weight mismatch")

if (!all(colnames(nmf_input) %in% colnames(module_abundance))) stop("Module-abundance profiles do not align")
module_abundance <- module_abundance[, colnames(nmf_input), drop = FALSE]

metadata$resolved_person_key <- metadata$donor_uid_v2
bad_donor <- is.na(metadata$resolved_person_key) |
  tolower(trimws(metadata$resolved_person_key)) %in% c("", "na", "nan", "none", "nr", "unknown", "__missing__", "__not_available__")
metadata$donor_resolution_G08 <- ifelse(bad_donor, "unresolved_placeholder", "resolved")
metadata$resolved_person_key[bad_donor] <- paste0(
  "UNRESOLVED::", metadata$study_canonical[bad_donor], "::", metadata$covarnet_sample_uid[bad_donor]
)

cross_study <- aggregate(study_canonical ~ resolved_person_key,
                         metadata[metadata$donor_resolution_G08 == "resolved", ],
                         function(x) length(unique(x)))
cross_study <- cross_study[cross_study$study_canonical > 1L, , drop = FALSE]
if (nrow(cross_study)) stop("Resolved donor keys cross studies; global joint block sampling is required")

study_donors <- split(metadata$resolved_person_key, metadata$study_canonical)
study_donors <- lapply(study_donors, unique)

make_draw <- function(draw) {
  set.seed(SEED + draw)
  selections <- list()
  columns <- character()
  counter <- 1L
  for (study in sort(names(study_donors))) {
    donors <- sort(study_donors[[study]])
    sampled <- sample(donors, length(donors), replace = TRUE)
    for (position in seq_along(sampled)) {
      donor <- sampled[position]
      profiles <- metadata$covarnet_sample_uid[
        metadata$study_canonical == study & metadata$resolved_person_key == donor
      ]
      selections[[counter]] <- data.frame(
        draw = draw, study_canonical = study, draw_position = position,
        selected_person_key = donor,
        donor_resolution = unique(metadata$donor_resolution_G08[
          metadata$study_canonical == study & metadata$resolved_person_key == donor
        ])[1L],
        n_profiles_retained = length(profiles),
        profile_ids = paste(profiles, collapse = ";"), stringsAsFactors = FALSE
      )
      columns <- c(columns, profiles)
      counter <- counter + 1L
    }
  }
  list(columns = columns, selection = do.call(rbind, selections))
}

source_audit <- data.frame(
  item = c("NMF matrix", "baseline NMF fit", "module weights", "raw proportions", "counts", "sample metadata"),
  file = required[c(1, 2, 3, 5, 6, 7)],
  sha256 = vapply(file.path(input_dir, required[c(1, 2, 3, 5, 6, 7)]), sha256, character(1)),
  rows = c(nrow(nmf_input), NA, nrow(weights), nrow(proportions), nrow(counts), nrow(metadata)),
  columns = c(ncol(nmf_input), NA, ncol(weights), ncol(proportions), ncol(counts), ncol(metadata)),
  stringsAsFactors = FALSE
)
write.table(source_audit, file.path(out_dir, "covarnet_source_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

if (MODE == "preflight") {
  timing <- list()
  comparison <- list()
  selection_rows <- list()
  for (draw in seq_len(N_PREFLIGHT)) {
    sampled <- make_draw(draw)
    selection_rows[[draw]] <- sampled$selection
    matrix <- nmf_input[, sampled$columns, drop = FALSE]
    colnames(matrix) <- make.unique(paste0("draw", draw, "::", sampled$columns))
    started <- proc.time()
    set.seed(SEED + 10000L + draw)
    fit <- nmf(matrix, rank = K, method = "nsNMF", nrun = N_RUN,
               seed = SEED + 10000L + draw, .options = paste0("p",workers))
    elapsed <- proc.time() - started
    timing[[draw]] <- data.frame(
      draw = draw, n_profiles = ncol(matrix), nrun = N_RUN,
      elapsed_seconds = elapsed[["elapsed"]], user_seconds = elapsed[["user.self"]],
      system_seconds = elapsed[["sys.self"]], stringsAsFactors = FALSE
    )
    comparison[[draw]] <- compare_basis(reference_basis, basis(fit), draw)
  }
  timing <- do.call(rbind, timing)
  timing$projected_999_wall_hours <- median(timing$elapsed_seconds) * 999 / 3600
  timing$projected_999_CPU_hours_workers <- timing$projected_999_wall_hours * workers
  write.table(timing, file.path(out_dir, "structural_refit_timing_preflight.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(do.call(rbind, comparison), file.path(out_dir, "structural_refit_preflight_matches.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(do.call(rbind, selection_rows), file.path(out_dir, "structural_refit_preflight_draw_index.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  status <- if (max(timing$projected_999_CPU_hours_workers) <= 24) "within_budget" else "paused_over_24_CPU_hour_budget"
  writeLines(c(
    paste0("status=", status),
    paste0("median_elapsed_seconds=", median(timing$elapsed_seconds)),
    paste0("projected_999_wall_hours=", median(timing$projected_999_wall_hours)),
    paste0("projected_999_CPU_hours_workers=", median(timing$projected_999_CPU_hours_workers)),
    paste0("projection_assumes_",workers,"_workers_and_is_not_a_runtime_guarantee")
  ), file.path(out_dir, "structural_refit_budget_decision.txt"))
  cat("CoVarNet structural preflight complete; status:", status, "\n")
  quit(save = "no", status = 0)
}

parent_table <- read.delim(file.path(input_dir,"CoVarNet_node_parent_map.tsv"),check.names=FALSE)
stopifnot(!anyDuplicated(parent_table$cell_state),setequal(parent_table$cell_state,rownames(proportions)))
parent <- parent_table$parent_lineage[match(rownames(proportions),parent_table$cell_state)]
stopifnot(!anyNA(parent),all(nzchar(parent)))
names(parent) <- rownames(proportions)
top_indices <- lapply(colnames(reference_basis), function(module) {
  order(reference_basis[, module], decreasing = TRUE)
})
names(top_indices) <- colnames(reference_basis)

selection_rows <- vector("list", N_BOOTSTRAP)
projection_rows <- vector("list", N_BOOTSTRAP)
structural_rows <- list()
structural_failures <- list()

for (draw in seq_len(N_BOOTSTRAP)) {
  sampled <- make_draw(draw)
  selection_rows[[draw]] <- sampled$selection
  cols <- sampled$columns
  raw <- proportions[, cols, drop = FALSE]
  boot_input <- nmf_input[, cols, drop = FALSE]
  boot_module <- module_abundance[, cols, drop = FALSE]
  studies <- metadata[cols, "study_canonical"]
  parent_adjusted <- rank_parent_residuals(raw, parent, studies)
  summaries <- list()
  row_counter <- 1L
  for (module in colnames(reference_basis)) {
    indices <- top_indices[[module]]
    for (top_n in TOP_VALUES) {
      selected <- indices[seq_len(top_n)]
      summaries[[row_counter]] <- data.frame(
        draw = draw, module = module, top_n = top_n,
        n_profiles = length(cols),
        mean_frozen_module_abundance = mean(boot_module[module, ], na.rm = TRUE),
        median_frozen_module_abundance = median(boot_module[module, ], na.rm = TRUE),
        raw_median_pairwise_spearman = median_pairwise_spearman(raw, selected),
        parent_adjusted_median_pairwise_spearman = median_pairwise_spearman(parent_adjusted, selected),
        stringsAsFactors = FALSE
      )
      row_counter <- row_counter + 1L
    }
  }
  projection_rows[[draw]] <- do.call(rbind, summaries)

  if (RUN_STRUCTURAL && draw <= STRUCTURAL_DRAWS) {
    matrix <- boot_input
    colnames(matrix) <- make.unique(paste0("draw", draw, "::", cols))
    fit <- tryCatch({
      set.seed(SEED + 10000L + draw)
      nmf(matrix, rank = K, method = "nsNMF", nrun = N_RUN,
          seed = SEED + 10000L + draw, .options = paste0("p",workers))
    }, error = function(e) e)
    if (inherits(fit, "error")) {
      structural_failures[[length(structural_failures) + 1L]] <- data.frame(
        draw = draw, reason = conditionMessage(fit), stringsAsFactors = FALSE
      )
    } else {
      structural_rows[[length(structural_rows) + 1L]] <- compare_basis(reference_basis, basis(fit), draw)
    }
  }
  if (draw == 1L || draw %% 100L == 0L || draw == N_BOOTSTRAP) {
    message("CoVarNet donor-block projection: ", draw, "/", N_BOOTSTRAP)
  }
}

selection <- do.call(rbind, selection_rows)
projection <- do.call(rbind, projection_rows)
selection_connection <- gzfile(file.path(out_dir, "donor_block_draw_index.tsv.gz"), "wt")
write.table(selection, selection_connection, sep = "\t", quote = FALSE, row.names = FALSE)
close(selection_connection)
projection_connection <- gzfile(file.path(out_dir, "frozen_module_donor_block_distribution.tsv.gz"), "wt")
write.table(projection, projection_connection, sep = "\t", quote = FALSE, row.names = FALSE)
close(projection_connection)

quantile_value <- function(x, probability) quantile(x, probability, na.rm = TRUE, names = FALSE)
summary_rows <- list()
counter <- 1L
for (module in unique(projection$module)) {
  for (top_n in unique(projection$top_n)) {
    x <- projection[projection$module == module & projection$top_n == top_n, ]
    for (metric in c(
      "mean_frozen_module_abundance", "median_frozen_module_abundance",
      "raw_median_pairwise_spearman", "parent_adjusted_median_pairwise_spearman"
    )) {
      values <- x[[metric]]
      summary_rows[[counter]] <- data.frame(
        module = module, top_n = top_n, metric = metric,
        n_attempted = N_BOOTSTRAP, n_valid = sum(is.finite(values)),
        median = median(values, na.rm = TRUE),
        q025 = quantile_value(values, 0.025), q975 = quantile_value(values, 0.975),
        q05 = quantile_value(values, 0.05), q95 = quantile_value(values, 0.95),
        stringsAsFactors = FALSE
      )
      counter <- counter + 1L
    }
  }
}
write.table(do.call(rbind, summary_rows), file.path(out_dir, "frozen_module_donor_block_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

if (length(structural_rows)) {
  structural <- do.call(rbind, structural_rows)
  structural_connection <- gzfile(file.path(out_dir, "structural_refit_donor_block_matches.tsv.gz"), "wt")
  write.table(structural, structural_connection, sep = "\t", quote = FALSE, row.names = FALSE)
  close(structural_connection)
  structural_summary <- aggregate(
    cbind(cosine_similarity, top10_jaccard) ~ reference_module,
    structural,
    function(x) c(n_valid = sum(is.finite(x)), median = median(x, na.rm = TRUE),
                  q025 = quantile_value(x, 0.025), q975 = quantile_value(x, 0.975),
                  fraction_lt_0.8 = mean(x < 0.8, na.rm = TRUE))
  )
  write.table(structural_summary, file.path(out_dir, "structural_refit_donor_block_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
if (length(structural_failures)) {
  write.table(do.call(rbind, structural_failures), file.path(out_dir, "structural_refit_failures.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

manifest <- list(
  status = "completed",
  analysis_id = "G08_COVARNET_DONOR_BLOCK_V1",
  mode = MODE,
  seed = SEED,
  n_bootstrap_attempts = N_BOOTSTRAP,
  n_nodes = nrow(nmf_input), n_profiles = ncol(nmf_input), rank = K,
  algorithm = "nsNMF", nrun = N_RUN,
  transformation = "row-wise min-max of all-atlas-cell proportions",
  transformation_max_abs_error = transformation_error,
  baseline_weight_max_abs_error = weight_error,
  n_resolved_people = length(unique(metadata$resolved_person_key[metadata$donor_resolution_G08 == "resolved"])),
  n_unresolved_profile_units = sum(metadata$donor_resolution_G08 != "resolved"),
  cross_study_resolved_people = nrow(cross_study),
  structural_refit_requested = RUN_STRUCTURAL,
  structural_draws_requested = if (RUN_STRUCTURAL) STRUCTURAL_DRAWS else 0L,
  structural_draws_valid = if (length(structural_rows)) length(unique(do.call(rbind, structural_rows)$draw)) else 0L,
  structural_draws_failed = length(structural_failures),
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  R = R.version.string, NMF = as.character(packageVersion("NMF"))
)
writeLines(paste(capture.output(dput(manifest)), collapse = "\n"), file.path(out_dir, "run_manifest_dput.txt"))
cat("CoVarNet donor-block full analysis complete\n")
