rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(NMF)
  library(ggplot2)
})

# =============================================================================
# CoVarNet follow-up validation
#
# This script starts from CoVarNet_final_v3. It never edits the atlas or the
# original CoVarNet outputs. Every new file is written below FOLLOWUP_ROOT.
# The workflow is resumable: expensive NMF fits are cached as RDS files.
# =============================================================================

ROOT <- "/home/anzhen_vhd/scRNAseq/codex/CoVarNet_final_v3"
PRIMARY_ROOT <- file.path(ROOT, "stable_NMF")
FOLLOWUP_ROOT <- file.path(ROOT, "followup_validation_v1")
CACHE_ROOT <- file.path(FOLLOWUP_ROOT, "cache_nmf_fits")
TABLE_ROOT <- file.path(FOLLOWUP_ROOT, "tables")
FIGURE_ROOT <- file.path(FOLLOWUP_ROOT, "figures")
GENETIC_SEARCH_ROOT <- "/home/anzhen_vhd/更新版MAGMA_scDRS"

dir.create(CACHE_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_ROOT, recursive = TRUE, showWarnings = FALSE)

SEED <- 20260807L
K_PRIMARY <- 14L
ADJACENT_RANKS <- c(12L, 13L, 15L, 16L)
ADJACENT_NRUN <- 100L
SENSITIVITY_NRUN <- 30L
N_SAMPLE_BOOTSTRAP <- 30L
N_STUDY_BOOTSTRAP <- 30L
N_DOWNSAMPLE <- 10L
DOWNSAMPLE_DEPTH <- 500L
TOP_NODES <- 10L
TOP_NODES_COHERENCE <- 5L
MIN_WITHIN_STUDY_N <- 5L
MIN_DISEASE_GROUP_N <- 2L
N_GENETIC_PERM <- 10000L

workers <- max(1L, min(8L, parallel::detectCores(logical = FALSE) - 1L))
nmf_options <- paste0("p", workers)
cat("NMF workers:", workers, "\n")
cat("Output:", FOLLOWUP_ROOT, "\n")

required_files <- c(
  file.path(ROOT, "CoVarNet_node_proportion_of_all_atlas_cells.tsv"),
  file.path(ROOT, "CoVarNet_node_proportion_among_stable_nodes_sensitivity.tsv"),
  file.path(ROOT, "CoVarNet_node_counts.tsv"),
  file.path(ROOT, "CoVarNet_sample_metadata.tsv"),
  file.path(PRIMARY_ROOT, sprintf("final_nmf_K%02d.rds", K_PRIMARY))
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}


# =============================================================================
# 1. Helpers
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

normalize_columns <- function(m) {
  z <- sweep(m, 2L, colSums(m), "/")
  z[!is.finite(z)] <- 0
  z
}

row_minmax <- function(m) {
  lo <- apply(m, 1L, min, na.rm = TRUE)
  hi <- apply(m, 1L, max, na.rm = TRUE)
  out <- sweep(m, 1L, lo, "-")
  out <- sweep(out, 1L, pmax(hi - lo, .Machine$double.eps), "/")
  out[!is.finite(out)] <- 0
  out
}

robust_row_minmax <- function(m, lower = 0.01, upper = 0.99) {
  lo <- apply(m, 1L, quantile, probs = lower, na.rm = TRUE, names = FALSE)
  hi <- apply(m, 1L, quantile, probs = upper, na.rm = TRUE, names = FALSE)
  out <- sweep(m, 1L, lo, "-")
  out <- sweep(out, 1L, pmax(hi - lo, .Machine$double.eps), "/")
  # Matrix subassignment preserves dimensions; pmin/pmax with a scalar first
  # argument can drop the dim attribute in some R versions.
  out[out < 0] <- 0
  out[out > 1] <- 1
  out[!is.finite(out)] <- 0
  out
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("_+", "_", x)
}

fit_nmf_cached <- function(m, rank, nrun, label, seed) {
  path <- file.path(CACHE_ROOT, paste0(safe_name(label), "_K", rank, ".rds"))
  if (file.exists(path)) {
    cat("Reuse:", basename(path), "\n")
    return(readRDS(path))
  }
  cat("Fit:", label, "K=", rank, "samples=", ncol(m), "nrun=", nrun, "\n")
  set.seed(seed)
  fit <- nmf(
    m, rank = rank, method = "nsNMF", nrun = nrun,
    seed = seed, .options = nmf_options
  )
  saveRDS(fit, path)
  fit
}

cosine_similarity <- function(A, B) {
  common <- intersect(rownames(A), rownames(B))
  A <- A[common, , drop = FALSE]
  B <- B[common, , drop = FALSE]
  A <- sweep(A, 2L, sqrt(colSums(A^2)), "/")
  B <- sweep(B, 2L, sqrt(colSums(B^2)), "/")
  A[!is.finite(A)] <- 0
  B[!is.finite(B)] <- 0
  t(A) %*% B
}

# Dependency-free Hungarian algorithm. This performs a globally optimal,
# one-to-one assignment instead of matching every reference module to its
# independently best candidate.
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
  for (jj in 2:(m + 1L)) {
    if (p[jj] > 0L) assignment[p[jj]] <- jj - 1L
  }
  assignment
}

max_weight_assignment <- function(score) {
  nr <- nrow(score)
  nc <- ncol(score)
  n <- max(nr, nc)
  padded <- matrix(0, n, n)
  padded[seq_len(nr), seq_len(nc)] <- score
  assignment <- hungarian_min(max(padded) - padded)
  assignment[seq_len(nr)]
}

compare_basis_one_to_one <- function(reference, candidate, run_type, run_id) {
  candidate <- candidate[rownames(reference), , drop = FALSE]
  sim <- cosine_similarity(reference, candidate)
  assignment <- max_weight_assignment(sim)
  ref_top <- lapply(seq_len(ncol(reference)), function(k) {
    rownames(reference)[order(reference[, k], decreasing = TRUE)[
      seq_len(min(TOP_NODES, nrow(reference)))
    ]]
  })

  rbindlist(lapply(seq_len(ncol(reference)), function(k) {
    j <- assignment[k]
    if (j > ncol(candidate)) {
      return(data.table(
        run_type = run_type, run_id = as.character(run_id),
        reference_module = colnames(reference)[k], matched_module = NA_character_,
        cosine_similarity = 0, top_node_jaccard = 0
      ))
    }
    cand_top <- rownames(candidate)[order(candidate[, j], decreasing = TRUE)[
      seq_len(min(TOP_NODES, nrow(candidate)))
    ]]
    data.table(
      run_type = run_type, run_id = as.character(run_id),
      reference_module = colnames(reference)[k],
      matched_module = if (!is.null(colnames(candidate))) colnames(candidate)[j] else as.character(j),
      cosine_similarity = sim[k, j],
      top_node_jaccard = length(intersect(ref_top[[k]], cand_top)) /
        length(union(ref_top[[k]], cand_top))
    )
  }))
}

random_effect_meta <- function(effect, variance) {
  ok <- is.finite(effect) & is.finite(variance) & variance > 0
  y <- effect[ok]
  v <- variance[ok]
  k <- length(y)
  if (k < 2L) {
    return(list(k = k, estimate = NA_real_, se = NA_real_, z = NA_real_,
                p = NA_real_, tau2 = NA_real_, I2 = NA_real_, Q = NA_real_))
  }
  w <- 1 / v
  fixed <- sum(w * y) / sum(w)
  Q <- sum(w * (y - fixed)^2)
  df <- k - 1L
  C <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - df) / pmax(C, .Machine$double.eps))
  wr <- 1 / (v + tau2)
  estimate <- sum(wr * y) / sum(wr)
  se <- sqrt(1 / sum(wr))
  z <- estimate / se
  p <- 2 * pnorm(-abs(z))
  I2 <- if (Q > 0) max(0, (Q - df) / Q) * 100 else 0
  list(k = k, estimate = estimate, se = se, z = z, p = p,
       tau2 = tau2, I2 = I2, Q = Q)
}

save_plot <- function(plot, filename, width, height) {
  pdf_path <- file.path(FIGURE_ROOT, paste0(filename, ".pdf"))
  tryCatch(
    ggsave(pdf_path, plot, width = width, height = height, device = cairo_pdf),
    error = function(e) {
      warning("cairo_pdf unavailable; using standard PDF device: ", conditionMessage(e))
      ggsave(pdf_path, plot, width = width, height = height, device = "pdf")
    }
  )
  ggsave(file.path(FIGURE_ROOT, paste0(filename, ".png")), plot,
         width = width, height = height, dpi = 300)
}

normalize_key <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("^.*\\|", "", x)
  gsub("[^a-z0-9]", "", x)
}


# =============================================================================
# 2. Input alignment and primary reference
# =============================================================================

prop <- read_matrix(file.path(ROOT, "CoVarNet_node_proportion_of_all_atlas_cells.tsv"))
prop_closed <- read_matrix(file.path(ROOT, "CoVarNet_node_proportion_among_stable_nodes_sensitivity.tsv"))
counts <- read_matrix(file.path(ROOT, "CoVarNet_node_counts.tsv"))
meta <- fread(file.path(ROOT, "CoVarNet_sample_metadata.tsv"))

common_samples <- Reduce(intersect, list(
  colnames(prop), colnames(prop_closed), colnames(counts), meta$covarnet_sample_uid
))
prop <- prop[, common_samples, drop = FALSE]
prop_closed <- prop_closed[rownames(prop), common_samples, drop = FALSE]
counts <- counts[rownames(prop), common_samples, drop = FALSE]
meta <- meta[match(common_samples, covarnet_sample_uid)]
stopifnot(all(meta$covarnet_sample_uid == common_samples))

primary_matrix <- row_minmax(prop)
primary_fit <- readRDS(file.path(PRIMARY_ROOT, sprintf("final_nmf_K%02d.rds", K_PRIMARY)))
W_reference <- normalize_columns(basis(primary_fit))
colnames(W_reference) <- sprintf("CM%02d", seq_len(ncol(W_reference)))
H_reference <- normalize_columns(coef(primary_fit))
rownames(H_reference) <- colnames(W_reference)

if (!identical(rownames(W_reference), rownames(primary_matrix))) {
  W_reference <- W_reference[rownames(primary_matrix), , drop = FALSE]
}

cat("Nodes:", nrow(primary_matrix), " Samples:", ncol(primary_matrix),
    " Studies:", uniqueN(meta$study_canonical),
    " Tissues:", uniqueN(meta$tissue_v2), "\n")


# =============================================================================
# 3. Adjacent-rank analysis with exact one-to-one matching
# =============================================================================

adjacent_rows <- list()
for (k in ADJACENT_RANKS) {
  fit <- fit_nmf_cached(
    primary_matrix, k, ADJACENT_NRUN,
    paste0("adjacent_rank_", k), SEED + 100L + k
  )
  W <- normalize_columns(basis(fit))
  colnames(W) <- sprintf("K%02d_M%02d", k, seq_len(ncol(W)))
  adjacent_rows[[as.character(k)]] <- compare_basis_one_to_one(
    W_reference, W, "adjacent_rank", paste0("K", k)
  )
}
adjacent <- rbindlist(adjacent_rows)
fwrite(adjacent, file.path(TABLE_ROOT, "01_adjacent_rank_one_to_one_matching.tsv"), sep = "\t")


# =============================================================================
# 4. Strict stability analyses
# =============================================================================

stability_rows <- list()

run_sensitivity_fit <- function(m, type, id, seed_offset) {
  if (is.null(dim(m)) || length(dim(m)) != 2L) {
    stop(
      "Sensitivity input is not a matrix: type=", type,
      "; id=", id, "; class=", paste(class(m), collapse = ",")
    )
  }
  if (ncol(m) < max(25L, 2L * K_PRIMARY)) return(NULL)
  fit <- tryCatch(
    fit_nmf_cached(
      m, K_PRIMARY, SENSITIVITY_NRUN,
      paste(type, id, sep = "__"), SEED + seed_offset
    ),
    error = function(e) {
      warning(type, " ", id, " failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(fit)) return(NULL)
  W <- normalize_columns(basis(fit))
  colnames(W) <- sprintf("candidate_%02d", seq_len(ncol(W)))
  compare_basis_one_to_one(W_reference, W, type, id)
}

# Alternative denominator and transformations.
transformations <- list(
  stable_node_denominator = row_minmax(prop_closed),
  robust_minmax_q01_q99 = robust_row_minmax(prop),
  hellinger_closed_composition = sqrt(sweep(
    prop, 2L, pmax(colSums(prop), .Machine$double.eps), "/"
  ))
)
for (nm in names(transformations)) {
  stability_rows[[paste0("transform_", nm)]] <- run_sensitivity_fit(
    transformations[[nm]], "transformation", nm, 1000L + match(nm, names(transformations))
  )
}

# Leave one study out.
studies <- unique(meta$study_canonical)
for (i in seq_along(studies)) {
  keep <- meta$study_canonical != studies[i]
  stability_rows[[paste0("LOS_", i)]] <- run_sensitivity_fit(
    primary_matrix[, keep, drop = FALSE], "leave_one_study", studies[i], 2000L + i
  )
}

# Leave one tissue out.
tissues <- unique(meta$tissue_v2)
for (i in seq_along(tissues)) {
  keep <- meta$tissue_v2 != tissues[i]
  stability_rows[[paste0("LOT_", i)]] <- run_sensitivity_fit(
    primary_matrix[, keep, drop = FALSE], "leave_one_tissue", tissues[i], 3000L + i
  )
}

# Fit healthy-only and disease-only atlases.
conditions <- intersect(c("healthy", "disease"), unique(meta$condition_v2))
for (i in seq_along(conditions)) {
  keep <- meta$condition_v2 == conditions[i]
  stability_rows[[paste0("LOC_", i)]] <- run_sensitivity_fit(
    primary_matrix[, keep, drop = FALSE], "condition_only", conditions[i], 4000L + i
  )
}

# Ordinary sample bootstrap, now with one-to-one module assignment.
set.seed(SEED + 5000L)
for (b in seq_len(N_SAMPLE_BOOTSTRAP)) {
  idx <- sample(seq_len(ncol(primary_matrix)), replace = TRUE)
  boot <- primary_matrix[, idx, drop = FALSE]
  colnames(boot) <- make.unique(colnames(boot))
  stability_rows[[paste0("sample_boot_", b)]] <- run_sensitivity_fit(
    boot, "sample_bootstrap", b, 5000L + b
  )
}

# Cluster bootstrap resamples studies, then samples within each selected study.
set.seed(SEED + 6000L)
for (b in seq_len(N_STUDY_BOOTSTRAP)) {
  drawn_studies <- sample(studies, length(studies), replace = TRUE)
  idx <- unlist(lapply(drawn_studies, function(s) {
    candidates <- which(meta$study_canonical == s)
    sample(candidates, length(candidates), replace = TRUE)
  }), use.names = FALSE)
  boot <- primary_matrix[, idx, drop = FALSE]
  colnames(boot) <- make.unique(colnames(boot))
  stability_rows[[paste0("study_boot_", b)]] <- run_sensitivity_fit(
    boot, "hierarchical_study_bootstrap", b, 6000L + b
  )
}

# Multinomial cell-depth sensitivity. Stable-node counts plus an "other cell"
# category are sampled to preserve the all-atlas denominator.
eligible_depth <- meta$n_total_atlas_cells >= DOWNSAMPLE_DEPTH
if (sum(eligible_depth) >= max(30L, 2L * K_PRIMARY)) {
  set.seed(SEED + 7000L)
  for (b in seq_len(N_DOWNSAMPLE)) {
    ds <- matrix(
      0, nrow = nrow(counts), ncol = sum(eligible_depth),
      dimnames = list(rownames(counts), colnames(counts)[eligible_depth])
    )
    eligible_indices <- which(eligible_depth)
    for (jj in seq_along(eligible_indices)) {
      j <- eligible_indices[jj]
      total <- meta$n_total_atlas_cells[j]
      other <- max(0, total - sum(counts[, j]))
      probability <- c(counts[, j], other)
      probability <- probability / sum(probability)
      draw <- as.vector(rmultinom(1L, DOWNSAMPLE_DEPTH, probability))
      ds[, jj] <- draw[seq_len(nrow(counts))] / DOWNSAMPLE_DEPTH
    }
    stability_rows[[paste0("downsample_", b)]] <- run_sensitivity_fit(
      row_minmax(ds), "cell_depth_downsample", paste0("depth", DOWNSAMPLE_DEPTH, "_rep", b),
      7000L + b
    )
  }
}

stability <- rbindlist(stability_rows, fill = TRUE)
fwrite(stability, file.path(TABLE_ROOT, "02_strict_stability_one_to_one.tsv"), sep = "\t")

stability_summary <- stability[, .(
  n_runs = .N,
  median_cosine = median(cosine_similarity, na.rm = TRUE),
  q05_cosine = as.numeric(quantile(cosine_similarity, 0.05, na.rm = TRUE)),
  minimum_cosine = min(cosine_similarity, na.rm = TRUE),
  median_top10_jaccard = median(top_node_jaccard, na.rm = TRUE),
  fraction_cosine_ge_0.80 = mean(cosine_similarity >= 0.80, na.rm = TRUE)
), by = .(run_type, reference_module)]
fwrite(stability_summary, file.path(TABLE_ROOT, "03_strict_stability_summary.tsv"), sep = "\t")


# =============================================================================
# 5. Study/tissue concentration of module abundance
# =============================================================================

sample_modules <- cbind(meta, as.data.table(t(H_reference)))
concentration <- rbindlist(lapply(colnames(W_reference), function(cm) {
  study <- sample_modules[, .(value = sum(get(cm))), by = study_canonical][order(-value)]
  study[, share := value / sum(value)]
  tissue <- sample_modules[, .(value = sum(get(cm))), by = tissue_v2][order(-value)]
  tissue[, share := value / sum(value)]
  data.table(
    module = cm,
    top_study = study$study_canonical[1L],
    top_study_share = study$share[1L],
    effective_n_studies = 1 / sum(study$share^2),
    n_studies_gt_1pct = sum(study$share > 0.01),
    top_tissue = tissue$tissue_v2[1L],
    top_tissue_share = tissue$share[1L],
    effective_n_tissues = 1 / sum(tissue$share^2),
    n_tissues_gt_1pct = sum(tissue$share > 0.01)
  )
}))
fwrite(concentration, file.path(TABLE_ROOT, "04_module_study_tissue_concentration.tsv"), sep = "\t")


# =============================================================================
# 6. Within-study top-node covariance meta-analysis
# =============================================================================

top_nodes <- rbindlist(lapply(seq_len(ncol(W_reference)), function(k) {
  ord <- order(W_reference[, k], decreasing = TRUE)
  data.table(
    module = colnames(W_reference)[k], rank = seq_len(TOP_NODES),
    cell_state = rownames(W_reference)[ord[seq_len(TOP_NODES)]],
    weight = W_reference[ord[seq_len(TOP_NODES)], k]
  )
}))
fwrite(top_nodes, file.path(TABLE_ROOT, "05_top_nodes_reference.tsv"), sep = "\t")

coherence_study_rows <- list()
pair_index <- 0L
for (cm in unique(top_nodes$module)) {
  nodes <- top_nodes[module == cm & rank <= TOP_NODES_COHERENCE, cell_state]
  pairs <- combn(nodes, 2L, simplify = FALSE)
  for (pair in pairs) {
    pair_index <- pair_index + 1L
    for (study in studies) {
      idx <- which(meta$study_canonical == study)
      if (length(idx) < MIN_WITHIN_STUDY_N) next
      x <- prop[pair[1L], idx]
      y <- prop[pair[2L], idx]
      if (sd(x) == 0 || sd(y) == 0) next
      rho <- suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
      if (!is.finite(rho)) next
      rho <- pmin(0.999999, pmax(-0.999999, rho))
      coherence_study_rows[[length(coherence_study_rows) + 1L]] <- data.table(
        module = cm, node1 = pair[1L], node2 = pair[2L],
        study_canonical = study, n_samples = length(idx), rho = rho,
        fisher_z = atanh(rho), variance = 1 / pmax(length(idx) - 3L, 1L)
      )
    }
  }
}
coherence_study <- rbindlist(coherence_study_rows, fill = TRUE)
fwrite(coherence_study, file.path(TABLE_ROOT, "06_within_study_node_correlations.tsv"), sep = "\t")

coherence_meta <- coherence_study[, {
  fit <- random_effect_meta(fisher_z, variance)
  list(
    n_studies = fit$k,
    pooled_fisher_z = fit$estimate,
    pooled_rho = tanh(fit$estimate),
    SE = fit$se, Z = fit$z, P = fit$p,
    tau2 = fit$tau2, I2_percent = fit$I2,
    direction_consistency = mean(rho > 0)
  )
}, by = .(module, node1, node2)]
coherence_meta[, FDR := p.adjust(P, method = "BH")]
fwrite(coherence_meta, file.path(TABLE_ROOT, "07_within_study_node_correlation_meta.tsv"), sep = "\t")

coherence_module <- coherence_meta[, .(
  n_top5_pairs_meta_analyzed = .N,
  median_pooled_rho = median(pooled_rho, na.rm = TRUE),
  minimum_pooled_rho = min(pooled_rho, na.rm = TRUE),
  positive_pair_fraction = mean(pooled_rho > 0, na.rm = TRUE),
  positive_FDR_pair_fraction = mean(pooled_rho > 0 & FDR < 0.05, na.rm = TRUE),
  median_direction_consistency = median(direction_consistency, na.rm = TRUE),
  median_I2_percent = median(I2_percent, na.rm = TRUE)
), by = module]
fwrite(coherence_module, file.path(TABLE_ROOT, "08_module_within_study_coherence_summary.tsv"), sep = "\t")


# =============================================================================
# 7. Study-specific disease effects and random-effects meta-analysis
# =============================================================================

sample_H <- t(H_reference)
sample_H <- sample_H / rowSums(sample_H)
log_H <- log(sample_H + 1e-6)
clr_H <- log_H - rowMeans(log_H)
colnames(clr_H) <- colnames(W_reference)
disease_data <- cbind(meta, as.data.table(clr_H))

disease_study_rows <- list()
for (cm in colnames(W_reference)) {
  for (study in studies) {
    z <- disease_data[
      study_canonical == study & condition_v2 %chin% c("healthy", "disease")
    ]
    nh <- sum(z$condition_v2 == "healthy")
    nd <- sum(z$condition_v2 == "disease")
    if (nh < MIN_DISEASE_GROUP_N || nd < MIN_DISEASE_GROUP_N) next

    paired <- z[, .(value = mean(get(cm))),
                by = .(covarnet_donor_key, tissue_v2, condition_v2)]
    paired <- dcast(
      paired, covarnet_donor_key + tissue_v2 ~ condition_v2,
      value.var = "value"
    )
    use_paired <- all(c("healthy", "disease") %in% names(paired)) &&
      nrow(paired[is.finite(healthy) & is.finite(disease)]) >= 3L

    if (use_paired) {
      paired <- paired[is.finite(healthy) & is.finite(disease)]
      delta <- paired$disease - paired$healthy
      effect <- mean(delta)
      variance <- var(delta) / length(delta)
      method <- "paired_donor"
      n_effective <- length(delta)
    } else {
      healthy <- z[condition_v2 == "healthy", get(cm)]
      disease <- z[condition_v2 == "disease", get(cm)]
      effect <- mean(disease) - mean(healthy)
      variance <- var(disease) / length(disease) + var(healthy) / length(healthy)
      method <- "independent_donors"
      n_effective <- nh + nd
    }
    if (!is.finite(effect) || !is.finite(variance) || variance <= 0) next
    disease_study_rows[[length(disease_study_rows) + 1L]] <- data.table(
      module = cm, study_canonical = study, method = method,
      n_healthy = nh, n_disease = nd, n_effective = n_effective,
      clr_mean_difference = effect, variance = variance,
      SE = sqrt(variance), Z = effect / sqrt(variance),
      P = 2 * pnorm(-abs(effect / sqrt(variance)))
    )
  }
}
disease_study <- rbindlist(disease_study_rows, fill = TRUE)
fwrite(disease_study, file.path(TABLE_ROOT, "09_disease_effects_within_study.tsv"), sep = "\t")

disease_meta <- disease_study[, {
  fit <- random_effect_meta(clr_mean_difference, variance)
  list(
    n_studies = fit$k,
    pooled_clr_difference = fit$estimate,
    SE = fit$se, Z = fit$z, P = fit$p,
    tau2 = fit$tau2, I2_percent = fit$I2,
    direction_consistency = mean(clr_mean_difference > 0),
    minimum_n_healthy = min(n_healthy),
    minimum_n_disease = min(n_disease)
  )
}, by = module]
disease_meta[, FDR := p.adjust(P, method = "BH")]
fwrite(disease_meta, file.path(TABLE_ROOT, "10_disease_effect_random_effects_meta.tsv"), sep = "\t")

disease_loso <- rbindlist(lapply(unique(disease_study$module), function(cm) {
  x <- disease_study[module == cm]
  if (nrow(x) < 3L) return(NULL)
  rbindlist(lapply(seq_len(nrow(x)), function(i) {
    fit <- random_effect_meta(x$clr_mean_difference[-i], x$variance[-i])
    data.table(
      module = cm, omitted_study = x$study_canonical[i],
      n_studies = fit$k, pooled_clr_difference = fit$estimate,
      SE = fit$se, P = fit$p, I2_percent = fit$I2
    )
  }))
}), fill = TRUE)
if (nrow(disease_loso)) {
  fwrite(disease_loso, file.path(TABLE_ROOT, "11_disease_meta_leave_one_study.tsv"), sep = "\t")
}


# =============================================================================
# 8. Optional all-module MAGMA/scDRS integration
# =============================================================================

genetic_audit <- data.table()
genetic_results <- data.table()

if (dir.exists(GENETIC_SEARCH_ROOT)) {
  magma_files <- list.files(
    GENETIC_SEARCH_ROOT, pattern = "MAGMA_Level2_with_FDR\\.tsv$",
    recursive = TRUE, full.names = TRUE
  )

  for (magma_file in magma_files) {
    run_root <- dirname(dirname(dirname(magma_file)))
    scdrs_files <- list.files(
      file.path(run_root, "02_scDRS"),
      pattern = "scdrs_group.*Level2|Level2.*scdrs_group",
      recursive = TRUE, full.names = TRUE
    )
    if (!length(scdrs_files)) next
    scdrs_file <- scdrs_files[1L]
    scenario <- basename(run_root)

    mg <- tryCatch(fread(magma_file), error = function(e) NULL)
    sd_tab <- tryCatch(fread(scdrs_file), error = function(e) NULL)
    if (is.null(mg) || is.null(sd_tab)) next

    mg_name <- intersect(c("cell_state", "FULL_NAME", "VARIABLE", "group"), names(mg))[1L]
    sd_name <- intersect(c("group", "cell_state", "FULL_NAME", "VARIABLE"), names(sd_tab))[1L]
    if (is.na(mg_name) || is.na(sd_name)) next
    mg_z_col <- if (all(c("BETA", "SE") %in% names(mg))) "BETA_over_SE" else
      intersect(c("Z", "ZSTAT", "BETA_STD"), names(mg))[1L]
    sd_z_col <- intersect(c("assoc_mcz", "assoc_z", "Z"), names(sd_tab))[1L]
    if (is.na(mg_z_col) || is.na(sd_z_col)) next

    mg[, key := normalize_key(get(mg_name))]
    sd_tab[, key := normalize_key(get(sd_name))]
    if (mg_z_col == "BETA_over_SE") {
      mg[, genetic_z := BETA / SE]
    } else {
      mg[, genetic_z := as.numeric(get(mg_z_col))]
    }
    sd_tab[, genetic_z := as.numeric(get(sd_z_col))]
    mg <- mg[is.finite(genetic_z), .(genetic_z = mean(genetic_z)), by = key]
    sd_tab <- sd_tab[is.finite(genetic_z), .(genetic_z = mean(genetic_z)), by = key]

    node_map <- data.table(
      cell_state = rownames(W_reference),
      Level1 = sub("\\s*\\|.*$", "", rownames(W_reference)),
      key = normalize_key(rownames(W_reference))
    )
    magma_node <- merge(node_map, mg[, .(key, genetic_z)], by = "key", all.x = TRUE)
    scdrs_node <- merge(node_map, sd_tab, by = "key", all.x = TRUE)

    for (method in c("MAGMA", "scDRS")) {
      node_score <- if (method == "MAGMA") magma_node else scdrs_node
      node_score <- node_score[match(rownames(W_reference), cell_state)]
      coverage <- mean(is.finite(node_score$genetic_z))
      genetic_audit <- rbind(
        genetic_audit,
        data.table(
          scenario = scenario, method = method,
          source_file = if (method == "MAGMA") magma_file else scdrs_file,
          n_nodes_matched = sum(is.finite(node_score$genetic_z)),
          n_nodes_total = nrow(node_score), coverage = coverage
        ), fill = TRUE
      )
      if (coverage < 0.70) next

      z <- node_score$genetic_z
      lineages <- node_score$Level1
      set.seed(SEED + match(method, c("MAGMA", "scDRS")) +
                 match(magma_file, magma_files) * 100L)

      for (cm in colnames(W_reference)) {
        weight <- W_reference[, cm]
        valid <- is.finite(z)
        weight_valid <- weight[valid]
        observed <- sum(weight_valid * z[valid]) / sum(weight_valid)
        null <- numeric(N_GENETIC_PERM)
        for (b in seq_len(N_GENETIC_PERM)) {
          z_perm <- z
          for (lin in unique(lineages[valid])) {
            idx <- which(valid & lineages == lin)
            if (length(idx) > 1L) z_perm[idx] <- sample(z_perm[idx])
          }
          null[b] <- sum(weight_valid * z_perm[valid]) / sum(weight_valid)
        }
        null_sd <- sd(null)
        genetic_results <- rbind(
          genetic_results,
          data.table(
            scenario = scenario, method = method, module = cm,
            matched_node_fraction = coverage,
            matched_weight_fraction = sum(weight[valid]),
            weighted_genetic_score = observed,
            permutation_z = if (null_sd > 0) (observed - mean(null)) / null_sd else NA_real_,
            permutation_P_one_sided = (1 + sum(null >= observed)) / (N_GENETIC_PERM + 1)
          ), fill = TRUE
        )
      }
    }
  }
}

if (nrow(genetic_audit)) {
  fwrite(genetic_audit, file.path(TABLE_ROOT, "12_genetic_input_mapping_audit.tsv"), sep = "\t")
}
if (nrow(genetic_results)) {
  genetic_results[, FDR := p.adjust(permutation_P_one_sided, method = "BH"),
                  by = .(scenario, method)]
  fwrite(genetic_results, file.path(TABLE_ROOT, "13_all_module_genetic_enrichment.tsv"), sep = "\t")
} else {
  writeLines(c(
    "No Level2 MAGMA/scDRS pair reached >=70% mapping coverage.",
    "The CoVarNet stability and disease analyses completed normally.",
    "Inspect 12_genetic_input_mapping_audit.tsv if it exists; the updated",
    "Level2 labels in MAGMA/scDRS must match the Level2 component after '|'."
  ), file.path(TABLE_ROOT, "13_genetic_integration_not_run.txt"))
}


# =============================================================================
# 9. Operational module tiers
# =============================================================================

adjacent_summary <- adjacent[, .(
  adjacent_rank_median_cosine = median(cosine_similarity),
  adjacent_rank_minimum_cosine = min(cosine_similarity),
  adjacent_rank_median_top10_jaccard = median(top_node_jaccard)
), by = .(module = reference_module)]

wide_stability <- dcast(
  stability_summary,
  reference_module ~ run_type,
  value.var = "q05_cosine"
)
setnames(wide_stability, "reference_module", "module")

module_tier <- Reduce(
  function(x, y) merge(x, y, by = "module", all = TRUE),
  list(adjacent_summary, wide_stability, concentration, coherence_module, disease_meta)
)

needed_stability <- c(
  "leave_one_study", "leave_one_tissue", "condition_only",
  "hierarchical_study_bootstrap", "transformation"
)
for (nm in needed_stability) {
  if (!nm %in% names(module_tier)) module_tier[, (nm) := NA_real_]
}

# These are operational reporting tiers, not inferential P-value thresholds.
# They are written to a text file so the criteria remain explicit.
module_tier[, stable_structure :=
  adjacent_rank_median_cosine >= 0.75 &
  leave_one_study >= 0.75 &
  leave_one_tissue >= 0.75 &
  condition_only >= 0.70 &
  hierarchical_study_bootstrap >= 0.65 &
  transformation >= 0.75]
module_tier[is.na(stable_structure), stable_structure := FALSE]

module_tier[, broad_context :=
  top_study_share < 0.50 & top_tissue_share < 0.70 &
  effective_n_studies >= 3 & effective_n_tissues >= 2]
module_tier[is.na(broad_context), broad_context := FALSE]

module_tier[, coherent_within_study :=
  median_pooled_rho >= 0.20 & positive_pair_fraction >= 0.70]
module_tier[is.na(coherent_within_study), coherent_within_study := FALSE]

module_tier[, module_tier := fifelse(
  stable_structure & broad_context & coherent_within_study,
  "high_confidence_pan_vascular",
  fifelse(
    stable_structure & coherent_within_study & !broad_context,
    "context_specific_reproducible",
    "candidate_or_rank_sensitive"
  )
)]
setorder(module_tier, module)
fwrite(module_tier, file.path(TABLE_ROOT, "14_final_module_evidence_and_tiers.tsv"), sep = "\t")

writeLines(c(
  "Operational module-tier criteria (predefined in this follow-up script)",
  "Stable structure: adjacent-rank median cosine >=0.75; leave-one-study q05 >=0.75; leave-one-tissue q05 >=0.75; condition-only q05 >=0.70; hierarchical study-bootstrap q05 >=0.65; transformation q05 >=0.75.",
  "Broad context: top-study share <0.50; top-tissue share <0.70; effective studies >=3; effective tissues >=2.",
  "Within-study coherence: median pooled top-five-node rho >=0.20 and >=70% of top-node pairs have positive pooled rho.",
  "These thresholds classify evidence strength; they do not replace biological interpretation."
), file.path(TABLE_ROOT, "14_module_tier_criteria.txt"))


# =============================================================================
# 10. Publication-oriented figures
# =============================================================================

adjacent_plot_data <- copy(adjacent)
adjacent_plot_data[, rank := factor(run_id, levels = paste0("K", ADJACENT_RANKS))]
p_adjacent <- ggplot(
  adjacent_plot_data,
  aes(x = rank, y = factor(reference_module, levels = rev(colnames(W_reference))),
      fill = cosine_similarity)
) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", cosine_similarity)), size = 2.6) +
  scale_fill_gradientn(colours = c("#F7FBFF", "#6BAED6", "#08306B"), limits = c(0, 1)) +
  theme_classic(base_size = 9) +
  theme(axis.title = element_blank(), axis.text = element_text(colour = "black")) +
  labs(fill = "Cosine", title = "Adjacent-rank preservation of K=14 modules")
save_plot(p_adjacent, "FigS1_adjacent_rank_one_to_one", 6.8, 5.8)

stability_plot_data <- stability_summary[
  run_type %chin% c(
    "leave_one_study", "leave_one_tissue", "condition_only",
    "sample_bootstrap", "hierarchical_study_bootstrap",
    "transformation", "cell_depth_downsample"
  )
]
stability_plot_data[, run_type := factor(
  run_type,
  levels = c(
    "transformation", "cell_depth_downsample", "condition_only",
    "leave_one_tissue", "leave_one_study", "sample_bootstrap",
    "hierarchical_study_bootstrap"
  )
)]
p_stability <- ggplot(
  stability_plot_data,
  aes(x = run_type, y = factor(reference_module, levels = rev(colnames(W_reference))),
      fill = q05_cosine)
) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", q05_cosine)), size = 2.2) +
  scale_fill_gradientn(colours = c("#FFF5F0", "#FB6A4A", "#67000D"), limits = c(0, 1)) +
  theme_classic(base_size = 8) +
  theme(
    axis.title = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1),
    axis.text = element_text(colour = "black")
  ) +
  labs(fill = "q05 cosine", title = "Strict one-to-one module stability")
save_plot(p_stability, "FigS2_strict_module_stability", 10.5, 6.2)

context_long <- melt(
  concentration,
  id.vars = "module",
  measure.vars = c("top_study_share", "top_tissue_share"),
  variable.name = "context", value.name = "share"
)
p_context <- ggplot(
  context_long,
  aes(x = share, y = factor(module, levels = rev(colnames(W_reference))), colour = context)
) +
  geom_vline(xintercept = c(0.5, 0.7), linetype = 3, colour = "grey60") +
  geom_point(size = 2.5) +
  scale_colour_manual(
    values = c(top_study_share = "#0072B2", top_tissue_share = "#D55E00"),
    labels = c(top_study_share = "Top study", top_tissue_share = "Top tissue")
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_classic(base_size = 9) +
  theme(axis.title.y = element_blank(), axis.text = element_text(colour = "black")) +
  labs(x = "Fraction of total module abundance", colour = NULL,
       title = "Study and tissue concentration")
save_plot(p_context, "FigS3_module_context_concentration", 7.2, 5.7)

p_coherence <- ggplot(
  coherence_module,
  aes(x = median_pooled_rho,
      y = factor(module, levels = rev(colnames(W_reference))),
      size = positive_FDR_pair_fraction, colour = median_I2_percent)
) +
  geom_vline(xintercept = 0.2, linetype = 3, colour = "grey60") +
  geom_point(alpha = 0.9) +
  scale_colour_gradient(low = "#009E73", high = "#CC79A7") +
  scale_size_continuous(range = c(2, 6), limits = c(0, 1)) +
  theme_classic(base_size = 9) +
  theme(axis.title.y = element_blank(), axis.text = element_text(colour = "black")) +
  labs(x = "Median meta-analytic rho among top-five node pairs",
       size = "Positive pairs\nFDR < 0.05", colour = "Median I2 (%)",
       title = "Within-study module coherence")
save_plot(p_coherence, "FigS4_within_study_module_coherence", 8.2, 5.7)

if (nrow(disease_meta)) {
  disease_meta[, lower := pooled_clr_difference - 1.96 * SE]
  disease_meta[, upper := pooled_clr_difference + 1.96 * SE]
  disease_meta[, module_factor := factor(module, levels = rev(colnames(W_reference)))]
  p_disease <- ggplot(
    disease_meta,
    aes(x = pooled_clr_difference, y = module_factor, colour = I2_percent)
  ) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_segment(aes(x = lower, xend = upper, yend = module_factor), linewidth = 0.5) +
    geom_point(size = 2.3) +
    scale_colour_gradient(low = "#0072B2", high = "#D55E00") +
    theme_classic(base_size = 9) +
    theme(axis.title.y = element_blank(), axis.text = element_text(colour = "black")) +
    labs(x = "Disease minus healthy module abundance (CLR)", colour = "I2 (%)",
         title = "Within-study random-effects disease meta-analysis")
  save_plot(p_disease, "FigS5_module_disease_meta_analysis", 8.2, 5.7)
}

if (nrow(genetic_results)) {
  p_genetic <- ggplot(
    genetic_results,
    aes(x = permutation_z,
        y = factor(module, levels = rev(colnames(W_reference))),
        colour = method, shape = FDR < 0.05)
  ) +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_point(size = 2.4) +
    facet_wrap(~ scenario, scales = "free_x") +
    scale_colour_manual(values = c(MAGMA = "#D55E00", scDRS = "#0072B2")) +
    theme_classic(base_size = 8) +
    theme(axis.title.y = element_blank(), axis.text = element_text(colour = "black")) +
    labs(x = "Lineage-stratified permutation Z", colour = NULL, shape = "FDR < 0.05",
         title = "Unbiased genetic enrichment across all CoVarNet modules")
  save_plot(p_genetic, "FigS6_all_module_genetic_enrichment", 11, 6.2)
}


# =============================================================================
# 11. Manifest and completion marker
# =============================================================================

manifest <- data.table(
  item = c(
    "primary_root", "followup_root", "K_primary", "adjacent_ranks",
    "n_nodes", "n_samples", "n_studies", "n_tissues",
    "sample_bootstrap_replicates", "study_bootstrap_replicates",
    "downsample_replicates", "downsample_depth", "genetic_permutations"
  ),
  value = c(
    PRIMARY_ROOT, FOLLOWUP_ROOT, K_PRIMARY, paste(ADJACENT_RANKS, collapse = ","),
    nrow(primary_matrix), ncol(primary_matrix), uniqueN(meta$study_canonical),
    uniqueN(meta$tissue_v2), N_SAMPLE_BOOTSTRAP, N_STUDY_BOOTSTRAP,
    N_DOWNSAMPLE, DOWNSAMPLE_DEPTH, N_GENETIC_PERM
  )
)
fwrite(manifest, file.path(FOLLOWUP_ROOT, "run_manifest.tsv"), sep = "\t")

writeLines(c(
  paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("Primary K:", K_PRIMARY),
  paste("Output:", FOLLOWUP_ROOT),
  "Review tables/14_final_module_evidence_and_tiers.tsv first."
), file.path(FOLLOWUP_ROOT, "ANALYSIS_COMPLETED.txt"))

cat("\n============================================================\n")
cat("CoVarNet follow-up validation completed.\n")
cat("Review:", file.path(TABLE_ROOT, "14_final_module_evidence_and_tiers.tsv"), "\n")
cat("Figures:", FIGURE_ROOT, "\n")
cat("============================================================\n")
