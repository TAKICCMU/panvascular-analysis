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

args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
ROOT <- normalizePath(args[1])
PRIMARY_ROOT <- file.path(ROOT, "stable_NMF")
FOLLOWUP_ROOT <- args[2]
CACHE_ROOT <- file.path(FOLLOWUP_ROOT, "cache_nmf_fits")
TABLE_ROOT <- file.path(FOLLOWUP_ROOT, "tables")
FIGURE_ROOT <- file.path(FOLLOWUP_ROOT, "figures")
GENETIC_SEARCH_ROOT <- "/path/to/magma_scdrs"

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

workers <- 1L
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
  binding <- list(input=m, rank=rank, nrun=nrun, seed=seed, method="nsNMF", adapter=tools::md5sum(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
  receipt <- paste0(path,".binding.rds")
  if (file.exists(path)) {
    if(!file.exists(receipt) || !identical(readRDS(receipt),binding)) stop("Content binding mismatch: ",path)
    return(readRDS(path))
  }
  cat("Fit:", label, "K=", rank, "samples=", ncol(m), "nrun=", nrun, "\n")
  set.seed(seed)
  fit <- nmf(
    m, rank = rank, method = "nsNMF", nrun = nrun,
    seed = seed, .options = nmf_options
  )
  saveRDS(fit, path)
  saveRDS(binding, receipt)
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
sidecar <- fread(Sys.getenv("COVARNET_SIDECAR"))
stopifnot(!anyDuplicated(sidecar$analysis_profile_uid))
mi <- match(meta$covarnet_sample_uid,sidecar$analysis_profile_uid);stopifnot(!anyNA(mi))
meta[,condition_v2:=sidecar$group_binary[mi]]
meta[,covarnet_donor_key:=sidecar$donor_block[mi]]
stopifnot(!any(meta$study_canonical=="GSE220774"), !any(grepl("GSM7025944|GSM7025945",meta$covarnet_sample_uid)))

stopifnot(!anyDuplicated(meta$covarnet_sample_uid),setequal(colnames(prop),meta$covarnet_sample_uid),identical(dimnames(prop),dimnames(counts)),identical(dimnames(prop),dimnames(prop_closed)))
common_samples <- Reduce(intersect, list(
  colnames(prop), colnames(prop_closed), colnames(counts), meta$covarnet_sample_uid
))
prop <- prop[, common_samples, drop = FALSE]
prop_closed <- prop_closed[rownames(prop), common_samples, drop = FALSE]
counts <- counts[rownames(prop), common_samples, drop = FALSE]
meta <- meta[match(common_samples, covarnet_sample_uid)]
stopifnot(all(meta$covarnet_sample_uid == common_samples))

primary_matrix <- read_matrix(file.path(PRIMARY_ROOT,"NMF_input_minmax.tsv"))
prop <- prop[rownames(primary_matrix),,drop=FALSE]
prop_closed <- prop_closed[rownames(primary_matrix),,drop=FALSE]
counts <- counts[rownames(primary_matrix),,drop=FALSE]
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
  if (ncol(m) < max(25L, 2L * K_PRIMARY)) {
    fwrite(data.table(type=type,id=id,n=ncol(m),status="NOT_ESTIMABLE_INSUFFICIENT_PROFILES"),file.path(TABLE_ROOT,paste0("skip_",safe_name(paste(type,id)),".tsv")),sep="\t"); return(NULL)
  }
  fit <- tryCatch(
    fit_nmf_cached(
      m, K_PRIMARY, SENSITIVITY_NRUN,
      paste(type, id, sep = "__"), SEED + seed_offset
    ),
    error = function(e) {
      stop(type, " ", id, " failed: ", conditionMessage(e))
    }
  )
  if (is.null(fit)) return(NULL)
  W <- normalize_columns(basis(fit))
  colnames(W) <- sprintf("candidate_%02d", seq_len(ncol(W)))
  compare_basis_one_to_one(W_reference, W, type, id)
}

# Alternative denominator and transformations.
transformations <- list(
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

studies <- unique(meta$study_canonical)
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

# Cluster bootstrap resamples studies, then samples within each selected study.
set.seed(SEED + 6000L)
for (b in seq_len(N_STUDY_BOOTSTRAP)) {
  drawn_studies <- sample(studies, length(studies), replace = TRUE)
  idx <- unlist(lapply(drawn_studies, function(s) {
    donors <- unique(meta[study_canonical==s,covarnet_donor_key])
    draws <- sample(donors,length(donors),replace=TRUE)
    unlist(lapply(draws,function(d) which(meta$study_canonical==s & meta$covarnet_donor_key==d)))
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
fwrite(data.table(covarnet_sample_uid=meta$covarnet_sample_uid,eligible_depth=eligible_depth,n_total=meta$n_total_atlas_cells),file.path(TABLE_ROOT,"downsample_eligibility.tsv"),sep="\t")
if(sum(eligible_depth)<max(30L,2L*K_PRIMARY)) stop("Insufficient profiles for required depth sensitivity")
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
      stopifnot(total >= sum(counts[, j]))
      other <- total - sum(counts[, j])
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
# 5. Study/tissue concentration of relative module coefficient
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



writeLines(capture.output(sessionInfo()),file.path(FOLLOWUP_ROOT,"sessionInfo.txt"))
writeLines("COMPLETED_REMAINING_REFITS",file.path(FOLLOWUP_ROOT,"REFITS_COMPLETED.txt"))
