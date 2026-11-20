####pseudobulk运行
file.exists("/home/anzhen_vhd/更新版MAGMA_scDRS/final_atlas_MAGMA_scDRS_v1/GWAS_gene_priority_master.tsv")

requireNamespace("edgeR", quietly = TRUE)
requireNamespace("metafor", quietly = TRUE)

rm(list = ls())

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(edgeR)
})


# =============================================================================
# 0. Configuration
# =============================================================================

ROOT <- "/home/anzhen_vhd/scRNAseq/codex/pseudobulk_target_discovery_v1"
GWAS_MASTER <- "/home/anzhen_vhd/更新版MAGMA_scDRS/final_atlas_MAGMA_scDRS_v1/GWAS_gene_priority_master.tsv"
LEVEL2_STABLE_MAP <- "/home/anzhen_vhd/更新版MAGMA_scDRS/final_atlas_MAGMA_scDRS_v1/00_common_input/MAGMA_specificity/Level2_variable_mapping.tsv"
OUT <- file.path(ROOT, "edgeR_meta_results_v4")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

MIN_CASE_SAMPLES <- 5L
MIN_CONTROL_SAMPLES <- 5L
MIN_DISEASE_CONTEXTS_FOR_META <- 2L
MIN_DISEASE_CONTEXTS_FOR_SHARED <- 2L
MIN_DISEASE_CONTEXTS_FOR_STRICT_SHARED <- 3L
MIN_DIRECTION_CONSISTENCY <- 0.75
MIN_COVARIATE_COMPLETENESS <- 0.90
RUN_LEAVE_ONE_STUDY <- TRUE
MIN_STUDIES_FOR_LOS <- 2L
detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (is.na(detected_cores) || detected_cores < 2L) detected_cores <- 2L
META_CORES <- max(1L, min(4L, detected_cores - 1L))

# These two contexts directly represent disease domains included in the
# Genomic SEM input. PAH is retained as an extended vascular-disease context,
# but it is not counted as a directly GWAS-aligned replication context.
GWAS_ALIGNED_CONTEXTS <- c("Aortic_dissection", "Atherosclerosis")

# Leave-one-study models are refitted genome-wide, but only GWAS-prioritised
# gene rows are retained. Keeping every gene from every omitted-study fit would
# create an unnecessarily large intermediate table without strengthening the
# target-prioritisation claim.
if (!file.exists(GWAS_MASTER)) {
  stop("GWAS gene-priority master table not found: ", GWAS_MASTER,
       "\nComplete 04_MAGMA_scDRS_final_atlas_Jupyter first.")
}
if (!file.exists(LEVEL2_STABLE_MAP)) {
  stop("Stable Level2 mapping table not found: ", LEVEL2_STABLE_MAP)
}

level2_stable_map <- fread(LEVEL2_STABLE_MAP)
if (!"group" %in% names(level2_stable_map)) {
  stop("Level2 mapping table must contain a 'group' column")
}
STABLE_LEVEL2_STATES <- unique(trimws(as.character(level2_stable_map$group)))
STABLE_LEVEL2_STATES <- STABLE_LEVEL2_STATES[!is.na(STABLE_LEVEL2_STATES) & STABLE_LEVEL2_STATES != ""]
if (length(STABLE_LEVEL2_STATES) != 71L) {
  warning("Expected 71 stable Level2 states, found ", length(STABLE_LEVEL2_STATES))
}

gwas_for_los <- fread(GWAS_MASTER)
required_gwas_columns <- c(
  "SYMBOL", "primary_FDR", "primary_Bonferroni",
  "retained_leave_one_CAD_nominal"
)
missing_gwas_columns <- setdiff(required_gwas_columns, names(gwas_for_los))
if (length(missing_gwas_columns)) {
  stop("GWAS master table is missing: ", paste(missing_gwas_columns, collapse = ", "))
}
gwas_for_los[, SYMBOL := toupper(trimws(as.character(SYMBOL)))]
GWAS_LOS_SYMBOLS <- unique(gwas_for_los[
  primary_FDR == TRUE | primary_Bonferroni == TRUE, SYMBOL
])
if (!length(GWAS_LOS_SYMBOLS)) stop("No primary MAGMA-prioritised genes found")

message(
  "GWAS master: ", nrow(gwas_for_los), " genes; ",
  sum(gwas_for_los$primary_Bonferroni %in% TRUE, na.rm = TRUE),
  " primary Bonferroni genes; ",
  sum(gwas_for_los$primary_FDR %in% TRUE, na.rm = TRUE),
  " primary FDR genes"
)
message("Stable Level2 states retained for formal analysis: ", length(STABLE_LEVEL2_STATES))


bh <- function(p) p.adjust(p, method = "BH")

normalise_condition <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[x %in% c("control", "normal", "non-disease", "nondisease")] <- "healthy"
  x[x %in% c("case", "affected")] <- "disease"
  x
}


# =============================================================================
# 1. Read Matrix Market pseudobulk exports
# =============================================================================

read_level <- function(level) {
  folder <- file.path(ROOT, level)
  message("Reading ", level)
  con <- gzfile(file.path(folder, "counts_genes_by_pseudobulk.mtx.gz"), "rt")
  on.exit(close(con), add = TRUE)
  counts <- readMM(con)
  counts <- as(counts, "dgCMatrix")
  features <- fread(file.path(folder, "features.tsv"))
  units <- fread(file.path(folder, "pseudobulk_unit_metadata.tsv"))
  stopifnot(nrow(counts) == nrow(features), ncol(counts) == nrow(units))
  stopifnot(!anyDuplicated(units$pseudobulk_id))

  if (level == "Level2") {
    keep_columns <- which(as.character(units$cell_state) %chin% STABLE_LEVEL2_STATES)
    excluded_states <- sort(setdiff(
      unique(as.character(units$cell_state)), STABLE_LEVEL2_STATES
    ))
    if (!length(keep_columns)) stop("No stable Level2 pseudobulk units matched the mapping table")
    if (length(excluded_states)) {
      message(
        "Excluding non-stable Level2 states from formal pseudobulk analysis: ",
        paste(excluded_states, collapse = "; ")
      )
    }
    counts <- counts[, keep_columns, drop = FALSE]
    units <- units[keep_columns]
  }

  rownames(counts) <- features$gene_id
  colnames(counts) <- units$pseudobulk_id
  units[, condition := normalise_condition(condition_v2)]
  units[, disease_context := trimws(as.character(disease_context))]
  list(level = level, counts = counts, features = features, units = units)
}


# =============================================================================
# 1b. Input audit before model fitting
# =============================================================================

input_audit <- rbindlist(lapply(c("Level1", "Level2"), function(level) {
  units <- fread(file.path(ROOT, level, "pseudobulk_unit_metadata.tsv"))
  if (level == "Level2") {
    units <- units[as.character(cell_state) %chin% STABLE_LEVEL2_STATES]
  }
  units[, condition := normalise_condition(condition_v2)]
  units[, disease_context := trimws(as.character(disease_context))]
  units[, annotation_level := rep(level, .N)]
  units[, .(
    n_pseudobulk_units = .N,
    n_case = sum(condition == "disease"),
    n_control = sum(condition == "healthy"),
    n_studies = uniqueN(study_canonical)
  ), by = c("annotation_level", "cell_state", "disease_context")]
}), fill = TRUE)
fwrite(input_audit, file.path(OUT, "00_pseudobulk_input_audit.tsv"), sep = "\t")
fwrite(
  data.table(Level2_cell_state = sort(STABLE_LEVEL2_STATES)),
  file.path(OUT, "00_stable_Level2_states_used.tsv"),
  sep = "\t"
)

present_contexts <- unique(input_audit$disease_context)
missing_aligned_contexts <- setdiff(GWAS_ALIGNED_CONTEXTS, present_contexts)
if (length(missing_aligned_contexts)) {
  stop(
    "GWAS-aligned disease contexts absent from pseudobulk input: ",
    paste(missing_aligned_contexts, collapse = ", ")
  )
}


# =============================================================================
# 2. Robust design construction and edgeR quasi-likelihood test
# =============================================================================

prepare_design <- function(meta) {
  meta <- copy(meta)
  meta[, condition := factor(condition, levels = c("healthy", "disease"))]
  meta[, study_canonical := factor(as.character(study_canonical))]
  meta[, sex_v2 := trimws(as.character(sex_v2))]
  meta[sex_v2 %chin% c("", "NA", "NaN", "nan", "None"), sex_v2 := NA_character_]
  meta[, sex_v2 := factor(sex_v2)]
  meta[, age_v2 := suppressWarnings(as.numeric(age_v2))]
  meta[, donor_tissue_key := paste(
    as.character(study_canonical), as.character(donor_uid_v2),
    as.character(tissue_v2), sep = "::"
  )]
  paired_donor_tissues <- meta[
    !is.na(condition), .(n_conditions = uniqueN(condition)), by = donor_tissue_key
  ][
    n_conditions > 1L, donor_tissue_key
  ]
  # A single stratification term avoids collinearity between study and donor:
  # paired observations receive a donor-specific block, whereas independent
  # observations share a study-specific block. The latter controls study-level
  # shifts while retaining case-control contrasts within that study.
  meta[, analysis_block := factor(fifelse(
    donor_tissue_key %in% paired_donor_tissues,
    donor_tissue_key,
    paste0("UNPAIRED_STUDY::", as.character(study_canonical))
  ))]
  
  keep <- !is.na(meta$condition)
  meta <- meta[keep]
  if (sum(meta$condition == "healthy") < MIN_CONTROL_SAMPLES ||
      sum(meta$condition == "disease") < MIN_CASE_SAMPLES) return(NULL)
  
  accepted <- character()
  # The study/donor stratification is the primary nuisance term. Sex and age are
  # added only when at least 90% complete and identifiable, preventing an entire
  # study being removed merely because demographic metadata are absent.
  for (term in c("analysis_block", "sex_v2", "age_v2")) {
    candidate <- c(accepted, term)
    cc <- complete.cases(meta[, c("condition", candidate), with = FALSE])
    m <- meta[cc]
    if (term %in% c("sex_v2", "age_v2") &&
        nrow(m) / nrow(meta) < MIN_COVARIATE_COMPLETENESS) next
    if (sum(m$condition == "healthy") < MIN_CONTROL_SAMPLES ||
        sum(m$condition == "disease") < MIN_CASE_SAMPLES) next
    if (term == "analysis_block" && uniqueN(m$analysis_block) < 2L) next
    if (term == "sex_v2" && uniqueN(m$sex_v2) < 2L) next
    if (term == "age_v2" && (uniqueN(m$age_v2) < 5L || sd(m$age_v2) == 0)) next
    f <- as.formula(paste("~", paste(c(candidate, "condition"), collapse = " + ")))
    mm <- model.matrix(f, data = m)
    if (qr(mm)$rank == ncol(mm) && "conditiondisease" %in% colnames(mm)) {
      accepted <- candidate
      meta <- m
    }
  }
  
  f <- as.formula(paste("~", paste(c(accepted, "condition"), collapse = " + ")))
  design <- model.matrix(f, data = meta)
  if (qr(design)$rank != ncol(design) || !"conditiondisease" %in% colnames(design)) {
    design <- model.matrix(~ condition, data = meta)
    accepted <- character()
  }
  list(
    meta = meta,
    design = design,
    covariates = paste(accepted, collapse = "+"),
    n_paired_donor_tissues = length(paired_donor_tissues),
    design_formula = paste("~", paste(c(accepted, "condition"), collapse = " + "))
  )
}


run_one_de <- function(counts, features, meta, level, state, context,
                       omitted_study = NA_character_, keep_symbols = NULL) {
  if (!is.na(omitted_study)) meta <- meta[study_canonical != omitted_study]
  prepared <- prepare_design(meta)
  if (is.null(prepared)) return(NULL)
  meta <- prepared$meta
  design <- prepared$design
  x <- counts[, meta$pseudobulk_id, drop = FALSE]
  
  gene_anno <- data.frame(
    gene_symbol = features$gene_symbol,
    row.names = features$gene_id,
    stringsAsFactors = FALSE
  )
  y <- DGEList(counts = x, genes = gene_anno)
  keep_gene <- filterByExpr(y, design = design, min.count = 10)
  if (sum(keep_gene) < 100L) return(NULL)
  y <- y[keep_gene, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  test <- glmQLFTest(fit, coef = which(colnames(design) == "conditiondisease"))
  tab <- as.data.table(topTags(test, n = Inf, sort.by = "none")$table, keep.rownames = "gene_id")
  tab[, FDR := bh(PValue)]
  tab[, `:=`(
    annotation_level = level,
    cell_state = state,
    disease_context = context,
    n_case = sum(meta$condition == "disease"),
    n_control = sum(meta$condition == "healthy"),
    n_studies = uniqueN(meta$study_canonical),
    n_paired_donor_tissues = prepared$n_paired_donor_tissues,
    covariates = prepared$covariates,
    design_formula = prepared$design_formula,
    omitted_study = omitted_study
  )]
  if (!is.null(keep_symbols)) {
    tab <- tab[toupper(as.character(gene_symbol)) %in% keep_symbols]
  }
  tab
}


# =============================================================================
# 3. Disease-context DE and leave-one-study sensitivity
# =============================================================================

full_results <- list()
los_results <- list()

for (level in c("Level1", "Level2")) {
  dat <- read_level(level)
  contexts <- sort(unique(dat$units$disease_context))
  contexts <- contexts[!is.na(contexts) & contexts != ""]
  states <- sort(unique(dat$units$cell_state))
  
  for (state in states) {
    state_meta <- dat$units[cell_state == state]
    for (context in contexts) {
      m <- state_meta[disease_context == context & condition %in% c("healthy", "disease")]
      if (nrow(m) == 0L) next
      # Primary inference is restricted to studies that contain both disease
      # and healthy donors for this cell state. This prevents a case-only study
      # being compared with an unrelated control-only atlas.
      paired_studies <- m[, .(
        n_case_in_study = sum(condition == "disease"),
        n_control_in_study = sum(condition == "healthy")
      ), by = study_canonical][
        n_case_in_study >= 2L & n_control_in_study >= 2L,
        study_canonical
      ]
      m <- m[study_canonical %in% paired_studies]
      if (nrow(m) == 0L) next
      key <- paste(level, state, context, sep = "__")
      message("DE: ", key)
      result <- tryCatch(
        run_one_de(dat$counts, dat$features, m, level, state, context),
        error = function(e) {
          message("  skipped: ", conditionMessage(e)); NULL
        }
      )
      if (is.null(result)) next
      full_results[[key]] <- result
      
      if (RUN_LEAVE_ONE_STUDY && length(GWAS_LOS_SYMBOLS) > 0L &&
          uniqueN(m$study_canonical) >= MIN_STUDIES_FOR_LOS) {
        for (study in unique(m$study_canonical)) {
          los <- tryCatch(
            run_one_de(
              dat$counts, dat$features, m, level, state, context, study,
              keep_symbols = GWAS_LOS_SYMBOLS
            ),
            error = function(e) NULL
          )
          if (!is.null(los)) {
            los_results[[paste(key, study, sep = "__")]] <- los
          }
        }
      }
    }
  }
}

if (length(full_results) == 0L) stop("No disease-context comparison passed sample thresholds")
de <- rbindlist(full_results, fill = TRUE)
rm(full_results)
if (exists("dat")) rm(dat)
gc()
fwrite(de, file.path(OUT, "all_disease_context_edgeR_results.tsv.gz"), sep = "\t")

design_audit <- unique(de[, .(
  annotation_level, cell_state, disease_context,
  n_case, n_control, n_studies, n_paired_donor_tissues,
  covariates, design_formula
)])
fwrite(design_audit, file.path(OUT, "edgeR_design_and_sample_audit.tsv"), sep = "\t")

# An exhaustive, model-free cross-context direction summary is retained for all
# tested genes. Random-effects fitting below is intentionally restricted to
# MAGMA-prioritised genes to keep the target-discovery analysis tractable and to
# avoid millions of unnecessary per-gene model fits.
de[, gene_symbol_match := toupper(trimws(as.character(gene_symbol)))]
all_gene_direction_summary <- de[
  is.finite(logFC),
  .(
    n_disease_contexts = uniqueN(disease_context),
    disease_contexts = paste(sort(unique(disease_context)), collapse = ";"),
    positive_fraction = mean(logFC > 0),
    direction_consistency = max(mean(logFC > 0), mean(logFC < 0)),
    median_context_logFC = median(logFC),
    n_context_FDR_lt_0.05 = sum(FDR < 0.05),
    minimum_context_FDR = min(FDR, na.rm = TRUE)
  ),
  by = .(annotation_level, cell_state, gene_id, gene_symbol, gene_symbol_match)
]
fwrite(
  all_gene_direction_summary,
  file.path(OUT, "all_genes_cross_context_direction_summary.tsv.gz"),
  sep = "\t"
)

if (length(los_results)) {
  los <- rbindlist(los_results, fill = TRUE)
  rm(los_results)
  gc()
  fwrite(los, file.path(OUT, "leave_one_study_edgeR_results.tsv.gz"), sep = "\t")
  los_summary <- los[, .(
    n_leave_one_runs = .N,
    positive_fraction_LOS = mean(logFC > 0, na.rm = TRUE),
    median_logFC_LOS = median(logFC, na.rm = TRUE),
    max_PValue_LOS = max(PValue, na.rm = TRUE)
  ), by = .(annotation_level, cell_state, disease_context, gene_id)]
  fwrite(los_summary, file.path(OUT, "leave_one_study_stability_summary.tsv.gz"), sep = "\t")
} else {
  los_summary <- data.table()
}


# =============================================================================
# 4. Cross-disease random-effects meta-analysis within each cell state
# =============================================================================

if (!requireNamespace("metafor", quietly = TRUE)) {
  stop("Install metafor in this R environment: install.packages('metafor')")
}

meta_one <- function(x) {
  x <- x[is.finite(logFC) & is.finite(F) & F > 0 & is.finite(PValue)]
  if (uniqueN(x$disease_context) < MIN_DISEASE_CONTEXTS_FOR_META) return(NULL)
  # For a one-degree-of-freedom QL F test, |logFC|/sqrt(F) is used as an
  # approximate standard error for cross-context effect-size synthesis.
  x[, se_approx := abs(logFC) / sqrt(F)]
  x <- x[is.finite(se_approx) & se_approx > 0]
  if (nrow(x) < 2L) return(NULL)
  fit <- tryCatch(
    metafor::rma.uni(yi = x$logFC, sei = x$se_approx, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  pos <- mean(x$logFC > 0)
  neg <- mean(x$logFC < 0)
  data.table(
    n_disease_contexts = uniqueN(x$disease_context),
    n_GWAS_aligned_contexts = uniqueN(
      x[disease_context %chin% GWAS_ALIGNED_CONTEXTS, disease_context]
    ),
    disease_contexts = paste(sort(unique(x$disease_context)), collapse = ";"),
    meta_logFC = as.numeric(fit$b),
    meta_SE = fit$se,
    meta_Z = fit$zval,
    meta_P = fit$pval,
    meta_CI_low = fit$ci.lb,
    meta_CI_high = fit$ci.ub,
    I2 = fit$I2,
    direction_consistency = max(pos, neg),
    consensus_direction = ifelse(pos >= neg, "up_in_disease", "down_in_disease"),
    n_context_FDR_lt_0.05 = sum(x$FDR < 0.05),
    median_context_logFC = median(x$logFC)
  )
}

de_priority <- de[gene_symbol_match %chin% GWAS_LOS_SYMBOLS]
eligible_meta_keys <- de_priority[, .(
  n_disease_contexts = uniqueN(disease_context)
), by = .(annotation_level, cell_state, gene_id)][
  n_disease_contexts >= MIN_DISEASE_CONTEXTS_FOR_META
]
de_priority <- de_priority[
  eligible_meta_keys,
  on = .(annotation_level, cell_state, gene_id),
  nomatch = 0L
]
if (!nrow(de_priority)) stop("No MAGMA-prioritised gene was testable in >=2 disease contexts")

meta_split <- split(
  de_priority,
  by = c("annotation_level", "cell_state", "gene_id"),
  keep.by = TRUE
)
meta_worker <- function(x) {
  ans <- meta_one(x)
  if (is.null(ans)) return(NULL)
  ans[, `:=`(
    annotation_level = x$annotation_level[1],
    cell_state = x$cell_state[1],
    gene_id = x$gene_id[1]
  )]
  ans
}
if (.Platform$OS.type == "unix" && META_CORES > 1L) {
  meta_list <- parallel::mclapply(meta_split, meta_worker, mc.cores = META_CORES)
} else {
  meta_list <- lapply(meta_split, meta_worker)
}
meta <- rbindlist(meta_list, fill = TRUE)
if (!nrow(meta)) stop("No cross-disease MAGMA-prioritised gene meta-analysis succeeded")
meta[, target_set_meta_FDR := bh(meta_P), by = .(annotation_level, cell_state)]

feature_map <- unique(de[, .(gene_id, gene_symbol)])
meta <- merge(meta, feature_map, by = "gene_id", all.x = TRUE)
setcolorder(meta, c(
  "annotation_level", "cell_state", "gene_id", "gene_symbol",
  setdiff(names(meta), c("annotation_level", "cell_state", "gene_id", "gene_symbol"))
))
fwrite(meta, file.path(OUT, "GWAS_prioritised_cross_disease_meta_analysis.tsv.gz"), sep = "\t")


# =============================================================================
# 5. Integrate disease-responsive genes with genomic prioritisation
# =============================================================================

gwas <- copy(gwas_for_los)
meta[, gene_symbol_match := toupper(trimws(as.character(gene_symbol)))]
candidates <- merge(meta, gwas, by.x = "gene_symbol_match", by.y = "SYMBOL", all.x = TRUE)

candidates[, strict_three_context_shared :=
             n_disease_contexts >= MIN_DISEASE_CONTEXTS_FOR_STRICT_SHARED &
             direction_consistency >= MIN_DIRECTION_CONSISTENCY]
candidates[, Tier := "GWAS_prioritised_without_shared_DE_support"]
candidates[
  primary_FDR == TRUE & target_set_meta_FDR < 0.05 &
    n_disease_contexts >= MIN_DISEASE_CONTEXTS_FOR_SHARED &
    direction_consistency >= MIN_DIRECTION_CONSISTENCY,
  Tier := "Tier2_GWAS_FDR_plus_multidisease_DE"
]
candidates[
  primary_Bonferroni == TRUE & target_set_meta_FDR < 0.05 &
    n_GWAS_aligned_contexts >= 2L &
    direction_consistency >= MIN_DIRECTION_CONSISTENCY &
    retained_leave_one_CAD_nominal == TRUE,
  Tier := "Tier1_GWAS_Bonferroni_GWAS-aligned_DE_LOT-CAD"
]
candidates[, tier_order := match(Tier, c(
  "Tier1_GWAS_Bonferroni_GWAS-aligned_DE_LOT-CAD",
  "Tier2_GWAS_FDR_plus_multidisease_DE",
  "GWAS_prioritised_without_shared_DE_support"
))]
setorder(candidates, tier_order, target_set_meta_FDR, meta_P)
fwrite(candidates, file.path(OUT, "GWAS_pseudobulk_integrated_candidates.tsv.gz"), sep = "\t")
fwrite(
  candidates[Tier %chin% c(
    "Tier1_GWAS_Bonferroni_GWAS-aligned_DE_LOT-CAD",
    "Tier2_GWAS_FDR_plus_multidisease_DE"
  )],
  file.path(OUT, "prioritised_target_candidates.tsv"), sep = "\t"
)

analysis_parameters <- data.table(
  parameter = c(
    "MIN_CASE_SAMPLES", "MIN_CONTROL_SAMPLES",
    "MIN_DISEASE_CONTEXTS_FOR_META", "MIN_DISEASE_CONTEXTS_FOR_SHARED",
    "MIN_DISEASE_CONTEXTS_FOR_STRICT_SHARED", "MIN_DIRECTION_CONSISTENCY",
    "MIN_COVARIATE_COMPLETENESS", "RUN_LEAVE_ONE_STUDY",
    "MIN_STUDIES_FOR_LOS", "META_CORES", "GWAS_ALIGNED_CONTEXTS"
  ),
  value = c(
    MIN_CASE_SAMPLES, MIN_CONTROL_SAMPLES,
    MIN_DISEASE_CONTEXTS_FOR_META, MIN_DISEASE_CONTEXTS_FOR_SHARED,
    MIN_DISEASE_CONTEXTS_FOR_STRICT_SHARED, MIN_DIRECTION_CONSISTENCY,
    MIN_COVARIATE_COMPLETENESS, RUN_LEAVE_ONE_STUDY,
    MIN_STUDIES_FOR_LOS, META_CORES, paste(GWAS_ALIGNED_CONTEXTS, collapse = ";")
  )
)
fwrite(analysis_parameters, file.path(OUT, "analysis_parameters.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))

cat("\nCompleted pseudobulk target-discovery analysis.\n")
cat("Results:", OUT, "\n")
