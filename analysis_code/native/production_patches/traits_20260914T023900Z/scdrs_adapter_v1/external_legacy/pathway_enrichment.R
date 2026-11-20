#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
})

set.seed(20260822)

out_root <- "/path/to/panvas_work/final_revision_work/pathway_enrichment"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

magma_file <- "/path/to/fuma/FUMA结果/FUMA_LDr20.1Q_bifactor/magma.genes.out"
pseudobulk_file <- "/path/to/downloads/pseudobulk_target_discovery_v1/edgeR_meta_results_v4/all_disease_context_edgeR_results.tsv.gz"

if (!file.exists(magma_file)) stop("Missing MAGMA input: ", magma_file)
if (!file.exists(pseudobulk_file)) stop("Missing pseudobulk input: ", pseudobulk_file)

collections <- list(
  GO_BP = c("C5", "GO:BP"),
  GO_CC = c("C5", "GO:CC"),
  GO_MF = c("C5", "GO:MF"),
  KEGG_LEGACY = c("C2", "CP:KEGG_LEGACY"),
  KEGG_MEDICUS = c("C2", "CP:KEGG_MEDICUS"),
  REACTOME = c("C2", "CP:REACTOME")
)

load_gene_sets <- function() {
  rows <- lapply(names(collections), function(nm) {
    spec <- collections[[nm]]
    x <- as.data.table(msigdbr(species = "Homo sapiens", collection = spec[1], subcollection = spec[2]))
    x[, database := nm]
    x[, .(database, pathway = gs_name, gene = gene_symbol, db_version)]
  })
  unique(rbindlist(rows))
}

gene_sets <- load_gene_sets()
fwrite(gene_sets, file.path(out_root, "MSigDB_2026.1_Hs_gene_sets_used.tsv.gz"), sep = "\t")

gsea_one <- function(stats, analysis_id, level, state) {
  stats <- stats[is.finite(stats) & !is.na(names(stats)) & names(stats) != ""]
  stats <- sort(tapply(stats, names(stats), max), decreasing = TRUE)
  ans <- lapply(names(collections), function(db) {
    gs <- gene_sets[database == db, split(gene, pathway)]
    z <- suppressWarnings(fgseaMultilevel(
      pathways = gs,
      stats = stats,
      minSize = 10,
      maxSize = 500,
      eps = 1e-50,
      scoreType = "std"
    ))
    if (!nrow(z)) return(NULL)
    z <- as.data.table(z)
    z[, `:=`(
      database = db,
      analysis_id = analysis_id,
      annotation_level = level,
      cell_state = state,
      n_ranked_genes = length(stats)
    )]
    z[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
    z
  })
  rbindlist(ans, fill = TRUE)
}

ora_one <- function(selected, universe, analysis_id, level, state, direction) {
  selected <- unique(selected[!is.na(selected) & selected != ""])
  universe <- unique(universe[!is.na(universe) & universe != ""])
  selected <- intersect(selected, universe)
  if (length(selected) < 3L) return(NULL)
  x <- gene_sets[gene %chin% universe, .(pathway_genes = list(unique(gene))), by = .(database, pathway)]
  if (!nrow(x)) return(NULL)
  n_sel <- length(selected)
  n_bg <- length(universe)
  x[, overlap_genes := lapply(pathway_genes, intersect, y = selected)]
  x[, overlap := lengths(overlap_genes)]
  x[, pathway_size := lengths(pathway_genes)]
  x <- x[overlap >= 3L]
  if (!nrow(x)) return(NULL)
  x[, pval := phyper(overlap - 1L, pathway_size, n_bg - pathway_size, n_sel, lower.tail = FALSE)]
  x[, padj := p.adjust(pval, method = "BH"), by = database]
  x[, odds_ratio := ((overlap + 0.5) * (n_bg - pathway_size - n_sel + overlap + 0.5)) /
       ((n_sel - overlap + 0.5) * (pathway_size - overlap + 0.5))]
  x[, overlap_genes := vapply(overlap_genes, paste, collapse = ";", FUN.VALUE = character(1))]
  x[, pathway_genes := NULL]
  x[, `:=`(
    analysis_id = analysis_id,
    annotation_level = level,
    cell_state = state,
    direction = direction,
    n_selected_genes = n_sel,
    n_background_genes = n_bg
  )]
  x
}

message("Running MAGMA enrichment")
magma <- fread(magma_file)
magma <- magma[!is.na(SYMBOL) & SYMBOL != ""]
setorder(magma, P, -ZSTAT)
magma <- magma[, .SD[1], by = SYMBOL]
magma[, bonferroni := P < 0.05 / .N]
magma[, fdr := p.adjust(P, method = "BH")]
fwrite(magma, file.path(out_root, "MAGMA_primary_gene_statistics.tsv.gz"), sep = "\t")

magma_stats <- setNames(magma$ZSTAT, magma$SYMBOL)
magma_gsea <- gsea_one(magma_stats, "MAGMA_primary_r2_0.10", "gene_level", "all_genes")
magma_ora <- rbindlist(list(
  ora_one(magma[bonferroni == TRUE, SYMBOL], magma$SYMBOL, "MAGMA_primary_r2_0.10", "gene_level", "all_genes", "Bonferroni"),
  ora_one(magma[fdr < 0.05, SYMBOL], magma$SYMBOL, "MAGMA_primary_r2_0.10", "gene_level", "all_genes", "FDR_lt_0.05")
), fill = TRUE)
fwrite(magma_gsea, file.path(out_root, "MAGMA_primary_GSEA_GO_KEGG_Reactome.tsv.gz"), sep = "\t")
fwrite(magma_ora, file.path(out_root, "MAGMA_primary_ORA_GO_KEGG_Reactome.tsv.gz"), sep = "\t")

message("Reconstructing cross-context pseudobulk random-effects meta-analysis")
pb <- fread(pseudobulk_file)
pb <- pb[is.finite(logFC) & is.finite(F) & F > 0 & !is.na(gene_symbol) & gene_symbol != ""]
pb[, SE := abs(logFC) / sqrt(F)]
pb <- pb[is.finite(SE) & SE > 0]

meta_group <- function(beta, se) {
  k <- length(beta)
  if (k < 2L) return(list(k = k, beta = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_, tau2 = NA_real_, I2 = NA_real_))
  v <- se^2
  w <- 1 / v
  fixed <- sum(w * beta) / sum(w)
  Q <- sum(w * (beta - fixed)^2)
  C <- sum(w) - sum(w^2) / sum(w)
  tau2 <- if (C > 0) max(0, (Q - (k - 1)) / C) else 0
  wr <- 1 / (v + tau2)
  b <- sum(wr * beta) / sum(wr)
  s <- sqrt(1 / sum(wr))
  z <- b / s
  p <- 2 * pnorm(-abs(z))
  I2 <- if (Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0
  list(k = k, beta = b, se = s, z = z, p = p, tau2 = tau2, I2 = I2)
}

pb_meta <- pb[, meta_group(logFC, SE), by = .(annotation_level, cell_state, gene_id, gene_symbol)]
pb_meta <- pb_meta[k >= 2L & is.finite(z)]
pb_meta[, fdr := p.adjust(p, method = "BH"), by = .(annotation_level, cell_state)]
pb_meta[, direction := fifelse(beta > 0, "higher_in_disease", "lower_in_disease")]
fwrite(pb_meta, file.path(out_root, "pseudobulk_full_cross_context_random_effects_meta.tsv.gz"), sep = "\t")

message("Running cell-state-resolved pseudobulk enrichment")
keys <- unique(pb_meta[, .(annotation_level, cell_state)])
gsea_rows <- vector("list", nrow(keys))
ora_rows <- vector("list", nrow(keys))
for (i in seq_len(nrow(keys))) {
  lev <- keys$annotation_level[i]
  st <- keys$cell_state[i]
  x <- pb_meta[annotation_level == lev & cell_state == st]
  aid <- paste("pseudobulk", lev, st, sep = "::")
  gsea_rows[[i]] <- gsea_one(setNames(x$z, x$gene_symbol), aid, lev, st)
  ora_rows[[i]] <- rbindlist(list(
    ora_one(x[fdr < 0.05 & beta > 0, gene_symbol], x$gene_symbol, aid, lev, st, "higher_in_disease"),
    ora_one(x[fdr < 0.05 & beta < 0, gene_symbol], x$gene_symbol, aid, lev, st, "lower_in_disease")
  ), fill = TRUE)
  if (i %% 10L == 0L) message("  completed ", i, "/", nrow(keys), " states")
}
pb_gsea <- rbindlist(gsea_rows, fill = TRUE)
pb_ora <- rbindlist(ora_rows, fill = TRUE)
fwrite(pb_gsea, file.path(out_root, "pseudobulk_cell_state_GSEA_GO_KEGG_Reactome.tsv.gz"), sep = "\t")
fwrite(pb_ora, file.path(out_root, "pseudobulk_cell_state_ORA_GO_KEGG_Reactome.tsv.gz"), sep = "\t")

gsea_summary <- pb_gsea[padj < 0.05, .(
  n_cell_states = uniqueN(paste(annotation_level, cell_state)),
  n_positive = sum(NES > 0),
  n_negative = sum(NES < 0),
  median_NES = median(NES),
  minimum_padj = min(padj)
), by = .(database, pathway)][order(minimum_padj, -n_cell_states)]
ora_summary <- pb_ora[padj < 0.05, .(
  n_cell_states = uniqueN(paste(annotation_level, cell_state)),
  n_higher = sum(direction == "higher_in_disease"),
  n_lower = sum(direction == "lower_in_disease"),
  minimum_padj = min(padj)
), by = .(database, pathway)][order(minimum_padj, -n_cell_states)]
fwrite(gsea_summary, file.path(out_root, "pseudobulk_recurrent_GSEA_pathways.tsv"), sep = "\t")
fwrite(ora_summary, file.path(out_root, "pseudobulk_recurrent_ORA_pathways.tsv"), sep = "\t")

parameters <- data.table(
  parameter = c("analysis_date", "random_seed", "MSigDB_version", "GSEA_engine", "GSEA_min_size", "GSEA_max_size", "ORA_test", "ORA_min_overlap", "ORA_background", "pseudobulk_meta_model"),
  value = c(as.character(Sys.Date()), "20260822", unique(gene_sets$db_version)[1], "fgseaMultilevel", "10", "500", "one-sided hypergeometric", "3", "all genes tested in the corresponding analysis", "DerSimonian-Laird random-effects model")
)
fwrite(parameters, file.path(out_root, "pathway_analysis_parameters.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))

summary <- data.table(
  analysis = c("MAGMA_GSEA", "MAGMA_ORA", "pseudobulk_meta_genes", "pseudobulk_GSEA", "pseudobulk_ORA"),
  n_rows = c(nrow(magma_gsea), nrow(magma_ora), nrow(pb_meta), nrow(pb_gsea), nrow(pb_ora)),
  n_FDR_lt_0.05 = c(
    sum(magma_gsea$padj < 0.05, na.rm = TRUE),
    sum(magma_ora$padj < 0.05, na.rm = TRUE),
    sum(pb_meta$fdr < 0.05, na.rm = TRUE),
    sum(pb_gsea$padj < 0.05, na.rm = TRUE),
    sum(pb_ora$padj < 0.05, na.rm = TRUE)
  )
)
fwrite(summary, file.path(out_root, "pathway_analysis_summary.tsv"), sep = "\t")
print(summary)
