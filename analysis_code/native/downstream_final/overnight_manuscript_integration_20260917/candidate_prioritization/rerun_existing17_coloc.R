#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  setDTthreads(1)
  library(coloc)
  library(ggplot2)
})

ROOT <- "/Users/taki/Desktop/Codex任务1/overnight_candidate_prioritization_20260917/coloc_existing17"
MANIFEST <- file.path(ROOT, "00_metadata", "qtl_download_manifest.tsv")
OUT <- file.path(ROOT, "02_results")
FIG <- file.path(ROOT, "03_figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

GWAS <- c(
  primary = "/Users/taki/Desktop/FUMA/full_bifactor_G_FUMA_all.tsv.gz",
  leave_one_CAD = "/Users/taki/Desktop/FUMA/LOT_CAD/leave_one_CAD_bifactor_G_FUMA_all.tsv.gz"
)
P1 <- 1e-4; P2 <- 1e-4; P12 <- 1e-5
P12_GRID <- c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4)
MIN_OVERLAP <- 100L

comp <- function(x) chartr("ACGT", "TGCA", toupper(x))
is_snv <- function(x) !is.na(x) & nchar(x) == 1L & x %chin% c("A", "C", "G", "T")
is_pal <- function(a, b) paste0(a, b) %chin% c("AT", "TA", "CG", "GC")
safe <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

manifest <- fread(MANIFEST)
manifest[, n_qtl_rows := as.integer(n_qtl_rows)]
manifest[, analysis_ready_source := !is.na(output) & output != "" &
           n_qtl_rows >= MIN_OVERLAP & file.exists(output)]

parse_qtl <- function(path) {
  x <- fread(path, showProgress = FALSE)
  x[, gene_id := sub("\\..*$", "", as.character(gene_id))]
  x[, `:=`(
    rsid = trimws(as.character(rsid)), ref = toupper(as.character(ref)),
    alt = toupper(as.character(alt)), beta = as.numeric(beta),
    se = as.numeric(se), maf = as.numeric(maf), an = as.numeric(an),
    pvalue = as.numeric(pvalue), position = as.integer(position)
  )]
  x <- x[grepl("^rs[0-9]+$", rsid) & is.finite(beta) & is.finite(se) & se > 0 &
           is.finite(maf) & maf > 0 & maf <= 0.5 & is.finite(an) & an > 0]
  setorder(x, rsid, pvalue)
  unique(x, by = "rsid")
}

qtl_rows <- manifest[analysis_ready_source == TRUE]
qtl_list <- lapply(qtl_rows$output, parse_qtl)
names(qtl_list) <- paste(qtl_rows$gene, qtl_rows$context, sep = "__")
all_rsids <- unique(unlist(lapply(qtl_list, `[[`, "rsid"), use.names = FALSE))
message("Analysis-ready QTL contexts: ", nrow(qtl_rows), "; unique rsIDs: ", length(all_rsids))

harmonise <- function(g, q) {
  z <- merge(g, q, by = "rsid", suffixes = c("_gwas", "_qtl"), all = FALSE)
  if (!nrow(z)) return(z)
  z[, `:=`(A1 = toupper(A1), A2 = toupper(A2), ref = toupper(ref), alt = toupper(alt),
            beta_qtl = beta, se_qtl = se)]
  z[, palindromic := is_snv(ref) & is_snv(alt) & is_pal(ref, alt)]
  z[, orientation := fifelse(A1 == alt & A2 == ref, "direct",
                      fifelse(A1 == ref & A2 == alt, "reverse",
                      fifelse(!palindromic & is_snv(A1) & is_snv(A2) & is_snv(ref) & is_snv(alt) &
                                A1 == comp(alt) & A2 == comp(ref), "strand_direct",
                      fifelse(!palindromic & is_snv(A1) & is_snv(A2) & is_snv(ref) & is_snv(alt) &
                                A1 == comp(ref) & A2 == comp(alt), "strand_reverse", NA_character_))))]
  z <- z[!is.na(orientation)]
  z <- z[!(palindromic & (MAF >= 0.42 | maf >= 0.42))]
  setorder(z, rsid, pvalue)
  z <- unique(z, by = "rsid")
  z[, beta_qtl_aligned := fifelse(orientation %chin% c("direct", "strand_direct"), beta_qtl, -beta_qtl)]
  z[, direction_interpretable := !palindromic]
  z[]
}

quiet_coloc <- function(...) {
  fit <- NULL
  invisible(capture.output(fit <- coloc.abf(...)))
  fit
}

classify <- function(h0, h1, h2, h3, h4, pg, pq) {
  cond <- ifelse(h3 + h4 > 0, h4 / (h3 + h4), NA_real_)
  if (!is.finite(pg) || pg >= 5e-8) return("no_regional_GWAS_signal")
  if (!is.finite(pq) || pq >= 1e-4) return("no_strong_regional_eQTL")
  if (h4 >= 0.80 && cond >= 0.80) return("strong_colocalisation")
  if (h4 >= 0.50 && cond >= 0.70) return("suggestive_colocalisation")
  if (h3 >= 0.80) return("distinct_signals_H3")
  "ambiguous_regional_evidence"
}

results <- list(); sensitivity <- list(); audits <- list()
ri <- 0L; si <- 0L; ai <- 0L

for (scenario in names(GWAS)) {
  message("Reading GWAS: ", scenario)
  g <- fread(GWAS[[scenario]], select = c("SNP", "CHR", "BP", "A1", "A2", "MAF", "BETA", "SE", "P", "N"))
  setnames(g, "SNP", "rsid")
  g <- g[rsid %chin% all_rsids]
  g[, `:=`(BETA = as.numeric(BETA), SE = as.numeric(SE), MAF = as.numeric(MAF),
            P = as.numeric(P), N = as.numeric(N))]
  g <- g[grepl("^rs[0-9]+$", rsid) & is.finite(BETA) & is.finite(SE) & SE > 0 &
           is.finite(MAF) & MAF > 0 & MAF <= 0.5 & is.finite(N) & N > 0]
  setorder(g, rsid, P); g <- unique(g, by = "rsid")

  for (j in seq_len(nrow(qtl_rows))) {
    meta <- qtl_rows[j]
    q <- qtl_list[[j]]
    h <- harmonise(g, q)
    ai <- ai + 1L
    audits[[ai]] <- data.table(
      scenario = scenario, gene = meta$gene, gene_id = meta$gene_id,
      context = meta$context, source = meta$source, nominated_lineage = meta$nominated_lineage,
      n_qtl_valid = nrow(q), n_harmonised = nrow(h),
      n_palindromic_retained = if (nrow(h)) sum(h$palindromic) else 0L,
      test_status = ifelse(nrow(h) >= MIN_OVERLAP, "tested", "insufficient_overlap")
    )
    if (nrow(h) < MIN_OVERLAP) next

    d1 <- list(beta = h$BETA, varbeta = h$SE^2, MAF = h$MAF,
               N = median(h$N), snp = h$rsid, type = "quant", sdY = 1)
    d2 <- list(beta = h$beta_qtl_aligned, varbeta = h$se_qtl^2, MAF = h$maf,
               N = as.numeric(meta$sample_size), snp = h$rsid, type = "quant", sdY = 1)
    fit <- quiet_coloc(d1, d2, p1 = P1, p2 = P2, p12 = P12)
    pp <- as.numeric(fit$summary[paste0("PP.H", 0:4, ".abf")])
    snp <- as.data.table(fit$results); setorder(snp, -SNP.PP.H4)
    top <- snp[1]; top_h <- h[match(top$snp, rsid)]
    min_pg <- min(h$P); min_pq <- min(h$pvalue)
    ri <- ri + 1L
    results[[ri]] <- data.table(
      scenario = scenario, gene = meta$gene, gene_id = meta$gene_id,
      nominated_lineage = meta$nominated_lineage, context = meta$context,
      qtl_source = meta$source, qtl_sample_size = meta$sample_size,
      n_overlap = nrow(h), gwas_effective_N = median(h$N),
      min_gwas_p = min_pg, min_eqtl_p = min_pq,
      PP.H0 = pp[1], PP.H1 = pp[2], PP.H2 = pp[3], PP.H3 = pp[4], PP.H4 = pp[5],
      H4_given_H3H4 = pp[5] / (pp[4] + pp[5]), top_shared_snp = top$snp,
      top_shared_SNP_PP.H4 = top$SNP.PP.H4, top_gwas_beta = top_h$BETA,
      top_eqtl_beta_aligned = top_h$beta_qtl_aligned,
      top_direction_interpretable = top_h$direction_interpretable,
      top_effect_direction = ifelse(top_h$direction_interpretable,
        ifelse(sign(top_h$BETA) == sign(top_h$beta_qtl_aligned), "same", "opposite"),
        "palindromic_ambiguous"),
      classification = classify(pp[1], pp[2], pp[3], pp[4], pp[5], min_pg, min_pq),
      primary_prior_p12 = P12, method = "coloc.abf_complete_cis_region"
    )
    for (p12 in P12_GRID) {
      fs <- quiet_coloc(d1, d2, p1 = P1, p2 = P2, p12 = p12)
      ps <- as.numeric(fs$summary[paste0("PP.H", 0:4, ".abf")])
      si <- si + 1L
      sensitivity[[si]] <- data.table(
        scenario = scenario, gene = meta$gene, context = meta$context,
        qtl_source = meta$source, p12 = p12, PP.H3 = ps[4], PP.H4 = ps[5],
        H4_given_H3H4 = ps[5] / (ps[4] + ps[5])
      )
    }
  }
  rm(g); invisible(gc())
}

audit <- rbindlist(audits, fill = TRUE)
res <- rbindlist(results, fill = TRUE)
sens <- rbindlist(sensitivity, fill = TRUE)
fwrite(audit, file.path(OUT, "input_harmonisation_audit.tsv"), sep = "\t")
fwrite(res, file.path(OUT, "Tier1_full_coloc_abf_results.tsv.gz"), sep = "\t")
fwrite(sens, file.path(OUT, "Tier1_coloc_prior_sensitivity.tsv.gz"), sep = "\t")

# Explicit all-gene coverage table, including genes/contexts that could not be tested.
coverage <- CJ(scenario = names(GWAS), gene = unique(manifest$gene), unique = TRUE)
tested_summary <- res[, .(
  n_contexts_tested = .N, n_strong = sum(classification == "strong_colocalisation"),
  n_suggestive = sum(classification == "suggestive_colocalisation"),
  n_distinct_H3 = sum(classification == "distinct_signals_H3"),
  maximum_PP.H4 = max(PP.H4), maximum_H4_given_H3H4 = max(H4_given_H3H4),
  best_context = context[which.max(PP.H4)], best_qtl_source = qtl_source[which.max(PP.H4)],
  best_classification = classification[which.max(PP.H4)]
), by = .(scenario, gene)]
source_summary <- manifest[, .(
  n_contexts_requested = .N, n_contexts_with_qtl = sum(n_qtl_rows >= MIN_OVERLAP, na.rm = TRUE),
  n_contexts_no_qtl = sum(n_qtl_rows < MIN_OVERLAP, na.rm = TRUE),
  n_remote_failures = sum(status == "TABIX_FAILED")
), by = gene]
coverage <- merge(coverage, source_summary, by = "gene", all.x = TRUE)
coverage <- merge(coverage, tested_summary, by = c("scenario", "gene"), all.x = TRUE)
coverage[is.na(n_contexts_tested), `:=`(n_contexts_tested = 0L, n_strong = 0L,
  n_suggestive = 0L, n_distinct_H3 = 0L, best_classification = "not_tested_no_complete_QTL_or_overlap")]
setorder(coverage, scenario, gene)
fwrite(coverage, file.path(OUT, "Tier1_gene_level_colocalisation_coverage.tsv"), sep = "\t")


writeLines(capture.output(sessionInfo()), file.path(ROOT, "00_metadata", "R_sessionInfo.txt"))
