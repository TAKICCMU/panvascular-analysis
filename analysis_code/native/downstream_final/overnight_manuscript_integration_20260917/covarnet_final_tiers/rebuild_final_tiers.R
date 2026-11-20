Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
library(data.table); setDTthreads(1)
WROOT <- "/data/anzhen_vhd/scRNAseq/codex/panvascular_downstream_storage_v1/expanded_parallel_20260916_v1"
ROOT <- file.path(WROOT,"covarnet"); PRIMARY_ROOT <- file.path(ROOT,"stable_NMF")
TABLE_ROOT <- file.path(WROOT,"overnight_manuscript_integration_20260917/covarnet_final_tiers")
read_matrix <- function(path) {d<-fread(path);r<-d[[1]];m<-as.matrix(d[,-1]);storage.mode(m)<-"double";rownames(m)<-r;m}
W_reference <- read_matrix(file.path(PRIMARY_ROOT,"module_node_weights_K14.tsv"))
H_reference <- t(read_matrix(file.path(PRIMARY_ROOT,"sample_module_relative_coefficient_K14.tsv")))
meta <- fread(file.path(ROOT,"CoVarNet_sample_metadata.tsv"))
meta <- meta[match(colnames(H_reference),covarnet_sample_uid)]
stopifnot(!anyNA(meta$covarnet_sample_uid),!anyDuplicated(meta$covarnet_sample_uid),max(abs(colSums(H_reference)-1))<1e-10)
latest<-fread(file.path(WROOT,"remaining_downstream_20260916_v1/covarnet/association_patch_v2/results/aligned_latest_metadata.tsv"))
latest<-latest[match(meta$covarnet_sample_uid,covarnet_sample_uid)]
present<-!is.na(latest$covarnet_sample_uid)
stopifnot(identical(meta$study_canonical[present],latest$study_canonical[present]),identical(meta$tissue_v2[present],latest$tissue_v2[present]))
fwrite(meta[!present,.(covarnet_sample_uid,study_canonical,tissue_v2)],file.path(TABLE_ROOT,"profiles_in_structure_not_association.tsv"),sep="\t")
prop <- read_matrix(file.path(ROOT,"CoVarNet_node_proportion_of_all_atlas_cells.tsv"))[,meta$covarnet_sample_uid,drop=FALSE]
stopifnot(identical(colnames(prop),colnames(H_reference)),all(rownames(W_reference)%in%rownames(prop)),ncol(prop)==223,ncol(W_reference)==14)
studies <- unique(meta$study_canonical)
TOP_NODES <- 10L; TOP_NODES_COHERENCE <- 5L; MIN_WITHIN_STUDY_N <- 5L
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



gate<-fread(file.path(WROOT,"remaining_downstream_20260916_v1/covarnet/stability_gate_audit_v1/module_stability_gate.tsv"))
setnames(gate,"reference_module","module")
module_tier<-Reduce(function(x,y)merge(x,y,by="module",all=TRUE),list(gate,concentration,coherence_module))
module_tier[, broad_context := top_study_share < .50 & top_tissue_share < .70 & effective_n_studies >=3 & effective_n_tissues >=2]
module_tier[is.na(broad_context),broad_context:=FALSE]
module_tier[,coherent_within_study := median_pooled_rho>=.20 & positive_pair_fraction>=.70]
module_tier[is.na(coherent_within_study),coherent_within_study:=FALSE]
module_tier[, module_tier:=fifelse(stable_structure & broad_context & coherent_within_study,"high_confidence_pan_vascular",fifelse(stable_structure & coherent_within_study & !broad_context,"context_specific_reproducible","candidate_or_rank_sensitive"))]
module_tier[,分级:=fifelse(module_tier=="high_confidence_pan_vascular","高置信泛血管模块",fifelse(module_tier=="context_specific_reproducible","情境特异可重复模块","候选或秩敏感模块"))]
setorder(module_tier,module)
fwrite(module_tier,file.path(TABLE_ROOT,"14_final_module_evidence_and_tiers.tsv"),sep="\t")
fwrite(module_tier[,.(模块=module,结构稳定=stable_structure,广泛覆盖=broad_context,研究内一致=coherent_within_study,最大研究占比=top_study_share,最大组织占比=top_tissue_share,有效研究数=effective_n_studies,有效组织数=effective_n_tissues,合并相关中位数=median_pooled_rho,正向节点对比例=positive_pair_fraction,未通过稳定门槛=failed_gates,分级)],file.path(TABLE_ROOT,"最终14模块中文分级.tsv"),sep="\t")
stopifnot(nrow(module_tier)==14,sum(module_tier$stable_structure)==1,module_tier[stable_structure==TRUE,module]=="CM09")
capture.output(sessionInfo(),file=file.path(TABLE_ROOT,"sessionInfo.txt"))
print(module_tier[,.(module,stable_structure,broad_context,coherent_within_study,median_pooled_rho,positive_pair_fraction,module_tier)])
