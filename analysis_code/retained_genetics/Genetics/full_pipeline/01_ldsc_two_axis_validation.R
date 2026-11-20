options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(psych)
})

ROOT <- "/path/to/panvas_analysis/Genomic_SEM"
OUT <- "/path/to/panvas_analysis/Genomic_SEM_validation/01_ldsc_two_axis"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(OUT)

AUDITED <- "/path/to/panvas_analysis/Genomic_SEM_validation/00_input_audit/LDSCoutput_audited.RData"
load(if (file.exists(AUDITED)) AUDITED else file.path(ROOT, "LDSCoutput.RData"))
stopifnot(exists("LDSCoutput"))

traits <- colnames(LDSCoutput$S)
expected <- c("AA", "Stroke", "AAA", "CAVS", "CAD", "TAAD", "PAD", "CA")
stopifnot(setequal(traits, expected))

S <- LDSCoutput$S
V <- LDSCoutput$V
R <- LDSCoutput$S_Stand
VR <- LDSCoutput$V_Stand

fwrite(as.data.table(S, keep.rownames = "trait"),
       file.path(OUT, "LDSC_genetic_covariance.tsv"), sep = "\t")
fwrite(as.data.table(R, keep.rownames = "trait"),
       file.path(OUT, "LDSC_genetic_correlation.tsv"), sep = "\t")

# The order used by GenomicSEM::ldsc is:
# (1,1), (1,2), ..., (1,p), (2,2), ..., (p,p).
pair_table <- function(nm) {
  rbindlist(lapply(seq_along(nm), function(i) {
    data.table(i = i, j = i:length(nm), trait1 = nm[i],
               trait2 = nm[i:length(nm)])
  }))
}

pairs <- pair_table(traits)
pairs[, pair_index := .I]
pairs[, rg := R[cbind(i, j)]]
pairs[, rg_se := sqrt(diag(VR))]
pairs[, rg_z := rg / rg_se]
pairs[, rg_p := 2 * pnorm(abs(rg_z), lower.tail = FALSE)]

aortic <- c("AA", "AAA", "TAAD")
athero <- c("CAD", "Stroke", "PAD", "CA", "CAVS")
pairs[, block := fifelse(
  trait1 %in% aortic & trait2 %in% aortic, "within_aortic",
  fifelse(
    trait1 %in% athero & trait2 %in% athero, "within_atherosclerosis_related",
    "between_axes"
  )
)]
pairs[i == j, block := "diagonal"]

fwrite(pairs, file.path(OUT, "pairwise_rg_with_axis_labels.tsv"), sep = "\t")

# Quantify whether within-axis correlations exceed between-axis correlations.
off <- pairs[i != j]
off_idx <- off$pair_index
y <- off$rg
Voff <- VR[off_idx, off_idx, drop = FALSE]

contrast_test <- function(group_a, group_b, label) {
  w <- rep(0, nrow(off))
  w[off$block == group_a] <- 1 / sum(off$block == group_a)
  w[off$block == group_b] <- -1 / sum(off$block == group_b)
  estimate <- sum(w * y)
  se <- sqrt(drop(t(w) %*% Voff %*% w))
  data.table(
    contrast = label,
    estimate = estimate,
    se = se,
    z = estimate / se,
    p = 2 * pnorm(abs(estimate / se), lower.tail = FALSE)
  )
}

block_summary <- off[, .(
  n_pairs = .N,
  mean_rg = mean(rg),
  median_rg = median(rg),
  min_rg = min(rg),
  max_rg = max(rg)
), by = block]

contrasts <- rbind(
  contrast_test("within_aortic", "between_axes",
                "within_aortic_minus_between"),
  contrast_test("within_atherosclerosis_related", "between_axes",
                "within_atherosclerosis_related_minus_between")
)

fwrite(block_summary, file.path(OUT, "axis_block_rg_summary.tsv"), sep = "\t")
fwrite(contrasts, file.path(OUT, "axis_block_rg_contrasts.tsv"), sep = "\t")

# EFA on the smoothed genetic correlation matrix.
R_smooth <- as.matrix(nearPD(R, corr = TRUE)$mat)
kmo <- psych::KMO(R_smooth)
efa1 <- factanal(covmat = R_smooth, factors = 1, rotation = "none")
efa2 <- factanal(covmat = R_smooth, factors = 2, rotation = "promax")

load1 <- as.data.table(unclass(efa1$loadings), keep.rownames = "trait")
load2 <- as.data.table(unclass(efa2$loadings), keep.rownames = "trait")
fwrite(load1, file.path(OUT, "EFA_one_factor_loadings.tsv"), sep = "\t")
fwrite(load2, file.path(OUT, "EFA_two_factor_promax_loadings.tsv"), sep = "\t")
fwrite(data.table(overall_KMO = kmo$MSA),
       file.path(OUT, "KMO_overall.tsv"), sep = "\t")

# Reproduce paLDSC only if the existing compiled helper is available.
pa_helper <- file.path(ROOT, "12.EFA_pro_lib.Rc")
if (file.exists(pa_helper)) {
  compiler::loadcmp(pa_helper, env = environment())
  if (exists("paLDSC")) {
    pdf(file.path(OUT, "parallel_analysis_full_LDSC.pdf"), width = 7, height = 5)
    paLDSC(
      S = LDSCoutput$S_Stand,
      V = LDSCoutput$V_Stand,
      r = 1000,
      p = 0.95,
      diag = FALSE,
      fa = TRUE,
      fm = "minres",
      nfactors = 1,
      save.pdf = FALSE
    )
    dev.off()
  }
}

print(block_summary)
print(contrasts)
print(load2)
writeLines(capture.output(sessionInfo()), file.path(OUT, "R_sessionInfo.txt"))
