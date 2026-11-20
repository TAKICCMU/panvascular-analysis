#!/usr/bin/env python3
"""Run a focused multi-signal colocalisation sensitivity analysis locally.

The analysis is deliberately restricted to primary ABF-positive Tier 1
gene-context pairs. It tests whether the reported positive regional signals
remain supported when more than one association signal is allowed in either
trait. It is not a new genome-wide or all-context discovery screen.
"""

from __future__ import annotations

import gzip
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path("/path/to/panvas_project")
COLOC_ROOT = ROOT / "Tier1_full_colocalisation_v1"
ABF_FILE = COLOC_ROOT / "02_results/Tier1_full_coloc_abf_results.tsv.gz"
QTL_ROOT = COLOC_ROOT / "01_qtl_regions_uniform"
PRIMARY_GWAS = Path("/path/to/FUMA/full_bifactor_G_FUMA_all.tsv.gz")
LOT_GWAS = Path(
    "/path/to/FUMA/LOT_CAD/"
    "leave_one_CAD_bifactor_G_FUMA_all.tsv.gz"
)
PLINK = Path("/path/to/FUMA/software/plink1.9/plink")
REFERENCE = ROOT / "data_maf0.01_rs_ref"
RSCRIPT = Path("/usr/local/bin/Rscript")
OUT = ROOT / "Tier1_positive_multisignal_coloc_local_v1"


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


for required in [
    ABF_FILE,
    PRIMARY_GWAS,
    LOT_GWAS,
    PLINK,
    RSCRIPT,
    Path(str(REFERENCE) + ".bed"),
    Path(str(REFERENCE) + ".bim"),
    Path(str(REFERENCE) + ".fam"),
]:
    require(required)

OUT.mkdir(parents=True, exist_ok=True)
ld_root = OUT / "01_ld_reference"
ld_root.mkdir(exist_ok=True)

abf = pd.read_csv(ABF_FILE, sep="\t")
selected = abf.loc[
    (abf["scenario"] == "primary")
    & abf["classification"].isin(
        ["strong_colocalisation", "suggestive_colocalisation"]
    )
].copy()
selected = selected.sort_values(["gene", "context"]).reset_index(drop=True)
if selected.shape[0] != 11:
    raise RuntimeError(
        f"Expected 11 primary ABF-positive contexts, observed {selected.shape[0]}"
    )

selected["qtl_file"] = selected.apply(
    lambda x: str(QTL_ROOT / f"{x['gene']}__{x['context']}__cis1Mb.tsv.gz"),
    axis=1,
)
for p in selected["qtl_file"]:
    require(Path(p))
selected.to_csv(OUT / "00_selected_ABF_positive_contexts.tsv", sep="\t", index=False)


def read_rsids(path: Path) -> set[str]:
    data = pd.read_csv(path, sep="\t", usecols=["rsid"])
    return set(data["rsid"].dropna().astype(str))


gene_ld_manifest = []
for gene, frame in selected.groupby("gene", sort=True):
    rsids: set[str] = set()
    for qtl_file in frame["qtl_file"]:
        rsids.update(read_rsids(Path(qtl_file)))
    extract_file = ld_root / f"{gene}.extract.txt"
    extract_file.write_text("\n".join(sorted(rsids)) + "\n", encoding="utf-8")
    prefix = ld_root / gene
    cmd_make = [
        str(PLINK),
        "--bfile",
        str(REFERENCE),
        "--extract",
        str(extract_file),
        "--keep-allele-order",
        "--make-bed",
        "--out",
        str(prefix),
    ]
    cmd_ld = [
        str(PLINK),
        "--bfile",
        str(prefix),
        "--keep-allele-order",
        "--r",
        "square",
        "gz",
        "--out",
        str(prefix),
    ]
    subprocess.run(cmd_make, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    subprocess.run(cmd_ld, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    bim = pd.read_csv(
        prefix.with_suffix(".bim"),
        sep=r"\s+",
        header=None,
        names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
    )
    gene_ld_manifest.append(
        {
            "gene": gene,
            "n_qtl_rsids_requested": len(rsids),
            "n_reference_variants": int(bim.shape[0]),
            "bfile_prefix": str(prefix),
            "ld_file": str(Path(str(prefix) + ".ld.gz")),
        }
    )

pd.DataFrame(gene_ld_manifest).to_csv(
    OUT / "01_LD_reference_manifest.tsv", sep="\t", index=False
)

r_script = OUT / "02_run_coloc_susie.R"
r_script.write_text(
    r'''
suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(susieR)
})

args <- commandArgs(trailingOnly = TRUE)
selected_file <- args[[1]]
primary_gwas <- args[[2]]
lot_gwas <- args[[3]]
ld_root <- args[[4]]
out_dir <- args[[5]]

selected <- fread(selected_file)
gwas_list <- list(primary = fread(primary_gwas), leave_one_CAD = fread(lot_gwas))

comp <- function(x) {
  x <- toupper(as.character(x))
  chartr("ACGT", "TGCA", x)
}

is_pal <- function(a, b) {
  paste0(toupper(a), toupper(b)) %chin% c("AT", "TA", "CG", "GC")
}

align_effect <- function(effect, other, beta, ref_a1, ref_a2) {
  effect <- toupper(as.character(effect))
  other <- toupper(as.character(other))
  ref_a1 <- toupper(as.character(ref_a1))
  ref_a2 <- toupper(as.character(ref_a2))
  out <- rep(NA_real_, length(beta))
  direct_same <- effect == ref_a1 & other == ref_a2
  direct_flip <- effect == ref_a2 & other == ref_a1
  strand_same <- comp(effect) == ref_a1 & comp(other) == ref_a2
  strand_flip <- comp(effect) == ref_a2 & comp(other) == ref_a1
  out[direct_same | strand_same] <- beta[direct_same | strand_same]
  out[direct_flip | strand_flip] <- -beta[direct_flip | strand_flip]
  out[is_pal(effect, other)] <- NA_real_
  out
}

result_rows <- list()
signal_rows <- list()
failure_rows <- list()

for (i in seq_len(nrow(selected))) {
  row <- selected[i]
  gene <- row$gene
  context <- row$context
  qtl <- fread(row$qtl_file)
  qtl <- qtl[!is.na(rsid) & !duplicated(rsid)]

  prefix <- file.path(ld_root, gene)
  bim <- fread(paste0(prefix, ".bim"), header = FALSE)
  setnames(bim, c("CHR", "SNP", "CM", "BP", "A1", "A2"))
  ld <- as.matrix(fread(paste0(prefix, ".ld.gz"), header = FALSE))
  if (nrow(ld) != nrow(bim) || ncol(ld) != nrow(bim)) {
    stop("LD/BIM dimension mismatch for ", gene)
  }
  rownames(ld) <- colnames(ld) <- bim$SNP

  for (scenario in names(gwas_list)) {
    gwas <- copy(gwas_list[[scenario]])
    gwas <- gwas[!is.na(SNP) & !duplicated(SNP)]
    common <- Reduce(intersect, list(bim$SNP, qtl$rsid, gwas$SNP))
    status <- "ok"
    error <- NA_character_
    n_gwas_cs <- NA_integer_
    n_qtl_cs <- NA_integer_
    max_pp_h4 <- NA_real_
    best_gwas_signal <- NA_integer_
    best_qtl_signal <- NA_integer_
    n_overlap <- length(common)

    tryCatch({
      if (n_overlap < 100) stop("Fewer than 100 overlapping variants")
      ref <- bim[match(common, SNP)]
      gg <- gwas[match(common, SNP)]
      qq <- qtl[match(common, rsid)]

      zg <- align_effect(gg$A1, gg$A2, gg$BETA / gg$SE, ref$A1, ref$A2)
      zq <- align_effect(qq$alt, qq$ref, qq$beta / qq$se, ref$A1, ref$A2)
      keep <- is.finite(zg) & is.finite(zq) & is.finite(gg$SE) & gg$SE > 0 &
        is.finite(qq$se) & qq$se > 0
      snps <- common[keep]
      if (length(snps) < 100) stop("Fewer than 100 allele-aligned variants")

      R <- ld[snps, snps, drop = FALSE]
      R <- (R + t(R)) / 2
      R <- 0.999 * R
      diag(R) <- 1
      rownames(R) <- colnames(R) <- snps

      gg <- gg[keep]
      qq <- qq[keep]
      N_gwas <- as.numeric(median(gg$N, na.rm = TRUE))
      N_qtl <- as.numeric(row$qtl_sample_size)
      d1 <- list(
        beta = align_effect(gg$A1, gg$A2, gg$BETA, ref$A1[keep], ref$A2[keep]),
        varbeta = gg$SE^2,
        snp = snps,
        N = N_gwas,
        type = "quant",
        sdY = 1,
        LD = R
      )
      d2 <- list(
        beta = align_effect(qq$alt, qq$ref, qq$beta, ref$A1[keep], ref$A2[keep]),
        varbeta = qq$se^2,
        snp = snps,
        N = N_qtl,
        type = "quant",
        sdY = 1,
        LD = R
      )
      s1 <- coloc::runsusie(
        d1, maxit = 1000, repeat_until_convergence = FALSE,
        L = 10, estimate_residual_variance = FALSE
      )
      s2 <- coloc::runsusie(
        d2, maxit = 1000, repeat_until_convergence = FALSE,
        L = 10, estimate_residual_variance = FALSE
      )
      n_gwas_cs <- if (is.null(s1$sets$cs)) 0L else length(s1$sets$cs)
      n_qtl_cs <- if (is.null(s2$sets$cs)) 0L else length(s2$sets$cs)
      cc <- coloc::coloc.susie(
        s1, s2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5
      )
      if (!is.null(cc$summary) && nrow(cc$summary) > 0) {
        sm <- as.data.table(cc$summary)
        pp_col <- grep("PP.H4", names(sm), value = TRUE)[1]
        if (is.na(pp_col)) stop("coloc.susie summary lacks PP.H4")
        best <- which.max(sm[[pp_col]])
        max_pp_h4 <- as.numeric(sm[[pp_col]][best])
        best_gwas_signal <- as.integer(sm$idx1[best])
        best_qtl_signal <- as.integer(sm$idx2[best])
        sm[, `:=`(
          scenario = scenario,
          gene = gene,
          context = context,
          n_variants = length(snps)
        )]
        signal_rows[[length(signal_rows) + 1L]] <- sm
      } else {
        status <- "no_signal_pair"
      }
      n_overlap <- length(snps)
    }, error = function(e) {
      status <<- "failed"
      error <<- conditionMessage(e)
      failure_rows[[length(failure_rows) + 1L]] <<- data.table(
        scenario = scenario, gene = gene, context = context,
        error = error
      )
    })

    result_rows[[length(result_rows) + 1L]] <- data.table(
      scenario = scenario,
      gene = gene,
      context = context,
      ABF_PP.H4_primary = row$PP.H4,
      ABF_classification_primary = row$classification,
      n_variants = n_overlap,
      n_GWAS_credible_sets = n_gwas_cs,
      n_QTL_credible_sets = n_qtl_cs,
      max_SuSiE_PP.H4 = max_pp_h4,
      best_GWAS_signal = best_gwas_signal,
      best_QTL_signal = best_qtl_signal,
      status = status,
      error = error
    )
  }
}

summary <- rbindlist(result_rows, fill = TRUE)
fwrite(summary, file.path(out_dir, "03_coloc_susie_summary.tsv"), sep = "\t")
if (length(signal_rows)) {
  fwrite(rbindlist(signal_rows, fill = TRUE),
         file.path(out_dir, "04_coloc_susie_signal_pairs.tsv.gz"), sep = "\t")
}
if (length(failure_rows)) {
  fwrite(rbindlist(failure_rows, fill = TRUE),
         file.path(out_dir, "05_coloc_susie_failures.tsv"), sep = "\t")
} else {
  fwrite(data.table(scenario=character(), gene=character(), context=character(), error=character()),
         file.path(out_dir, "05_coloc_susie_failures.tsv"), sep = "\t")
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo.txt"))
''',
    encoding="utf-8",
)

subprocess.run(
    [
        str(RSCRIPT),
        str(r_script),
        str(OUT / "00_selected_ABF_positive_contexts.tsv"),
        str(PRIMARY_GWAS),
        str(LOT_GWAS),
        str(ld_root),
        str(OUT),
    ],
    check=True,
)

summary = pd.read_csv(OUT / "03_coloc_susie_summary.tsv", sep="\t")
manifest = {
    "analysis": "Focused multi-signal colocalisation sensitivity analysis",
    "scope": "Primary ABF strong/suggestive Tier 1 gene-context pairs only",
    "n_selected_contexts": int(selected.shape[0]),
    "n_scenarios": 2,
    "n_tests": int(summary.shape[0]),
    "n_ok": int((summary["status"] == "ok").sum()),
    "n_no_signal_pair": int((summary["status"] == "no_signal_pair").sum()),
    "n_failed": int((summary["status"] == "failed").sum()),
    "thresholds": {
        "strong": "maximum signal-pair PP.H4 >= 0.80",
        "suggestive": "0.50 <= maximum signal-pair PP.H4 < 0.80",
    },
    "important_interpretation": (
        "This sensitivity analysis assesses the positive ABF shortlist under a "
        "multi-signal model. It is not an independent discovery analysis and "
        "does not classify untested contexts as negative."
    ),
    "generated_at": datetime.now(timezone.utc).isoformat(),
}
(OUT / "run_manifest.json").write_text(
    json.dumps(manifest, indent=2), encoding="utf-8"
)
print(json.dumps(manifest, indent=2))
print(summary.to_string(index=False))
