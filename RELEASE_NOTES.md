# Release notes

## 22 November 2026 — Submission documentation and metadata

The current manuscript contains four main figures, ten Extended Data figures, 24 Supplementary Tables and six Supplementary Data tables. The consolidated table and panel maps in `06_Maps` remain the current navigation points. Older file and map names are retained for source correspondence.

### Disease scope

The primary pseudobulk comparisons use six studies and 53 sample records across three backgrounds:

| Disease background | Studies |
| --- | --- |
| Atherosclerosis | GSE159677; GSE248610 |
| Aortic dissection | GSE213740; GSE222318; GSE254132 |
| Pulmonary arterial hypertension | GSE228644 |

REML pools available estimates for a gene and cell state across these backgrounds. The candidate selection criterion requires results in both atherosclerosis and aortic dissection, the two backgrounds aligned with the component GWAS phenotypes. Pulmonary arterial hypertension contributes to pooling where an estimate exists but does not count towards that criterion. The 27 selected fibroblast pairs use the first two backgrounds; COL4A1 in endothelial cells and ARHGEF12 in macrophages also include pulmonary arterial hypertension. The candidate set and effect estimates are unchanged.

### GSE131778 donor metadata

The eight GEO records represent four patients:

| Source donor | GEO sample records |
| --- | --- |
| rca1 | GSM3819856; GSM3819857 |
| rca2 | GSM3819858; GSM3819859 |
| rca3 | GSM3819860; GSM3819861; GSM3819862 |
| rca4 | GSM3819863 |

The submission workbooks now distinguish these donor identities from the keys used in the saved analyses. The existing CoVarNet input contains eight grouping keys for this study; those original keys remain available for tracing the existing results. This metadata correction does not constitute a model refit. The existing Figure 1f summary and recorded donor-based summaries have not been recomputed.

### Reporting changes

The manuscript names the diseases contributing to the primary expression results and adds the ethics statement for secondary analysis of public, deidentified data. The Nature Reporting Summary follows the current manuscript and tables. The auxiliary `Rng` field was removed from the displayed NMF rank-selection block and its dictionary; random seeds required to reproduce the analysis remain documented.

This repository update changes the introduction and release documentation. Scientific source code, original execution records, stored results and model fits are unchanged.
