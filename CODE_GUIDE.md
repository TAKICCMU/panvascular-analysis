# Code navigation

1. `06_Maps/MODULE_CODE_RESULT_CROSSWALK.tsv` preserves the original 32 workflow stage identifiers. Stages24–25 are marked removed from the current manuscript; the remaining stages are linked to current table numbers.
2. `analysis_code/` contains the actual retained computational scripts, local wrappers, configurations and execution provenance. `recovered_sources/` adds the original QC notebook cells, full-count builder, scIB benchmark code and plotting path map.
3. `03_Plot_Code/ADOPTED_CODE_ORDER.md` gives the adopted plotting/editing chain. The authoritative 14 SVG masters are `04_Current_Figures/Editable_SVG`. Final plotting patches include Figure1a body/vessel alignment and the merged Figure4. Earlier builders alone do not reproduce the final visual edits.
4. `table_code/20261117` contains the 60-table cleaning/renumbering code; `table_code/final_20261118` contains the latest Table8/24/59/60 and corresponding data-format corrections. The current workbook itself is in the data archive.
5. `tests/review_20260922/` contains existing test sources and original logs, including the original 76-test synthetic run and QLF/paired-block/REML/DL checks. These were not rerun during packaging. Historical result references inside logs retain the table/figure numbering at execution time; current numbering is mapped separately.

Unpack the code and data archives into the same parent directory; their top-level folders complement each other. `06_Maps/CURRENT_PANEL_DATA_CODE_MAP.tsv` uses paths relative to that combined root. Author notebooks, manual annotation decisions, several external count/LD resources and native model environments remain required for a full reproduction. This is an auditable source archive, not a claimed one-command portable workflow. Do not automatically run historical queue launchers or overwrite accepted outputs.

Environments and seeds are retained in the original parameter/receipt files and current Supplementary Table60. The original 76-test log includes shared package functions beyond the current manuscript scope; a test pass is not evidence of scientific validity. No new test result has been generated here.
