# Original spatial module runner

Recovered on 2026-11-21 from the original compressed project work records. The source at `code/module_spatial_refresh_20260917.py` is byte identical to the script identified by the existing `SOURCE_SHA256.json`: **9,106 bytes, 76 lines**, SHA256 `856541b63a13b5dbfedd8a49914807a17749d1ad330382f62e8cd59a6ef3bb02`.

The original server temporary file was absent. Historical source archives retained output receipts but not this runner. The original creation text and its two subsequent literal edits were recovered from the 2026-09-17 work record (records 771, 791 and 835). No historical shell command was executed, and no algorithm was reconstructed from the manuscript. The final matching hash distinguishes this strict89 script from the earlier 128-state null-pool draft.

## Original invocation and dependencies

The archived invocation used `/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python /tmp/module_spatial_refresh_20260917.py`. Historical absolute paths are preserved in the recovered source. The script writes to the original `publication_visualization_modules_20260917_strict89` result directory and has no output-directory argument. Do not invoke it against accepted outputs as a packaging check.

The script uses NumPy, pandas, SciPy, scikit-learn and threadpoolctl. It loads eight original helper definitions from `../overnight_manuscript_integration_20260917/spatial_joint_evidence_v1/code/legacy_spatial.py`: `zscore`, `matrix_correlations`, `upper_mean`, `spatial_cross_matrix`, `module_scores`, `lineage_from_state`, `bh_fdr` and `strip_cell2location_prefix`. Its own `neighbour_lag` supplies the coordinate-based six-neighbour calculation. The helper's original SHA256 is `c01db4d1116f4b42f35d9d91a0b850686a4c862bd6be744bfc742d3d1a1c400c`.

All 34 non-runner entries in `SOURCE_SHA256.json` were read and hashed on the server: the helper, 28 named/nuisance abundance files, module weights, node map, spot source, section audit and reference-factor map. Every hash matches the original receipt. This recovery did not reinstall or execute the original environment; it does not infer historical package versions from the currently installed software.

## Source-to-result correspondence

| Implementation in recovered source | Existing result |
| --- | --- |
| Lines 15–18: 89 nodes, 14 modules, 49 sections, 45,348 spots, 128 mapped factors, 6 stored programmes; seed 20260917 | Original completion receipt and result dimensions |
| Lines 19–34: 1,000 common draws per module and null method; strict-parent draws from the eligible 89-node pool | `null_draw_identity.tsv.gz`, 28,000 rows |
| Lines 36–50: q05 and mean abundance; section coherence and programme associations | `module_section_stats.tsv`, 1,372 rows |
| Line 51: medians from sections to existing units, then datasets | `module_biological_unit_stats.tsv`, 980 rows; `module_dataset_stats.tsv`, 196 rows |
| Lines 53–67: common-draw aggregation and BH across 14 modules × 2 metrics within each posterior/null method | `module_coherence_meta.tsv`, 112 rows; current Supplementary Table 21 |
| Lines 68–74: original positive-direction test and BH across 14 modules × 6 programmes within each posterior | `module_program_direction_meta.tsv`, 168 rows; current Supplementary Table 23B |

The current `Supplementary_Data.xlsx` sheets `Module_section_statistics` and `Module_dataset_statistics` retain the 1,372 and 196 records. All cells in these sheets and Supplementary Tables 21 and 23B were compared with the original result TSVs. Text, identities, row order and missing values match; numeric values agree within spreadsheet serialization precision (largest absolute difference below 6 × 10⁻¹⁶). No values were rewritten.

The six result files and `COMPLETE.json` match the current manuscript package's module source data byte for byte. The three module files used by the retained Figure4 plotting chain also match the server outputs. Historical `Figure5` source paths refer to provenance, not a fifth current main figure. `SOURCE_RECOVERY.json` contains the input checks, output hashes, worksheet mapping and recovery record references.

The source was syntax-parsed only. No analysis, null draws, significance tests, models, plots or workbooks were regenerated.
