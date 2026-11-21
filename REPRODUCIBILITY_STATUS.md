# Source and reproducibility status

Static inspection on 2026-11-21 found that all 825 Python scripts and 214 R scripts parse. The 52 source lines in the existing model-definition index match the archived files. These checks do not execute the scientific workflow.

## Recovered final spatial module entry point

On 2026-11-21, retained project source text restored `/tmp/module_spatial_refresh_20260917.py` exactly: 9,106 bytes, 76 lines, SHA256 `856541b63a13b5dbfedd8a49914807a17749d1ad330382f62e8cd59a6ef3bb02`. This equals the hash already present in the original server `SOURCE_SHA256.json`.

The source is now retained at `analysis_code/native/downstream_final/publication_visualization_modules_20260917_strict89/code/module_spatial_refresh_20260917.py`. The original `legacy_spatial.py` helper remains unchanged. All 34 other files named by the original receipt match their recorded hashes on the server. The six result files and completion receipt match the current local source data byte for byte. Supplementary Table 21, Table 23B, and the module section/dataset detail sheets agree with those result files, allowing spreadsheet numeric serialization precision.

This resolves the missing stage29 entry source. `RECOVERY.md` and `SOURCE_RECOVERY.json` beside the retained receipts record the original invocation, helper dependencies, source lines, input hashes and output correspondence. The recovered script was parsed but not imported or run. No analysis, manuscript, figure, workbook or scientific result was changed; full reproduction still requires the original inputs and environment described below.

## Existing dependencies

Historical absolute paths occur in 290 source/configuration files. External count matrices, LD resources, fitted mapping objects, runtime environments and some manual annotation steps remain necessary. Their paths and identities are retained rather than silently rewritten. Earlier test logs remain historical execution records. No scientific analysis or model fitting was run during this inspection.

Public preparation utilities and metadata copies omit operational correspondence. `PUBLIC_SOURCE_NOTES.md` lists these non-computational publication edits. Original computational paths, model settings and archived result hashes are retained.
