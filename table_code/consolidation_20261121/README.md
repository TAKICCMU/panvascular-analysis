# Supplementary table consolidation

These scripts reorganise the existing 60 tables into 24 Supplementary Tables and six Supplementary Data tables. They do not refit scientific models. The original 60-table inputs remain linked in `06_Maps/ALL60_TABLE_CODE_MAP.tsv`, and the current numbering is in `06_Maps/Supplementary_Table_Crosswalk.tsv`.

The scripts retain their original filesystem dependencies. `prepare_consolidation.py` reads the original workbook and generates intermediate JSON. `build_workbooks.mjs` uses the artifact spreadsheet runtime to write the reorganised workbooks. `verify_workbooks.py` compares their values. The two JSON reports are the original execution records, not new test results.

Review the original paths and supply the stated workbook/runtime dependencies before executing. The current workbooks are maintained in the manuscript package. The 2026-11-20 data archive predates this consolidation and retains its original 60-table layout.
