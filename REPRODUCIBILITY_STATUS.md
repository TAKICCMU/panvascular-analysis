# Source and reproducibility status

Static inspection on 2026-11-21 found that all 825 Python scripts and 214 R scripts parse. The 52 source lines in the existing model-definition index match the archived files. These checks do not execute the scientific workflow.

## Missing final module entry point

The source receipt for the final 14-module, 89-node spatial analysis records `/tmp/module_spatial_refresh_20260917.py` with SHA256 `856541b63a13b5dbfedd8a49914807a17749d1ad330382f62e8cd59a6ef3bb02`. That file is not included in the current source archive and was not found at its original server location or in the searched project code archives. The receipt and `COMPLETE.json` are retained under `analysis_code/native/downstream_final/publication_visualization_modules_20260917_strict89/`.

The available `legacy_spatial.py` supplies helper functions and is not a substitute for that final runner. The final coherence results remain unchanged, but stage29 cannot yet be reproduced from this repository alone. Recover the original runner and confirm its recorded hash before marking this source gap resolved.

## Existing dependencies

Historical absolute paths occur in 290 source/configuration files. External count matrices, LD resources, fitted mapping objects, runtime environments and some manual annotation steps remain necessary. Their paths and identities are retained rather than silently rewritten. Earlier test logs remain historical execution records. No scientific analysis or model fitting was run during this inspection.
