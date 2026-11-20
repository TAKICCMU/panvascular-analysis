#!/bin/bash
set -euo pipefail
TASK_PYTHON='/Users/taki/.cache/panvascular-panel-plotting/bin/python'
TASK_CODE_DIR="$(cd "$(dirname "$0")" && pwd)"
for TASK_SCRIPT in export_native_panels.py redraw_data_panels.py redraw_evidence_panels.py restore_standalone_context.py reflow_native_maps.py redraw_genetic_matrices.py upgrade_frozen_umap_raster.py finalize_delivery.py; do
  "$TASK_PYTHON" "$TASK_CODE_DIR/$TASK_SCRIPT"
done
