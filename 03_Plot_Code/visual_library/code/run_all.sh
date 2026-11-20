#!/bin/bash
set -euo pipefail
TASK_PYTHON="${PANVASCULAR_PLOT_PYTHON:-/Users/taki/.cache/panvascular-panel-plotting/bin/python}"
TASK_CODE_DIR="$(cd "$(dirname "$0")" && pwd)"
"$TASK_PYTHON" "$TASK_CODE_DIR/redraw_main.py"
"$TASK_PYTHON" "$TASK_CODE_DIR/refine_native_and_ed.py"
"$TASK_PYTHON" "$TASK_CODE_DIR/package_review.py"
