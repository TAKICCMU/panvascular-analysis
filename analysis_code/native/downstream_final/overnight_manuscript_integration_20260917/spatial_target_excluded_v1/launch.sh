#!/usr/bin/env bash
set -euo pipefail
JOB=/data/anzhen_vhd/scRNAseq/codex/panvascular_downstream_storage_v1/expanded_parallel_20260916_v1/overnight_manuscript_integration_20260917/spatial_target_excluded_v1
export CUDA_VISIBLE_DEVICES=GPU-ee506680-afc9-11fd-3ed0-99b220571975
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMBA_NUM_THREADS=4
export MPLBACKEND=Agg PYTHONUNBUFFERED=1
PY=/home/anzhen_vhd/miniconda3/envs/cell2location_GPU/bin/python
exec 9>"$JOB/.launch.lock"
flock -n 9 || { echo 'This independent branch is already running'; exit 73; }
# Only this branch lock; dispatcher owns shared resource/queue admission.
if [[ ! -f "$JOB/smoke/PASS.json" ]]; then
 "$PY" "$JOB/code/spatial_target_excluded.py" smoke >>"$JOB/smoke.log" 2>&1
fi
exec "$PY" "$JOB/code/spatial_target_excluded.py" production >>"$JOB/production.log" 2>&1
