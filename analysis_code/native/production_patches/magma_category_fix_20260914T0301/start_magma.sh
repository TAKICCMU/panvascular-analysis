#!/usr/bin/env bash
set -euo pipefail
JOB=/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908/production_launches/magma_category_fix_20260914T0301
exec 9>"$JOB/.magma_launch.lock"
flock -n 9 || exit 73
[[ ! -e "$JOB/MAGMA_EXIT_STATUS.json" ]] || exit 74
export CUDA_VISIBLE_DEVICES="" PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
export OMP_NUM_THREADS=12 OPENBLAS_NUM_THREADS=12 MKL_NUM_THREADS=12 NUMBA_NUM_THREADS=12 NUMEXPR_NUM_THREADS=12
/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python -B "$JOB/scripts/resource_gate.py"
date -Is > "$JOB/MAGMA_STARTED.txt"
set +e
/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python -B /home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908/production_launches/magma_category_fix_20260914T0301/magma_adapter_v1/scripts/run_downstream.py --config /home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908/production_launches/magma_category_fix_20260914T0301/configs/magma_config_APPROVED.yaml --execute --modules cohort magma
rc=$?
set -e
printf '{"exit_code":%s,"completed_at":"%s"}\n' "$rc" "$(date -Is)" > "$JOB/MAGMA_EXIT_STATUS.json"
exit "$rc"
