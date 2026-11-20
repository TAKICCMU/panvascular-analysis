#!/usr/bin/env bash
set -euo pipefail
JOB=/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908/production_launches/traits_20260914T023900Z
exec 9>"$JOB/.magma_queue.lock"
flock -n 9 || exit 73
if [[ -e "$JOB/MAGMA_EXIT_STATUS.json" ]]; then exit 74; fi
export CUDA_VISIBLE_DEVICES=""
export PYTHONDONTWRITEBYTECODE=1
export MPLBACKEND=Agg
export OMP_NUM_THREADS=12 OPENBLAS_NUM_THREADS=12 MKL_NUM_THREADS=12 NUMBA_NUM_THREADS=12 NUMEXPR_NUM_THREADS=12
/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python -B "$JOB/scripts/wait_for_phase1_r2.py"
date -Is > "$JOB/MAGMA_STARTED.txt"
set +e
/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python -B /home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908/downstream_CPU_L1L2_v2_1_release_compat_20260914/scripts/run_downstream.py --config /home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908/production_launches/traits_20260914T023900Z/configs/magma_config_APPROVED_r2.yaml --execute --modules cohort magma
rc=$?
set -e
printf '{"exit_code":%s,"completed_at":"%s"}\n' "$rc" "$(date -Is)" > "$JOB/MAGMA_EXIT_STATUS.json"
exit "$rc"
