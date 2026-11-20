#!/usr/bin/env bash
set -euo pipefail
JOB=/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908/production_launches/phase1_/20260913T183720_920345Z
exec 9>"$JOB/.cpu_execution.lock"
flock -n 9 || { echo 'This job is already running'; exit 73; }
if [[ -e "$JOB/CPU_EXIT_STATUS.json" ]]; then
  echo 'This job has already executed; inspect outputs instead of rerunning.'; exit 74
fi
export CUDA_VISIBLE_DEVICES=""
export PYTHONDONTWRITEBYTECODE=1
export MPLBACKEND=Agg
cd /home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908/downstream_CPU_L1L2_v2_1_release_compat_20260914
date -Is > "$JOB/STARTED.txt"
set +e
/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python -B /home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908/downstream_CPU_L1L2_v2_1_release_compat_20260914/scripts/run_downstream.py --config /home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908/production_launches/phase1_/20260913T183720_920345Z/configs/cpu_config_APPROVED.yaml --execute --modules cohort counts markers figures dotplot
rc=$?
set -e
printf '{"exit_code":%s,"completed_at":"%s"}\n' "$rc" "$(date -Is)" > "$JOB/CPU_EXIT_STATUS.json"
exit "$rc"
