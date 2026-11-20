# Cell 2: 启动其他13个谱系的独立CPU批处理，不改动巨噬细胞对象/结果。
from pathlib import Path
import importlib.util

OTHER_L2_ROOT = Path(
    "/home/anzhen_vhd/scRNAseq/codex/"
    "panvascular_atlas_rebuild_v2_20260907/"
    "06_reanalysis_takeover_20260908"
)
OTHER_L2_CODE = OTHER_L2_ROOT / "Other_L2_CPU_overnight_code"
_other_batch_spec = importlib.util.spec_from_file_location(
    "other_l2_batch_driver", OTHER_L2_CODE / "other_lineages_cpu_batch.py"
)
other_l2_batch = importlib.util.module_from_spec(_other_batch_spec)
_other_batch_spec.loader.exec_module(other_l2_batch)

OTHER_L2_RUN, OTHER_L2_PID = other_l2_batch.launch_other_lineages(
    source=OTHER_L2_ROOT / "panvascular_CPU_theta2_level1_frozen_v1.h5ad",
    max_workers=2,
    total_threads=40,
    resume_run=globals().get("OTHER_L2_RUN"),
    macrophage_run=globals().get("MACROPHAGE_L2_RUN"),
)
print("其他谱系输出：", OTHER_L2_RUN)
print("调度器日志：", OTHER_L2_RUN / "scheduler.log")
print("当前巨噬细胞变量未清除；无需保持本cell一直显示输出。")
