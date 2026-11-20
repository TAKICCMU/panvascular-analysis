# Cell 1: 只读导出 study/sample/临床信息；不清空kernel，不读取表达矩阵。
from pathlib import Path
import importlib.util

OTHER_L2_ROOT = Path(
    "/home/anzhen_vhd/scRNAseq/codex/"
    "panvascular_atlas_rebuild_v2_20260907/"
    "06_reanalysis_takeover_20260908"
)
OTHER_L2_CODE = OTHER_L2_ROOT / "Other_L2_CPU_overnight_code"
OTHER_L2_SOURCE = OTHER_L2_ROOT / "panvascular_CPU_theta2_level1_frozen_v1.h5ad"

_other_meta_spec = importlib.util.spec_from_file_location(
    "other_l2_metadata_export", OTHER_L2_CODE / "sample_metadata_export.py"
)
other_l2_metadata = importlib.util.module_from_spec(_other_meta_spec)
_other_meta_spec.loader.exec_module(other_l2_metadata)

SAMPLE_METADATA_EXPORT = other_l2_metadata.export_sample_metadata(
    source=OTHER_L2_SOURCE
)
print("主CSV：", SAMPLE_METADATA_EXPORT / "sample_study_clinical_metadata.csv")
