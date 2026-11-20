# 在原Macrophage Notebook下面运行；不清空或覆盖已有adata/Macrophage_L2。
from pathlib import Path
import importlib.util

L2_UPDATE_ROOT = Path(
    '/home/anzhen_vhd/scRNAseq/codex/'
    'panvascular_atlas_rebuild_v2_20260907/'
    '06_reanalysis_takeover_20260908'
)
L2_UPDATE_PACKAGE = L2_UPDATE_ROOT / 'L2_followup_metadata_v4_code'
L2_UPDATE_SCRIPT = L2_UPDATE_PACKAGE / 'run_L2_followup_and_update.py'

# 使用已明确冻结的全量母图谱。若另有更新的完整版本，请明确填入真实路径，
# 不按修改日期猜文件。程序会核验其矩阵与冻结基线相同，且不覆盖更新的冲突标签。
L2_SOURCE_FULL_H5AD = L2_UPDATE_ROOT / 'panvascular_CPU_theta2_level1_frozen_v1.h5ad'

if not L2_UPDATE_SCRIPT.is_file():
    raise FileNotFoundError(f'请先上传并解压整个代码包到：{L2_UPDATE_PACKAGE}')
_l2_update_spec = importlib.util.spec_from_file_location('l2_targeted_update_job', L2_UPDATE_SCRIPT)
l2_targeted_update_job = importlib.util.module_from_spec(_l2_update_spec)
_l2_update_spec.loader.exec_module(l2_targeted_update_job)

L2_TARGETED_RUN, L2_UPDATED_H5AD = l2_targeted_update_job.run_from_notebook(
    root=L2_UPDATE_ROOT,
    source=L2_SOURCE_FULL_H5AD,
    macrophage_run=globals().get('MACROPHAGE_L2_RUN'),
    other_run=globals().get('OTHER_L2_RUN'),
    total_threads=40,
    max_workers=2,
    resume_run=globals().get('L2_TARGETED_RUN'),
    run_targeted=True,
)
print('新的全量h5ad：', L2_UPDATED_H5AD)
print('待讨论结果：', L2_TARGETED_RUN / 'GLOBAL_TARGETED_REVIEW.html')
print('分谱系轻量审核包：', L2_TARGETED_RUN / 'reviews')
