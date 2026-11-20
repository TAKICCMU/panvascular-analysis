# 在原Jupyter Notebook的一个cell中运行；不要重启kernel，不覆盖adata变量。
from pathlib import Path
import importlib.util

L2_FINAL_ROOT = Path(
    '/home/anzhen_vhd/scRNAseq/codex/'
    'panvascular_atlas_rebuild_v2_20260907/'
    '06_reanalysis_takeover_20260908'
)
L2_FINAL_PACKAGE = L2_FINAL_ROOT / 'L2_final_marker_parent_freeze'
L2_FINAL_SCRIPT = L2_FINAL_PACKAGE / 'finalize_level2.py'

if not L2_FINAL_SCRIPT.is_file():
    raise FileNotFoundError(f'请先上传并解压完整代码包：\n{L2_FINAL_PACKAGE}')
if 'L2_UPDATED_H5AD' not in globals():
    raise RuntimeError(
        "请先明确设置上一轮全量L2工作文件：\n"
        "L2_UPDATED_H5AD = Path('/实际完整路径/上一轮全量L2工作文件.h5ad')\n"
        '不要使用冻结L1-only文件，也不要按修改时间猜测最新版。'
    )
L2_FINAL_SOURCE = Path(L2_UPDATED_H5AD).resolve()
print('本次将先备份、再覆盖的上一版全量工作文件：', L2_FINAL_SOURCE)

# 普通内存对象可以保留；若同一待覆盖文件以backed方式打开，先关闭句柄。
for _object_name, _object in list(globals().items()):
    if _object_name.startswith('_'):
        continue
    if getattr(_object, 'isbacked', False) is True:
        _filename = getattr(_object, 'filename', None)
        if _filename and Path(_filename).resolve() == L2_FINAL_SOURCE:
            raise RuntimeError(
                f'{_object_name}以backed模式打开了待替换文件。'
                f'请先执行 {_object_name}.file.close()，并关闭/删除这个backed对象后再运行。'
            )

_final_l2_spec = importlib.util.spec_from_file_location('l2_final_freeze_job', L2_FINAL_SCRIPT)
l2_final_freeze_job = importlib.util.module_from_spec(_final_l2_spec)
_final_l2_spec.loader.exec_module(l2_final_freeze_job)

L2_FINAL_RELEASE_RUN, L2_FINAL_H5AD = l2_final_freeze_job.run_from_notebook(
    target_h5ad=L2_FINAL_SOURCE,
    root=L2_FINAL_ROOT,
    review_run=globals().get('L2_NEXT_RUN'),
    threads=4,  # 本轮只核验及保存，不运行新科学分析；40核无需全部占用。
)
L2_UPDATED_H5AD = L2_FINAL_H5AD

print('\n最终正式全量h5ad：', L2_FINAL_H5AD)
print('备份、审批与核验记录：', L2_FINAL_RELEASE_RUN)
print('正式ontology：', L2_FINAL_ROOT / 'atlas_release/L1L2_final_v1/level2_ontology.tsv')
print('当前内存的旧adata没有自动刷新；请勿再用旧对象覆盖已冻结文件。')
