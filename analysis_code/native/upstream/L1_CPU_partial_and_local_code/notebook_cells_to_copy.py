# Copy sections into your EXISTING notebook.
# %% 1
from pathlib import Path
import importlib.util

if "adata" not in globals() or "markers" not in globals():
    raise RuntimeError("请在仍有全量adata和markers的原Notebook中执行；不要使用轻量DotPlot对象。")

ROOT = Path(
    "/home/anzhen_vhd/scRNAseq/codex/"
    "panvascular_atlas_rebuild_v2_20260907/"
    "06_reanalysis_takeover_20260908"
)
MODULE_CANDIDATES = [
    Path.cwd() / "l1_cpu_local_refinement.py",
    ROOT / "L1_CPU_partial_and_local_code/l1_cpu_local_refinement.py",
]
MODULE_FILE = next((p for p in MODULE_CANDIDATES if p.is_file()), None)
if MODULE_FILE is None:
    raise FileNotFoundError("请先把l1_cpu_local_refinement.py上传到当前目录或ROOT/L1_CPU_partial_and_local_code/。")

spec = importlib.util.spec_from_file_location("l1job", MODULE_FILE)
l1job = importlib.util.module_from_spec(spec)
spec.loader.exec_module(l1job)

# 使用母图谱登记的工作文件夹，不依赖可能被绘图代码覆盖的 OUT 变量
WORK_PARENT = Path(adata.uns["output_directory"])
L1_RUN = l1job.prepare_partial_L1(
    adata,
    markers=markers,
    work_parent=WORK_PARENT,
    patch_dir=ROOT / "annotation/HVG_clip_review_20260910",
)
print("当前L1数量：")
print(adata.obs["L1"].value_counts())
print("本次全部结果目录：", L1_RUN)

# %% 2
# 阻塞式运行：CPU-only，2个分支并行×每分支20线程
# 正式注释不会自动回填；当前adata保持第一步的初步L1。
LOCAL_STATUS = l1job.run_local_batch(
    L1_RUN,
    max_workers=2,
    total_threads=40,
)
LOCAL_STATUS

# %% 3
# 查看文件位置；不需要同时把13个分支全载入内存
print("母图谱：", L1_RUN / "CPU_theta2_L1_partial_full.h5ad")
print("离线图集：", L1_RUN / "GLOBAL_LOCAL_L1_REVIEW.html")
print("审核包：", L1_RUN / "Local_L1_review_materials.zip")

# 需要检查某个分支时，取消下一行的注释：
# C0 = l1job.load_local(L1_RUN, "C0")
# C5 = l1job.load_local(L1_RUN, "C5")
# print(C5)
# print(C5.obs["local_leiden_r0p5"].value_counts())

# 重启kernel之后，可以明确记录好的目录重用/恢复：
# L1_RUN = Path("上面打印出的本次目录")
# adata = sc.read_h5ad(L1_RUN / "CPU_theta2_L1_partial_full.h5ad")
# 不要先重新执行prepare_partial_L1创建另一个相同初步版本。
