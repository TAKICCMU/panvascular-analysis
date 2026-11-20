# CPU 从 raw counts 开始的 OmicVerse 工作本
入口：CPU_from_raw_counts_omicverse.ipynb。先重启 omicverse kernel，再顺序执行。
代码已做语法检查和已安装接口核对；没有执行全量重算。
现有 full_atlas_CPU_review.h5ad 不包含 counts；从 uns 的 raw_counts_connection 加载 counts。
沿用项目已经验证的 Scanpy 1.11.5 逐样本 clip 补丁（依赖服务器 PATCH_DIR），不修改 site-packages。
OmicVerse 负责 CPU 模式、scale、PCA、Harmony 和 UMAP。基础 normalize/log1p、显式CPU邻居和igraph Leiden使用Scanpy接口，保持既往流程。
输出到新建时间戳目录。raw counts 保留外部连接，全基因 log1p 保存供 DotPlot；工作 AnnData 的 X 只含 2000 HVG。
线程上限40（原CPU基线BLAS32、neighbors48），因此不保证按位复现；生物学预处理及显式算法参数保持原CPU基线。
