# CPU-Harmony 补充分支

入口：Harmony_CPU_vs_GPU_report.md；离线图集 review.html；全指标 benchmark_all_metrics_full_precision.tsv。输入均为只读历史路径链接。代码与适配差异位于 code/，全部审计位于 audit/。仅新增 CPU Harmony；原 benchmark 不重命名、不覆盖。

指标邻居搜索与 UMAP 为 CPU；Leiden 未执行，原指标使用 KMeans。评分 JAX 保持历史 GPU 后端以保持一致。
