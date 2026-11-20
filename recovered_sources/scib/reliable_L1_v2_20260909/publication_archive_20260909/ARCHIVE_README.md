# 论文与作图冻结存档

来源：reliable_L1_v2_20260909。benchmark_all_metrics_full_precision.csv/tsv 为统一宽表，保留源文件完整数值精度，按主分析总分排序。请用此表作图，避免从聊天表格抄录四舍五入值。原始指标、评分、参数、标签摘要、embedding 路径和图件附后；不覆盖历史结果。

全部 8 种方法、10 个官方指标。指标为 scib-metrics 官方返回的方向统一分数（越高越好），没有跨方法 min-max；cLISI/iLISI 并非未转换的原始 LISI 指数。Biology 为五项 biology 指标均值，Batch 为五项 batch 指标均值，Overall=0.6 Biology+0.4 Batch。补充 min-max 表不能混用为主分析。

Primary 共315949细胞，统一评价182968个可信L1细胞；132981个Unknown不作为评价细胞类型。scANVI与scPoli使用同一可信L1监督，L1 biology不是独立验证。PCA为未整合基线；CellHint未纳入本轮统一embedding排名。K10/20/30未运行；采用官方指标各自默认设置。完整设计见README.md和evaluation_manifest.json。
