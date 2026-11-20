# 冻结 benchmark 输入核验

PASS：315,949 个 Primary 细胞、247 个 sample_uid、2000 个 HVG 及其顺序、30 维原 PCA、182,968 个可信细胞评价集合均与历史冻结记录一致。

逐细胞比对 cell_id、sample_uid、原标签、reliable_L1、supervision_label、trusted 与独立保存的标签和评价清单；文件 SHA256 与历史 manifest 一致。表达源为冻结 primary.h5ad；本轮仅引用，不从全量 atlas 重建。

整合 seed=0，评分 seed=0；官方指标专用邻居数 15/50/90。主分数不做跨方法 min-max，biology/batch 权重 0.6/0.4。

旧 HVG/PCA 原样保留，包括此前报告的 38 个基因差异；本轮不实施 clip 更新。

详见 audit/INPUT_VALIDATION.json 和 historical_hashes.json。
