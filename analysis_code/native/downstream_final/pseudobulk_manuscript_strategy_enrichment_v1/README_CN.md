# Pseudobulk主分析选择与富集

用户2026-09-16指定：pseudobulk_manuscript_strategy_20260916_v1为主分析；三大背景扩展版为补充。保留纠正后的样本/供者信息，不恢复旧metadata错误。

复用79项全基因edgeR结果，不重算DE。原稿通路主线为至少2个疾病背景的DL随机效应汇总，SE沿用abs(logFC)/sqrt(F)近似，不称为精确QL标准误；原候选靶点REML结果不修改。用全基因meta Z做GSEA，另做79个背景内支持性富集。ORA按全基因FDR<0.05及方向选基因，全部可检验非空通路参与BH后才显示overlap>=3。重复基因symbol整体排除并审计。fgsea原生无法估计项保留NA。

使用冻结MSigDB 2026.1.Hs的GO BP/CC/MF、KEGG LEGACY/MEDICUS、REACTOME。所有已有结果只读；新结果均保存在本目录。ANALYSES.json提供状态名称与结果目录映射，PROGRESS.json/COMPLETE.json分别表示进度与经过自动核验的完成。
