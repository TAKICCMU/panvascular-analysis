# Neutrophil L2 CPU Round1 — 候选审核
输入：冻结 level1_annotation == Neutrophil; 18,509 个细胞，全部基因保留。
母图谱和巨噬细胞分支未修改；未应用正式level2_annotation。
各分辨率簇数：{0.4: 7, 0.8: 9, 1.2: 13}；基线0.8仅为审核锚点。

1. 每个分辨率均使用自身FindMarkers与自身sample分层marker表；不得混用编号。
2. 命名风格：特征marker + 有证据支持的状态描述 + 细胞类型。名称/亚群数不预设。
3. 经典marker和FindMarkers共同支持人工审核，不以最高module score自动分配。
4. 经审核确认属于其他谱系或身份不明时，只把新的level2_annotation记为Unknown；冻结L1不动。
5. 单个标记、HLA、IFN、应激、细胞毒或ECM共享程序本身不等于另一细胞类型。
6. 少于2细胞的sample只不参与HVG方差估计，细胞仍保留在后续分析。
7. 未重新运行Scrublet，无细胞删除，无疾病推断、CoVarNet、CellChat或轨迹。

流程：raw工作副本→Pearson residual HVG 3000，逐sample clip=sqrt(n)
→全基因CP10k+log1p→临时HVG矩阵scale（max=10）→PCA 50
→CPU Harmony(theta=2)→neighbors 30→Leiden 0.4/0.8/1.2→UMAP。
FindMarkers：全基因中至少3细胞检出的基因，Scanpy Wilcoxon，layer=log1p，use_raw=False。
FDR/效应筛选仅用于画特征marker图；不改聚类或L1，不用来自动删小群。
DotPlot：含零细胞的平均log1p；点大小=raw counts>0；相对图是每基因内0–1缩放。
临床字段、donor缺失、NR及PROVISIONAL均保留；覆盖样本数不等于独立患者数。
L1筛选后同一细胞的全基因log1p与母图谱本应可一致，本程序按同一公式在工作副本重算。
单study状态不自动删除，也不能宣称跨研究稳定；采样/疾病不能指导聚类。
所有status均为待审核，不是已验证的正式L2亚型。
冻结输入SHA256: 8d36e05a8eab608d330a14531f7651a28b526817b5723a6da04aae98a09cda1b
