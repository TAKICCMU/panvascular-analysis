# Macrophage L2 CPU Round1 — 人工审核材料

状态：PAUSED_FOR_MACROPHAGE_L2_ROUND1_REVIEW
输入：冻结 level1_annotation == Macrophage；210,967细胞；全基因raw counts完整保留。
母图谱只读，未修改level1_annotation、全局坐标、benchmark或历史PROJECT_STATE。
本轮各分辨率计算簇数：{0.4: 8, 0.8: 17, 1.2: 21}。基线r=0.8仅为审核锚点，不是已确认最佳分辨率。
所有FindMarkers分辨率均使用各自的分组，不把r0.8的marker套给另一个分区。

## 注释规则
1. 经典marker + FindMarkers + 多基因组合 + 来源分布，审核后采用marker + 状态 + macrophage命名。
2. 不预设必须出现多少状态；也不因两个簇都表达同一宽程序就自动合并。
3. 若审核确认是其他谱系或身份不能可靠判定，后续仅将level2_annotation设为Unknown；level1_annotation保持Macrophage。
4. HLA、IFN、细胞毒、吞噬、ECM等共享程序或单个基因不能单独证明其他谱系。边界表是描述性结果，不是自动Unknown分类器。
5. 本轮不创建正式level2_annotation。所有候选表仍为PENDING，不自动命名，不将全部细胞标Unknown。
6. 不删除细胞，不重跑Scrublet，不回填L1，不做CoVarNet、疾病关联、轨迹或细胞通讯。

## 输入及方法
- X：本谱系全部基因raw counts；layers[log1p]：从独立counts工作副本重新CP10k + log1p。
- Pearson residual HVG=3000，theta=100，沿用审计补丁，每sample clip=sqrt(n_sample)。
- 本谱系中只有1个细胞的sample仅不参与HVG方差估计；这些细胞保留在PCA/Harmony/所有后续结果。
- Scale max_value=10；CPU PCA=50；Harmony theta=2.0，明确device=cpu；neighbors=30；Leiden igraph；UMAP min_dist=0.3。
- 仅sample_uid用于Harmony；age/sex/tissue/disease_context等只用于描述，不指导分群。
- FindMarkers：Scanpy 1.11.5 Wilcoxon（tie_correct=True），全基因中至少3个细胞检出的基因，one-versus-rest；不是仅HVG或预设marker。
- logfoldchanges沿用Scanpy基于平均log表达的近似值。P值来自聚类后同批细胞，不是独立验证或donor级疾病推断。
- 所有基因测试结果保留；仅marker候选图使用logFC>0.25、pct_in>=0.10、pct差>=0.05、FDR<0.05，并避开MT/RPL/RPS等泛表达锚点。这些条件不参与聚类/自动注释。
- DotPlot颜色：全部细胞平均log1p；点：raw counts>0比例。相对图0–1为每基因内缩放，不是绝对表达。
- 高类别数来源的UMAP仅展示前25类+Other，完整来源保留于表格。
- 暂未推断donor缺失值；来源集中、少数样本均为审核线索，不是自动删除条件。

## 需要共同回答
哪些簇有清楚的经典程序？哪些FindMarkers支持保留新的细粒度状态？哪些仅为应激或来源差异？
哪些存在另一谱系的完整程序而应设L2 Unknown？共享程序是否被误当成另一细胞类型？
不要从UMAP间距或marker数量直接推断细胞来源、转分化、功能或准确率。

原始冻结母文件SHA256：8d36e05a8eab608d330a14531f7651a28b526817b5723a6da04aae98a09cda1b
