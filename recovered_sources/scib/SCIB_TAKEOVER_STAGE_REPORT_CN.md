# scIB 接管阶段报告

本报告为阶段状态，**尚未完成九方法比较，不能据此锁定最终方法**。

## 已接管的分析口径

- 独立 study-wise marker 预注释；`batch_key=sample_uid`，不使用 HCA generic sample_id。
- 主集 315,949 cells；pilot 114,833 cells；247 samples，29 source study files、26 canonical studies。
- Primary 每 sample 最多1,500，pilot最多500；ASW/PCR计算子集也只按sample抽样。
- K=10/20/30等权；Leiden resolution=0.2/0.5/1.0/1.5，分别报告固定网格最佳NMI与ARI。
- 全方法共同按K缩放，0.4 batch + 0.6 biology；同时报告共同图指标与原流程metric分组敏感性。
- scANVI/scPoli/CellHint使用输入标签，统一表中保留该标志。用于训练/协调的同一标签也参与生物学评分，不能将高分解释为独立标签预测验证。

## 已完成与仍运行

8方法既有raw metrics完整；新的统一评分脚本对其重算已完成，结果在 `primary_provisional_unified/`。CellHint尚缺，评分与排名均为暂定。

| Method | Unified applicable score | Common graph score |
|---|---:|---:|
| scANVI | 0.8511 | 0.8778 |
| scPoli | 0.7272 | 0.8386 |
| Harmony theta1 | 0.7261 | 0.8035 |
| scVI | 0.7254 | 0.7749 |
| PCA | 0.6277 | 0.6767 |
| SCALEX | 0.6040 | 0.5407 |
| Scanorama | 0.3613 | 0.2974 |
| LIGER | 0.1808 | 0.2188 |

CellHint原pilot meta3/meta5继续运行，接管时为207/247与206/247；不终止、不修改、不重复启动。独立primary共享harmonization已提前启动，因为meta-neighbor只影响随后integrate，两个候选配置的harmonization完全相同。

共享harmonization局部缓存同一sample的15次相同confusion-matrix请求，返回深拷贝。360-cell/4-sample小样本比较证实：原native和缓存版最终relation、逐细胞reannotation、每步选择的两个阈值完全相同；不改grid、不改顺序、不改filter。新任务限制4核，与现有两个10-thread pilot合计不超过24核预算。

控制器待pilot graph自然完成后依次：冻结副本→原始邻居审计→共同scIB评估→选择meta3/meta5→等待共享primary harmonization→仅运行所选primary integrate→评估→合并九方法统一/共同图指标评分。日志、状态均写隔离目录。

## 查出并处理的问题

1. 旧评分分轨缩放，不能跨轨比较；接管版所有方法在同一K共同缩放。
2. 旧稿流程曾将已经scaled的cLISI再次反转；当前scib_metrics `clisi_knn(scale=True)`直接使用，higher-is-better。
3. 原CellHint导出遇邻居不足会补self；接管版不允许补齐。原pilot导出须先通过无self补位、无重复、有限非负距离、有效索引审计，否则明确失败并停在该依赖点。
4. 现有embedding `NeighborsResults`输入含self，而旧CellHint导出先去self。独立评估适配层在已通过上述无补位审计的真实邻居前明确插入一个标准self，再取K-1真实邻居，使K含义一致。这与邻居不足时用重复self凑数不同。
5. 新评分进一步禁止重复method/config/K行、缺K配置参与缩放，并验证每个NMI/ARI最佳值确实来自固定resolution网格。运行失败不能由skipna掩盖；CellHint仅ASW/PCR三项作为结构性NA。
6. CellHint原生n_neighbors未显式设置，247samples下默认741。接管版保留这一算法参数，保存native graph，再取距​​离排序的K子集评分，不中途为加速更改邻居定义。
7. 旧metric cache记录输入cell-order、PCA、metadata hash及embedding size，却未普遍记录embedding内容hash。可核验既有签名并冻结当前内容hash，但不能声称具备过去未保存的内容哈希证据。

## 与原稿流程的区别

沿用原稿多K、统一缩放和40/60思想；本次母图谱已去除重复来源，标签由旧Harmony注释改为独立marker预注释，成员抽样改为sample-only，因此旧稿分数与本轮分数不能直接当作同一试验的升降比较。

当前稿权威来源为 `/Users/taki/Desktop/泛血管最终投稿版与展示/G09_release_patch/clean_editorial/Supplementary_Tables.xlsx` 的 Supplementary Table 27（提取为 `audit/original_article/SuppTable27_original.tsv`），包含6方法、3个K。已逐项验证其18行×34列内容与历史 `Panvas_Nature_V1.9.2_submission_package/Source_Data/ExtendedData2_scIB_full_benchmark.csv` 数值完全一致（容差1e-12；见 `audit/G09_VS_V1_9_2_TABLE27_EQUIVALENCE.json`），因此上述比较适用于G09现稿。G09报告跨K均值：Harmony0.7797、scANVI0.6154、scVI0.5607、SCALEX0.4927、PCA0.4676、Scanorama0.2075。历史V1.5 Figure1另曾为各方法挑选最佳K；该历史展示方式不冒充G09当前ED2。接管流程使用三个K等权。

G09发布代码notebook另含0.5/0.5汇总语句，与G09 Table27数值实际0.6biology+0.4batch不一致。本接管以当前成稿/表27的40/60为准，并记录代码与成稿不一致，不能称原代码全程一致。

对**同一旧数据表只修复cLISI方向**，其余指标、分组、权重、K保持原样，跨K分数分别变为Harmony0.8026、scANVI0.7354、scVI0.6059、PCA0.5767、SCALEX0.3727、Scanorama0.1857。该计算仅证明旧评分方向错误的影响，不能挽救旧数据的重复身份及标签来源问题。完整旧/新对照见 `OLD_MANUSCRIPT_VS_CURRENT_BENCHMARK.tsv`，每行均明确标记跨run分数不可直接比较。

现统一主分组把graph connectivity归入batch（与常用scIB分组一致），并额外保留原稿将其归入biology的历史分组敏感性评分；共同图分数只使用各方法均可计算的图指标。最终九方法结果还可能因CellHint加入而改变统一缩放边界，当前八方法分数不是最终九方法分数。

## 路径

- 服务器主输出：`/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908/scib`
- 服务器状态：`TAKEOVER_STATUS.json`
- 日志：`logs/controller.log`、`logs/cellhint_primary_shared.log`
- 原主线输入与结果仅只读复用；所有新写入均在本接管目录。


## 既有 Harmony 收敛核查（2026-09-08）

Primary 315,949-cell theta1 的原日志明确在第15轮达到 OmicVerse 原生收敛判据（阈值0.0001，cap20，sample_uid，seed0）；pilot114,833-cell theta1在第17轮达到。因此当前8方法的所选Harmony配置并非cap20截断。此处是原生日志收敛报告，既有cache未保存objective完整轨迹，不能替代绝对改变量审计；full117万细胞需31轮与这两个subset不同，不应外推。未为此补跑或修改冻结cache。完整逐配置证据见 audit/harmony_convergence_existing/EXISTING_HARMONY_CONVERGENCE_AUDIT.json。
