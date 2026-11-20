# 其余13谱系L2：CPU独立批处理 + 样本临床元数据CSV导出

## 本任务不会做的事情
- 不重跑、覆盖或重命名巨噬细胞数据；不读取Notebook内的巨噬细胞AnnData。
- 不清空kernel、不修改当前adata、Macrophage_L2、MACROPHAGE_L2_RUN等变量。
- 不改冻结母h5ad、原始counts、L1、全局图/坐标、旧PROJECT_STATE或既有benchmark。
- 排除Macrophage（已完成首轮）、Schwann cell（用户要求不做L2）、L1 Unknown。
- 不自动贴正式L2标签；若后续人工审核确认其他谱系/身份不明，只在L2设Unknown，不回改L1。
- 共享程序、单个marker、HLA、细胞毒/IFN/ECM等不自动触发Unknown。

## 如何运行
完整解压此文件夹至：
`/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908/Other_L2_CPU_overnight_code`

继续原omicverse Notebook即可，不要为了本代码重启kernel。
将 `RUN_IN_CURRENT_NOTEBOOK.ipynb` 中的前两个cell依次复制到当前巨噬细胞代码后运行。
第一个cell快速导出元数据；第二个cell启动独立的服务器CPU调度进程后返回。
第三个cell仅用于查看进度，可重复执行，不重启任务。

执行目录：`ROOT/L2_CPU_others/Other_lineages_round1_时间戳/`。
运行时第一次会校验母文件SHA256及HVG补丁/依赖，可能需等待。初步启动不代表计算已完成。

## 覆盖范围
13个分支合计903,398个细胞：
Fibroblast 251371；SMC 243346；Endothelial 104415；T 156349；Fibromyocyte 29347；
Monocyte 31046；NK 25719；Dendritic 10556；B 10496；Plasma 4802；Mast 12929；
Neutrophil 18509；Pericyte 4513。

各自：全基因raw counts→Pearson residual HVG 3000（每sample clip=sqrt(n)）→全基因CP10k+log1p
→scale(max10)→PCA50→CPU Harmony(theta2)→neighbors30→Leiden0.4/0.8/1.2→UMAP。
每个分支只运行一次HVG/PCA/Harmony/邻居图；三个Leiden分辨率共用。
每档分辨率分别做自己的全表达基因Wilcoxon FindMarkers、经典marker图、特征marker图和sample分层marker表。
单细胞样本仅不参与HVG方差估计，不从其他步骤删细胞。无Harmony单batch则跳过并记录。
FindMarkers至少3细胞检出过滤只作用于被检验基因，完整raw矩阵保留所有基因。
不把gene-wise相对DotPlot值写成真实均值；同时保存原始均值和raw检出率。

## 计算和继续运行
- 默认2个谱系并行，每个20线程，总预算40线程；CPU-only（显式device=cpu，子进程CUDA隐藏）。
- 输入通过一次完整SHA256校验后，子任务检查文件状态和验证票据；批次结束再完整校验一次。
- 每个分支独立目录/UID、日志和断点。一个失败不取消其他分支。
- 同run重新启动：已完成且h5ad校验通过的分支不重算；未完成分支复用其raw/embedding/FindMarkers缓存。
- 默认记录精确ACTIVE_BATCH指针，误重复运行启动cell不会自动创建第二套相同任务。
- 不保证一夜完成：Wilcoxon、HDF5压缩、磁盘速度、服务器负载都会影响耗时。
- 关闭浏览器一般不影响已脱离Notebook输出管道的进程；服务器重启、作业系统清理/超时或强制终止仍会中断。
- 当前kernel中的巨噬细胞变量仅在kernel存活期间保留；不要点Restart/Shut Down。
- 当可用内存小于48GiB或磁盘小于12GiB时，不再启动新的分支，已有任务继续。释放空间后可继续排队。
- 大谱系保存raw、embedding检查点和最终审核h5ad，可能需要数十到上百GiB额外磁盘；不是13个1MB小文件。

## 结果位置
主表：`all_other_L2_cluster_candidates.tsv`（全为PENDING，无正式标签）。
进度：`batch_progress.csv`、`BATCH_STATUS.json`、`scheduler.log`。
统一入口：`OTHER_L2_REVIEW_INDEX.html`（服务器目录内打开，保留相对目录）。
统一包：`Other_L2_all_review_materials.zip`。
该统一ZIP含13个独立已压缩谱系审核包；再解压需要的谱系ZIP即可打开其自带图片HTML。
单分支：`lineages/<slug>/LINEAGE_L2_ROUND1_OFFLINE.html` 和 `lineage_L2_round1_review_materials.zip`。
表达：`lineages/<slug>/data/lineage_L2_round1_review.h5ad`。

局部h5ad中：X仍是该谱系全部基因raw；layers[log1p]为新算log1p；
obs[level1_annotation]保持父类，正式level2_annotation尚未写入。
obs[lineage_L2_r0p4/r0p8/r1p2]为计算簇，完整UID含谱系前缀。
obsm[X_pca_lineage_L2/X_harmony_lineage_L2/X_umap_lineage_L2]为局部结果。
每档分辨率目录分别有cluster_sample_marker_mean_log1p.tsv.gz和检出率、FindMarkers以及来源分布。

## 元数据CSV
第一个cell仅用h5py/anndata read_elem读取冻结h5ad的obs，不读取X/log1p/图，也不修改h5ad。
主文件：`sample_metadata_review_时间戳/sample_study_clinical_metadata.csv`，每个sample_uid一行（预期247行）。
- 保留原study、donor、age、sex、tissue、condition、disease、disease_context及来源字段。
- 同sample字段冲突显示[CONFLICT]，不取第一条、不投票、不平均年龄、不猜sex代码。
- 全部原值、出现细胞数、缺失和PROVISIONAL保存在values_long/wide_detail与冲突表。
- 联合组合表保存原study-donor-age-sex-背景组合，不把不同donor的属性交叉配对。
- disease/condition是记录中的粗分组；细病种仍看disease_context，不从study名字推断。
- 新信息请写verified_*和verification_reference/note。人工编辑CSV不会自动回填h5ad。
- 核对CSV期间请不要覆盖正在被分析的冻结h5ad。之后另行建立metadata更新版本。

## 科学解释边界
该轮生成marker和候选状态证据，不追求固定状态数，也不预先合并所有共享宽程序的簇。
经典marker面板不是分类器；FindMarkers P值为聚类后细胞层面描述，不是donor层面疾病检验。
来源覆盖和donor键数量不是独立验证；PROVISIONAL/缺失如实保留并单列。
配套marker清单承接前面讨论的方向，但不是声称本数据中必须出现的亚群。

## 测试范围
已完成语法、13谱系/903398计数、独立UID、marker清单、稀疏均值/检出率、metadata冲突保留测试。
使用假数据worker做了真实子进程调度测试：13分支、一次失败、其余继续、断点跳过12个成功任务、来源不变和打包。
本环境缺少OmicVerse/Scanpy/AnnData，未执行服务器全量科学计算，也未对真实h5ad运行导出；实际以日志核验。
