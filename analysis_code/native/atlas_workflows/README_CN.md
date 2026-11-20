# metadata v4 donor补全修订（2026-09-13）

本包与此前L2_followup_metadata_v4_code使用同一个目录名和主程序。
科学分析代码未改；只更新8行metadata的donor_uid、donor_id、biological_sample_uid，
以及对应缺失统计、外部修改记录和软件包校验清单。
用户明确确认“GSE166676和GSE184073内的样本都来自不同donor”。
内部donor ID由study+GSM命名，表示本次确认的身份分组，不声称是新发现的公共Individual accession。
主CSV仍为247行、30列，没有新增备注或不确定性字段。

## 上传位置
ROOT = /home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908
请先备份服务器旧 L2_followup_metadata_v4_code 文件夹，再上传本包同名文件夹。
不要在旧任务仍运行时替换其代码或输入。
CSV必须与本包 BUNDLE_CHECKSUMS.json、inputs/expected_summary.json 配套；只替换CSV会触发原程序完整性检查。

## 重新执行
原Notebook入口继续可用，但本次新metadata必须使用 resume_run=None。
原程序resume_run指向旧目录时会直接返回旧h5ad，不会应用新CSV。
若之前定点计算已完成，本次可用 run_targeted=False，只创建新全量h5ad而不重复做pairwise。
若定点计算尚需执行，则保留 run_targeted=True。
既有run目录、旧输入h5ad及历史SHA256不得修改。
这次没有在用户服务器执行，也没有修改任何h5ad；原L1/L2标签和sample分组未改。


---

# Macrophage＋其余13谱系：L2部分回填、定点复核、metadata v4

## 本轮执行范围

已上传的巨噬细胞和其余13谱系的表格是本轮依据，不能读取用户服务器kernel的实时状态。
这次代码从已保存的完整母h5ad与各L2首轮h5ad读取，不清空/替换已有Notebook对象。

用户已请求应用清楚的名称与明确Unknown；本包将上一轮 PRIORITY_NAMING（Macrophage）
和 NAME_CANDIDATE（其他谱系）作为本轮工作命名审批，共75项、604,732细胞。
明确建议的Unknown共21项、39,172细胞。它们是当前工作版本，非最终L2冻结或稳定网络节点。
所有 NAME_REVIEW / STATE_PENDING / LINEAGE_REVIEW / PENDING_STATE 等保留 Pending，共470,461细胞。
L1 Unknown的57,918细胞在L2仍为Unknown；Schwann 2,832细胞为Not_assessed。
因此新全量L2列：604,732有工作名称；97,090为Unknown；470,461 Pending；2,832 Not_assessed。
不把尚未命名当成身份错误，不把共享基因程序自动当成另一谱系。

Macrophage按r1.2，其他谱系按r0.8执行这次已有审批。完整UID与成员哈希逐项校验。
已有r1.2的Fibroblast/Endothelial等候选另作复核；不同分辨率的标签不混合回填。
即使某个工作命名仍参加pairwise比较，也不会根据新P值自动删除或改名。

## metadata处理边界

1. 每个sample_uid始终保留为一条分析样本记录，不删除或合并任何细胞。
2. GSE237230按用户最后的配对表及本次GEO核对：
   GSM7597522/23=patient4；24/25=patient1；26/27=patient2；28/29=patient3。
   8条样本、4个donor，不将技术重复改写成不同人。
3. GSE210248 rep1/rep2维持两条样本记录、同一原始donor。
4. GSE213128沿用第一份GSA登记证据的5个明确Individual accession与年龄性别。
5. 年龄和性别保持既有证据：76个样本缺年龄、56缺性别。2026-09-13用户明确确认GSE166676的6个样本及GSE184073的2个样本各来自不同供者，本修订据此填入8个独立内部donor键。
   已核实同一donor的行之间也没有可传播的新年龄或性别值（本次核查为0）。
6. 年龄/性别/供者的来源与未解析说明放在外部audit CSV，不广播到obs作为标签或特征。
   obs中缺失为真正NA，不把NR/PROVISIONAL虚构成一个已知人，也不伪造独立donor。
7. vascular_bed/vascular_bed_group由现有tissue_detail的明确解剖名称标准化，不推断未报告的分支位置。
8. disease_context是研究/比较背景；source_diagnosis继承来源说明；donor_disease_context与
   tissue_disease_status只有直接证据才填写，不用研究主题补造对照者的个人诊断。
9. 新版本主要是清理与结构化，不声称已独立外部复核全部247个样本。
10. 与GSE237230/其他技术重复相关的细胞独立性风险仍在外部审计中保存。
    本轮不做丰度网络/pseudobulk推断，保留样本记录不等于证明独立生物重复。

## 使用

将整个 L2_followup_metadata_v4_code 文件夹上传到项目ROOT。
在现有Notebook里运行 RUN_ONE_CELL.py 的全部代码（Notebook版本也只有一个cell）。
默认读取ROOT/panvascular_CPU_theta2_level1_frozen_v1.h5ad。
确有较新完整文件时显式修改launcher的L2_SOURCE_FULL_H5AD，程序核验其表达/坐标与冻结基线一致。
已有变量 MACROPHAGE_L2_RUN / OTHER_L2_RUN 用来定位各分支；不存在时使用本次已审核RESULTS记录的
确切路径。路径迁移则显式传入，不搜索“日期最新”并默默套用错误对象。

先完成完整样本metadata与部分L2映射，另存一个新h5ad到ROOT根目录，随后运行CPU定点分析。
不重写冻结文件，不改level1_annotation，不重跑HVG/PCA/Harmony/Leiden/UMAP。
新文件命名：panvascular_CPU_level1_L2_partial_metadata_v4_<timestamp>.h5ad。
X/layers/var/varm/varp/obsm/obsp从输入完整复制，仅更新obs和必要uns；逐数据集流式哈希核验。
文件含全量raw counts和log1p，不是外部链接或轻量表达连接对象。
原Notebook的adata和Macrophage_L2不改变；返回L2_UPDATED_H5AD是新权威文件路径。

## 定点后续分析

- 复用已有分群/坐标，经典marker、候选marker与来源图以本次metadata重新汇总。
- 所有待定群选择一个同谱系最近的已命名中心作为表达对照（只选对照，不赋予身份）。
- 加上此前明确列出的126项pairwise比较；自动去重。另有已有高分辨率的定点比较。
- 每项直接比较采用全基因log1p表达的Wilcoxon、tie correction、BH校正，测试至少3细胞检出的基因。
- 输出两侧效应/检出比例、top100与候选图；原始全基因结果留服务器，不塞入轻量包。
- 同sample两侧各>=20细胞才进入描述性方向一致性表；无患者层面的显著性检验。
- 检查raw共同检出与多基因检出数；不按这些阈值自动判定Unknown。
- 汇总assay、文库复杂度、同sample深度分位和来源，以检查scRNA/snRNA及深度影响。
- 增加其他已冻结L1类别的marker均值参照；它们不进入局部聚类，也不是外部验证。
- 后续生成的是证据，不自动产生或应用新的名称。L2仍未冻结。

## 运行资源与恢复

2个CPU进程，每个20线程。启动子进程前隐藏GPU；不依赖当前kernel是否初始化CUDA。
使用同一Python环境，沿用已审核Scanpy 1.11.5。
至少预留一个完整h5ad副本空间及额外结果空间。大谱系pairwise可能耗时数小时，不能承诺短时完成。
完整母文件在计算前先保存成功，后续某分支失败不丢失已批准注释与metadata。
每个分支/每个pair都有日志与DONE记录。重新运行并传resume_run可复用通过校验的完成结果。
如果进程中断而L2_TARGETED_RUN未返回，从ROOT/L2_followup_request_*.result.json找确切run/h5ad路径。
无活动进程的异常中断可能留下SCHEDULER.lock；确认对应PID已退出后才移除该锁，不自动夺锁。
程序不创建新的后台定时任务、不清理旧文件或其他分析进程。

## 审核文件

run/GLOBAL_TARGETED_REVIEW.html 是入口；run/reviews下分谱系输出独立ZIP，目标35MiB/份。
原始PDF、全基因DE、逐细胞联合检出长表保留在run/branches，不全部打进上传包。
同一谱系有多个编号包时一并提交或解压到同一目录才能查看完整HTML。
输出状态为PAUSED_FOR_GLOBAL_L2_TARGETED_REVIEW或PARTIAL_FAILURE_REVIEW_AVAILABLE。

## 测试与限制

通过：语法、真实1,175,115细胞的成员与UID映射、75命名/21Unknown计数、冲突输入拒绝、
metadata保留缺失与4对donor、稀疏表达均值/检出率不修改raw、流式HDF5内容哈希变更检测。
本地没有Scanpy/AnnData；尝试安装因网络不可用失败。因此未声称测试真实Scanpy差异分析、
DotPlot或AnnData写回，科学计算与完整新h5ad需要服务器runtime校验。详见tests/TEST_RESULTS.json。
