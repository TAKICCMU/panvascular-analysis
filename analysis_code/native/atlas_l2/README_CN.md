# 最终L2：固定映射、marker+父类命名、保护raw、备份后替换

## 本轮决定

不再进行HVG/PCA/Harmony/Leiden、FindMarkers或模型拟合。现有多轮审核足以支持这次采取的保守注释粒度。

- 应用75个明确的待定单元决定，全部以实际cell_id和当前Pending残余成员为依据。
- 父身份有支持但缺乏独立细功能证据的24个单元（222,383细胞），归入12个`marker+父类`标签；不使用`unspecified`。
- 8个身份混杂/不可靠的单元，共21,130细胞，设为L2 Unknown。3个是上轮明确提案，另5个是按本次“实在混杂即Unknown”原则作出的保守决定。
- 原有97,090个Unknown不回收，L1不改。Schwann 2,832细胞保留Not_assessed。
- MacL2_R1_r1p2_13的14,189细胞更名为CCL3 activated resident-like macrophage。
- SPP1 lipid-associated macrophage 7,757细胞保持不变，CYP27A1群仍独立。
- SMC6/12共同归入PRKG1 contractile SMC，不强行声明两个已验证功能亚型。
- 其他此前已批准名称保留，合并只在明确决定表范围内执行。

**marker+父类不是一个额外功能状态发现。** CD68 macrophage、DCN fibroblast等标签表示父类表达支持及细分不足，不要求该marker在每个细胞中检出，也不声称其独占该群。外部ontology保留label_scope解释；不会给母obs增加新的置信度或审核状态列。116个命名标签不等于116个经过独立验证的生物学亚型，更不自动等于116个CoVarNet合格节点。

最终预计：
- 1,054,063个命名L2细胞。
- 118,220个Unknown（包括原L1 Unknown 57,918个）。
- 2,832个Not_assessed。
- 0个Pending/unspecified。
- 116个非特殊L2名称，特殊值Unknown/Not_assessed另外统计。

代码中没有依照module最高分自动注释，也没有按目标比例回收细胞。

## 输入与运行

本包整个文件夹上传至：

`/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908/L2_final_marker_parent_freeze/`

在原Notebook中执行RUN_ONE_CELL.py的全部内容；RUN_ONE_CELL.ipynb也只有这一个cell。

必须保留或明确设置：
- `L2_UPDATED_H5AD`：上一轮更新后、临床metadata已导入的全量工作文件。
- `L2_NEXT_RUN`：上一轮14谱系有限复核目录。可从当前h5ad的`uns['L2_last_update_run']`解析，不需要在内存中加载巨噬细胞。

输入必须已是690,772 named /384,421 Pending /97,090 Unknown /2,832 Not_assessed的已审核版本。程序核对真实成员，不只看数量。若你后来另行修改过标签，会停止而非覆盖。

不要按修改时间自动搜最新h5ad，不可把冻结L1-only/raw counts源文件当覆盖目标。源文件必须位于指定ROOT。

同一项目40线程上限内，本轮只使用4线程预算；科学运算没有新增，主要耗时是大型文件复制和多轮哈希。前台等待至COMMITTED_FINAL_APPROVED，不关闭kernel，不同时运行修改这个母文件的其他脚本。普通内存对象无需清空；待覆盖文件的backed句柄必须关闭。现有局部Macrophage等对象不会被修改。

## 保存策略：同时满足“覆盖上一版”和固定最终入口

固定正式入口：

`ROOT/panvascular_CPU_theta2_L1L2_final_v1.h5ad`

1. 保存上一版标签及批准记录，校验当前样本、L1、原分区及75个残余review_unit。
2. 读取当前母文件obs/var及必要uns；所有raw/log1p/图只流式校验，不整体读入内存。
3. 运行服务器安装环境下的3细胞AnnData原生读写烟雾检查。
4. 复制完整母h5ad到同一磁盘的临时文件，只更改最终L2及必要冻结uns。
5. 校验全部X/layers/var/varm/varp/obsm/obsp（以及若存在raw）与旧文件逻辑哈希一致，临床和L1逐项不变。
6. 保留旧完整h5ad备份，再替换明确指定的原工作路径。
7. 固定最终文件名与更新后的工作文件名指向同一个完整HDF5文件（Linux硬链接）。不是表达外部链接，不重复保存第二份新表达数据；任一名字均可独立读取完整h5ad。
8. 两个名字都设只读，生成manifest/SHA256/PROJECT_STATE和外部ontology。

至少需要一份完整新h5ad大小加2GB的空闲磁盘；运行前确认有余量。旧版备份是旧inode，与新的正式文件不同。读取程序可以使用固定正式路径；不要用旧的内存AnnData重新写回这两个名字，否则会产生版本混淆。

没有自动启动任何下游科学模型。数据ready_for_downstream仅表示输入发布完成；统计设计和执行授权仍需审核。

## 最终结构

- X：1,175,115 x37,336全基因raw counts，不变。
- layers['log1p']：全基因CP10k+自然log1p，不变。
- obs['level1_annotation']：原冻结L1，不变。
- obs['level2_annotation']：最终L2；旧L2_annotation同义列被移除。
- 所有临床、donor、sample、study、assay及QC字段维持原值和缺失。
- 全局PCA/Harmony/UMAP/图、基因与PCA loading不变；局部图仍在各分支原h5ad。
- 外部level2_ontology.tsv提供稳定state ID、所属L1、marker、人数与名称语义，不新增obs注释标签。

附表位于`ROOT/atlas_release/L1L2_final_v1/`。
备份和逐细胞变更位于`ROOT/L2_finalization_<timestamp>/`。
阶段日志位于`ROOT/L2_final_freeze_launch_<timestamp>/`。
最终标志为COMMITTED_FINAL_APPROVED。

## 中断和重复运行

- 同一项目的最终化进程有排他锁，上一轮focus_scheduler仍运行时拒绝发布。
- 标签或分区变化、数据哈希变化、存在未知Pending等都会停止。
- 更名时备份旧文件，安装失败会尝试回滚旧工作文件。
- 已提交完整h5ad但附表写入中断时，再运行同一入口可根据已验证manifest完成附表；需校验新文件、旧备份和所有附表摘要。
- 已成功发布后再运行只复核，不重新应用分类。
- 不自动删除失败留下的partial文件。发生错误时根据FREEZE_STATUS和日志处理，不能把partial文件当正式发布文件。

## 测试范围

见tests/TEST_RESULTS.json。已使用真实上传的1,175,115细胞ID/L1登记、14个局部成员表及170个审核汇总进行映射测试；raw表达未提供，不能在本地验证服务器完整raw矩阵。全量计划测试中其他15类中的未进入L2成员仅核对其真实ID/L1及应保留的特殊标签，不伪造其样本信息。

本地已测试稀疏HDF5逻辑哈希、元数据不变性、备份/覆盖/回滚、硬链接、附表中断恢复、临床缺失和冲突处理、重复审批及错误输入拒绝。

本地没有AnnData/Scanpy，因此**没有执行本地原生AnnData最终写回或服务器百万细胞I/O**。尝试安装AnnData未成功。代码会在服务器用真实metadata类型先做小型原生测试，再写全量并读回核验。不要把准备代码说成已经完成服务器冻结。

## 主文件

- finalize_level2.py：主程序。
- inputs/final_pending_decisions.tsv：75个逐单元最终决定。
- inputs/final_new_Unknown_units.tsv：8个新增Unknown单元。
- inputs/marker_parent_label_summary.tsv：12个marker+父类分组及24个来源单元。
- inputs/projected_final_ontology.tsv：基于实际上传成员重建的预计116个名称，服务器发布时重新核验。
- evidence/final_anchor_expression_evidence.tsv：本轮新增标签锚点在真实残余单元中的均值/检出率。
- CHECKSUMS.json：整个交付包的静态内容校验，不能只改映射表绕过审核。
