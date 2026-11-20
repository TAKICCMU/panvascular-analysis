# 已批准三大疾病背景 pseudobulk 与 CKB 重点审核

本轮用户授权：“开始分析吧，帮我看一下CKB的变化”。只在新版本目录写入结果，旧25个单研究DE、冻结h5ad和标签保持不变。

## 输入与分组

使用已审核的 disease_background_3groups_review_v2/sample_background_sidecar_v2.tsv，SHA256：852c08353c53a1b4d32b0b5928c9d1ac6a0278e365748912adfaf2ff30cf2ce4。

- G1：动脉粥样硬化及钙化性主动脉瓣病，80疾病+51健康候选样本。
- G2：主动脉夹层及动脉瘤，46疾病+41健康候选样本。
- G3：肺高压，11疾病+12健康候选样本。
- 6个无对应部位病例的肱/髂动脉不进入主分析。
- 保留6个PRJNA562645_2健康瓣膜，包括儿童。健康为用户定义的参照组，不抹除原始诊断。

复用已验证的技术样本级raw pseudobulk（sample_L1 / sample_L2），不复用旧生物样本合并结果。原始输入包含完整1,175,115细胞与37,336特征；按原稿映射门槛选36,601个唯一Ensembl特征。CKB唯一对应ENSG00000166165。

## 方法

- 按更新后的供者×组织×条件，在各背景和细胞群内合并原始计数。低计数的技术分片先合并，再检查每单位/细胞群至少20细胞。
- 每个背景/细胞群至少5独立供者/臂；同一人多个标本不计作多个独立人。
- 按原稿排除Schwann/Epithelial/Unknown的正式靶点DE，L2使用冻结MAGMA资格表103个父谱系限定标签，不使用89个CoVarNet节点代替。每个细胞群最终是否可检验单独记录。
- 用户已批准额外健康对照，因此取消旧的study必须同时有足量两臂的自动纳入门槛；不能借此取消study调整。
- edgeR 4.0.16：filterByExpr(min.count=10)、TMM、estimateDisp(robust=TRUE)、glmQLFit(robust=TRUE, legacy=TRUE)、glmQLFTest。
- 强制保留study/研究分支、重复供者、组织、assay的混杂调整空间。重复供者在同一背景中使用供者固定块，单一单位供者按study块；这是对旧脚本只处理同供者同组织配对的补充，避免多组织同供者被当成完全独立重复。
- 仅删除线性冗余设计列，不删除其代表的调整空间。若疾病效应被该空间完全解释则报告不可识别，不自动退回~group。
- 性别/年龄沿用原90%完整度与可识别条件，完整病例排除逐模型保存；年龄至少5个不同值。纳入人口学项造成样本变化、或未调整人口学项都可追溯。
- 全基因BH在每个细胞群×背景内独立计算；CKB只提取，不为CKB重新调参或只在显著基因中校正。
- 多背景的共有效应是工作假设；研究来源与疾病强耦合、儿童对照和组织差异仍须在解释时保留。只含健康的study虽然纳入拟合，对研究调整后的疾病效应未必提供直接信息；逐study设计信息量单独输出。

方法依据：原投稿包05a/05b及edgeR官方用户手册（https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf）。不将这次三大背景重分组说成完全复刻旧对比。

## 已实际测试

原生合成测试覆盖：配对供者块、设计满秩、完全混杂拒绝、额外健康/儿童对照保留、封装与直接edgeR结果一致、CKB变化方向、全基因BH范围。
实际输入：原技术级计数/基因/单位文件SHA256与原manifest匹配；247个样本的细胞数与新sidecar匹配；重聚合按基因计数守恒检查；唯一CKB基因身份。
完整PREPARATION.json与内容绑定RECEIPT.json在全部输入导出后生成。提前执行的L1模型各自有EARLY_RECEIPT，不等待L2输入全部就绪；完成后经哈希核对由主队列复用。

## 运行安排

共享最多4个CPU重任务和40实际线程预算。其他CoVarNet、scDRS、通讯、GPU任务继续。为优先获得CKB，early_ckb.py暂存本轮新输入导出进程1218024，按序运行Macrophage/Fibroblast/SMC三个背景的完整全基因DE；finally恢复导出，主队列复用已完成模型。不要把该短暂SIGSTOP当成旧任务故障，也不要单独SIGCONT导致超预算。

实际入口：

```bash
/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python -u code/prepare.py
/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python -u code/early_ckb.py 1218024
/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python -u code/await_and_dispatch.py
```

入口在新版本目录运行，读取自身脚本路径确定目录。不要无条件重启，否则会触发防覆盖检查。PID可复用，检查进程完整命令与启动时间后再操作。

进度文件：preparation.log、EARLY_PROGRESS.json、early_ckb.log、LAUNCHED.json、PROGRESS.json、dispatch.log。
结果：results/<model>/DE.tsv、CKB_result.tsv、CKB_unit_expression.tsv、MODEL.json。
最终汇总：tables/CKB_review_results.tsv、CKB_REVIEW_CN.html及figures/CKB_L1_logFC.*。

未完成的模型不能称为已完成；CKB未通过表达过滤不能称为无变化。当前不自动运行新的跨背景meta，也不把QL-F导出的近似SE当成精确标准误。旧疾病通讯入口1191582仍暂停，本次pseudobulk授权不直接恢复旧细分背景代码。
