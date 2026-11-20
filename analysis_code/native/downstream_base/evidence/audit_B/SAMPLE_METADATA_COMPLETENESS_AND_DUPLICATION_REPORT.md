# 单细胞图谱样本元数据完整性与重复审计

## 审计范围

- 输入文件：`sample_study_clinical_metadata.csv`，247 行。原文件未被覆盖。
- 校验图谱：更正后的全基因计数对象为 1,175,115 个细胞 x 37,336 个基因；冻结 HVG 对象为 1,175,115 个细胞 x 2,000 个基因。
- 核对层级：技术样本、生物学样本、供者、组织、疾病背景和测序类型分开审计。
- 来源：GEO family SOFT/BioSample/SRA 记录、原始论文及补充材料、HCA Table S13、PRJNA562645 归档元数据、作者公开 H5AD（Zenodo 10.5281/zenodo.17520898）和现有元数据 crosswalk。

## 总体结果

- 新版文件与当前更正图谱的 247 个 `sample_uid` 完全一致；每个样本的细胞数和总细胞数均完全一致。
- 当前元数据包含 226 个生物学样本键。242/247 行可映射到 163 个已解析供者键。
- 年龄从 139/247（56.3%）提高到 171/247（69.2%）；性别从 161/247（65.2%）提高到 191/247（77.3%）。
- 疾病或参照分析背景从 53/247（21.5%）提高到 247/247（100.0%）；`source_diagnosis`、`tissue_detail`和`assay_type`也已全部填写。
- 测序类型为 235 个 scRNA-seq 技术样本和 12 个 snRNA-seq 技术样本；后者全部来自 GSE255895。
- 自动验收共 41 项，41 项通过，0 项失败。

## 主要修正

1. 将原 `sample_id` 中的 study-level 值替换为真实的 GSM、SRX 或来源样本键，并新增 `atlas_source_id`、`source_sample_accession`、`biological_sample_uid` 和 `donor_uid_curated`。
2. 修正了供者层级：HCA 为 38 个技术样本、37 个生物学样本和 6 位供者；GSE224559 为 2 位供者；GSE224587 为 5 位供者；GSE248609 为 4 位供者；GSE247238/GSE294291 在当前图谱中为 14 位供者。
3. 将技术性 `library_size_raw_cpu` 审计从临床元数据冲突中分离，因此原文件中 247 行全部标为临床冲突的情况已修正。
4. 区分了“参照组织”与“健康供者”。例如 GSE155468 的参照为非动脉瘤疾病对照，GSE213128 来自终末期心力衰竭移植受者的非病变瓣膜，GSE248609 的 NC 为糖尿病患者的非钙化区域，GSE255895 的 IMA 来自接受 CABG 的 CAD 患者。
5. PRJNA562645 的两个导入来源已分开：第一个包含 CAVD 与非 CAVD 对照，第二个是非 CAVD 瓣膜的年龄比较队列。
6. GSE255895 已按原始研究更正为 snRNA-seq，不再将整个图谱简化描述为仅含 scRNA-seq。

## 重复与异常重叠审计

### 未发现的问题

- `sample_uid` 重复数为 0，H5AD `obs_names` 重复数为 0。
- 全部 247 个样本两两比较的 2,000-HVG pseudobulk 指纹中，未发现“生物学样本键不同但表达指纹近乎完全相同”的样本对。
- 已知 GSE198750/PMID33626908 跨来源重复表示已处理：当前 H5AD 中没有 GSE198750 输入表示，保留了 PMID33626908 表示。新文件中 `study` 保留与 H5AD 对齐的历史分组名，实际保留来源应以 `atlas_source_id=PMID33626908` 判定。

### 需要在图谱冻结前解决的一项风险

- `GSE210248__GSM6424655_rep1` 与 `GSE210248__GSM6424655_rep2` 指向同一 GSM、BioSample 和 SRX，HVG pseudobulk Pearson r=0.995979。
- 两者分别含 9,632 和 9,269 个细胞，共享 112 个 barcode 字符串。这些共享 barcode 中没有完全相同的 2,000-HVG 或 37,336-基因表达行。
- 该结果排除了简单的整行重复，但不能区分“同一 GEM well 的多次测序被分别计数”与“同一供者的两个独立捕获”。GEO 仅提供一个合并处理矩阵，因此建议追溯原处理脚本，或用 GEO 的单个合并矩阵重建该样本后再冻结图谱。在查明前不应仅凭 barcode 或高相关性删除任意一组。

### 合理的重复设计

- GSE131778、GSE237230 和 HCA 包含同一生物学样本的多个技术捕获或技术重复。
- GSE159677、GSE224559、GSE224587、GSE248609、GSE247238/GSE294291、GSE220774 和 HCA 包含同一供者的不同区域、层次或组织。这些样本应保留，但不能作为独立供者重复计数。
- GSE220774 有 15 个 GEO 处理样本记录，论文报告 13 个独立 cDNA 文库。供者分组已核实为 5 人，但在将 `sample_uid` 解释为独立文库前，仍需追溯 15 个处理记录与 13 个文库的对应关系。

## 仍不能安全补齐的信息

- 年龄仍缺失 76/247，性别仍缺失 56/247。原因是公开记录未报告个体值，或论文只给出组水平汇总且无法与 GSM 一一配对。本审计未将组均值写入个体行。
- GSE213128 的 5 个瓣膜样本有独立 BioSample、年龄和性别，但无公开跨瓣膜供者对照表。`donor_uid_curated` 因此保持缺失，而不是用 `NR` 把 5 个样本错误合并。
- 没有基因型或供者身份 crosswalk 时，表达指纹审计不能证明不同项目之间绝对没有匿名供者重叠。当前结论是“未发现异常证据”，不是“已证明完全独立”。

## 下游使用规则

- 去批次和技术可视化可用 `analysis_sample_uid`；这不代表生物学独立性。
- pseudobulk、细胞比例、细胞通讯和其他样本水平推断应先合并同一 `biological_sample_uid` 的技术重复，并在需要时以 `donor_uid_curated` 建模或分层。
- 对同一供者的配对区域或多组织数据，应使用配对、供者随机效应或供者聚类稳健方差，不应将每个组织当作独立人。
- 建模时保留 `assay_type`，因为当前图谱包含 scRNA-seq 和 snRNA-seq。

## 主要公开来源

- GEO accession records: https://www.ncbi.nlm.nih.gov/geo/
- HCA arterial atlas: https://doi.org/10.1016/j.xgen.2025.101034
- HCA full text and supplementary context: https://pmc.ncbi.nlm.nih.gov/articles/PMC12802643/
- GSE159677 source paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC9556750/
- GSE220774 source paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC11675841/
- GSE224559 source paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC11211506/
- GSE248609 source paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC11471023/
- GSE255895 source paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC11551873/
- GSE271435 source paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC12912402/
- GSE254132 source paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC11686284/
- GSE247238/GSE294291 author object: https://doi.org/10.5281/zenodo.17520898
