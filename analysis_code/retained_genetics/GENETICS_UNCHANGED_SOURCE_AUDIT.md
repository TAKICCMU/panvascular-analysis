# 未变遗传背景：归档源表数值核验

本附录只读核对V1.9.3稿件、Supplementary Tables 04/06/07/09及10/11、Main Tables 1–3、Figure2原始归档导出，共12项关键数值组。未重算EFA/CFA/GWAS/LD/精细定位，也未在线刷新novelty。每项源路径、字段、原值与状态见同名TSV；逐artifact SHA256在GENETICS_UNCHANGED_SOURCE_HASHES.tsv。**归档数值一致≠新图谱重分析完成，也不等于对源GWAS研究重新审计。**

|审计|主张/核对结果|判断|
|---|---|---|
|G01|8研究表中cases合计263,790、controls合计3,962,345|与正文一致；只是贡献研究样本数之和，不能写独立去重参与者|
|G02|LDSC 28唯一trait pair有26正相关，范围−0.118189至0.716072|支持widespread positive；不等于26对均显著|
|G03|EFA aortic F2 0.734819–0.917585、F1 −0.023373–0.063221；其他5trait F1 0.515428–0.815885、F2 −0.091781–0.093734|与正文三位小数范围一致|
|G04|two-factor χ²18.224075/df20/P.572650/CFI1/SRMR.062836；bifactor χ²10.728656/df13/P.633537/CFI1/SRMR.046201|正文及Main Table1舍入一致|
|G05|ECV_general=.4091519279，即40.9%|一致|
|G06|full axis correlation=.32085342，leave-CAD=.260160896（P=1.98777e−7）|支持留一CAD后跨轴相关仍存在；此次未重新核chr split|
|G07|FUMA QC输入5,448,367 variants；primary finemapping归档Neff=1,271,027|符合约5.45m和正文Neff；Neff是下游输入旁证，未重新估计|
|G08|FUMA清理输入GWS=9,434；删除3,552，保留5,882；原始GWAS plot export GWS=9,437|正文9,434→5,882正确，但应注明FUMA清理输入口径，不能与原始GWAS输出混写|
|G09|直接计数归档GenomicRiskLoci.txt：primary155、r².6 sensitivity171、leave-CAD67|与正文及Main Table2一致|
|G10|归档interval-union intersection：full r².1与leave-CAD r².1重叠25|一致；未重新计算LD|
|G11|442个方法行对应222唯一scenario-region；analysis_ready=True为182唯一region；z-only182、median-Nhat180|**原文“222 analysis-ready regional tasks”不一致，应改**|
|G12|Supplementary Table10与Main Table3均为5个同名候选：rs7525267、rs75135846、rs78051210、rs333929、rs7236320|名单一致，novelty只能作为历史审计结果|

建议精细定位措辞：**“Mismatch-aware SuSiE-RSS analyses were archived for 222 regional tasks, of which 182 met the analysis-ready criteria in at least one prespecified implementation.”** 若要描述最终主方法，需要指定z-only（182）或median-Nhat（180），不能把两种方法442行当作442区域。40个区域未在任一方法获analysis-ready标记；不宜称全部222通过。这里检查的是归档布尔标记，不是重新评估所有LD/可信集质量。

Novelty时间边界：Supplementary Table11的GWAS Catalog release为**2026-07-20**；Table10明确PhenoScanner于**2026-08-11超时、独立确认待完成**。此处未重新访问搜索服务，不能升级为当前仍然新颖或多资源全部独立确认。可保留“历史多资源审计下的component-disease-novel候选”，同时保留非全球新颖、非因果的限定。

本附录独立于正在运行的新eQTL下载/ABF控制器，不修改其状态、现行主报告或交付hash；源审计hash单独提供。
