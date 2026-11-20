from pathlib import Path
import csv, hashlib, json
import fitz
O=Path('/Users/taki/Desktop/Panvascular_Spatial_Evidence_Assessment_20260920')
B=Path('/Users/taki/Desktop/Panvascular_Marker_Spatial_Review_20260919')
(O/'一页结论.md').write_text('''# 空间证据盘点：结论与主图用途

本轮仅盘点原29候选、7研究、49切片，共1,421个gene–section组合；读取原结果并制作描述性对照，未运行新统计。下述总体结果为Figure5c原q05/all_sections/cell_abundance结果。

|证据问题|LIPA|STAT3|
|---|---|---|
|当前示例有表达吗？|CA2：686/1,254 spots阳性（54.70%），目标总UMI 2,943|LCA：1,043/2,072 spots阳性（50.34%），目标总UMI 1,702|
|当前示例有已核实的表达空间自相关吗？|未提供相应检验结果|未提供相应检验结果|
|总体与模型推断Fib背景有关吗？|原β=0.0615，95%CI −0.0409–0.1639，FDR=0.452：当前总体结果不支持|原β=0.0473，95%CI 0.0217–0.0729，FDR=0.00280：支持正向条件关联|
|展示切片自身是否支持上述正关联？|CA2的原β=−0.6903，名义P=1.33×10⁻¹⁶，为负向条件关联|LCA的原β=−0.0176，名义P=0.465，未支持正向关联|
|区域／层次富集经过组织学确认了吗？|没有定位到独立组织分区及对应检验；清亮ROI仍未确定|没有定位到独立组织分区及对应检验；边界spots保留|

**如何理解模型。** 因变量为spot目标原始UMI，offset为完整文库UMI的log；预测项为标准化Fib推断丰度、fibrous-cap表达程序及二者交互。Fib是所有Fibroblast状态丰度之和，q05是各状态边际q05之和。β表示程序处于本切片均值时的条件丰度关联，不是“成纤维细胞内表达”，程序名称也不是病理分区。

**现有候选的证据边界。** 找到GSE243179两切片的56条原生表达Moran记录，但spot集合、变换及邻接定义不全；其中LIPA/STAT3的原FDR约0.989–0.995，不能移到当前两张示例。另有13候选、497个当前背景匹配的历史残差Moran记录，是旧模型残差诊断，不能代替原表达空间自相关。没有定位到可直接支持原29候选层次富集的独立逐spot病理证据。缺项均为未提供／未检验，不是0或无支持。

**依赖与不确定性。** 原代码确实在reference和mapping训练前排除了29目标；LIPA/STAT3本来不在其背景程序中。但这仍不是独立组织染色验证。原汇总先合并同一记录单位切片，再直接对单位做随机效应汇总，并非研究层meta；35记录单位不能全部当作已确认独立供者，串行切片协方差及丰度估计不确定性未传播。LIPA/STAT3总体I²分别91.2%/64.9%，切片结果存在方向差异。

**推荐主图功能。** Figure5a只承担“固定切片中的表达检出与组织背景展示”；Figure5c承担“模型推断细胞背景的条件关联”。不能合称组织学验证或层特异富集。七个LIPA疑点仍保留，真实RNA读数不等于细胞来源已确定；未标记spots也不代表逐spot已获病理确认。正式Figure5和其余图件未修改，ORA继续PENDING。

状态：PAUSED_FOR_SPATIAL_EVIDENCE_ASSESSMENT
''')
(O/'PROPOSED_NOT_RUN_最小新增分析计划.md').write_text('''# 最小新增分析计划（PROPOSED_NOT_RUN）

本文件是待作者授权的科学诊断方案，不是绘图校验；本轮未执行任何以下检验。先固定问题与可检验范围，再注册具体参数和比较家族，不能以是否显著决定调整。

## A. 是否存在非随机空间结构

- 范围保持原29候选×49切片；原spot纳入集合、七个LIPA疑点及STAT3边界spots全部保留。基因未测量或低信息量时记录不可评价。建议沿用已有≥100 spots且目标≥20阳性spots的可估门槛，授权前确认并锁定，不能事后改变。
- 根据真实Visium网格／物理坐标确定直接邻接；物理标定缺失时使用原阵列直接邻接而不声称毫米距离。组织片段边界须独立于目标表达确定，切断跨空腔／分离组织片段的边；不使用无距离上限的kNN补足邻居。不删除孤立spot，只报告其不可用于边统计的状态。
- 主问题与调整后问题分开：固定log1p(CP10k)表达的Moran I回答描述性结构；深度、检测能力和组成影响需另行预注册模型。不能把残差Moran当成原表达Moran，也不能假定任意置换满足可交换性。授权前明确置换单位、连通分量和零模型；固定一次算法、seed和次数。
- 报告I及其零模型期望、效应差、P/FDR和可评价spot／边数。建议将所有预先符合门槛的gene–section检验作为一个明确BH家族，完整列出不可评价组合；原已发表／历史P另列，不并入重校正。不得调邻居数、截断或方法来获得显著性。

## B. 是否与原细胞背景相关

- 先完成GSE314851文库到真实供者关系核验；没有身份闭环时，不能按独立供者泛化。当前结果保留为已完成的描述性模型，不因补诊断自动重拟合。
- 若批准进一步稳健性检验，保留目标排除、原结局与offset定义，预先处理同供者多切片依赖、研究层异质性和推断丰度不确定性。优先使用已保存posterior，不为此自动重做去卷积；若现有posterior不足则标记输入阻断。具体模型需在读取新检验结果前批准。
- 明确区分无条件丰度关联、给定程序的条件关联与交互，不把三者互换。报告β、CI、研究／供者方向分布及有效重复数。沿用原29×3项的87项主问题家族定义；新增诊断独立登记，不替换原FDR。组成和程序共线性、文库offset与目标表达耦合应作为解释限制，不自动移除协变量以追求显著。

## C. 是否存在组织区域／层次富集

- 当前阻断：缺少独立逐spot组织分区。先由病理审核者在不看目标表达的条件下，依据现有H&E及可核实原标注固定区域与边界；LIPA清亮ROI不得直接命名为管腔。没有可信分区则不做区域富集检验。
- 注册估计捕获圆盘跨边界spot的统一处理规则；不按表达决定区域，也不临时扩大／缩小ROI。原图和所有spot身份保留，混合／不确定区域另列，不能暗中并入背景。
- 以供者为生物学重复，区域为供者内配对，串行切片不能作独立供者。采用预先选定的区域汇总计数或可显式处理空间依赖的模型；原始全库深度作offset。检测率、细胞组成可解释层次差异，需要分别声明“总体区域差异”与“组成调整后的差异”，不能临时互换。
- 固定原候选、事先认可的区域对比及多重比较家族；报告区域效应、CI、供者覆盖和不确定性。重复不足时只交描述性组织对照，不做泛化显著性声明。

建议顺序：先确认供者关系与独立病理分区可得性；获批后做A的固定结构诊断。B仅在依赖问题可处理且获批后执行；C必须等独立分区。任何环节允许结果不显著或数据不足，不以换基因、换切片或反复调参弥补。
''')
(O/'README.md').write_text('''# 限定空间证据盘点

入口：`一页结论.md`；状态为`PAUSED_FOR_SPATIAL_EVIDENCE_ASSESSMENT`。本轮未改正式图、正文、科学源数据或统计结果，ORA仍PENDING。

- `tables/gene_section_spatial_evidence.tsv`：29基因×49切片完整1,421行；原始UMI、阳性spot数、单位身份、原模型及空间结构来源。
- `tables/LIPA_STAT3_当前展示切片证据.tsv`、`LIPA_STAT3_全部98个gene_section证据.tsv`：当前示例与全部其余切片，避免只看显著切片。
- `tables/candidate_evidence_level_summary.tsv`：29候选证据可得性汇总。检出率不能用于排列证据强弱。
- `review/same_section_diagnostic_2pages.pdf`：两张当前切片的完整视野，对照H&E、QC、目标表达、固定谱系marker、已有Fib丰度及已有程序。对应PNG和逐spot表同目录／tables。
- `notes/model_audit.md`：实际模型、target-excluded代码行为、聚合和多重校正来源；模型表保留全部后验、敏感性、单位、切片及未估项。
- `notes/native_structure_and_histology.md`：原生表达Moran、历史残差Moran及病理标注的适用范围。
- `PROPOSED_NOT_RUN_最小新增分析计划.md`：仅方案，未执行。

## 表字段读法

`section_*`仅为该切片已有估计，原切片表未提供CI/FDR时明确标注，不推导或填入总体值。`overall_record_unit_meta_*`是总体单位汇总，为追溯而重复联接，绝非切片显著性。`section_pathology_*`沿用原代码字段名，实际上为表达程序项，不代表病理分区。

`native_*`为作者原生空间结构结果；其适用spot集合／变换／邻接未完全提供。`historical_*`为旧locked71模型的残差Moran，不能与原表达Moran混称。`metadata_status`保留供者身份强弱；recorded unit不必然是独立供者。

19个gene–section原始来源未测量／未提供，联合基因矩阵内的填充零没有解释为生物学零。缺项用“未提供”或状态字段说明。原表真实0值保留；原P数值0也不解释为精确零概率。

## 本轮实际完成的检查

只读原始计数提取与已有表达比较，最大log1p(CP10k)数值误差8.88e−16；spot身份及原文库分母逐项一致。原模型24项LIPA/STAT3汇总与采用Source Data逐值一致。诊断图实际看图与脚本检查分列于notes；这些检查不证明组织学正确。

图中固定参照为DCN、COL1A1、LST1、ACTA2，未新增程序评分。H&E保持原图，所有视图同坐标／完整视野；零、NA及既有疑点保持区别。脚本、来源路径与SHA随包保留；没有启动新模型或检验。
''')
rows=list(csv.DictReader((B/'图件来源与尺寸.tsv').open(),delimiter='\t'))
checks=[]
for r in rows:
 for kind in ('pdf','svg'):
  p=Path(r[kind]); sha=hashlib.sha256(p.read_bytes()).hexdigest()
  checks.append({'figure':r['figure'],'kind':kind,'path':str(p),'sha256':sha,'matches_baseline_manifest':sha==r[kind+'_sha256']})
assert len(rows)==15 and all(r['matches_baseline_manifest'] for r in checks)
assert len(fitz.open(O/'review/same_section_diagnostic_2pages.pdf'))==2
(O/'notes/final_scope_checks.json').write_text(json.dumps({'baseline_figures':15,'baseline_pdf_svg_hashes':checks,'diagnostic_pdf_pages':2,'scientific_tests_run':False,'formal_figure_replaced':False,'ORA':'PENDING','histology':'UNCONFIRMED','status':'PAUSED_FOR_SPATIAL_EVIDENCE_ASSESSMENT'},ensure_ascii=False,indent=2))
print('Completed notes; all 30 baseline PDF/SVG hashes match. Diagnostic PDF: 2 pages.')
