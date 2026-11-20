from pathlib import Path
import json,csv,re,hashlib,shutil,datetime
import pandas as pd
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');O=P/'图件编辑资料/全流程代码审核_20260922';B=O/'Panvascular_code_review';A=B/'audit'
base='analysis_code/native/downstream_final/';over=base+'overnight_manuscript_integration_20260917/'
rows=[]
def add(stage,code,data,out,unit,status,limit):
 rows.append(dict(stage=stage,adopted_code=code,input_evidence=data,result_and_table=out,statistical_unit_or_encoding=unit,audit_scope=status,remaining_limit=limit,code_path_exists=(B/code).exists()))
add('01 原始QC/Scrublet','recovered_sources/original_QC_notebook_code.json','原始研究矩阵、QC metadata；本轮恢复原notebook代码及execution_count','Table01/02/28；1,175,115保留细胞','按原sample批次QC、原counts双细胞识别','SOURCE_RECOVERED_NOT_RERUN','不是247个样本全部执行日志；不能将单个notebook例子当作所有样本逐一验证')
add('02 纠正重复来源及全基因counts','recovered_sources/raw_count_build.py','原样本counts和去重映射','最终全量counts；发布manifest','cell_id/原始矩阵；不重复计入被排除来源','SOURCE_RECOVERED_NOT_RERUN','原始counts外置，不随审核包上传')
add('03 Pearson residuals样本clip','analysis_code/native/atlas_cpu/patched_function.py','原counts；sample_uid；theta100；每样本sqrt(n)clip','2000HVG；sample_clips.tsv','官方排名/合并语义；局部clip修正','SOURCE_AND_PARAMETERS','本轮不重新计算残差；原通过记录另存')
add('04 全量CPU theta2基线','analysis_code/native/upstream/CPU_theta2_author_executed_notebook.py','pipeline_parameters.json和HVG2000.tsv','PCA30、CPU Harmony theta2、neighbors30、Leiden1、UMAP；全量冻结坐标','全细胞；CP10k/log1p、HVG scale','SOURCE_PARAMETER_RESULT_LINK','混合notebook含theta1历史单元及未执行单元；不可直接逐行跑整个.py')
add('05 scIB同输入比较/CPU补充分支','recovered_sources/scib/benchmark_CPU_Harmony_additional_v1/code','冻结Primary315949；评价182968；247样本；2000HVG30PC','Table27、ED1；CPU-only出版表','相同评价细胞/标签；biology与batch分别保留','SOURCE_RECOVERED_AND_SCORES_MATCH','benchmark的冻结输入不是最终全量theta2图谱；嵌入/H5AD外置')
add('06 L1局部与人工映射','analysis_code/native/upstream/L1_CPU_partial_and_local_code','保留原cell_id、raw连接、候选及批准映射','最终16类L1；Figure1/ED2','局部证据+人工审批，不重新推断标签','SOURCE_MANIFEST_COUNTS_MATCH','未重新执行人工判定；不同历史L1版本不可互换')
add('07 L2分支','analysis_code/native/upstream/L2_CPU_others','冻结L1；分谱系原counts、局部参数','116命名状态；128父类限定完整分组','状态程序；含Unknown/Not_assessed','SOURCE_RETAINED','不是所有命名状态均经独立生物学验证')
add('08 L2冻结/发布契约','analysis_code/native/atlas_finalization/finalize_level2.py','决策表、原始cell_id、metadata v4','最终L1L2 h5ad manifest；Table03/29；Figure1','X原counts；layer log1p；正式level1/2字段','MANIFEST_SOURCE_DATA_COUNTS_MATCH','完整h5ad按用户要求不打包；发布时ready与执行授权属于不同记录')
add('09 样本/供者/疾病sidecar','analysis_code/native/analysis_sidecar','各分析实际采用config/receipt，而非最新mtime','Table02/61；原稿策略病例对照配置','technical capture/生物样本/供者区别','SOURCE_AND_CONFIG_RETAINED','不同模块有不同批准范围；不以一个新sidecar覆盖全部历史输入')
add('10 描述统计/marker汇总','analysis_code/native/downstream_base/src/expression.py','冻结X和log1p；本类型分母','Table26/28/30/31/36；Figure1g、ED2–5','mean含零；raw>0检出；0与NA分开','SYNTHETIC_PASS_AND_STORED_MATRIX_MATCH','缺少cell/sample×gene测量掩膜的问题不因图好看而消除')
add('11 遗传输入/LDSC/Genomic SEM','analysis_code/retained_genetics/Genetics/full_pipeline','原GWAS与LD参考、性状定义','Figure2；Table04/06–11/25/64','遗传相关、双因子、primary与leave-CAD','SOURCE_RETAINED_NOT_NATIVE_RERUN','基因型输入和LD参考外置；不把atlas更换当成GWAS重估')
add('12 MAGMA细胞类型富集','analysis_code/native/production_patches/magma_recovery_roundtrip_fix_20260914T0930','3场景×L1/L2×2权重；正式冻结labels','Table05/42/52；Figure3a/b','Average保留于回归、结果过滤后mapping一对一再BH','FINAL_PATCH_IDENTIFIED','历史v1/v2修复不是生产入口；本轮未重新调用MAGMA二进制')
add('13 scDRS','analysis_code/native/production_patches/scdrs_repair_20260914T0910','最终基因映射、3遗传场景、covariates','Table05/52；Figure3a/b','原trait score及群体汇总','SOURCE_RETAINED','scDRS原生运行/大矩阵不在本轮重跑；不以FDR单项代替双项支持')
add('14 主pseudobulk','analysis_code/native/downstream_final/pseudobulk_manuscript_strategy_20260916_v1/code/native.R','原策略样本分组、原counts聚合及设计；每组至少5','原稿策略主DE；Table43','edgeR TMM/robust QLF；设计块、可用sex/age；REML候选meta','NATIVE_SYNTHETIC_PASS','SE=abs(logFC)/sqrt(F)为继承近似，不是精确QL SE；三背景版不当主版')
add('15 Pseudobulk通路','analysis_code/native/downstream_final/pseudobulk_manuscript_strategy_enrichment_v1/code','全基因DE；DL signed Z；原gene sets','Table20/37完整结果；Figure3c右','DL通路排名，与REML候选meta区分','NATIVE_DL_TEST_PASS_AND_DISPLAY_MATCH','继承SE近似；本轮不重跑GSEA/ORA')
add('16 MAGMA通路','analysis_code/native/downstream_base/external_legacy/pathway_enrichment.R','遗传rank与gene sets；Table19原结果','Table19；Figure3c左17885基因9通路','preranked GSEA；历史ORA单独保留','DISPLAY_TO_COMPLETE_RESULT_MATCH','此源文件为保留原策略；完整原生重跑未执行；ORA按作者接受口径保持')
add('17 候选联合筛选','analysis_code/native/downstream_final/overnight_manuscript_integration_20260917/candidate_prioritization/build_candidates.py','2199遗传优先基因；REMLgene-state结果','29候选对；Table13；Figure2a/3d-e','五项联合规则；非空间/模块新增筛选','SOURCE_RULE_AND_29_ROW_VALUES_MATCH','输入全名单按源路径追溯；不可从Venn并集1658反推2199')
add('18 ABF共定位/QTL','analysis_code/native/downstream_final/overnight_manuscript_integration_20260917/candidate_prioritization','已有及新增QTL对齐、两GWAS场景','Table12/40/41/46；Figure3e/f','ABF后验；top shared SNP按SNP.PP.H4最大','SOURCE_RETAINED','源QTL/LD外置；未重跑QTL拉取或coloc')
add('19 SuSiE-RSS/QSNP','analysis_code/retained_genetics/Genetics/09e_susieRSS_mismatch_aware_rerun_Jupyter一键运行版.py','Z+LD；mismatch-aware有限参考N路径','Figure2、ED遗传；region/credible-set source','z_only不等于没有LD；原接受规则','SOURCE_RETAINED','原参考LD与GWAS外置；不把历史失败运行当最终输出')
add('20 CoVarNet主网络','analysis_code/native/downstream_final/code/run_covarnet.py','已批准sample/donor/group；89资格节点','Table14/15/44/61；Figure4','样本组成、NMF、原rank评价','SOURCE_RESULT_TABLE_RETAINED','K14由本轮结果记录；不硬设新候选数；程序中保留已选K引用')
add('21 CoVarNet敏感性/疾病关联','analysis_code/native/downstream_final/remaining_downstream_20260916_v1/covarnet','多种扰动、分层重采样、疾病模型','Table45/55/63；ED9','study/donor分层；保留不可估计','SOURCE_RETAINED','最终association_patch_v2优先于被修复脚本；未真实重拟合')
add('22 CoVarNet最终分级','analysis_code/native/downstream_final/overnight_manuscript_integration_20260917/covarnet_final_tiers/rebuild_final_tiers.R','结构稳定性+来源覆盖+研究内一致性','14模块等级；Figure4b；Table14/16','三维证据；不将结构通过等同全面高置信','SOURCE_AND_FINAL_TABLE_RETAINED','分类值不因可视化更换而改变')
add('23 细胞比例疾病关联','analysis_code/native/downstream_final/remaining_downstream_20260916_v1/composition','当前采用样本/donor；n_state/n_parent','Table35/60','beta-binomial GLMM；合成正向控制','SYNTHETIC_BACKEND_PASS','不能把基类测试通过称真实最终模型全部重算')
add('24 L1/L2逐样本通讯','analysis_code/native/downstream_final/remaining_downstream_20260916_v1/communication','批准样本信息、LIANA/CellChat门槛','Table38/39/49；通讯补充图','每生物样本推断；跨方法键对齐','SOURCE_RETAINED','采用thread_limited_r5与续跑记录；不启动旧队列')
add('25 通讯疾病/供者敏感性','analysis_code/native/downstream_final/overnight_manuscript_integration_20260917/review_communication.py','逐样本结果、疾病/背景、供者','Table38/39/49/62','按7键配对；same direction与both supported分开','SOURCE_RULE_RETAINED','Table39是记录计数非比例；全部未测量/不支持分开')
add('26 空间reference/mapping','analysis_code/native/downstream_final/spatial_gpu_resume_v2','sc参考与Visium原counts/image/scalefactors','Table21/22/54；原posterior文件','cell2location；q05和mean不同','SOURCE_AND_RUN_RECEIPTS','模型权重/原始空间h5ad外置；不重复训练')
add('27 target-excluded mapping','analysis_code/native/downstream_final/overnight_manuscript_integration_20260917/spatial_target_excluded_v1/code/spatial_target_excluded.py','identity/alias排除列表、原参考/映射gene集合','Table53；最终排除版posterior','reference及mapping基因身份排除','SOURCE_AND_COMPLETION_RECEIPT','不能仅凭分支名推断无循环性，保留exclusion.py实际检查')
add('28 空间联合模型/汇总','analysis_code/native/downstream_final/overnight_manuscript_integration_20260917/spatial_joint_evidence_v1/code/run_joint_mapping_adapter_v2.py','target-excluded mapping-root；49切片→35批准单位','Table32/48/57–59；Figure5b/c/d','Poisson loglink+library offset、丰度/程序/交互、block cluster SE；87BH','SOURCE_AND87_ROWS_MATCH','切片→单位固定效应→单位meta；不是按study再次meta；条件关联≠细胞内表达或组织学层次验证')
add('29 空间模块证据','analysis_code/native/downstream_final/overnight_manuscript_integration_20260917/spatial_joint_evidence_v1/code/legacy_spatial.py','模块权重、丰度、程序、原物理坐标','Table17/18/33/34；Figure5e/f','原section/dataset层描述；五程序完整保留','SOURCE_RESULTS_RETAINED','不把model丰度称为独立染色；原空间结构定义不改变')
add('30 正式绘图/编辑链','figure_code/native_plot_code','Source Data＋原脚本＋编辑补丁＋semantic_colors_current.json','15个当前PDF/SVG','编辑层与科学数值层分开','MASTER_HASH_AND_RENDER_TEST','多阶段编辑链；旧绝对路径和缺PATH_MAP导致不能跨机直接从数据一次重画')
add('31 15幅SVG导出','figure_code/render_current_svg.py','frozen_svg中逐图哈希匹配采用版本','当前15个母版PDF等尺寸导出','只做SVG→PDF，不是数据→图验证','ACTUAL_EXPORT_TEST','需要Cairo原生库、CairoSVG及PyMuPDF；字体/版本影响需对照')
add('32 补充表字段/数据适配','evidence/supplementary_tables/read_supplementary_table.py','最新xlsx+column_name_map+reviewed_field_definitions','64表、五项字典修正与8个计数header','机器字段与显示名明确映射','ALL64_WORKBOOK_DATA_EQUAL','旧脚本不能按美化后的显示名猜原变量；不要运行旧表头推断来覆盖已核查字典')
pd.DataFrame(rows).to_csv(A/'MODULE_CODE_RESULT_CROSSWALK.tsv',sep='\t',index=False)
# Fixed map to all 64 numbered tables
by={}
for ns,st in [([1,2,28],['01','02','09']),([3,26,29,30,31,36],['06','07','08','10']),([4,6,7,8,9,10,11,25,64],['11','19']),([5,42,52],['12','13']),([12,40,41,46],['18']),([13,43],['14','17']),([14,15,16,44,45,55,61,63],['20','21','22']),([17,18,33,34],['29']),([19],['16']),([20,37],['15']),([21,22,23,32,48,53,54,56,57,58,59],['26','27','28']),([24,47,50,51],['08','32']),([27],['05']),([35,60],['23']),([38,39,49,62],['24','25'])]:
 for n in ns:by[n]=st
trs=[]
for n in range(1,65):
 f=B/f'evidence/supplementary_tables/Table {n:02d}.tsv.gz';title=pd.read_csv(f,sep='\t',nrows=0).columns[0]
 trs.append(dict(table=n,title=title,modules=';'.join(by[n]),data=str(f.relative_to(B)),scope='全表数据/列名一致性已核对；科学模型按模块清单所列核查深度'))
pd.DataFrame(trs).to_csv(A/'ALL64_TABLE_CODE_MAP.tsv',sep='\t',index=False)
# Figure roles plus last edit files, no claim one script alone recreates complete history
fm=[]
for r in pd.read_csv(A/'FIGURE_MASTER_BINDINGS.tsv',sep='\t').to_dict('records'):
 name=r['figure'];n=int(re.search(r'(\d+)$',name).group(1));main=name.startswith('Figure')
 if main:
  mods={1:'03;04;06;07;08;10',2:'11;17;19',3:'12;13;14;15;16;17;18',4:'20;21;22',5:'26;27;28;29'}[n]
  native={1:'figure_code/native_plot_code/current/adopt_Figure1g_B.py',2:'figure_code/native_plot_code/current/fig12_visual_0921.py',3:'figure_code/native_plot_code/current/figure3_alignment_0921.py',4:'figure_code/native_plot_code/current/fig45_visual_0921.py',5:'figure_code/native_plot_code/点线可见性收尾_20260921/visibility.py'}[n]
  src='evidence/source_data/Figure'+str(n)
 else:
  mods={1:'01;05',2:'10',3:'07;10',4:'07;10',5:'07;10',6:'23',7:'11',8:'11;19',9:'20;21;22',10:'18;28;29'}[n]
  native='figure_code/native_dependencies' if n not in [2,3,4,5] else 'figure_code/native_plot_code/点线可见性收尾_20260921/visibility.py'
  src='evidence/source_data/ExtendedData'+('3_4_5' if n in [3,4,5] else str(n))
  if not (B/src).exists():src='evidence/source_data/Figure'+('2' if n==7 else '4')
 fm.append({**r,'modules':mods,'last_edit_or_native_collection':native,'source_data_collection':src,'current_override_index':'evidence/source_data/Current_Display_20260920/CURRENT_PANEL_SOURCE_INDEX.tsv','recreation_level':'SVG export tested; full multi-stage source-data rebuilding NOT tested'})
pd.DataFrame(fm).to_csv(A/'FIGURE_CODE_SOURCE_MAP.tsv',sep='\t',index=False)
# code-reading excerpts for exact statistical choices, without imports/jobs
queries=[(rows[13]['adopted_code'],r'allowed|glmQL|meta_one|p.adjust'),(base+'pseudobulk_manuscript_strategy_enrichment_v1/code/meta_DL.R',r'DL|SE_approx|meta_group'),(over+'candidate_prioritization/build_candidates.py',r'if b\(|strict_three'),(over+'spatial_joint_evidence_v1/code/run_joint.py',r'eligible=|groupby\(|fixed_effect|random_effects|87|target_removed'),(over+'spatial_joint_evidence_v1/code/legacy_spatial.py',r'GLM|offset=|cov_type|lineage_abundance|Poisson'),('analysis_code/native/downstream_base/external_legacy/pathway_enrichment.R',r'overlap|p.adjust|phyper')]
ex=[]
for rel,pattern in queries:
 for i,line in enumerate((B/rel).read_text().splitlines(),1):
  if re.search(pattern,line):ex.append({'code':rel,'line':i,'source_line':line})
pd.DataFrame(ex).to_csv(A/'MODEL_DEFINITION_SOURCE_LINES.tsv',sep='\t',index=False)
(B/'01_全流程代码与结果导航.md').write_text('# 全流程导航\n\n这是本轮重新核对后的代码审核入口。原件代码保留，采用版本按下表定位；不要按文件mtime选择入口。\n\n'+ '\n\n'.join(f"## {r['stage']}\n\n- 代码：[{Path(r['adopted_code']).name}]({r['adopted_code']})\n- 输入：{r['input_evidence']}\n- 对应结果：{r['result_and_table']}\n- 单位/定义：{r['statistical_unit_or_encoding']}\n- 本轮核查：{r['audit_scope']}\n- 边界：{r['remaining_limit']}" for r in rows),encoding='utf-8')
(B/'02_流程依赖图.md').write_text('''# 数据与代码依赖

```mermaid
flowchart TD
  Q[原始研究矩阵 / QC / Scrublet] --> C[纠正重复来源 / 全基因counts]
  C --> H[按样本Pearson HVG2000 / PCA30 / CPU Harmony theta2]
  H --> L[人工审阅L1 / 分谱系L2 / 冻结契约]
  B[冻结Primary scIB输入] --> BC[多方法benchmark / CPU Harmony补充分支]
  BC -.方法比较证据；不是全量输入.-> H
  L --> D[样本供者sidecar / 描述统计 / marker]
  L --> PB[原稿策略pseudobulk QLF]
  PB --> R[REML候选meta]
  PB --> DL[全基因DL signed Z / GSEA]
  G[GWAS / LDSC / Genomic SEM / MAGMA gene statistics] --> MG[MAGMA细胞富集 / scDRS / 遗传GSEA]
  L --> MG
  MG --> R
  R --> T[五项联合规则 / 29候选]
  T --> EQ[QTL / ABF / SuSiE证据]
  D --> CN[89节点CoVarNet / 稳定性 / 分级]
  D --> CC[逐样本通讯 / 疾病差异 / 供者敏感性]
  D --> CP[疾病细胞比例]
  L --> SP[cell2location reference / mapping / target exclusion]
  T --> SP
  SP --> JM[切片模型 / 批准单位汇总 / 87检验家族]
  CN --> JM
  JM --> SD[固定Source Data / 64表 / 15图]
  DL --> SD
  EQ --> SD
  CC --> SD
  CP --> SD
```

箭头表示输入依赖，不代表独立验证。L1/L2人工决定不可由历史代码自动重新批准。
''')
(B/'03_当前决定与科学口径.md').write_text('''# 当前采用口径

- 主输入为最终 CPU theta2 L1/L2 图谱：1,175,115细胞、37,336基因、247 sample_uid；16个L1、116命名L2，128个父类限定完整分组。完整h5ad不在审核包内，发布manifest及hash保留。
- pseudobulk主分析是原稿策略复刻版。候选整合为REML，通路signed Z为DL；继承的abs(logFC)/sqrt(F)是近似SE，不能改称原生精确QL SE。
- Figure1g是14×28既有均值的逐基因min–max B；点面积沿用原counts>0比例。逐基因颜色不允许跨基因比较绝对均值；未新增注释或测量覆盖验证。
- Figure3c左为MAGMA preranked GSEA，右为pseudobulk DL GSEA，两个输入和通路集合不相同；历史ORA保留Table19，不能再写成主图的配对GSEA/ORA。
- ORA按作者现决定保持：overlap≥3在先、数据库/输入定义内BH在后。原存档README的PENDING文字为旧记录；本包不重新开启此待办，也未重算FDR。
- 空间49切片→35分析单位的原分组已经作者确认。保留实际单位汇总，不将旧34个展示ID改写成49或60位供者。代码原记录metadata_status是当时状态，不应覆盖作者新确认。
- Figure5a保留CA1 1,779 spots与T1077 1,363 spots，H&E/推断Fib/FN1描述背景；Fib为各状态边际q05之和。Figure5b为60个gene–unit估计、34个记录ID；Figure5c为29个候选的条件丰度关联。
- 空间条件关联、表达检出、空间结构与组织学层次是不同问题。本包没有新建病理验证。
- 本包内历史分支、测试fixture、旧输出名（包括CKB）仅为代码血缘，不自动成为当前正文采用结果。正式候选和图件以当前Source Data与采用索引为准。
''')
(B/'04_已发现的复现限制.md').write_text('''# 复现限制与本轮处理

1. **旧包漏收上游材料**：本轮从服务器只读恢复原QC notebook的代码/执行计数、纠正counts代码，以及257个scIB代码/配置/小结果文件。每项记录源路径与SHA256；没有下载大矩阵，没有执行其生产入口。
2. **混合notebook不是一键流水线**：CPU原作者notebook内含历史theta1和未执行单元。正式theta2参数以pipeline_parameters.json及最终manifest连接。保留所有单元，不伪称它们全部执行过。
3. **绘图路径清单遗漏**：继承的current/resolve_paths.py期待native_plot_code/merge_records/PATH_MAP.json，但旧包没有这个文件。本轮把原始PATH_MAP恢复到recovered_sources，未擅自改写所有历史路径。完整从Source Data跨机重画仍需显式重定位和顺序验证。
4. **绘图有两种不同复现级别**：本轮实际测试15张采用SVG→PDF导出，并与原PDF同尺寸渲染比较。它不证明全部历史数据→panel→拼图链已重跑；后者本轮未执行，也不能自动覆盖稿件。
5. **大输入与外部依赖**：原始counts、最终h5ad、LD参考/QTL、GWAS大文件及cell2location模型外置，按原manifest取用。代码中有服务器绝对路径、R环境与原生MAGMA/CUDA依赖，需先对照receipt再运行，不能仅pip install后声称全流程可执行。
6. **旧说明失效**：原存档README中ORA和供者未核实措辞保留为原件，当前解释以03_当前决定与科学口径.md为准。未修改原科学结果或重新开启已确认问题。
7. **本轮测试覆盖的含义**：语法检查不验证算法；合成小测试验证已选函数/入口，不替代真实全量复算；64表和关键图表数值对照验证存储结果契合，不代表新的生物学验证。其余模块已收代码、参数与结果索引，未逐模型重新拟合。
''')
(B/'05_安全审核与运行说明.md').write_text('''# 审核运行入口

先阅读00_README_CN.md、01导航及audit/MODULE_CODE_RESULT_CROSSWALK.tsv。

可在解压目录运行的只读证据检查（只写本包audit输出）：

```bash
python audit/code/audit_review.py
python audit/code/check_workbook.py
```

依赖：Python、numpy、pandas、PyMuPDF；第二项另需openpyxl。check_workbook从真实xlsx逐格比对机器TSV，保留NA/0。

母版导出（不是从数据重新绘图）：

```bash
python figure_code/render_current_svg.py --out /path/to/new/review_export
```

需要CairoSVG、PyMuPDF及Cairo原生库。Mac本次测试补充了DYLD_FALLBACK_LIBRARY_PATH指向已存在的Cairo库；没有修改原分析环境。只写新输出目录，不覆盖正式图。

合成测试执行记录见audit/logs与audit/TEST_EXECUTION.md。不要运行tests/run_test_suite.py中的生产配置入口，也不要执行controller、queue、launch、dispatch脚本来“看看能不能跑”。这些脚本会启动真实工作，不属于本次审核。

`analysis_code/`与`figure_code/`是源代码原件，保留原路径和历史常量；`recovered_sources/`是本轮补收原件。没有全图谱h5ad，也没有凭空创造一个统一run_all假装已端到端验证。正式复算应先指定模块、输入hash和新输出路径。
''')
# R/Python source audit language for code review users
(B/'audit/TEST_EXECUTION.md').write_text('''# 本轮实际测试

2026-09-22，从本投稿包原件复制到服务器临时目录`/tmp/panvas_code_review_20260922`；CUDA_VISIBLE_DEVICES为空，OMP/MKL/OpenBLAS最多2线程；不读取真实atlas。

```bash
cd /tmp/panvas_code_review_20260922/downstream_base
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 CUDA_VISIBLE_DEVICES="" /home/anzhen_vhd/miniconda3/envs/omicverse/bin/python -B -m unittest tests.test_core tests.test_release_compat tests.test_v2_1_projection tests.test_models -v
/usr/bin/Rscript /tmp/panvas_code_review_20260922/primary/code/native.R /tmp/panvas_code_review_20260922/primary test
/usr/bin/Rscript /tmp/panvas_code_review_20260922/enrichment/code/meta_DL.R /tmp/panvas_code_review_20260922/enrichment test
```

第一组76测试PASS（含counts流式聚合、缺失分母、annotation契约、映射/发布完整性、模拟edgeR方向与奇异设计、meta、比例模型）。另外原稿主策略自己的原生QLF/配对块/REML数值小测试PASS，DL与metafor对照PASS。

后两者只复制各自code及evidence/original文件到临时目录，使用脚本明确的`test`分支；没有执行fit/run。测试中名为CKB的人工阳性控制不是本项目CKB的新分析。

Python AST / R parse / shell bash -n仅语法级。所有实际结果详见日志，不复用历史PASS冒充本轮测试。MAGMA/scDRS、整合、真实DE、CoVarNet、通讯、空间模型均未重跑。
''')
print('modules',len(rows),'missingcode',[x['stage'] for x in rows if not x['code_path_exists']])
