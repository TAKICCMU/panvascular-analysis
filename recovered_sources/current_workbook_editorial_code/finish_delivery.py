from pathlib import Path
import csv,json,gzip,shutil,hashlib
from openpyxl import load_workbook
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919'); O=P/'图件编辑资料/补充表与方法同步_20260921';B=O/'轻量备份';S=O/'待替换'
def backup(p):
 q=B/p.relative_to(P);q.parent.mkdir(parents=True,exist_ok=True)
 if p.exists() and not q.exists():shutil.copy2(p,q)
def update(p,text):backup(p);p.write_text(text)
def tsv(p,rows,fields):
 with p.open('w') as f:w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)
# Internal evidence ledger distinguishes executed evidence from availability gaps.
rows=[]
def add(parameter,value,source,scope,used,limitation=''):
 rows.append(dict(parameter=parameter,adopted_value=value,source=source,scope=scope,used_in_existing_results=used,limitation=limitation))
qc='/home/anzhen_vhd/scRNAseq/数据预处理.ipynb (cell 3; executed PRJNA562645_1 log)'
old='/Users/taki/Desktop/Codex任务1/final_submission_update_20260826_work/originals/Manuscript.docx (QC Methods)'
for a,v in [('detected genes','500 <= n_genes_by_counts <= 8000'),('raw UMI','1000 <= total_counts <= 100000'),('mitochondrial','pct_counts_mt < 20; MT- prefix'),('ribosomal','pct_counts_ribo <= 30; RPS/RPL prefixes'),('minimum technical sample','>=50 cells after hard filters, before doublet filtering'),('doublets','OmicVerse Scrublet per batch; automatic thresholds; filter_doublets=True')]:add(a,v,qc+'; '+old,'upstream QC','workflow source plus retained metadata evidence','All retained cells satisfy stored hard-QC thresholds; complete 247-sample historical execution logs/platform exceptions not recovered. Automatic threshold confirmed in available executed example, not individually revalidated for every sample.')
add('current reconstruction','post-QC mother -> remove identified excluded runs; 1,175,115 retained cell IDs','panvascular_atlas_rebuild_v2_20260907/code/05_build_corrected_full_gene_counts.py; 01_corrected_counts/audit/corrected_counts_build_status.json','current atlas','Yes; saved build provenance','No QC rerun in this editing task')
base='Public_Release/Analysis_Code/retained_genetics/Genetics/'
prep=base+'09c_susieRSS_finemapping_Jupyter一键运行版.py';fit=base+'09e_susieRSS_mismatch_aware_rerun_Jupyter一键运行版.py'
for a,v,s in [('region/reference','FUMA risk loci +/-100 kb; GRCh37 EUR 1000G, 502 individuals; reference MAF 0.01',prep),('variant limits','at least 20 aligned polymorphic variants; preparation stops above 20,000',prep+'; '+fit),('alignment','lowest-P duplicate retained; counted-allele Z alignment; palindromic/unresolved/reference-absent excluded; mean dosage imputation; zero variance excluded',prep+'; '+fit),('LD','standardized dosage cross-product/(N-1); symmetric; diagonal=1; finite Z/LD subset',fit),('version/function','susieR 0.16.6; susie_rss',fit+'; source version output'),('primary implementation','finite502_EB_z_only uses Z AND reference LD; n=NULL (no GWAS N), not LD-free',fit),('model parameters','L=min(10,p); R_finite=502; R_mismatch=eb; estimate_residual_variance=False; max_iter=200; tol=0.001; coverage=0.95; min_abs_corr=0.5',fit),('refine','TRUE attempted; FALSE only after fit error',fit),('region acceptance','converged and R_reliability_flag=False',fit),('credible set acceptance','accepted region plus within-set min_abs_corr>=0.5',fit),('recorded but not extra filters','Q_art/artifact_flag/R_sensitivity; PIP>=0.5 or >=0.9 display counts are not acceptance gates',fit),('sensitivity scope','median-N and single-effect diagnostics separate; not substituted for primary CS acceptance',fit)]:add(a,v,s,'primary / leave-CAD','Yes; actual production source')
add('credible-set deduplication','no additional project-level set deduplication gate found',fit,'credible sets','No additional gate added','Internal historical susieR default for deduplication not independently recovered; do not infer from current software defaults.')
add('adopted counts','primary 155 regions, 136 with accepted sets, 156 accepted sets; leave-CAD 67 attempted, 35 with accepted sets','Source_Data_Figures/Figure2/all_region_summary_with_acceptance.tsv; native genetics_assets/credible_sets.tsv','read-only result recount','Yes')
add('ED8g PIP rows','Primary: >=0.5 99, >=0.9 37; leave-CAD: >=0.5 45, >=0.9 22','Public_Release/Figure_Code/native_dependencies/Panvascular_FigureYa_20260919/code/genetics_assets/PIP_counts.tsv','variant-region rows, native display family','Yes','Not claimed as globally unique SNPs or used as CS gates')
tsv(O/'QC与Fine_mapping参数来源.tsv',rows,list(rows[0]))
# Only authorized files, original same-name destinations.
for name in ['Supplementary_Tables.xlsx','Manuscript.docx','Supplementary_Information.docx','Figure_Legends.docx']:
 dest=P/name if name!='Figure_Legends.docx' else P/'Nature_Submission'/name
 backup(dest);shutil.copy2(S/name,dest)
for name in ['Manuscript','Supplementary_Information','Figure_Legends']:
 dest=P/'Nature_Submission'/f'{name}.pdf';backup(dest);shutil.copy2(O/'PDF预览'/f'{name}.pdf',dest)
backup(P/'column_name_map.tsv');shutil.copy2(O/'column_name_map.tsv',P/'column_name_map.tsv')
# TSV compatibility: preserve original machine header row; update titles and approved semantic display values.
R=P/'Public_Release/Supplementary_Tables_TSV';maps=list(csv.DictReader((O/'column_name_map.tsv').open(),delimiter='\t')); by={}
for r in maps:by.setdefault(r['Worksheet'],[]).append(r['Original column'])
w=load_workbook(S/'Supplementary_Tables.xlsx',read_only=True)
for sheet in w:
 out=R/(sheet.title+'.tsv.gz');backup(out)
 with gzip.open(out,'wt',newline='') as f:
  writer=csv.writer(f,delimiter='\t')
  for i,row in enumerate(sheet.values,1):writer.writerow(by[sheet.title] if i==3 and sheet.title in by else row)
w.close();shutil.copy2(O/'column_name_map.tsv',R/'column_name_map.tsv');shutil.copy2(O/'formatted_structure.json',R/'table_structure.json')
(R/'README.md').write_text('''# Supplementary table machine-readable entry

Numbered TSV files preserve original machine headers on row 3 and the original column/record order. Rows 1–2 contain the current title and description. Scientific numbers, missingness and fit states are retained. Approved genetic-domain display labels are Aortic/Atherosclerotic; code IDs remain traceable through column_name_map.tsv.

The submission Excel uses readable headers; use the adapter below to restore original field names. table_structure.json gives exact data start/end rows and excludes trailing notes. Information and Data dictionary are navigation aids, not additional numbered tables. Empty cells are not zero; Unknown remains a cell label.
''')
(R/'read_supplementary_table.py').write_text('''from pathlib import Path
import csv, json
from openpyxl import load_workbook

def read_table(workbook, number, metadata_dir=None):
    """Return original machine headers and typed rows without changing scientific values."""
    metadata_dir = Path(metadata_dir or Path(__file__).parent)
    sheet = f"Table {number:02d}"
    with (metadata_dir / "column_name_map.tsv").open() as f:
        mapping = [r for r in csv.DictReader(f, delimiter="\\t") if r["Worksheet"] == sheet]
    structure = next(r for r in json.loads((metadata_dir / "table_structure.json").read_text()) if r["sheet"] == sheet)
    wb = load_workbook(workbook, read_only=True, data_only=False)
    ws = wb[sheet]
    actual = list(next(ws.iter_rows(min_row=3, max_row=3, values_only=True)))
    assert actual == [r["Display column"] for r in mapping], "Header mismatch; stop instead of silently remapping"
    rows = list(ws.iter_rows(min_row=structure["data_start"], max_row=structure["data_end"], values_only=True))
    wb.close()
    assert len(rows) == structure["rows"]
    return [r["Original column"] for r in mapping], rows
''')
# Update current explanatory entry points, not archived scientific code/results.
f=P/'Nature_Submission/AUTHOR_ACTIONS_CN.md';t=f.read_text();t=t.replace('ORA统计有效性、未核实供者身份和测量覆盖限制仍为PENDING；未运行新统计。作者需决定投稿前解决还是保留为明确限制，不因打包称为已解决。','作者已接受现有ORA方案和空间分组，本轮不再列为待修复事项。保留实际分析单位、单细胞metadata缺失及测量覆盖限制；正式投稿声明仍需作者确认。');update(f,t)
f=P/'Source_Data_Figures/Current_Display_20260920/本轮修改说明.md';t=f.read_text();i=t.index('**仍待审核：**');t=t[:i]+'**当前状态（2026-09-21）：**以上为绘图来源记录。Figure1g当前采用逐基因min–max B显示，原均值仍保留。作者已接受现有ORA方案和空间分组；49切片、35分析单位及60项估计/34个ID的计算口径不变。模型丰度不等于组织学染色，图件采用不增加组织学验证。本轮未修改图件。\n';update(f,t)
f=P/'README_当前状态.txt';t=f.read_text().replace('31页文字＋15页图','30页文字＋15页图').replace('64表数据单元未改。','64表科学数据保留；标题、列名、导航和显示格式已整理，629个数值文本转为数值存储。').replace('ORA、供者身份、测量覆盖等仍保留PENDING，不因组装自动解决。','作者已接受现有ORA方案和空间分析分组；不再列为待修复事项。单细胞缺失信息、测量覆盖和模型解释边界按实际方法保留。');t+='\n本轮补充表/Methods更新见 补充表与方法修改清单.md；完整列名映射见column_name_map.tsv。内部来源与轻量备份在图件编辑资料/补充表与方法同步_20260921。\n';update(f,t)
log='''# 本轮补充表与方法同步

- 已在同一投稿包替换Excel、正文、Supplementary Information和独立Figure Legends；未建立另一套桌面投稿包。
- 64表编号、顺序和完整记录保留。新增可点击Information目录及非编号Data dictionary；列名映射见column_name_map.tsv。Table28按真实内容标为post-QC细胞数，Table49准确描述队列贡献。Table18的P值为描述性中位数，未写成合并检验P。
- 逐项比较2,394,224个数据单元，科学数值、方向、缺失及逻辑状态保留；629个数值文本转为可排序数值，转换记录可追溯。未改变任何P/FDR或重跑科学分析。
- ORA保留真实的overlap≥3在先、数据库/输入集内BH在后；清除正式文件中过期待审任务句。作者已确认空间分组，保留49切片→35分析单位及60项估计/34个ID，删除分组待修复判断。
- QC补充500–8,000基因、1,000–100,000 counts（含端点），mt<20%、ribo≤30%，doublet前每技术样本≥50细胞及逐样本Scrublet。当前保留细胞的存储QC值均符合阈值。完整247样本历史Scrublet日志/平台例外未找齐；不声称本轮逐样本重新验证QC。
- Fine-mapping补充输入对齐、LD来源、实际SuSiE参数及收敛/reliability/min_abs_corr接受条件；保留155/136/156及67/35原计数。内部z_only名称仅指未提供GWAS N，仍使用LD。未新增PIP或诊断筛选条件。
- 遗传领域展示名统一Aortic/Atherosclerotic；ED9b/c改为节点对散点和median菱形，保留真实研究纳入分母。
- 64表全部表头和首/中/尾/最长文本抽样已检查；表格预览为脚本生成的字段对照，非原生Excel截图。Word已实际渲染看图；格式、引用字段和86条参考文献XML保留。全部15幅图未修改。

修改明细、数值核对、转换记录及QC/fine-mapping来源表均在“图件编辑资料/补充表与方法同步_20260921”。所有编辑供作者最终审阅，未上传期刊或公开发布。
''';update(P/'补充表与方法修改清单.md',log)
status=json.loads((P/'投稿核查/FINAL_PACKAGE_STATUS.json').read_text());status['pending']=['Official reporting form and author declarations','Corrected release deposition DOI/access','Documented single-cell metadata missingness and measurement-coverage limitations'];status['author_accepted']=['Existing ORA workflow/results','Existing spatial analysis-unit grouping'];status['supplementary_edit_receipt']='图件编辑资料/补充表与方法同步_20260921/workbook_validation.json';update(P/'投稿核查/FINAL_PACKAGE_STATUS.json',json.dumps(status,ensure_ascii=False,indent=2))
print('Replaced current documents; preserved original-header TSV entry; updated current status only.',flush=True)
