from pathlib import Path
import csv,hashlib,io,json,os,zipfile,time
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919')
OUT=P/'投稿包_供Agent审核_500MB以内.zip'
TMP=OUT.with_suffix('.zip.tmp')
LIMIT=500_000_000
root_files=['Manuscript.docx','Cover_Letter.docx','Supplementary_Information.docx','Supplementary_Tables.xlsx','全部图件_15页审阅.pdf','图件采用版本.tsv','README_当前状态.txt','SUBMISSION_FILE_LIST.tsv','UPLOAD_MANIFEST.tsv','column_name_map.tsv','补充表与方法修改清单.md']
roots=['Figures','Extended_Data','Nature_Submission','Source_Data_Figures','Supplementary_Table20_full_results','Public_Release','投稿核查']
files=[P/n for n in root_files];excluded=[]
for folder in roots:
 for p in (P/folder).rglob('*'):
  if not p.is_file():continue
  rel=p.relative_to(P)
  reason=None
  if p.name=='.DS_Store' or '__pycache__' in rel.parts:reason='Operating-system/cache file'
  elif 'document_previews' in rel.parts or 'render_test' in rel.parts:reason='Redundant local rendering preview'
  elif folder=='投稿核查' and p.name in ['manuscript_before_text.txt','manuscript_changes.json']:reason='Historical manuscript snapshot; current manuscript and correction summary included'
  elif p.suffix.lower() in ['.h5ad','.h5','.hdf5']:reason='Full expression/matrix object excluded'
  elif rel.parts[:2]==('Public_Release','Derived_GWAS') and p.suffix=='.gz':reason='Large full-variant derived GWAS input; schema, source SHA256, manifests and figure/table results included'
  elif p==Path(__file__):reason='Review-package builder; not scientific analysis code'
  if reason:excluded.append({'path':str(rel),'bytes':p.stat().st_size,'reason':reason})
  else:files.append(p)
files=sorted(set(files))
def sha_file(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
manifest=[]
for p in files:
 manifest.append({'path':str(p.relative_to(P)),'bytes':p.stat().st_size,'sha256':sha_file(p)})
readme='''# 请先读：当前投稿包的独立Agent审核版

本包是当前采用投稿包的精简审核副本，不是新的科学结果版本。原目录：/Users/taki/Desktop/投稿包_最新版_20260919。压缩文件严格小于500,000,000字节；本轮已同步正文、Methods、图注与补充表；没有重新计算科学结果或修改图件。

## 建议阅读顺序
1. Nature_Submission/Manuscript_with_15_Figures.pdf：当前正文＋完整15幅图，保留原生尺寸。
2. Manuscript.docx、Nature_Submission/Figure_Legends.docx、Supplementary_Information.docx、Supplementary_Tables.xlsx。
3. 投稿核查/核对与完成说明_CN.md 和 Nature_Submission/AUTHOR_ACTIONS_CN.md：已完成核验及尚未解决事项。
4. Source_Data_Figures/README.md 和 SOURCE_FILE_INDEX.tsv：数值与当前panel对应；完整未筛选富集结果在Supplementary_Table20_full_results/。
5. Public_Release/Analysis_Code/README.md、MODULE_INDEX.tsv：实际原生分析代码与环境/参数/回执；Figure_Code/README.md：当前绘图代码、可编辑SVG与已测试导出入口。

## 本轮审核任务
请依据本包实际文件独立审阅，不把历史“PASS”“采用”或本助手的完成记录当作科学正确性的证明。
- 核对正文Results/Methods、图表引用顺序、图注、Table/Supplementary内容和当前Source Data是否一致。
- 核对统计单位、分母、样本/供者独立性、效应与CI、P/FDR家族、数据资格规则和科学措辞；区分表达、条件丰度关联、程序关联、交互与组织学验证。
- 审阅代码的实际实现与所写Methods是否等价。说明代码检查、直接数值检查和无法在当前包验证的事项，不把语法通过当作分析验证。
- 在原生页面尺寸下检查图件可读性；不要以聊天缩略图或仅PDF文本解析判定视觉错误。
- 输出有文件/页码/panel/表格行或代码行证据的发现，区分：必须修正、建议修正、需要补充数据才能判断、未发现问题。提供具体修正建议；不要为了找问题重新定义候选或统计口径。
- 不自动修改文件、不启动全量模型、不重新注释或重跑富集；需要新科学检验时提出最小方案并标记PROPOSED_NOT_RUN。

## 当前应保留的口径
- Figure1g：14类×28基因，已存逐基因min–max色值，面积为原检出率；不是新注释验证。A原均值展示仍在Supporting_Visuals。
- Figure3a：Fibromyocyte与Smooth muscle不同；实心为两方法原FDR均<0.05，空心不表示两项均不显著。
- Figure3c：MAGMA遗传GSEA＋原稿策略pseudobulk GSEA；不是旧七通路GSEA＋ORA。
- Figure5：FN1两张切片的描述性三视图；FN1/SCARB1共60 gene–unit estimates、34 recorded unit IDs，不能把60项估计当作60个单位；空间分组按作者确认的35分析单位方案保留。Fib图为各状态边际q05之和。单位β±SE与总体95%CI不可混称。
- ED4为EC＋T/NK；ED5为髓系。pseudobulk主分析为原稿策略复刻版，三背景分支是补充。
- 作者已接受现有ORA方案与空间分组；实际统计单位和方法保留。单细胞metadata缺失、测量覆盖与组织学解释边界按实际证据保留。不得把编辑组装称为新增统计验证、已投稿或已公开。

## 已包含与有意排除
包含当前正文/投稿信、15张PDF与SVG、64张补充表、完整Table20富集记录、图件Source Data、代码/配置/环境与运行回执、轻量冻结注释和坐标、当前核对记录。
不包含全量图谱h5ad、原始表达矩阵、大型模型、历史投稿包、重复的六个上传ZIP，也不包含三个大型全变异GWAS统计文件。GWAS的真实来源、schema及SHA256仍保留，相关图表结果也保留。详见REVIEW_EXCLUDED_FILES.tsv。
因此，本包适用于文稿/图表/代码/证据一致性审核，不足以离线从原始counts重跑全部科学分析。若某项结论必须依赖被排除的原始输入才能验证，请明确报告，不能假定已验证。

代码保留原始绝对路径以追溯。解压后请优先使用本包相对路径；不要自动访问他人的路径或启动历史生产launcher。历史脚本/回执是来源证据，不是另一版当前投稿稿件。

## 完整性
REVIEW_CONTENT_MANIFEST.tsv是本审核ZIP中实际纳入文件的SHA256清单。包内原SUBMISSION_FILE_LIST/UPLOAD_MANIFEST是完整版来源记录，不代表这里包含其列出的全部大型上传包。额外说明文件为本次审核入口，不改科学数据。
'''
def tsv(rows,fields):
 s=io.StringIO();w=csv.DictWriter(s,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows);return s.getvalue()
extras={
'00_START_HERE_请先阅读.md':readme,
'REVIEW_CONTENT_MANIFEST.tsv':tsv(manifest,['path','bytes','sha256']),
'REVIEW_EXCLUDED_FILES.tsv':tsv(excluded,[ 'path','bytes','reason']),
'REVIEW_SCOPE.json':json.dumps({'created_local':'2026-09-21','base_package':str(P),'purpose':'independent manuscript, figure, data and code audit','max_bytes':LIMIT,'scientific_reanalysis':False,'full_atlas_included':False,'full_variant_GWAS_included':False,'n_original_files':len(files)},ensure_ascii=False,indent=2)}
with zipfile.ZipFile(TMP,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True,strict_timestamps=False) as z:
 for p in files:z.write(p,str(p.relative_to(P)))
 for name,data in extras.items():z.writestr(name,data)
size=TMP.stat().st_size
assert size<LIMIT,(size,LIMIT)
with zipfile.ZipFile(TMP) as z:
 assert z.testzip() is None
 for row in manifest:
  assert hashlib.sha256(z.read(row['path'])).hexdigest()==row['sha256'],row['path']
 assert len(z.namelist())==len(manifest)+len(extras)
 assert all(f'Figures/Figure{i}.pdf' in z.namelist() for i in range(1,6))
 assert all(f'Extended_Data/ExtendedData_Figure{i}.pdf' in z.namelist() for i in range(1,11))
os.replace(TMP,OUT)
result={'file':str(OUT),'bytes':size,'MB_decimal':round(size/1e6,2),'MiB':round(size/1048576,2),'files':len(files)+len(extras),'CRC':'PASS','all_original_file_SHA256':'PASS','sha256':sha_file(OUT)}
(P/'投稿核查/AGENT_REVIEW_PACKAGE_VALIDATION.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
print(json.dumps(result,ensure_ascii=False),flush=True)
