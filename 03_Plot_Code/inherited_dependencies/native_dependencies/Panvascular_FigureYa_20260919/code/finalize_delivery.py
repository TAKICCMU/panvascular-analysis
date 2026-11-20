from pathlib import Path
import fitz,json,csv,hashlib,shutil,numpy as np
from PIL import Image,ImageDraw
O=Path(__file__).resolve().parents[1];P=Path('/Users/taki/Desktop/投稿包_最新版_20260919')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
seq=[('Figure'+str(i),'main','Figures') for i in range(1,6)]+[('ExtendedData_Figure'+str(i),'extended_data','Extended_Data') for i in range(1,11)]
merged=fitz.open();compare=fitz.open();rows=[]
for i,(name,folder,basefolder) in enumerate(seq,1):
 pdf=O/folder/(name+'.pdf');svg=O/'editable_svg'/(name+'.svg');old=P/basefolder/(name+'.pdf')
 assert pdf.exists() and svg.exists()
 d=fitz.open(pdf);before=fitz.open(old);assert len(d)==1
 merged.insert_pdf(d)
 # The before/after comparison preserves real mother-PDF size, with 1:1 page inserts.
 page=compare.new_page(width=before[0].rect.width+d[0].rect.width+24,height=max(before[0].rect.height,d[0].rect.height)+25)
 page.insert_text((8,14),'Before | current submission',fontsize=8)
 page.insert_text((before[0].rect.width+20,14),'After | FigureYa adaptation',fontsize=8)
 page.show_pdf_page(fitz.Rect(0,25,before[0].rect.width,before[0].rect.height+25),before,0)
 page.show_pdf_page(fitz.Rect(before[0].rect.width+24,25,before[0].rect.width+24+d[0].rect.width,d[0].rect.height+25),d,0)
 spans=[s for b in d[0].get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans'] if s['text'].strip()]
 out=[s['text'] for s in spans if s['bbox'][0]<-.2 or s['bbox'][1]<-.2 or s['bbox'][2]>d[0].rect.width+.2 or s['bbox'][3]>d[0].rect.height+.2]
 assert not out,(name,out)
 rows.append({'page':i,'figure':name,'pdf':str(pdf),'svg':str(svg),'width_mm':round(d[0].rect.width*25.4/72,3),'height_mm':round(d[0].rect.height*25.4/72,3),'pdf_sha256':sha(pdf),'svg_sha256':sha(svg),'baseline_pdf_sha256':sha(old),'minimum_extracted_font_pt':round(min(s['size'] for s in spans),2),'font_note':'4.2pt mathematical superscript on 6pt base text' if name=='ExtendedData_Figure1' else '5.5pt GWAS abbreviations retained in approved panel a; ordinary labels >=6pt' if name=='Figure2' else 'ordinary labels >=6pt','data_integrity':'CHECKED_SOURCE_AND_GEOMETRY','layout_legibility':'CHECKED_AT_NATIVE_SIZE','visual_review':'AGENT_AND_ROOT_VIEWED; AUTHOR_REVIEW_PENDING'})
mp=O/'全部图件_FigureYa_15页审阅.pdf';merged.save(mp,garbage=4,deflate=True);merged.close();compare.save(O/'review/全15图_同尺寸前后对照.pdf',garbage=4,deflate=True);compare.close()
merged=fitz.open(mp);assert len(merged)==15
render_deltas=[]
for i,r in enumerate(rows):
 d=fitz.open(r['pdf']);assert merged[i].rect==d[0].rect
 a=np.frombuffer(merged[i].get_pixmap().samples,dtype=np.uint8).astype(np.int16);b=np.frombuffer(d[0].get_pixmap().samples,dtype=np.uint8).astype(np.int16);delta=np.abs(a-b)
 assert delta.max()<=1,(r['figure'],int(delta.max()))
 render_deltas.append(dict(figure=r['figure'],maximum_8bit_channel_difference=int(delta.max()),changed_channels=int(np.count_nonzero(delta)),page_rectangle_exact=True))
with (O/'图件版本与尺寸.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
# Verify all original submission-file hashes, not just scientific numerical values.
original_checks=[]
for r in csv.DictReader((P/'文件来源与校验.tsv').open(),delimiter='\t'):
 p=P/r['file'];assert p.exists(),p
 assert sha(p)==r['sha256'],p
 original_checks.append(r['file'])
(O/'review/final_checks.json').write_text(json.dumps({'pages':15,'merged_page_geometry_equals_individual':True,'merged_render_check':render_deltas,'render_note':'Identical page bounds; some pages differ only by at most 1/255 in raster channel rounding, not text/shape displacement. No scaling or clipping applied.','baseline_manifest_files_unchanged':len(original_checks),'baseline':str(P),'template_commit':'f627917b79f28558779fb2e3ea2014cede41b89e','science_recomputed':False,'manuscript_changed':False,'author_acceptance':False,'status':'PAUSED_FOR_FIGUREYA_VISUAL_REVIEW'},ensure_ascii=False,indent=2))
# Exact downloaded template provenance and licence; do not claim manuscript data are CC licensed.
t=[]
for p in sorted((O/'templates').rglob('*')):
 if p.is_file() and p.suffix.lower() in ['.r','.rmd']:
  rel=p.relative_to(O/'templates');t.append({'template_file':str(rel),'sha256':sha(p),'url':'https://github.com/ying-ge/FigureYa/blob/f627917b79f28558779fb2e3ea2014cede41b89e/'+str(rel)})
with (O/'templates/TEMPLATE_PROVENANCE.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['template_file','sha256','url'],delimiter='\t');w.writeheader();w.writerows(t)
(O/'templates/ATTRIBUTION.md').write_text('FigureYa — Ying Ge and module authors. Lu X et al. (2025), FigureYa: A Standardized Visualization Framework for Enhancing Biomedical Data Interpretation and Research Efficiency. iMetaMed. DOI:10.1002/imm3.70005.\nOfficial source: https://github.com/ying-ge/FigureYa\nPinned commit: f627917b79f28558779fb2e3ea2014cede41b89e.\nDownloaded template code is CC BY-NC-SA 4.0; see LICENSE. Adapted plotting code retains this attribution and applicable upstream licence. Changes are documented in notes and individual scripts. Scientific datasets, tissue images, and unmodified manuscript components retain their original rights and are not relicensed by this notice. No implied endorsement by FigureYa authors.\n')
(O/'README_先看这里.md').write_text('''# FigureYa模板适配审阅版

请先看 **全部图件_FigureYa_15页审阅.pdf**；主图PDF在main，Extended Data在extended_data，同名可编辑SVG在editable_svg。原投稿包未覆盖，正文未改。本轮状态：PAUSED_FOR_FIGUREYA_VISUAL_REVIEW。

## 主要变化
- ED2–5：表达均值用完整色块，检出率用空心圆面积，保留全部4400个marker图元；深底用白色轮廓、浅底用深灰轮廓。固定状态色与细胞坐标不变。
- Figure3、ED6、ED10：CI点及表格对齐；Figure3/ED10的29候选矩阵矢量化，共享行不变。
- Figure4、ED9：模块权重改为细线＋端点，减少大色块；模块色、等级、权重和顺序保留。
- Figure5c：实际用R运行FigureYa森林图模板的适配代码，将基因、细胞状态、CI和原有FDR分栏对齐；空间示例及其他panel保持。
- 遗传图：采用保序热图、集合算法核对、紧凑图例和统计表格样式。已批准的候选流程、P/E/M定义和全部真实计数保留。
- 不适合模板替换的UMAP、组织图和已批准流程图继续沿用，不能将本轮称为全部panel由原版Rmd重跑。

## 实现与检查
Figure5c为原生R forestplot适配；其余主要为FigureYa绘图语法移植到现有Python/SVG代码。只使用绘图部分，未执行模板中的聚类、富集、统计检验或重新排序。具体来源/行号/源数据路径与哈希见notes，原始模板及许可在templates。

全部5张主图183×170mm、10张ED180×170mm。普通标签以6–7pt为主；Figure2a保留批准版的5.5pt GWAS缩写；ED1对数刻度的4.2pt为6pt基字的上标，非普通正文。各页实际尺寸见图件版本与尺寸.tsv。

本地脚本检查：源数据/冻结几何、矩阵单元数、CI端点、29行绝对中心、原投稿包哈希、合并稿与单图逐页几何一致；部分页面8位渲染通道有最大1/255舍入差异，无几何位移、缩放或裁切。实际看图检查：全部15幅图及必要局部；已修复深色底圆圈对比和遗传图旧代码导致的Venn/色条标签回退。三类检查分开记录；不代表作者已批准。

## 对照及后续
review/全15图_同尺寸前后对照.pdf：左右各为母版真实尺寸，未压缩母版高度。review还保存局部预览；不是正式投稿附图。必要图注建议见notes/图注建议汇总.md，未写回DOCX。FigureYa引用和CC BY-NC-SA 4.0源码许可详见templates/ATTRIBUTION.md及LICENSE。

不重算P/FDR，不改候选、表达、模块等级或统计定义。等待作者视觉审核后，才考虑替换当前投稿图件。
''')
(O/'notes/图注建议汇总.md').write_text('''# 本轮显示编码的必要图注建议（未写回DOCX）

1. ED2b、ED3b、ED4b/d、ED5b：Colour tiles indicate the original mean-expression statistic; circle area indicates the proportion of cells with detected expression. Outline brightness is adjusted only for contrast and carries no additional biological meaning. 各panel仍使用其原色标（ED2b为Within-lineage scaled mean，其余Mean log1p(CP10k)），不可混称。
2. Figure3d、ED6a/b：菱形表示原效应估计，横线表示原95%CI。ED10c中菱形为All sections、方形为Exclude QC-fallback；其余定义保持。
3. Figure4c/d、ED9a：点表示原NMF权重，从零起的细线仅辅助读数，组内rank与跨组共同横轴不变。
4. Figure5c：分列显示基因与细胞状态，并列出原BH FDR。Fib=Fibroblast，EC=Endothelial cell，Mac=Macrophage；置信区间、原模型及87检验家族不变。
5. 建议在作图方法/代码声明中引用FigureYa（Lu et al., iMetaMed 2025, doi:10.1002/imm3.70005），注明采用并适配其可视化模板；不能写成使用FigureYa重新完成科学分析。
''')
print(mp)
