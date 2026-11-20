from pathlib import Path
import json,hashlib,html,re,shutil,sys
import pandas as pd,numpy as np
from lxml import etree as E
import pymupdf as fitz
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');O=P/'独立小图_最新版'
files=['native_panel_manifest.json','redraw_manifest.json','evidence_redraw_manifest.json','context_manifest.json','map_reflow_manifest.json','genetic_matrix_manifest.json','highres_manifest.json'];rows={}
for fn in files:
 for r in json.loads((O/'review'/fn).read_text()):rows[(r['figure'],r['panel'])]=r
assert len(rows)==72
sourcehash={}
for f in (O/'review').glob('*source_hashes.json'):sourcehash.update(json.loads(f.read_text()))
checks=[]
for f,h in sourcehash.items():checks.append({'check':'source_read_only','item':f,'pass':Path(f).is_file() and hashlib.sha256(Path(f).read_bytes()).hexdigest()==h})
reg=pd.read_csv(P/'图件采用版本.tsv',sep='\t')
for _,r in reg.iterrows():
 for ext in ['pdf','svg']:
  path=Path(r['package_'+ext]);checks.append({'check':'adopted_master_unchanged','item':str(path),'pass':hashlib.sha256(path.read_bytes()).hexdigest()==r['package_'+ext+'_sha256']})
index=[];outside=[]
for (figure,panel),r in rows.items():
 name=figure+panel;folder=O/figure;pdf=folder/(name+'.pdf');svg=folder/(name+'.svg');png=folder/(name+'.png');assert all(p.exists() for p in [pdf,svg,png])
 doc=fitz.open(pdf);assert len(doc)==1;doc.set_metadata({});doc.del_xml_metadata();tmp=pdf.with_suffix('.tmp.pdf');doc.save(tmp,garbage=4,deflate=True);doc.close();tmp.replace(pdf);doc=fitz.open(pdf);page=doc[0]
 spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans'] if s['text'].strip()];font_sizes=[s['size'] for s in spans]
 for s in spans:
  if not page.rect.contains(fitz.Rect(s['bbox'])):outside.append({'panel':name,'text':s['text'],'bbox':s['bbox']})
 sr=E.parse(str(svg));assert sr.xpath('//*[local-name()="text"]'),name
 external=[v for el in sr.iter() for k,v in el.attrib.items() if k.endswith('href') and not v.startswith(('#','data:'))]
 assert not external,(name,external)
 sourcefiles=r.get('source_files',[r.get('source_svg','')]);sourcefiles=[x for x in sourcefiles if x]
 index.append({'figure':figure,'panel':panel,'panel_name':name,'method':r['mode'],'width_mm':round(page.rect.width/72*25.4,2),'height_mm':round(page.rect.height/72*25.4,2),'font_p10_pt':round(float(np.percentile(font_sizes,10)),2) if font_sizes else None,'pdf':str(pdf.relative_to(O)),'svg':str(svg.relative_to(O)),'preview':str(png.relative_to(O)),'source_files':';'.join(sourcefiles),'source_objects':';'.join(r.get('source_objects',[])),'change':r['change'],'pdf_sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),'svg_sha256':hashlib.sha256(svg.read_bytes()).hexdigest(),'data_integrity':'existing values/identities retained; no scientific analysis rerun','layout_legibility':'standalone labels, legends, geometry reviewed','visual_review':'agent render review completed; author approval pending'})
index.sort(key=lambda r:(r['figure'].startswith('Extended'),int(re.search(r'Figure(\d+)',r['figure'])[1]),r['panel']))
assert all(r['pass'] for r in checks)
pd.DataFrame(index).to_csv(O/'逐panel来源与修改清单.tsv',sep='\t',index=False)
pd.DataFrame(checks).to_csv(O/'review/READ_ONLY_CHECKS.tsv',sep='\t',index=False)
(O/'review/TEXT_BOUNDS.json').write_text(json.dumps(outside,ensure_ascii=False,indent=2))
# Convenient, separate state identity tables for maps formerly read alongside marker matrices.
M=pd.read_csv(P/'Public_Release/Figure_Code/native_dependencies/Panvascular_Marker_Communication_Enrichment_20260920/source_data/display_marker_values.tsv',sep='\t')
for fig,nr in [('ExtendedData_Figure3',41),('ExtendedData_Figure4',9),('ExtendedData_Figure5',38)]:
 sizes=M.groupby('vmax').state_id.nunique();caps=sizes[sizes.eq(nr)].index.tolist();sub=M[M.vmax.isin(caps)][['state_id','display_label']].drop_duplicates()
 if fig=='ExtendedData_Figure4':sub=pd.concat([sub,M[M.vmax.isin(sizes[sizes.eq(16)].index)][['state_id','display_label']].drop_duplicates().query('state_id != display_label')]).drop_duplicates()
 sub.to_csv(O/fig/'地图编号与完整状态名称.tsv',sep='\t',index=False)

# Single-panel-at-a-time browser: navigation only, not a combined figure.
payload=json.dumps(index,ensure_ascii=False)
page='''<!doctype html><html lang="zh"><meta charset="utf-8"><title>泛血管图谱 · 独立小图审阅</title>
<style>body{margin:0;font:15px -apple-system,BlinkMacSystemFont,Arial;background:#edf0f3;color:#203039}header{background:#fff;padding:15px 25px;border-bottom:1px solid #ccd3d9;position:sticky;top:0}h1{font-size:20px;margin:0 0 10px}select,button,a{font:inherit;margin-right:12px}button,select{padding:7px;border:1px solid #bdc8d0;border-radius:5px;background:white}main{padding:20px;max-width:1300px;margin:auto}#stage{padding:16px;background:white;overflow:auto;text-align:center}#plot{max-width:100%;height:auto}#desc{line-height:1.6;font-size:13px;color:#465964}a{color:#245a80}small{color:#61737f}</style>
<header><h1>泛血管图谱 · 72个独立小图</h1><select id="fig"></select><select id="pan"></select><button id="prev">上一张</button><button id="next">下一张</button><a id="pdf">PDF</a><a id="svg">可编辑SVG</a><small>本页仅导航，未拼成整图；每张图等待作者审阅。</small></header><main><p id="desc"></p><div id="stage"><img id="plot"></div></main><script>
const rows=PAYLOAD;let k=0;const fs=[...new Set(rows.map(x=>x.figure))];const f=document.querySelector('#fig'),p=document.querySelector('#pan');fs.forEach(x=>f.add(new Option(x,x)));
function show(){let r=rows[k];f.value=r.figure;p.innerHTML='';rows.forEach((x,i)=>{if(x.figure==r.figure)p.add(new Option(x.panel_name,i))});p.value=k;document.querySelector('#plot').src=r.preview;document.querySelector('#plot').alt=r.panel_name;document.querySelector('#pdf').href=r.pdf;document.querySelector('#svg').href=r.svg;document.querySelector('#desc').textContent=r.panel_name+' · '+r.width_mm+' × '+r.height_mm+' mm · '+r.method+' — '+r.change;}
f.onchange=()=>{k=rows.findIndex(x=>x.figure===f.value);show()};p.onchange=()=>{k=+p.value;show()};document.querySelector('#prev').onclick=()=>{k=(k-1+rows.length)%rows.length;show()};document.querySelector('#next').onclick=()=>{k=(k+1)%rows.length;show()};show();</script></html>'''.replace('PAYLOAD',payload)
(O/'打开逐图审阅.html').write_text(page)
mode=pd.Series([r['method'] for r in index]).value_counts().to_dict()
(O/'README_交付说明.md').write_text(f'''# 独立小图审阅版

基线：现有投稿包《图件采用版本.tsv》登记的5幅主图与10幅Extended Data。15幅原PDF/SVG的SHA256全部吻合，且完成后再次核对未改变。

## 文件

- 共15个Figure文件夹、72个小图。每个小图均为一页PDF、可编辑SVG、PNG预览。
- 打开 `打开逐图审阅.html` 可逐张查看；本轮没有制作整图或15页合并PDF。
- `逐panel来源与修改清单.tsv` 记录每张小图的来源、处理方式、实际尺寸和文件哈希。
- `code/` 为实际运行的绘图代码。地图编号的完整身份表放在对应ED文件夹。

## 实际修改范围

- 29个panel按现有汇总数据重新绘图：marker点图、完整森林图、富集、类别/研究数、模块证据、QC/benchmark及遗传矩阵等。
- 4个panel重新安排冻结地图对象；Figure1b另用完全相同坐标、绘图顺序、Unknown显示掩码及固定颜色导出600dpi点层，左/右分别由146×122和187×156提高到870×725及1118×932像素。
- 其余38个panel保留采用版矢量/图像对象，按独立panel整理边界、可读字号和必要上下文；其中EC地图对象名、ED10b共享行名已补齐。**这些没有被称为数据重绘。**
- marker点面积与原检测比例线性对应；零检出为零面积；颜色沿用当前蓝色顺序色带。Figure1g直接读取已存储B版0–1颜色值；ED marker使用各自原mean log1p及原上限，未重新标准化。
- Figure3c仍为MAGMA GSEA与pseudobulk GSEA，两侧原9条通路与右侧63单元不变。Figure3a实心/空心、Figure3f情景/支持分类保留原定义。
- Figure5的基因、两张切片、空间图值、色限、总体模型与29候选保持。未进行空间新检验或组织学确认。

## 检查与边界

本地脚本核对：72套文件齐备且PDF均单页、SVG保留文本、没有外链图片；采用母版和读取科学表的哈希未改变；逐表保留统计单位与映射。marker检查覆盖392个主图单元及3357个ED单元。CI绘图区覆盖全部已存储端点。

实际看图：逐组查看全部72个panel的渲染；对marker、遗传矩阵、森林图、地图身份、空间图及拆图后图例进行局部检查，并修复检查中发现的裁切/遗漏。看图通过仅指本轮agent排版检查，不代表作者最终采用或科学模型有效性验证。

这些是**独立编辑与审阅的小图**，尺寸由内容决定，未承诺未来重新拼版后仍满足整幅Nature制作尺寸。后续排版时应按清单中的真实毫米尺寸检查字号，不能任意拉伸。未改正文、原15图、统计结果、P/FDR或分析代码。

正式投稿母版保持原位置。本目录为唯一当前独立小图输出，没有复制新投稿包或增加桌面历史版本。
''')
(O/'参考样式与代码.md').write_text('''# 使用的参考与实际采纳点

1. [Nature figure guide: Building and exporting figure panels](https://research-figure-guide.nature.com/figures/building-and-exporting-figure-panels/)。本轮采用：常规可编辑字体、矢量文字、紧凑但按信息密度分配空间、颜色承载数据含义、图例邻近数据、避免装饰图标。指南中的89/183 mm和整幅170 mm高度是后续整图制作要求；本轮交付独立panel，未拼整图。
2. [An organotypic atlas of human vascular cells, Nature Medicine (2024)](https://www.nature.com/articles/s41591-024-03376-x)及[Teichlab公开代码](https://github.com/Teichlab/vascular_atlas)。已阅读其marker绘图notebook的绘图单元：状态分组、独立名称区、显式分类颜色、点图表达与检出比例分工。仅学习可视化结构，不复制该研究marker或数据，也未执行其分析流程。
3. [A single-cell atlas of the human brain vasculature across development, adulthood and disease, Nature (2024)](https://www.nature.com/articles/s41586-024-07493-y)。参考资源总览、细胞身份与表达证据分层呈现；不借用其组织图。
4. [FigureYa公开绘图模板](https://github.com/ying-ge/FigureYa)作为类型导航参考；本轮实际代码为本目录Python脚本，未运行模板附带统计分析。

保留本项目既定L1/L2、模块和支持方向颜色，不为“期刊风格”重新定义科学编码。不能将“借鉴顶刊”当成期刊已认可图件或研究结论。
''')
(O/'review/FINAL_VALIDATION.json').write_text(json.dumps({'panels':72,'figures':15,'modes':mode,'input_hash_checks':len(checks),'all_read_only_checks_passed':all(r['pass'] for r in checks),'text_bounds_flags':outside,'scientific_analysis_rerun':False,'author_approval':'PENDING'},ensure_ascii=False,indent=2))
print(json.dumps({'panels':len(index),'modes':mode,'read_only_checks':len(checks),'text_outside_flags':len(outside)},ensure_ascii=False))
