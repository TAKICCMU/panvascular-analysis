from common_main import *
from resolve_paths import resolve_legacy_path
from PIL import Image
BASE=Path('/Users/taki/Desktop/Panvascular_Main_Figure_Refinement_20260920/reproducibility_dependencies/Panvascular_Figure3_Group_Spacing_20260920')
v=pd.read_csv(BASE/'FIGURE_VERSIONS.tsv',sep='\t');v['source_pdf']=v['source_pdf'].map(resolve_legacy_path);v['source_svg']=v['source_svg'].map(resolve_legacy_path);changed=['Figure1','Figure3','Figure4','Figure5'];checks=[]
for fig in changed:
 pdf=O/'figures'/f'{fig}.pdf';svg=O/'editable_svg'/f'{fig}.svg';d=fitz.open(pdf);p=d[0];assert len(d)==1
 assert abs(p.rect.width/mm-183)<.001 and abs(p.rect.height/mm-170)<.001
 bad=[]
 for b in p.get_text('dict')['blocks']:
  for ln in b.get('lines',[]):
   for s in ln['spans']:
    z=fitz.Rect(s['bbox'])
    if z.x0<-.2 or z.y0<-.2 or z.x1>p.rect.width+.2 or z.y1>p.rect.height+.2:bad.append({'text':s['text'],'bbox':s['bbox']})
 checks.append({'figure':fig,'width_mm':p.rect.width/mm,'height_mm':p.rect.height/mm,'text_outside_page':bad})
 assert not bad,(fig,bad)
 row=v.figure==fig;v.loc[row,'source_pdf']=str(pdf);v.loc[row,'source_svg']=str(svg);v.loc[row,'version']=O.name;v.loc[row,'status']='AUTHOR_REVIEW_PENDING';v.loc[row,'sha256']=sha(pdf);v.loc[row,'svg_sha256']=sha(svg)
 p.get_pixmap(matrix=fitz.Matrix(2,2)).save(O/'review'/f'{fig}_final.png')
 # Same-size before/after per affected figure, preview only; mother PDFs never scaled.
 prior=Path(resolve_legacy_path(pd.read_csv(BASE/'FIGURE_VERSIONS.tsv',sep='\t').set_index('figure').loc[fig,'source_pdf']))
 old=fitz.open(prior);pi=old[0].get_pixmap(matrix=fitz.Matrix(1.5,1.5));pn=p.get_pixmap(matrix=fitz.Matrix(1.5,1.5))
 im0=Image.frombytes('RGB',[pi.width,pi.height],pi.samples);im1=Image.frombytes('RGB',[pn.width,pn.height],pn.samples)
 im=Image.new('RGB',(im0.width+im1.width+16,max(im0.height,im1.height)),'white');im.paste(im0,(0,0));im.paste(im1,(im0.width+16,0));im.save(O/'review'/f'{fig}_before_after.png')
merged=fitz.open()
for _,row in v.iterrows():
 d=fitz.open(row.source_pdf);merged.insert_pdf(d)
 typ,grp=d.xref_get_key(d[0].xref,'Group')
 if typ=='dict':
  assert not re.search(r'\d+ \d+ R',grp)
  merged.xref_set_key(merged[-1].xref,'Group',grp)
book=O/'全部图件_15页_主图调整版.pdf';merged.save(book,garbage=4,deflate=True);merged.close();m=fitz.open(book);assert len(m)==15
for i,(_,row) in enumerate(v.iterrows()):
 d=fitz.open(row.source_pdf);assert d[0].rect==m[i].rect;assert d[0].get_pixmap().samples==m[i].get_pixmap().samples
v.to_csv(O/'FIGURE_VERSIONS.tsv',sep='\t',index=False)
# The 29 evidence y centres remain numerically identical to the corresponding forest rows.
y=pd.read_csv(O/'notes/Figure3_shared_row_check.tsv',sep='\t');assert y.delta_pt.abs().max()==0
checks.append({'merge_pages':15,'all_native_page_dimensions_and_render_equal':True,'unchanged_figures':11,'scientific_models_rerun':False,'ORA':'PENDING','donor_identity':'PENDING_AS_BEFORE','spatial_histology':'NOT_INDEPENDENTLY_CONFIRMED'})
(O/'notes/DELIVERY_CHECKS.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2))
(O/'本轮修改说明.md').write_text('''# 本轮主图调整（待作者审核）

- Figure1：16类名称分布在UMAP四周；引线交叉为0。右图改为6类血管床的Tissue着色，不放图面图例；完整配色映射存Source Data。新增g为14主要谱系×28 canonical markers（排除Unknown/Schwann），圆面积表示检出比例、颜色为原mean log1p(CP10k)，392项统计值与原表一致。marker图展示注释的表达依据，并非独立证明所有注释正确。
- Figure3c：每个富集区三条分组竖线统一左对齐，名称成列；全部通路及原数值保留。3d减少左边空白、扩大森林图区；d/e共用29个原页面行中心。3e删除Not evaluated图例（现有矩阵没有该状态），其他三类一行排列。
- Figure3f：保留原五个基因—组织上下文和两个场景；新增PP.H3及共同SNP数两列，与PP.H4互补。仅提取既有结果，圆/方形及实心/空心沿用原场景与证据分类；不是增加新的显著结论。
- Figure4c/d：原NMF权重改为水平柱形，长度、零基线、顺序与颜色不变。
- Figure5a：两行完整组织三视图：胫前动脉CA1（1,779 spots）和冠脉T1077（1,363 spots），均为H&E / Inferred Fib abundance / FN1。对象标签居中。5b中FN1图例在零线左、SCARB1在右，pooled行FDR/I²对齐；60个gene–unit estimates、34个recorded unit IDs不变。5c–f保留原对象及布局。放大空间配套页在companion。

所有主图母版183×170 mm，原生PDF/SVG可编辑；其余11幅不重绘。15页合并稿逐页尺寸及渲染与单图一致。已实际看图核对名称/引线、分组竖线、共享行、图例和空间三视图；数值与版本检查另见notes。SVG中Figure3仅新增/替换明确指定的块，Figure4仅替换10个权重符号。

**仍待审核：**本轮图件采用；空间组织学来源仍未独立确认，T1077是recorded library ID，不能称已确认供者；ORA统计有效性及既有供者身份待办继续PENDING。未重跑科学模型、未更改正文或原投稿包。各图图注建议见notes，本轮新增Figure1g及Figure3f的两列需在统一文稿合并时对应说明。
''')
# Crosswalk into unchanged scientific source data; no duplicate scientific matrices necessary.
idx=[{'figure':'Figure1','source':'source_data/Figure1g_major_lineage_marker_values.tsv','meaning':'14x28 original mean/detection values; new display panel g'},
 {'figure':'Figure1','source':'source_data/Figure1b_tissue_sample_mapping.tsv','meaning':'Read-only sample-to-major-vascular-bed display mapping'},
 {'figure':'Figure3','source':str(BASE.parent/'Panvascular_Figure3_Genetic_Transcriptomic_20260920/source_data'),'meaning':'Original9 MAGMA and9x7 pseudobulk pathways, unchanged'},
 {'figure':'Figure3','source':'source_data/Figure3d_all29_unchanged.tsv','meaning':'29 original candidates and original effects/CI'},
 {'figure':'Figure3','source':'source_data/Figure3f_same5_contexts_both_scenarios.tsv','meaning':'5 frozen gene-tissue contexts,2 original scenarios; PP.H3/4+n_overlap'},
 {'figure':'Figure4','source':str(S/'Figure4/05_top_nodes_reference.tsv'),'meaning':'Original module weight values; display geometry only'},
 {'figure':'Figure5','source':'notes/Figure5_Source_Data_index.tsv','meaning':'Two complete sections, unchanged meta/unit estimates; detailed source index'}]
pd.DataFrame(idx).to_csv(O/'SOURCE_DATA_INDEX.tsv',sep='\t',index=False)
print('PASS: 4 revised figures; 11 inherited; 15 pages same native dimensions and rendering; no text beyond page.')
