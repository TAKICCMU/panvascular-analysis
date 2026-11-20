from pathlib import Path
from zipfile import ZipFile
import xml.etree.ElementTree as E
import json,hashlib,shutil
import numpy as np,pandas as pd,fitz
O=Path(__file__).resolve().parents[1];B=Path('/Users/taki/Desktop/Panvascular_Figure5_Evidence_Focus_20260920');D=O/'Source_Data/Figure5';MSS=Path('/Users/taki/Desktop/投稿包_最新版_20260919');XF=MSS/'Supplementary_Tables.xlsx';ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
zz=ZipFile(XF);ss=[]
if 'xl/sharedStrings.xml' in zz.namelist():ss=[''.join(x.itertext()) for x in E.fromstring(zz.read('xl/sharedStrings.xml')).findall('m:si',ns)]
# Resolve table name through workbook relationships, rather than assume sheet numbering.
r=E.fromstring(zz.read('xl/workbook.xml'));sheet=next(x for x in r.findall('m:sheets/m:sheet',ns) if x.get('name')=='Table 32');rel=sheet.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id');rels=E.fromstring(zz.read('xl/_rels/workbook.xml.rels'));target=next(x.get('Target') for x in rels if x.get('Id')==rel);part=target.lstrip('/') if target.startswith('/') else 'xl/'+target
rows=[]
for row in E.fromstring(zz.read(part)).findall('m:sheetData/m:row',ns):
 vals={}
 for c in row:
  v=c.find('m:v',ns);st=c.find('m:is',ns);val=ss[int(v.text)] if c.get('t')=='s' and v is not None else ''.join(st.itertext()) if st is not None else v.text if v is not None else ''
  vals[''.join(x for x in c.get('r') if x.isalpha())]=val
 if int(row.get('r'))>=4 and vals.get('A'):
  rows.append({'gene':vals['A'],'term':vals['L'],'workbook_sheet':'Table 32','excel_row':int(row.get('r')),'excel_gene_cell':'A'+row.get('r'),'excel_FDR_cell':'Q'+row.get('r'),'excel_beta':float(vals['D']),'excel_FDR':float(vals['Q'])})
x=pd.DataFrame(rows);raw=pd.read_csv(D/'all29_original_meta.tsv',sep='\t',float_precision='round_trip');p=raw[raw.posterior.eq('q05')&raw.sensitivity.eq('all_sections')].copy();p['original_TSV_line']=p.index+2;p=p.merge(x,on=['gene','term'],validate='one_to_one');assert len(p)==87;assert np.allclose(p.beta,p.excel_beta,rtol=1e-13,atol=1e-14);assert np.allclose(p.FDR_BH_all_29_pairs_x_3_terms,p.excel_FDR,rtol=1e-13,atol=1e-16)
p['already_significant_original_FDR_lt_005']=p.FDR_BH_all_29_pairs_x_3_terms.lt(.05);p.to_csv(D/'P86_Table32_all87_original_entries.tsv',sep='\t',index=False)
hit=p[p.already_significant_original_FDR_lt_005].copy();hit.to_csv(D/'P86_Table32_10_8_1_significant_entries.tsv',sep='\t',index=False)
assert hit.groupby('term').size().to_dict()=={'cell_abundance':10,'interaction':1,'pathology':8}
u=pd.read_csv(D/'Figure5b_60_gene_unit_estimates.tsv',sep='\t',float_precision='round_trip');meta=p[p.term.eq('cell_abundance')&p.gene.isin(['FN1','SCARB1'])];meta.to_csv(D/'Figure5b_pooled_FDR_I2_original.tsv',sep='\t',index=False)
c=pd.read_csv(D/'all29_original_section_models.tsv',sep='\t',float_precision='round_trip');s=c[c.gene.eq('SCARB1')&c.posterior.eq('q05')&c.term.eq('cell_abundance')&c.section_id.eq('GSE248608::GSM7917894_CA1')].iloc[0];p1=u[u.gene.eq('SCARB1')&u.biological_unit_id.eq('GSE248608::P1')].iloc[0]
lines=[]
for term in ['cell_abundance','pathology','interaction']:
 t=hit[hit.term.eq(term)].sort_values('excel_row');items='；'.join(f"{r.gene}（第{r.excel_row}行，{r.excel_FDR_cell}）" for _,r in t.iterrows());lines.append(f'- `{term}`：{len(t)}项。'+items)
(O/'图注与引用建议_待采用.md').write_text('''# Figure5a/b图注与引用建议（未写回正文）

## Figure5a/b正式图注建议

a, Original H&E, model-inferred Fibroblast abundance and FN1 log1p(CP10k) expression in GSE248608/GSM7917894_CA1, shown separately with the same orientation, coordinates and full-image window. All 1,779 originally included spots are retained. The Fib image shows the sum of marginal q05 abundances across Fibroblast states, not the q05 of the total-abundance distribution. Observed expression zeros are grey; positive expression retains the previously fixed colour limits, and abundance retains its own scale. These views provide descriptive spatial context, not evidence of FN1 downregulation within fibroblasts, histologically confirmed layer enrichment or independent validation. The CA1 SCARB1 view and both genes' full-field displays remain in the companion material and Source Data.

b, All 60 gene–unit effect estimates for FN1 (34) and SCARB1 (26), representing 34 distinct recorded unit IDs, from the original primary q05/all-section conditional Fib-abundance model. Colour and point shape identify genes; direction is indicated by the position relative to zero. Study rows group estimates by source and are not study-level meta-analyses. All opposite-direction and nonsignificant estimates are retained. Diamonds and horizontal intervals show the original pooled estimates and 95% confidence intervals; adjacent fields report their original FDR and I². These estimates are conditional on the section-centred expression programme. Recorded IDs are not all verified independent donors. The detailed companion retains each unit identity and β ± its original SE; those bars are not confidence intervals. The former detection-fraction matrix is retained in the Figure5 Source Data, without a new Extended Data number.

c–f文字中的统计定义和科学内容不变；本轮不重新安排它们的版面。

## 配套CA1 SCARB1图的必要说明

The CA1 section-level SCARB1 coefficient is +0.053896 (nominal P=0.581553; rounded P=0.582). The combined GSE248608::P1 unit, which includes CA1 and NC1, has coefficient −0.024183. Neither the overall FDR nor the overall direction is assigned to this section image. The image remains descriptive and is not a section-level confirmation of the overall positive association.

中文：CA1图单切片P=0.582；P1合并单位方向为负；总体正向结果仍保留，但不得用于给CA1图贴显著标签。配套两页PDF保持原文件，以上说明随附而不改其数据或误差条。

## P86引用修正（只建议，不写回）

原句：In the primary q05/all-section analysis, 10 abundance terms, 8 regional-programme terms and 1 interaction met FDR<0.05 within the joint 87-test family (Fig. 5b-d and Supplementary Table 32).

建议：In the primary q05/all-section analysis, 10 conditional cell-abundance terms met FDR<0.05 (Fig. 5c; Supplementary Table 32, cell_abundance entries). Eight expression-programme terms and one abundance-by-programme interaction also met FDR<0.05 within the same 87-test family (Supplementary Table 32, pathology and interaction entries, respectively). The FN1 and SCARB1 unit-level distributions are shown in Fig. 5b, and abundance-term posterior/QC sensitivity is summarised in Fig. 5d.

原因：5b是两基因单位分布，5c是29个丰度项，5d是丰度项四设置；它们没有绘制全部程序项与交互项。原字段`pathology`在实际模型中是表达程序，不能称独立病理分区。下列行号来自当前Supplementary_Tables.xlsx的Table 32工作表，已与原TSV逐项核对；未重算P或BH。

'''+ '\n'.join(lines)+'''

完整87项位置、gene/term、原始效应与FDR见Source_Data/Figure5/P86_Table32_all87_original_entries.tsv。19项显著结果另列10/8/1索引；这只是筛出既有FDR达标条目，不是重新检验或新的候选筛选。

ORA继续PENDING，供者/文库身份待办继续保留。没有因为完成图件而将其标为解决。
''')
# Promote matrix and required spatial input assets into a formal Source Data subtree.
params=json.loads((B/'notes/display_parameters.json').read_text());image=Path(params['image']);shutil.copyfile(image,D/'CA1_original_H_and_E.png');(D/'CA1_display_parameters_unchanged.json').write_text(json.dumps(params,indent=2))
idx=pd.read_csv(O/'Source_Data_index.tsv',sep='\t');extra=[]
for fname,panel,meaning,src in [
 ('Figure5b_60_gene_unit_estimates.tsv','5b','60 gene-unit estimates, 34 unique recorded unit IDs',D/'all29_original_unit_models.tsv'),
 ('Figure5b_34_recorded_unit_IDs.tsv','5b','Identity crosswalk: 26 IDs overlap between genes; missing is not zero',D/'all29_original_unit_models.tsv'),
 ('Figure5b_point_mapping.tsv','5b','Exact original beta plus display-only vertical offset',D/'all29_original_unit_models.tsv'),
 ('Figure5b_pooled_FDR_I2_original.tsv','5b','Original overall effect/95%CI/FDR/I2, not section statistics',D/'all29_original_meta.tsv'),
 ('P86_Table32_all87_original_entries.tsv','P86 reference','Exact Table32 and original TSV entry mapping, all 87 tests',XF),
 ('P86_Table32_10_8_1_significant_entries.tsv','P86 reference','Existing significant 10/8/1 terms, unchanged original FDR',XF),
 ('CA1_original_H_and_E.png','5a','Byte-identical original image, full field',image),
 ('CA1_display_parameters_unchanged.json','5a; companion','Inherited normalization, scales, metadata; CA1 no new crop',B/'notes/display_parameters.json')]:
 extra.append({'panel':panel,'file':'Source_Data/Figure5/'+fname,'source':str(src),'source_sha256':sha(src),'output_sha256':sha(D/fname),'meaning':meaning})
idx=pd.concat([idx,pd.DataFrame(extra)],ignore_index=True);idx.to_csv(O/'Source_Data_index.tsv',sep='\t',index=False)
(D/'README.md').write_text('''# Figure5正式Source Data增量交付

本目录独立保存，未覆盖原投稿包。Source_Data_index.tsv位于交付根目录，列出panel、原始来源及SHA256。

- CA1_existing_spot_values.tsv：全部1,779 spots，原始FN1/SCARB1表达与Fib丰度均保留；正式5a只绘FN1，SCARB1未被删除。
- CA1_original_H_and_E.png与CA1_display_parameters_unchanged.json：原图、固定色限、full-field方向、原scale信息。
- Figure5b_60_gene_unit_estimates.tsv / Figure5b_34_recorded_unit_IDs.tsv：60项估计，34个不同ID；不是60个单位或34名已确认独立供者。SCARB1的26个ID都在FN1集合内。
- Figure5b_pooled_FDR_I2_original.tsv：两个pooled行新增字段的原值；CI为总体原CI，不是单位SE。
- Former_5b_detection_fraction.tsv及Former_5b_detection_matrix.pdf/SVG：原检出率矩阵的正式可追溯保留位置，含全部7研究×29基因；不新增ED编号，不只留内部review。
- all29_original_*：原全部meta、单位与切片结果；Fig5c_all29.tsv：完整29候选森林图输入。均复制原文件，不重算。
- P86_Table32_all87_original_entries.tsv与10_8_1索引：P86真实表格条目位置。
- module_*：5e/f原结果，布局与数值均未改变。

所有SE仍是SE；未知文库身份、ORA和组织学未确认状态保持原有待办。配套两页PDF及可编辑SVG在companion目录中，按字节沿用基线。
''')
(O/'修改与计数说明.md').write_text('''# Figure5a/b收尾交付

- **5a**：仅保留当前CA1的H&E／Inferred Fib abundance／FN1三视图。沿用原完整视野，不另裁组织，全部1,779 spots、方向、色限、原图及数值保持。零表达图例有独立位置。
- **5b**：60个gene–unit estimates来自**34个去重recorded unit IDs**。FN1为34项，SCARB1为26项；共享26个ID。不是60个不同单位，也不是34名已确认独立供者。来源行只分组，不重做研究层meta。
- pooled行新增原始FDR和I²：FN1为1.79×10⁻¹⁸／97.1%，SCARB1为0.0229／88.9%。点的颜色／形状表示基因，零线两侧表示效应方向。
- **5c–f完全沿用**：29行森林图、四设置、14模块和五项模块程序均未压缩或改动。对应SVG元素和下部PDF渲染与基线一致。
- **配套两页原文件保留**：全部单位身份、β±原SE和总体95%CI未变。CA1的SCARB1图仍在配套与Source Data；其单切片P=0.582、P1合并单位系数−0.024183明确保留。
- 原检出率矩阵已进入正式Source_Data/Figure5，含数值表及PDF/SVG，不仅存在review。
- P86的10丰度／8程序／1交互分别定位到Supplementary Table32真实行号及原TSV；b–d不再被建议承担未画出的全部模型项。图注和引用只作为待采用文字，未写回正文。

Figure5：183×170 mm。完成脚本计数/来源核验及实际看图；作者采用仍待确认。未修改其他14图或历史文件，未运行新模型、推导单位CI或调整P/FDR。ORA和供者身份待办继续PENDING。

PAUSED_FOR_FIGURE5_ADOPTION_REVIEW
''')
# Checks cover this task only, not new scientific approval.
check=json.loads((O/'notes/checks.json').read_text());doc=fitz.open(O/'figures/Figure5.pdf');pg=doc[0];spans=[s for b in pg.get_text('dict')['blocks'] if 'lines'in b for l in b['lines'] for s in l['spans']];check.update({'actual_visual_review':'Figure5_after.png inspected; labels, pooled fields and dedicated legend readable; author confirmation pending','dimensions_mm':[pg.rect.width*25.4/72,pg.rect.height*25.4/72],'minimum_font_pt':min(s['size'] for s in spans),'text_outside_page':sum(not pg.rect.contains(fitz.Rect(s['bbox'])) for s in spans),'Table32_87_entries_matched':True,'primary_significant_terms':hit.groupby('term').size().to_dict(),'CA1_SCARB1_original_section_P':float(s.p),'P1_SCARB1_original_unit_beta':float(p1.beta),'formal_matrix_promoted':True})
assert check['text_outside_page']==0 and check['minimum_font_pt']>=5.99
for p,h in check['source_hashes'].items():assert sha(p)==h
(O/'notes/checks.json').write_text(json.dumps(check,indent=2))
print('Table32 matched; source index and captions written; all source hashes unchanged.')
