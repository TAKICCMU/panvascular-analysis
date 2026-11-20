from common_main import *
import shutil
PK=Path('/Users/taki/Desktop/投稿包_最新版_20260919');BK=PK/'图件编辑资料/轻量备份/Figure1g_B采用前_20260921';SD=PK/'Source_Data_Figures/Current_Display_20260920';NOUT=O/'notes';E=O/'reproducibility_dependencies/Panvascular_Figure3_Genetic_Transcriptomic_20260920'
# Existing figures and scientific tables are inputs, not candidates for recomputation.
protected={str(p):sha(p) for p in PK.rglob('*') if p.is_file() and p.suffix.lower() in ['.docx','.xlsx']}
science={str(p):sha(p) for p in (O/'marker_audit_20260921/tables').glob('*.tsv')}
for n in ['MAGMA_display9.tsv','MAGMA_complete_GSEA.tsv','pseudobulk_display9_all7_lineages.tsv','pseudobulk_display9_selection.tsv','pseudobulk_input_manifest.tsv','pathway_display_labels.tsv']:science[str(E/'source_data'/n)]=sha(E/'source_data'/n)
read=lambda p:pd.read_csv(p,sep='\t',float_precision='round_trip')
m=read(E/'source_data/MAGMA_display9.tsv');full=read(E/'source_data/MAGMA_complete_GSEA.tsv');j=m.merge(full,on=['database','pathway'],suffixes=['_display','_source'],validate='one_to_one');assert len(j)==9
for c in ['NES','pval','padj','size']:assert np.allclose(j[c+'_display'],j[c+'_source'],rtol=0,atol=1e-12)
pb=read(E/'source_data/pseudobulk_display9_all7_lineages.tsv');manifest=read(E/'source_data/pseudobulk_input_manifest.tsv');assert len(pb)==63
for _,r in manifest.iterrows():
 assert sha(r.source)==r.sha256
 raw=read(r.source);a=pb[pb.analysis_id.eq(r.analysis_id)].merge(raw,on=['database','pathway'],suffixes=['_display','_source'],validate='one_to_one')
 assert len(a)==9
 for c in ['NES','pval','padj','size']:assert np.allclose(a[c+'_display'],a[c+'_source'],rtol=0,atol=1e-12,equal_nan=True)
# Copy the already selected display tables, not new pathway selections or test results.
(SD/'Figure3c_current').mkdir(exist_ok=True)
for n in ['MAGMA_display9.tsv','MAGMA_complete_GSEA.tsv','pseudobulk_display9_all7_lineages.tsv','pseudobulk_display9_selection.tsv','pseudobulk_input_manifest.tsv','pathway_display_labels.tsv']:shutil.copy2(E/'source_data'/n,SD/'Figure3c_current'/n)
sp=read(O/'notes/Figure5_Source_Data_index.tsv');sm_path=Path(sp.loc[sp.panel.eq('5a tissue identity'),'source_file'].iloc[0]);sm=read(sm_path);par=json.loads((O/'notes/Figure5_checks.json').read_text());records=[]
for row,key,short,n in [(1,'original','CA1',1779),(2,'added','T1077',1363)]:
 spec=par[key];v=read(sp.loc[sp.panel.eq('5a '+short),'source_file'].iloc[0]);assert len(v)==n and v.section_id.nunique()==1 and v.section_id.iloc[0]==spec['section'] and v.spot_id.nunique()==n
 assert v.in_tissue.eq(1).all() and v.library_size_original.ge(500).all()
 rr=sm[sm.section_id.eq(spec['section'])];assert len(rr)==1;rr=rr.iloc[0]
 records.append({'row':row,'study':spec['section'].split('::')[0],'section':spec['section'].split('::')[1],'full_section_id':spec['section'],'recorded_unit_id':v.biological_unit_id.iloc[0] if 'biological_unit_id' in v else rr.biological_unit_id,'n_existing_spots':n,'FN1_colour_min':spec['limits']['FN1'][0],'FN1_colour_max':spec['limits']['FN1'][1],'Fib_colour_min':spec['limits']['Fib_abundance'][0],'Fib_colour_max':spec['limits']['Fib_abundance'][1],'spot_source':sp.loc[sp.panel.eq('5a '+short),'source_file'].iloc[0],'spot_source_sha256':sha(sp.loc[sp.panel.eq('5a '+short),'source_file'].iloc[0]),'image':spec['image'],'image_sha256':sha(spec['image']),'identity_source':str(sm_path),'selection_record':str(O/'notes/Figure5_checks.json'),'selection_review':'Existing recorded choice verified; not blinded; no new selection','donor_claim':'CA1 paired P1 as recorded' if row==1 else 'T1077 library ID; not independently confirmed donor'})
section_table=pd.DataFrame(records);section_table.to_csv(NOUT/'Figure5a_current_sections.tsv',sep='\t',index=False);shutil.copy2(NOUT/'Figure5a_current_sections.tsv',SD/'Figure5a_current_sections.tsv')
text_cn='''# 本轮待采用图注与来源（未写回DOCX）

本文件是Figure1g、Figure3c和Figure5a的当前图注入口；历史“同七通路GSEA＋ORA”和单CA1图注不再用于当前图件。其历史结果保留，不删除科学数据。

## Figure1g — B显示已采用

**待采用英文：** Canonical marker expression across 14 displayed L1 classes: Endothelial cell, Fibroblast, Smooth muscle cell, Fibromyocyte, Pericyte, Macrophage, Monocyte, Neutrophil, Dendritic cell, T cell, Natural killer cell, B cell, Plasma cell and Mast cell. Unknown and Schwann cell are omitted from this display. Circle area is proportional to the fraction of cells with a positive raw count within each class; zero detection has zero area. Colour shows the previously audited mean log1p(CP10k) expression after per-gene min–max scaling across the 14 displayed classes, on a 0–1 scale. The unscaled means, detection fractions, stable gene identifiers and fixed row/column order are provided in Source Data; the unscaled A display is retained as companion material.

说明：保留全部28基因及所有跨类型单元。TCF21不是本数据中Fibromyocyte的特异marker；显示缩放不是新的注释验证。此前审计复现了存储矩阵统计，但独立逐sample/cell测量覆盖掩码仍未提供，不能据本次采用称其获得新增验证。

## Figure3c — MAGMA遗传通路＋原稿策略pseudobulk GSEA

**待采用英文：** Left, nine pathways from the existing preranked GSEA of MAGMA gene statistics (17,885 ranked genes; Supplementary Table 19). Within each of GO BP, KEGG LEGACY and Reactome, the three displayed pathways were selected from the existing adjusted-P < 0.05 results by increasing adjusted P, with pathway identifier used to resolve ties. Blue bar lengths indicate the original normalised enrichment score (NES); the adjacent numbers report the original adjusted P values. Right, existing preranked GSEA results from the manuscript-strategy pseudobulk primary DerSimonian–Laird meta-analysis, ranked by its signed Z statistics, for seven L1 populations (EC, Fib, SMC, Mac, Mo, T and NK; complete results associated with Supplementary Table 20). For each database, three pathways were previously selected by the number of populations with adjusted P < 0.05, followed by median absolute NES and pathway identifier. All seven population results are retained for each displayed pathway. Red and blue encode positive and negative NES, respectively, on the existing symmetric scale; dots denote original adjusted P < 0.05. Light versus dark dot colours are used only for contrast. The two blocks summarise distinct existing genetic and transcriptomic enrichment analyses; they are not paired GSEA/ORA results or an independent validation experiment.

来源与口径：左侧Table 19中Analysis=preranked_GSEA的记录，当前9条来源及数值已与MAGMA_complete_GSEA逐行核对。右侧9×7单元已与7个原生fgsea_all.tsv.gz核对，分析ID及SHA256在pseudobulk_input_manifest.tsv。两侧均各含GO BP/KEGG LEGACY/Reactome各3条，但不是同一套通路。右侧颜色范围为原±{vmax:.14g}（图面四舍五入为±3.39）。未重新排序选材、未重跑富集、未改P/FDR。ORA原统计有效性待办继续PENDING；移出当前3c不表示该待办已解决。

## Figure5a — 两行空间三视图

**待采用英文：** Original section views from GSE248608, GSM7917894_CA1 (top; anterior tibial artery, 1,779 included spots), and GSE314851, GSM9414346_T1077 (bottom; coronary artery, 1,363 included spots). In each row, the panels show the original H&E image, inferred Fib abundance and observed FN1 log1p(CP10k) expression using the same field, orientation and included spot coordinates. Inferred Fib abundance is the sum of the marginal q05 abundances of Fibroblast states, not the q05 of their summed abundance distribution. Both rows retain the original in-tissue and total-UMI ≥ 500 inclusion rule. FN1 uses the previously stored section-specific positive-expression colour limits (2nd–98th percentiles); zero is shown separately in grey and missing values are treated separately. The Fib-abundance colour limits are the stored section-specific minimum and maximum. Colour intensity should not be compared directly across the separately scaled sections. These displays provide descriptive spatial context, not evidence of cell-intrinsic downregulation, pathology-confirmed layer enrichment or independent histological validation. The sections were chosen during results presentation after review of existing results and H&E. T1077 is a recorded library identifier, not a confirmed independent donor.

来源已同步：两行身份均可追溯至原sample_manifest.tsv、Figure5_checks.json选材记录、原spot表、原H&E及scalefactors。上行配对单位为原记录GSE248608::P1，下行为GSE314851::T1077；不合并改称两名已确认供者。T1077此前选材记录为冠脉组织补充、H&E主体完整、无QC fallback及原FN1条件系数方向与总体一致；本轮不重新挑图、未将总体FDR标成单切片显著性。逐行精确色限与来源hash见Figure5a_current_sections.tsv。图面已补GSE/section短标识，完整GSM标识在本图注与索引。旧SCARB1 CA1描述图及不显著说明仍保留在配套材料。
'''.format(vmax=pb.NES.abs().max())
(NOUT/'CURRENT_PENDING_CAPTIONS.md').write_text(text_cn)
# Authoritative mapping to current content; originals and A remain separately traceable.
index=[{'panel':'Figure1g B','display':'14 x 28, approved per-gene 0-1 colour; original area fractions','source':str(SD/'marker_audit_20260921/per_gene_scaled_mean.tsv'),'supporting_source':str(SD/'marker_audit_20260921/original_mean_log1p.tsv'),'caption':'CURRENT_PENDING_CAPTIONS.md / Figure1g'},
 {'panel':'Figure1g A companion','display':'Unscaled mean; not formal g; no new ED','source':str(PK/'图件编辑资料/marker_audit_20260921/figures/Figure1g_A_original.pdf'),'supporting_source':str(SD/'marker_audit_20260921/original_detection_fraction.tsv'),'caption':'CURRENT_PENDING_CAPTIONS.md / Figure1g'},
 {'panel':'Figure3c left','display':'9 existing MAGMA preranked GSEA pathways; not ORA','source':str(SD/'Figure3c_current/MAGMA_display9.tsv'),'supporting_source':str(PK/'Supplementary_Tables.xlsx')+' :: Table 19 / preranked_GSEA','caption':'CURRENT_PENDING_CAPTIONS.md / Figure3c'},
 {'panel':'Figure3c right','display':'9 pathways x 7 L1 original-strategy pseudobulk GSEA','source':str(SD/'Figure3c_current/pseudobulk_display9_all7_lineages.tsv'),'supporting_source':str(SD/'Figure3c_current/pseudobulk_input_manifest.tsv'),'caption':'CURRENT_PENDING_CAPTIONS.md / Figure3c'},
 {'panel':'Figure5a rows 1-2','display':'CA1 1779 spots; T1077 1363 spots; descriptive context','source':str(SD/'Figure5a_current_sections.tsv'),'supporting_source':str(sm_path),'caption':'CURRENT_PENDING_CAPTIONS.md / Figure5a'}]
pd.DataFrame(index).to_csv(NOUT/'CURRENT_PANEL_SOURCE_INDEX.tsv',sep='\t',index=False)
for row in index:assert Path(row['source']).exists(),row['source']
# Restricted merge: update only Figure1 and Figure5; verify all 13 other native files unchanged.
v=read(O/'FIGURE_VERSIONS.tsv').sort_values('page');pv=read(PK/'图件采用版本.tsv').set_index('figure');changed=['Figure1','Figure5'];checks=[]
for i,row in v.iterrows():
 pdf=Path(row.source_pdf);svg=Path(row.source_svg);doc=fitz.open(pdf);p=doc[0]
 if row.figure in changed:
  bad=[s['text'] for b in p.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans'] if s['bbox'][0]<-.2 or s['bbox'][1]<-.2 or s['bbox'][2]>p.rect.width+.2 or s['bbox'][3]>p.rect.height+.2];assert not bad,(row.figure,bad)
  v.loc[i,'sha256']=sha(pdf);v.loc[i,'svg_sha256']=sha(svg);v.loc[i,'version']='Current in-place: Figure1g B adopted / Figure5 section IDs; 20260921'
  v.loc[i,'status']='FIGURE1G_B_AUTHOR_APPROVED' if row.figure=='Figure1' else 'CURRENT_DISPLAY_IDENTIFIERS_UPDATED'
  shutil.copy2(pdf,Path(pv.loc[row.figure,'package_pdf']));shutil.copy2(svg,Path(pv.loc[row.figure,'package_svg']))
 else:assert sha(pdf)==row.sha256 and sha(svg)==row.svg_sha256
 checks.append({'figure':row.figure,'width_mm':p.rect.width/mm,'height_mm':p.rect.height/mm,'changed_this_round':row.figure in changed})
v.to_csv(O/'FIGURE_VERSIONS.tsv',sep='\t',index=False)
merged=fitz.open()
for _,row in v.iterrows():
 d=fitz.open(row.source_pdf);merged.insert_pdf(d);typ,val=d.xref_get_key(d[0].xref,'Group')
 if typ=='dict':assert not re.search(r'\d+ \d+ R',val);merged.xref_set_key(merged[-1].xref,'Group',val)
book=O/'全部图件_15页_主图调整版.pdf';tmp=book.with_suffix('.tmp.pdf');merged.save(tmp,garbage=4,deflate=True);merged.close();tmp.replace(book);shutil.copy2(book,PK/'全部图件_15页审阅.pdf')
mmrg=fitz.open(book);assert len(mmrg)==15
for i,(_,row) in enumerate(v.iterrows()):
 p=fitz.open(row.source_pdf)[0];assert p.rect==mmrg[i].rect and p.get_pixmap().samples==mmrg[i].get_pixmap().samples
 for c in v.columns:
  if c!='figure':pv.loc[row.figure,c]=row[c]
 pv.loc[row.figure,'package_pdf_sha256']=sha(pv.loc[row.figure,'package_pdf']);pv.loc[row.figure,'package_svg_sha256']=sha(pv.loc[row.figure,'package_svg'])
pv.reset_index().sort_values('page').to_csv(PK/'图件采用版本.tsv',sep='\t',index=False)
for name in ['CURRENT_PENDING_CAPTIONS.md','CURRENT_PANEL_SOURCE_INDEX.tsv','Figure5a_current_sections.tsv','Figure1g_B_adoption_checks.json']:
 shutil.copy2(NOUT/name,PK/'图件编辑资料/本轮说明_20260920'/name)
shutil.copy2(NOUT/'CURRENT_PANEL_SOURCE_INDEX.tsv',PK/'图件编辑资料/Source_Data_INDEX_current.tsv');shutil.copy2(NOUT/'CURRENT_PANEL_SOURCE_INDEX.tsv',SD/'CURRENT_PANEL_SOURCE_INDEX.tsv')
for f in ['adopt_Figure1g_B.py','finalize_B_adoption.py']:shutil.copy2(O/'code'/f,PK/'图件编辑资料/code/current_20260920'/f)
assert protected=={str(p):sha(p) for p in PK.rglob('*') if p.is_file() and p.suffix.lower() in ['.docx','.xlsx']}
assert all(sha(p)==h for p,h in science.items())
check={'Figure1g_B_adopted':True,'stored_statistics_reused_without_recalculation':True,'annotation_newly_validated':False,'measurement_coverage_newly_validated':False,'Figure3c_captions_and_current_sources_synced':True,'Figure5a_both_rows_sources_and_captions_synced':True,'Figure5a_blocked_rows':[],'ORA':'PENDING','donor_limitations':'UNCHANGED','unchanged_figure_count':13,'all_10_ED_unchanged':True,'15_native_pages_render_equal':True,'docx_xlsx_unchanged':True,'statistics_tables_unchanged':True,'visual_review':'Figure1g and Figure5 row labels visually inspected, no clipped labels; author-approved B display adopted','checks':checks}
(NOUT/'B_ADOPTION_DELIVERY_CHECKS.json').write_text(json.dumps(check,indent=2));shutil.copy2(NOUT/'B_ADOPTION_DELIVERY_CHECKS.json',PK/'图件编辑资料/VERSION_EQUALITY.json')
summary='''# 本轮完成说明

- **B显示已采用**：正式Figure1g使用原14类×28基因已生成的逐基因0–1颜色表，点面积仍为原检出率；按整图目标位置原生重导，未拉伸独立A/B页面。A版保留在原配套材料及Source Data，不新增ED。
- **计算沿用已核查存储统计**：本轮不重算表达、检出率或缩放颜色值。原均值、检出率、稳定ID及固定顺序不变。
- **未新增验证**：TCF21仍保留，不称其为Fibromyocyte特异marker；注释边界、逐sample/cell测量覆盖与ORA待办未因采用B而解决。
- **Figure3c说明已同步**：当前为9条MAGMA排序GSEA及9×7的pseudobulk GSEA，原数值逐项对照通过；当前图注/索引不再使用“同七通路GSEA＋ORA”。
- **Figure5a说明已同步**：已补GSE248608/CA1（1779 spots）及GSE314851/T1077（1363 spots）。两行来源与选材记录可追溯，无受阻行。T1077仅为文库ID，不称独立供者；各自色限和原证据边界保留。

仅覆盖现有Figure1、Figure5 PDF/SVG及15页合并稿；其余13幅（含10幅ED）原文件不变，正文未写回。待采用文字与逐panel来源以`图件编辑资料/本轮说明_20260920/CURRENT_PENDING_CAPTIONS.md`和`CURRENT_PANEL_SOURCE_INDEX.tsv`为准。轻量备份仅在现有图件编辑资料内；没有新增桌面版本或整套投稿包。
'''
(O/'本轮修改说明.md').write_text(summary);(PK/'图件编辑资料/本轮修改说明.md').write_text(summary)
# Update entrypoints; historical audit records remain explicitly historical.
for p in [PK/'README_当前状态.txt',PK/'图件编辑资料/README_CN.txt',O/'README_唯一最新版.md']:
 s=p.read_text();s=s.replace('正式Figure1g仍保留原均值颜色版。','正式Figure1g已采用B版逐基因min–max颜色；原均值A版仍保留配套及Source Data。').replace('不自动替换正式1g。','本轮已按作者授权采用B版，原A仍保留配套。').replace('marker_audit_20260921是正式Figure1g之外的A/B与计算审查，不自动覆盖原值；','marker_audit_20260921保留A/B与既有计算审查；当前正式Figure1g已按作者授权采用B显示，');s+='\n当前Figure1g、3c和5a待采用图注：图件编辑资料/本轮说明_20260920/CURRENT_PENDING_CAPTIONS.md；正文未改。\n';p.write_text(s)
# Existing checksum inventory refreshed after content-bound synchronization.
idxpath=PK/'文件来源与校验.tsv';rows=read(idxpath).to_dict('records');known={x['file'] for x in rows}
for row in rows:
 p=PK/row['file']
 if p.is_file():row['sha256']=sha(p)
for p in PK.rglob('*'):
 if p.is_file() and p!=idxpath and str(p.relative_to(PK)) not in known:rows.append({'file':str(p.relative_to(PK)),'kind':'B adoption or current caption/source crosswalk','source':str(O),'sha256':sha(p)})
pd.DataFrame(rows).to_csv(idxpath,sep='\t',index=False)
print('PASS: B adopted; 2 PDFs/SVGs replaced; 13 files incl all ED unchanged; 15 native pages match; captions/source indices updated; DOCX unchanged')
