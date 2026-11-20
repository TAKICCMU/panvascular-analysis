from pathlib import Path
import csv,json,hashlib,shutil
import fitz
O=Path(__file__).resolve().parents[1]
B=Path('/Users/taki/Desktop/Panvascular_Marker_Spatial_Review_20260919')
F=Path('/Users/taki/Desktop/Panvascular_Figure5_AB_Finalization_20260920')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
rows=list(csv.DictReader((B/'图件来源与尺寸.tsv').open(),delimiter='\t'))
assert len(rows)==15
for name in ['Figure5.pdf','Figure5.svg']:shutil.copyfile(F/'figures'/name,O/'Figure5'/name)
for folder in ['Source_Data','companion']:shutil.copytree(F/folder,O/folder,dirs_exist_ok=True)
shutil.copyfile(F/'Source_Data_index.tsv',O/'Source_Data_index.tsv')
shutil.copyfile(F/'图注与引用建议_待采用.md',O/'待统一文稿合并/Figure5ab图注及P86引用_待采用.md')
shutil.copyfile(F/'修改与计数说明.md',O/'Figure5/修改与计数说明_采用来源.md')
# Preserve the recorded baseline lineage rather than infer versions from mtimes.
shutil.copyfile(B/'图件来源与尺寸.tsv',O/'原14幅图版本登记来源.tsv')
merged=fitz.open();adopted=[];sources=[]
for i,r in enumerate(rows):
 name=r['figure'];is5=name=='Figure5'
 srcpdf=F/'figures/Figure5.pdf' if is5 else Path(r['pdf'])
 srcsvg=F/'figures/Figure5.svg' if is5 else Path(r['svg'])
 if not is5:
  assert sha(srcpdf)==r['pdf_sha256'] and sha(srcsvg)==r['svg_sha256']
 else:
  assert sha(srcpdf)==sha(O/'Figure5/Figure5.pdf') and sha(srcsvg)==sha(O/'Figure5/Figure5.svg')
 doc=fitz.open(srcpdf);assert len(doc)==1
 pg=doc[0];merged.insert_pdf(doc,from_page=0,to_page=0)
 # insert_pdf omits the page-level transparency group; preserve the original direct dictionary.
 group_type,group_value=doc.xref_get_key(pg.xref,'Group')
 if group_type=='dict':
  import re
  assert not re.search(r'\d+\s+\d+\s+R',group_value), 'Resolve referenced Group objects before copying'
  merged.xref_set_key(merged[-1].xref,'Group',group_value)
 else:assert group_type=='null',group_type
 adopted.append({'page':i+1,'figure':name,'adopted_version':F.name if is5 else B.name,'adoption_scope':'FIGURE5_ADOPTED_THIS_ROUND' if is5 else 'INHERITED_RECORDED_VERSION_NO_NEW_REVIEW','source_pdf':str(srcpdf),'source_svg':str(srcsvg),'delivery_pdf':str(O/'Figure5/Figure5.pdf') if is5 else str(srcpdf),'delivery_svg':str(O/'Figure5/Figure5.svg') if is5 else str(srcsvg),'width_mm':round(pg.rect.width*25.4/72,4),'height_mm':round(pg.rect.height*25.4/72,4),'pdf_sha256':sha(srcpdf),'svg_sha256':sha(srcsvg),'baseline_manifest':str(B/'图件来源与尺寸.tsv') if not is5 else str(F/'notes/checks.json')})
 sources.append(srcpdf)
merged.set_toc([[1,r['figure'],r['page']] for r in adopted]);dest=O/'全部图件_15页_Figure5采用版.pdf';merged.save(dest,garbage=0,deflate=False);merged.close()
combined=fitz.open(dest);assert len(combined)==15;checks=[]
for i,src in enumerate(sources):
 doc=fitz.open(src);a=doc[0];b=combined[i]
 geom=all(tuple(getattr(a,k))==tuple(getattr(b,k)) for k in ['rect','mediabox','cropbox']) and a.rotation==b.rotation
 text=a.get_text()==b.get_text()
 ap=a.get_pixmap(matrix=fitz.Matrix(1.5,1.5),alpha=False);bp=b.get_pixmap(matrix=fitz.Matrix(1.5,1.5),alpha=False)
 render=(ap.width,ap.height,ap.samples)==(bp.width,bp.height,bp.samples)
 assert geom and text and render,(i,geom,text,render)
 checks.append({'page':i+1,'figure':adopted[i]['figure'],'geometry_rotation_identical':geom,'text_identical':text,'render_pixels_identical':render})
with (O/'ADOPTED_FIGURE_VERSIONS.tsv').open('w') as out:
 w=csv.DictWriter(out,fieldnames=list(adopted[0]),delimiter='\t');w.writeheader();w.writerows(adopted)
md=['# 当前15页采用版本清单','',f'Figure5采用 `{F.name}`；其他14幅图沿用 `{B.name}/图件来源与尺寸.tsv` 中已登记文件，未重新绘制、未在本轮重新批准。','', '|页|图|本次来源版本|原生尺寸（mm）|','|---|---|---|---|']
for r in adopted:md.append(f"|{r['page']}|{r['figure']}|{'Figure5_AB_Finalization_20260920（本轮采用）' if r['figure']=='Figure5' else 'Marker_Spatial_Review_20260919（沿用登记）'}|{r['width_mm']:g} × {r['height_mm']:g}|")
md+=['','逐文件路径及SHA256见ADOPTED_FIGURE_VERSIONS.tsv。合并采用原生PDF页面插入，15页的尺寸、CropBox、旋转、文本及渲染像素与各自母版一致；没有缩放、裁切或重新绘图。','','Figure5口径保持：60个gene–unit estimates，34个不同recorded unit IDs；不是60个独立单位，也不是34名已确认供者。原检出率矩阵在Source_Data/Figure5，完整两页配套图在companion。CA1的SCARB1描述性图及单切片P=0.582、P1合并方向为负的说明均保留。','','Figure5a/b图注及P86拆分引用建议位于“待统一文稿合并”，尚未写回DOCX。ORA与供者身份问题仍PENDING；本轮采用不代表这些问题解决，也不代表其他14幅图获得新的最终审核。']
(O/'逐图版本清单.md').write_text('\n'.join(md))
# Verify every preserved Figure5 package asset, not merely that copies exist.
copies=[]
for folder in ['Source_Data','companion']:
 for p in (F/folder).rglob('*'):
  if p.is_file():
   q=O/folder/p.relative_to(F/folder);assert sha(p)==sha(q);copies.append(str(q.relative_to(O)))
assert sha(F/'Source_Data_index.tsv')==sha(O/'Source_Data_index.tsv')
for r in adopted:assert sha(r['source_pdf'])==r['pdf_sha256'] and sha(r['source_svg'])==r['svg_sha256']
(O/'MERGE_CHECK.json').write_text(json.dumps({'merged_file':str(dest),'merged_sha256':sha(dest),'pages':checks,'Figure5_adopted':True,'other14_new_author_approval':False,'native_page_insert_only':True,'science_rerun':False,'DOCX_modified':False,'ORA':'PENDING','donor_identity':'PENDING','counts':{'gene_unit_estimates':60,'distinct_recorded_unit_IDs':34},'preserved_assets_byte_identical':copies,'status':'MERGE_COMPLETE_FIGURE5_ADOPTED_NO_FURTHER_REDESIGN'},ensure_ascii=False,indent=2))
print('Merged 15 pages: original page geometry, text and render pixels match all 15 source PDFs.')
print(dest)
