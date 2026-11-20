from pathlib import Path
import csv,json,hashlib
import fitz
O=Path(__file__).resolve().parents[1];B=Path('/Users/taki/Desktop/Panvascular_Adopted_Figure5_Merged_20260920')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
rows=list(csv.DictReader((B/'ADOPTED_FIGURE_VERSIONS.tsv').open(),delimiter='\t'));out=fitz.open();sources=[]
for i,r in enumerate(rows):
 p=O/'figures/Figure1.pdf' if i==0 else Path(r['source_pdf']);s=fitz.open(p);out.insert_pdf(s);kind,val=s.xref_get_key(s[0].xref,'Group')
 if kind=='dict':out.xref_set_key(out[-1].xref,'Group',val)
 sources.append(p)
 if i==0:r.update(source_pdf=str(p),source_svg=str(O/'figures/Figure1.svg'),delivery_pdf=str(p),delivery_svg=str(O/'figures/Figure1.svg'),adopted_version=O.name,adoption_scope='LABEL_COMPLETION_AUTHOR_REVIEW_PENDING',pdf_sha256=sha(p),svg_sha256=sha(O/'figures/Figure1.svg'))
 else:assert sha(p)==r['pdf_sha256']
out.set_toc([[1,r['figure'],i+1] for i,r in enumerate(rows)]);out.save(O/'全部图件_15页_Figure1b标注补齐.pdf');out.close();merged=fitz.open(O/'全部图件_15页_Figure1b标注补齐.pdf')
for i,p in enumerate(sources):
 s=fitz.open(p);assert s[0].rect==merged[i].rect;assert s[0].get_pixmap(matrix=fitz.Matrix(1.5,1.5)).samples==merged[i].get_pixmap(matrix=fitz.Matrix(1.5,1.5)).samples
with (O/'FIGURE_VERSIONS.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
old=fitz.open('/Users/taki/Desktop/Panvascular_Marker_Spatial_Review_20260919/main/Figure1.pdf');new=fitz.open(O/'figures/Figure1.pdf');checks=[]
for mmrect in [(0,0,183,37),(0,86,183,170),(94,37,183,86)]:
 rect=fitz.Rect(*[v*72/25.4 for v in mmrect]);same=old[0].get_pixmap(matrix=fitz.Matrix(2,2),clip=rect).samples==new[0].get_pixmap(matrix=fitz.Matrix(2,2),clip=rect).samples;assert same;checks.append({'unchanged_region_mm':mmrect,'pixels_identical':same})
x=json.loads((O/'notes/checks.json').read_text());x.update({'visual_review':'ACTUALLY_VIEWED_Figure1_after.png; 16 numbers and full names readable; author review pending','unchanged_other_regions':checks,'merged_pages':15,'all_merged_pages_native_identical':True,'Figure5_adopted_version_preserved':rows[4]['source_pdf']});(O/'notes/checks.json').write_text(json.dumps(x,indent=2))
(O/'修改说明.md').write_text('''# Figure1b标注补齐

原图图例包含全部16类，但UMAP图内仅标注Macrophage、Fibroblast、Smooth muscle、Endothelial、T和Natural killer共6类，缺少另外10类的图内对应。

本轮补齐16类编号和引线；图例给出编号对应的完整细胞类型名称。编号仅用于读图，不是重新分群编号。锚点选自该L1类别内的真实细胞，代表定位，不暗示该类别所有细胞都集中在这一处；完整cell ID及冻结坐标见notes/Figure1b_16_label_mapping.tsv。

全部细胞、冻结坐标、颜色、图中嵌入的UMAP像素均未改变。其他Figure1面板及右侧study图保持，采用版Figure5和其余图未修改。新15页合并稿逐页保留原生尺寸，渲染与各自单图一致。

图注可追加（未写回DOCX）：Numbers identify the 16 major cell categories listed in the adjacent legend; leader lines indicate representative positions on the unchanged global embedding.

仅补图形标注，未进行新注释或科学分析。等待作者审核。
''')
print('Figure1b 16 labels complete; 15-page merge checked; all other figure regions unchanged.')
