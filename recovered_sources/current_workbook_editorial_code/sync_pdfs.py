from pathlib import Path
import fitz,shutil,json
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');O=P/'图件编辑资料/补充表与方法同步_20260921';N=P/'Nature_Submission'
d=fitz.open(N/'Manuscript.pdf');m=len(d);f=fitz.open(P/'全部图件_15页审阅.pdf');assert len(f)==15
d.insert_pdf(f);out=N/'Manuscript_with_15_Figures.pdf';b=O/'轻量备份/Nature_Submission'/out.name
if not b.exists():shutil.copy2(out,b)
d.save(str(out)+'.tmp');d.close();Path(str(out)+'.tmp').replace(out)
z=fitz.open(out)
for i,page in enumerate(f):
 assert z[m+i].rect==page.rect
 assert z[m+i].get_pixmap().samples==page.get_pixmap().samples
(O/'pdf_validation.json').write_text(json.dumps({'manuscript_pages':m,'figure_pages':15,'combined_pages':len(z),'original_figure_page_dimensions_and_render_identical':True,'revised_docx_render_visual_review':['QC Methods','fine-mapping Methods','ED9 caption','Supplementary table index','supplementary fine-mapping details']},indent=2))
print('PDF merge verified',len(z))
