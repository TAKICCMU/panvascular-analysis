"""Install display-only Figure4a refinement; preserve all protected content."""
from pathlib import Path
import shutil,json,hashlib
import pymupdf as fitz
import numpy as np
P=Path('/Users/taki/Desktop/泛血管稿件_批注修订_20261116')
W=Path(__file__).resolve().parents[1]
E=P/'图件编辑资料/本轮定稿核对_20261117'
D=E/'Figure4a'
B=D/'backup_before_model_simplification'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest=E/'DELIVERY_SHA256.json'; orig=json.loads(manifest.read_text())
protected={p:sha(P/p) for p in orig if not p.startswith('Figures/Figure4') and p!='Figures_and_ExtendedData_审阅.pdf'}
old=fitz.open(B/'Figure4.pdf');new=fitz.open(W/'review/Figure4.pdf')
assert old[0].rect==new[0].rect
clip=fitz.Rect(0,25*72/25.4,old[0].rect.width,old[0].rect.height)
r0=old[0].get_pixmap(matrix=fitz.Matrix(2,2),clip=clip,alpha=False)
r1=new[0].get_pixmap(matrix=fitz.Matrix(2,2),clip=clip,alpha=False)
assert r0.samples==r1.samples, 'Protected b-d pixels changed'
for letter in 'bcd':
 for ext in ['pdf','svg','png']:
  assert (W/'panels'/f'Figure4{letter}.{ext}').read_bytes()==(P/'图件编辑资料/Figure4/panels'/f'Figure4{letter}.{ext}').read_bytes()
for sub in ['panels','code','review']:
 (D/sub).mkdir(parents=True,exist_ok=True)
 for f in (W/sub).iterdir():
  if f.is_file():shutil.copy2(f,D/sub/f.name)
for ext in ['pdf','svg']:shutil.copy2(W/'review'/f'Figure4.{ext}',P/'Figures'/f'Figure4.{ext}')
# Replace only page 4. Other 13 pages are copied from the actual current combined PDF.
m=fitz.open(P/'Figures_and_ExtendedData_审阅.pdf');out=fitz.open()
assert len(m)==14
for i in range(14):out.insert_pdf(new if i==3 else m,from_page=0 if i==3 else i,to_page=0 if i==3 else i)
tmp=P/'Figures_and_ExtendedData_审阅.tmp.pdf';out.save(tmp,garbage=4,deflate=True);out.close()
t=fitz.open(tmp)
for i in range(14):
 source=new[0] if i==3 else m[i]
 assert source.rect==t[i].rect
 assert source.get_text('text')==t[i].get_text('text')
 assert source.get_pixmap(alpha=False).samples==t[i].get_pixmap(alpha=False).samples
m.close();t.close();tmp.replace(P/'Figures_and_ExtendedData_审阅.pdf')
for p,h in protected.items():assert sha(P/p)==h
for p in ['Figures/Figure4.pdf','Figures/Figure4.svg','Figures_and_ExtendedData_审阅.pdf']:orig[p]=sha(P/p)
manifest.write_text(json.dumps(orig,ensure_ascii=False,indent=2)+'\n')
new[0].get_pixmap(matrix=fitz.Matrix(96/72,96/72),alpha=False).save(D/'review/Figure4_actual_size_96dpi.png')
report={'panel_a_mm':[183,25],'full_figure_mm':[183,170],'minimum_panel_font_pt':6.0,'bcd_pixels_unchanged':True,'bcd_source_files_byte_identical':True,'other_13_merged_pages_render_identical':True,'documents_and_tables_unchanged':True,'schematic_only_no_science_rerun':True,'biorender_status':'PENDING_CURRENT_CANVAS_SYNC','new_figure4_pdf_sha256':sha(P/'Figures/Figure4.pdf')}
(D/'review/detail_refinement_checks.json').write_text(json.dumps(report,indent=2)+'\n')
(D/'README.md').write_text('''# Current Figure 4a\n\n183 × 25 mm vector schematic; assembled Figure 4 remains 183 × 170 mm.\n\nThe latest refinement clarifies the cell-proportion matrix, NMF factors and module groups; the spatial branch now displays candidate expression and inferred cell abundance as an adjusted association within each section, followed by pooling. Model terms and the depth offset remain in Methods; no scientific model was changed. Matrices, tissue outlines, nodes and aggregation glyphs are schematic and do not encode fabricated measurements or inferred edges. Numerical counts are unchanged. No scientific analysis was rerun.\n\nPanels b–d retain the same files, position, dimensions and rendered pixels. Panel letter a is added only in the assembly / as a separate BioRender text object.\n\nCurrent source: code/redraw_figure4a.py.\nCurrent standalone artwork: panels/Figure4a.svg and .pdf.\nFull preview: review/Figure4.png.\nBioRender canvas: https://app.biorender.com/illustrations/6af982741be2d4f7c66d0f97\n\nThe previous affected artwork and code are retained under the affected-file backup subdirectories.\n''')
print(json.dumps(report,indent=2))
