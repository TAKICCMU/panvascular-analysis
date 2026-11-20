from pathlib import Path
import copy,xml.etree.ElementTree as E,subprocess,os,json
import fitz
B=Path('/Users/taki/Desktop/Panvascular_Content_Alignment_ST_Audit_20260919');O=Path('/Users/taki/Desktop/Panvascular_Marker_Spatial_Review_20260919');ns='{http://www.w3.org/2000/svg}'
p=B/'editable_svg/ExtendedData_Figure10.svg';tree=E.parse(p);r=tree.getroot();targets=[e for e in r.iter(ns+'text') if e.text=='Fraction'];assert len(targets)==1;t=targets[0];t.text='Pair fraction';parent=next(e for e in r.iter() if t in list(e));new=copy.deepcopy(t);new.text='(ρ > 0; FDR < 0.05)';new.set('y',str(float(t.get('y'))+7.5));new.set('transform',f'rotate(0 {new.get("x")} {new.get("y")})');parent.append(new)
out=O/'editable_svg/ExtendedData_Figure10.svg';tree.write(out,encoding='utf-8',xml_declaration=True)
env=os.environ.copy();env['DYLD_FALLBACK_LIBRARY_PATH']='/Users/taki/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/lib'
subprocess.run(['/Users/taki/Desktop/Codex任务1/.venv_figure/bin/python','/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Local_Touchup_20260918/code/export_svg_pdf.py',str(out),str(O/'review/ED10_label_render.pdf')],env=env,check=True)
old=fitz.open(B/'extended_data/ExtendedData_Figure10.pdf');patch=fitz.open(O/'review/ED10_label_render.pdf');doc=fitz.open(B/'extended_data/ExtendedData_Figure10.pdf');page=doc[0];box=fitz.Rect(115,208.6,190,225);w,h=page.rect.width,page.rect.height
page.add_redact_annot(box,fill=(1,1,1));page.apply_redactions(images=0,graphics=0)
page.show_pdf_page(box,patch,0,clip=box);doc.save(O/'extended_data/ExtendedData_Figure10.pdf',garbage=4,deflate=True);page.get_pixmap(matrix=fitz.Matrix(2,2)).save(O/'review/ED10_axis_label.png')
# All other content is source PDF reused at original coordinates; verify shared c/e panel region exactly.
q=fitz.Rect(69*72/25.4,0,w,h);assert page.get_pixmap(matrix=fitz.Matrix(2,2),clip=q).samples==old[0].get_pixmap(matrix=fitz.Matrix(2,2),clip=q).samples
(O/'notes/ED10_label_change.json').write_text(json.dumps({'change':'Fraction -> Pair fraction / (rho > 0; FDR < 0.05)','source_field':'positive_FDR_pair_fraction','denominator':'meta-analyzed pairs of the five highest-weight nodes in each module; n_top5_pairs_meta_analyzed (10 in plotted modules)','numerator':'pooled_rho > 0 AND BH FDR < 0.05','other_panel_scope':'PDF outside xlabel patch reused exactly; c/e shared row region pixel-identical','plot_values_unchanged':True},indent=2))
