from pathlib import Path
import json,csv,hashlib,shutil,copy,io
import fitz,numpy as np,xml.etree.ElementTree as E
from PIL import Image,ImageDraw
O=Path(__file__).resolve().parents[1];B=Path('/Users/taki/Desktop/Panvascular_Content_Alignment_ST_Audit_20260919');sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();N={'s':'http://www.w3.org/2000/svg'}
seq=[('Figure'+str(i),'main') for i in range(1,6)]+[('ExtendedData_Figure'+str(i),'extended_data') for i in range(1,11)];changed={'Figure3','Figure5','ExtendedData_Figure2','ExtendedData_Figure3','ExtendedData_Figure4','ExtendedData_Figure5','ExtendedData_Figure10'}
for name,folder in seq:
 if name not in changed:
  shutil.copy2(B/folder/(name+'.pdf'),O/folder/(name+'.pdf'));shutil.copy2(B/'editable_svg'/(name+'.svg'),O/'editable_svg'/(name+'.svg'))
# Verify baseline unchanged against its pre-existing frozen output manifest.
for row in csv.DictReader((B/'图件版本与尺寸.tsv').open(),delimiter='\t'):
 for kind in ['pdf','svg']:assert sha(row[kind])==row[kind+'_sha256'],row[kind]
# Figure5b-f and panel letters are literally retained SVG subtrees.
b=E.parse(B/'editable_svg/Figure5.svg').getroot();n=E.parse(O/'editable_svg/Figure5.svg').getroot();fgb=b.find("s:g[@id='figure_1']",N);fgn=n.find("s:g[@id='figure_1']",N)
removed={'axes_1','axes_2','axes_3','axes_4','text_133','text_134','text_135','legend_1'}
for el in fgb:
 if el.get('id') not in removed:
  after=fgn.find("s:g[@id='%s']"%el.get('id'),N)
  assert after is not None and E.tostring(el)==E.tostring(after),el.get('id')
# All existing spatial observations, limits and flags unchanged.
old=json.load(open(B/'audit/Figure5_display_and_alignment.json'))['spatial_display'];new=json.load(open(O/'audit/Figure5a_display.json'))['data_regions']
for a,z in zip(old,new):
 for key in ['gene','n_spots','n_positive','n_zero','vmin','vmax','xlim','ylim','flagged_barcodes']:assert a[key]==z[key],key
 assert z['same_H_and_E_spot_scale']
assert sum(len(z['flagged_barcodes']) for z in new)==7
# Source data hashes checked by marker layer audit; no statistical recalculation.
merged=fitz.open();rows=[]
for i,(name,folder) in enumerate(seq,1):
 pdf=O/folder/(name+'.pdf');svg=O/'editable_svg'/(name+'.svg');d=fitz.open(pdf);assert len(d)==1;pg=d[0]
 target_width=183 if folder=='main' else 180;assert abs(pg.rect.width/72*25.4-target_width)<.01 and abs(pg.rect.height/72*25.4-170)<.01
 spans=[s for b in pg.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans'] if s['text'].strip()]
 outside=[s['text'] for s in spans if not fitz.Rect(-.2,-.2,pg.rect.width+.2,pg.rect.height+.2).contains(fitz.Rect(s['bbox']))];assert not outside,(name,outside)
 if name not in changed:assert sha(pdf)==sha(B/folder/(name+'.pdf')) and sha(svg)==sha(B/'editable_svg'/(name+'.svg'))
 merged.insert_pdf(d);rows.append({'page':i,'figure':name,'action':'LIMITED_MODIFICATION' if name in changed else 'BASELINE_BYTE_IDENTICAL_COPY','pdf':str(pdf),'svg':str(svg),'width_mm':target_width,'height_mm':170,'baseline_pdf_sha256':sha(B/folder/(name+'.pdf')),'pdf_sha256':sha(pdf),'svg_sha256':sha(svg),'minimum_extracted_font_pt':round(min(s['size'] for s in spans),2)})
mp=O/'全部图件_Marker简化与空间审阅_15页.pdf';merged.save(mp,garbage=4,deflate=True);merged.close();merged=fitz.open(mp);assert len(merged)==15
mergechecks=[]
for i,row in enumerate(rows):
 d=fitz.open(row['pdf']);assert d[0].rect==merged[i].rect
 a=np.frombuffer(d[0].get_pixmap().samples,np.uint8).astype(int);z=np.frombuffer(merged[i].get_pixmap().samples,np.uint8).astype(int);mx=int(np.abs(a-z).max());assert mx<=1,(row['figure'],mx)
 mergechecks.append({'figure':row['figure'],'geometry_exact':True,'max_8bit_render_channel_difference':mx})
with (O/'图件来源与尺寸.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(rows)
# True same-size marker comparisons using page insertion without rescaling.
comp=fitz.open()
for i in range(2,6):
 name=f'ExtendedData_Figure{i}';a=fitz.open(B/'extended_data'/(name+'.pdf'));z=fitz.open(O/'extended_data'/(name+'.pdf'));w=a[0].rect.width;h=a[0].rect.height
 pg=comp.new_page(width=2*w+20,height=h+22);pg.insert_text((5,13),name+' BEFORE',fontsize=8);pg.insert_text((w+25,13),name+' AFTER',fontsize=8)
 pg.show_pdf_page(fitz.Rect(0,22,w,h+22),a,0);pg.show_pdf_page(fitz.Rect(w+20,22,2*w+20,h+22),z,0)
comp.save(O/'review/marker_同尺寸前后对照.pdf',garbage=4,deflate=True)
# Figure5 a same-page crop comparison; not a separate redesigned figure.
comp=fitz.open();ims=[]
for file in [B/'main/Figure5.pdf',O/'main/Figure5.pdf']:
 d=fitz.open(file);clip=fitz.Rect(0,0,100*72/25.4,60*72/25.4);pm=d[0].get_pixmap(matrix=fitz.Matrix(3,3),clip=clip);ims.append(Image.frombytes('RGB',(pm.width,pm.height),pm.samples))
out=Image.new('RGB',(ims[0].width*2+20,ims[0].height+25),'white');dr=ImageDraw.Draw(out)
for i,(im,lab) in enumerate(zip(ims,['BEFORE','AFTER'])):out.paste(im,(i*(im.width+20),25));dr.text((i*(im.width+20)+5,5),lab,fill='black')
out.save(O/'review/Figure5a_前后对照.png')
# Keep attribution without copying historical datasets or redundant builds.
for file in ['FigureYa_ATTRIBUTION.md','FigureYa_LICENSE']:
 if (B/'notes'/file).exists():shutil.copy2(B/'notes'/file,O/'notes'/file)
# Record final a content boundary separately from panel allocation and actual data regions.
d=fitz.open(O/'main/Figure5.pdf');p=d[0];clip=fitz.Rect(0,0,100*72/25.4,60*72/25.4);pm=p.get_pixmap(matrix=fitz.Matrix(3,3),clip=clip);ar=np.asarray(Image.frombytes('RGB',(pm.width,pm.height),pm.samples));ys,xs=np.where(ar.min(2)<245)
content_mm=[float(xs.min())/3/72*25.4,float(ys.min())/3/72*25.4,float(xs.max()+1)/3/72*25.4,float(ys.max()+1)/3/72*25.4]
s=json.load(open(O/'audit/Figure5a_display.json'));s['full_content_boundary_including_letter_mm_xyxy']=content_mm;(O/'audit/Figure5a_display.json').write_text(json.dumps(s,indent=2))
hist=fitz.open(O/'review/Histology_ROI_three_view_review.pdf');assert len(hist)==5
status={'status':'PAUSED_FOR_MARKER_SIMPLIFICATION_AND_SPATIAL_REVIEW','layout':'COMPLETED_SCOPED_CHECKS_AUTHOR_REVIEW_PENDING','numerical_preservation':'CONFIRMED_WITHIN_SCOPED_SOURCE_AND_GEOMETRY_CHECKS','histology':'UNRESOLVED_AUTHOR_PATHOLOGY_REVIEW_REQUIRED','ORA_statistical_validity':'PENDING','scientific_reanalysis':False,'manuscript_modified':False,'unchanged_figures_copied_without_redraw':8,'modified_figures':sorted(changed),'Figure5_b_f_and_letters_identical_SVG_subtrees':True,'spatial_all_counts_coords_flags_and_limits_preserved':True,'histology_review_pages':5,'merged_pages':15,'merge_checks':mergechecks}
(O/'audit/final_checks.json').write_text(json.dumps(status,ensure_ascii=False,indent=2));print(json.dumps({k:v for k,v in status.items() if k!='merge_checks'},ensure_ascii=False,indent=2))
