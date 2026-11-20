from common import *
from PIL import Image,ImageDraw
BASE=Path('/Users/taki/Desktop/Panvascular_Marker_Communication_Enrichment_20260920');old=ET.parse(BASE/'editable_svg/Figure3.svg').getroot();new=ET.parse(O/'editable_svg/Figure3.svg').getroot()
for key in ['axes_1','figure3_b_content','axes_5','axes_6','axes_7','figure3_f_original_fill_key']:
 a=next(e for e in old.iter() if e.get('id')==key);b=next(e for e in new.iter() if e.get('id')==key);assert ET.tostring(a)==ET.tostring(b),key
assert not any(e.get('id')=='enrichc_figure_1' for e in new.iter())
# Exact retained scientific values: no transformation of NES, P, FDR or set membership.
m=pd.read_csv(O/'source_data/MAGMA_display9.tsv',sep='\t');full=pd.read_csv(O/'source_data/MAGMA_complete_GSEA.tsv',sep='\t');j=m.merge(full,on=['database','pathway'],suffixes=['_plot','_full'],validate='one_to_one')
for col in ['pval','padj','NES','size']:assert np.allclose(j[col+'_plot'],j[col+'_full'],rtol=0,atol=1e-12,equal_nan=True)
p=pd.read_csv(O/'source_data/pseudobulk_display9_all7_lineages.tsv',sep='\t');manifest=pd.read_csv(O/'source_data/pseudobulk_input_manifest.tsv',sep='\t');n=0
for _,row in manifest.iterrows():
 src=pd.read_csv(row.source,sep='\t');z=p[p.cell_state.eq(row.state)];qj=z.merge(src,on=['pathway','database'],suffixes=['_plot','_full'],validate='one_to_one');assert len(qj)==len(z)
 for col in ['pval','padj','NES','size']:assert np.allclose(qj[col+'_plot'],qj[col+'_full'],rtol=0,atol=1e-12,equal_nan=True)
 n+=len(qj)
assert n==63
versions=pd.read_csv(BASE/'FIGURE_VERSIONS.tsv',sep='\t');rows=[];book=fitz.open()
for _,r in versions.iterrows():
 changed=r.figure=='Figure3';pdf=O/'figures/Figure3.pdf' if changed else Path(r.delivered_pdf);svg=O/'editable_svg/Figure3.svg' if changed else Path(r.delivered_svg);doc=fitz.open(pdf);assert len(doc)==1;book.insert_pdf(doc);t,v=doc.xref_get_key(doc[0].xref,'Group')
 if t=='dict':assert not re.search(r'\d+ \d+ R',v);book.xref_set_key(book[-1].xref,'Group',v)
 rows.append({'page':int(r.page),'figure':r.figure,'source_pdf':str(pdf),'source_svg':str(svg),'version':O.name if changed else r.source_version,'status':'AUTHOR_REVIEW_PENDING' if changed else r.status,'width_mm':doc[0].rect.width/mm,'height_mm':doc[0].rect.height/mm,'sha256':sha(pdf),'svg_sha256':sha(svg)})
book.save(O/'全部图件_15页_遗传与转录通路更新.pdf',garbage=4,deflate=True);book.close();book=fitz.open(O/'全部图件_15页_遗传与转录通路更新.pdf');assert len(book)==15
for row,page in zip(rows,book):
 source=fitz.open(row['source_pdf'])[0];assert page.rect==source.rect;assert page.get_pixmap().samples==source.get_pixmap().samples,row['figure']
pd.DataFrame(rows).to_csv(O/'FIGURE_VERSIONS.tsv',sep='\t',index=False)
# Same-size, native-coordinate local comparisons.
rect=fitz.Rect(0,48*mm,183*mm,79.5*mm);ims=[]
for name,path in [('before',BASE/'figures/Figure3.pdf'),('after',O/'figures/Figure3.pdf')]:
 page=fitz.open(path)[0];pix=page.get_pixmap(matrix=fitz.Matrix(3,3),clip=rect);im=Image.frombytes('RGB',(pix.width,pix.height),pix.samples);im.save(O/'review'/f'Figure3c_{name}.png');ims.append(im)
canvas=Image.new('RGB',(max(i.width for i in ims),sum(i.height for i in ims)+50),'white');dd=ImageDraw.Draw(canvas);y=0
for label,im in zip(['BEFORE','AFTER'],ims):dd.text((5,y+3),label,fill='black');canvas.paste(im,(0,y+23));y+=im.height+25
canvas.save(O/'review/Figure3c_before_after.png')
# At 1:1 native scale, detect text collisions within separate c blocks (ignore dots/colorbar glyphs).
page=fitz.open(O/'figures/Figure3.pdf')[0];words=page.get_text('words');cwords=[w for w in words if 50*mm<w[1]<79*mm];pairs=[]
for i,w in enumerate(cwords):
 for v in cwords[i+1:]:
  r=fitz.Rect(w[:4])&fitz.Rect(v[:4])
  if not r.is_empty and r.width*r.height>.2:pairs.append([w[4],v[4],r.width*r.height])
# Document script results separately from visual review.
(O/'notes/VALIDATION.json').write_text(json.dumps({'exact_original_MAGMA_records':9,'exact_original_pseudobulk_records':n,'protected_Figure3_a_b_d_e_f_objects':'BYTE_IDENTICAL','other14_pages':'REUSED_UNCHANGED','all15_merge_native_geometry_and_render':'IDENTICAL','text_overlap_candidates':pairs,'visual_review':'Actual Figure3 and c detail reviewed; author approval pending','ORA':'PENDING; archived comparison retained; absent from new main c','no_scientific_models_or_tests_run':True},indent=2))
pd.DataFrame([{'file':str(f.relative_to(O)),'sha256':sha(f)} for f in (O/'source_data').iterdir()]).to_csv(O/'SOURCE_DATA_INDEX.tsv',sep='\t',index=False)
print('PASS: 9 MAGMA + 63 pseudobulk records; 15 pages; overlaps',pairs)
