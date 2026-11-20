from pathlib import Path
import fitz,json,csv,hashlib,shutil,base64,html,collections
import numpy as np,pandas as pd,xml.etree.ElementTree as E
from PIL import Image,ImageDraw
O=Path(__file__).resolve().parents[1];B=Path('/Users/taki/Desktop/Panvascular_FigureYa_20260919');P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
N={'s':'http://www.w3.org/2000/svg'}
seq=[('Figure'+str(i),'main') for i in range(1,6)]+[('ExtendedData_Figure'+str(i),'extended_data') for i in range(1,11)]
merged=fitz.open();rows=[]
for i,(name,folder) in enumerate(seq,1):
 pdf=O/folder/(name+'.pdf');svg=O/'editable_svg'/(name+'.svg');doc=fitz.open(pdf);assert len(doc)==1
 merged.insert_pdf(doc);pg=doc[0];spans=[s for b in pg.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans'] if s['text'].strip()]
 outside=[s['text'] for s in spans if not fitz.Rect(-.2,-.2,pg.rect.width+.2,pg.rect.height+.2).contains(fitz.Rect(s['bbox']))];assert not outside,(name,outside)
 rows.append({'page':i,'figure':name,'pdf':str(pdf),'svg':str(svg),'width_mm':round(pg.rect.width/72*25.4,3),'height_mm':round(pg.rect.height/72*25.4,3),'pdf_sha256':sha(pdf),'svg_sha256':sha(svg),'baseline_pdf_sha256':sha(B/folder/(name+'.pdf')),'minimum_extracted_font_pt':round(min(s['size'] for s in spans),2),'data_integrity':'geometry/source checked; display flags in Figure5a only','layout_legibility':'bbox/legend/gutter checks; no page overflow','visual_review':'actual PDF previews reviewed; author confirmation pending'})
mp=O/'全部图件_内容对齐与空间核查_15页审阅.pdf';merged.save(mp,garbage=4,deflate=True);merged.close();merged=fitz.open(mp)
assert len(merged)==15
checks=[]
for i,r in enumerate(rows):
 d=fitz.open(r['pdf']);assert d[0].rect==merged[i].rect
 a=np.frombuffer(d[0].get_pixmap().samples,dtype=np.uint8).astype(int);b=np.frombuffer(merged[i].get_pixmap().samples,dtype=np.uint8).astype(int);delta=np.abs(a-b);assert delta.max()<=1,(r['figure'],delta.max())
 checks.append({'figure':r['figure'],'page_geometry_identical':True,'max_8bit_channel_delta':int(delta.max()),'different_channels':int(np.count_nonzero(delta))})
with (O/'图件版本与尺寸.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(rows)
# Baseline and original submission must remain byte-identical to their existing manifests.
sourcechecks=[]
for r in csv.DictReader((B/'图件版本与尺寸.tsv').open(),delimiter='\t'):
 for key in ['pdf','svg']:assert sha(r[key])==r[key+'_sha256'];sourcechecks.append(r[key])
for r in csv.DictReader((P/'文件来源与校验.tsv').open(),delimiter='\t'):assert sha(P/r['file'])==r['sha256'];sourcechecks.append(str(P/r['file']))
for p,h in json.load(open(O/'audit/spatial_counts_input_sha256.json')).items():assert sha(p)==h;sourcechecks.append(p)
# Figure5b-f scientific paths, expressions and categorical colours must be unchanged.
old=E.parse(B/'editable_svg/Figure5.svg').getroot();new=E.parse(O/'editable_svg/Figure5.svg').getroot();groups=['axes_5','axes_6','FigureYa6_native_forest','axes_8','axes_9','axes_10']
def signature(el):
 return {'paths':[(p.get('d'),p.get('style')) for p in el.findall('.//s:path',N)],'images':[p.get('{http://www.w3.org/1999/xlink}href') for p in el.findall('.//s:image',N)],'text':[(p.text,p.get('style'),p.get('x'),p.get('y')) for p in el.findall('.//s:text',N)]}
for g in groups:assert signature(old.find(f".//s:g[@id='{g}']",N))==signature(new.find(f".//s:g[@id='{g}']",N)),g
f5=fitz.open(O/'main/Figure5.pdf');p=f5[0];matrix=[im for im in p.get_images() if im[2:4]==(29,7)];assert len(matrix)==1 and f5.xref_get_key(matrix[0][0],'Interpolate')[1]=='false'
words=[s for b in p.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans']]
def box(t):return fitz.Rect(next(s['bbox'] for s in words if s['text']==t))
fgap=box('CM01 | fibrous cap').y0-max(box('Dataset').y1,box('Median').y1);assert fgap>1.5,fgap
# True Figure5a before/after crop of the exported mother PDFs (same page scale).
ims=[]
for path in [B/'main/Figure5.pdf',O/'main/Figure5.pdf']:
 d=fitz.open(path);pix=d[0].get_pixmap(matrix=fitz.Matrix(3,3),clip=fitz.Rect(0,0,100*72/25.4,60*72/25.4));ims.append(Image.frombytes('RGB',(pix.width,pix.height),pix.samples))
out=Image.new('RGB',(ims[0].width*2+20,ims[0].height+25),'white');dr=ImageDraw.Draw(out)
for i,(im,lab) in enumerate(zip(ims,['BEFORE','AFTER'])):out.paste(im,(i*(im.width+20),25));dr.text((i*(im.width+20)+5,5),lab,fill='black')
out.save(O/'review/Figure5a_actual_before_after.png')
# Figure3 focused comparison as an easy-review, editable-PDF-preserving pair of pages.
comp=fitz.open()
for region,rect in [('b_c',fitz.Rect(0,0,518.74,226)),('d_f',fitz.Rect(0,225,518.74,481.89))]:
 before=fitz.open(B/'main/Figure3.pdf');after=fitz.open(O/'main/Figure3.pdf');pg=comp.new_page(width=rect.width*2+24,height=rect.height+24)
 pg.insert_text((4,13),'BEFORE',fontsize=8);pg.insert_text((rect.width+20,13),'AFTER',fontsize=8)
 pg.show_pdf_page(fitz.Rect(0,24,rect.width,24+rect.height),before,0,clip=rect)
 pg.show_pdf_page(fitz.Rect(rect.width+24,24,rect.width*2+24,24+rect.height),after,0,clip=rect)
comp.save(O/'review/Figure3_局部前后对照.pdf',garbage=4,deflate=True)
# Standalone offline HTML with embedded audit images and literal report text.
report=(O/'notes/Figure5a_核查记录.md').read_text();htmlparts=['<!doctype html><meta charset="utf-8"><title>Figure5a核查</title><style>body{font:16px sans-serif;color:#222;max-width:1300px;margin:32px auto;padding:0 24px}pre{white-space:pre-wrap;font:15px/1.6 sans-serif}img{max-width:100%;height:auto}table{border-collapse:collapse;font-size:13px}td,th{border:1px solid #ddd;padding:6px}</style><h1>Figure5a组织学疑点与绘图核查</h1><pre>'+html.escape(report)+'</pre>']
for path in [O/'review/Figure5a_actual_before_after.png']+sorted((O/'review').glob('spatial_ROI_*.png')):
 htmlparts.append('<h2>'+html.escape(path.stem)+'</h2><img src="data:image/png;base64,'+base64.b64encode(path.read_bytes()).decode()+'">')
for pth in sorted((O/'audit').glob('spatial_counts_roi_*_center_included.tsv')):
 if 'adjacent' in pth.name:continue
 z=pd.read_csv(pth,sep='\t');z['spot_number_in_ROI_image']=np.arange(1,len(z)+1)
 cols=['spot_number_in_ROI_image','barcode','raw_gene_count','matrix_total_UMI_all_features','raw_gene_CP10k','raw_gene_log1p_CP10k','in_tissue','in_plot','entire_physical_spot_disk_inside_roi'];htmlparts+=['<h2>'+html.escape(pth.stem)+'</h2>',z[cols].to_html(index=False,float_format=lambda v:f'{v:.6g}')]
(O/'Figure5a_核查记录_含图.html').write_text('\n'.join(htmlparts))
for fn in ['ATTRIBUTION.md','LICENSE']:shutil.copy2(B/'templates'/fn,O/'notes'/('FigureYa_'+fn))
checksout={'status':'PAUSED_FOR_CONTENT_ALIGNMENT_AND_FIGURE5A_REVIEW','pages':15,'source_files_verified_unchanged':len(sourcechecks),'science_recomputed':False,'raw_count_validation_only':True,'manuscript_changed':False,'Figure5_b_f_geometry_and_values_unchanged':True,'Figure5b_interpolation':False,'Figure5f_legend_to_first_label_gap_pt':fgap,'merge_checks':checks,'data_integrity':'PASS scoped source and geometry checks','layout_legibility':'script checks plus actual agent/root rendered-PDF inspection','visual_review':'AUTHOR_REVIEW_PENDING; LIPA histological source remains unresolved','source_checks':sourcechecks}
(O/'audit/final_checks.json').write_text(json.dumps(checksout,ensure_ascii=False,indent=2));print(json.dumps({k:v for k,v in checksout.items() if k not in ['source_checks','merge_checks']},ensure_ascii=False,indent=2))
