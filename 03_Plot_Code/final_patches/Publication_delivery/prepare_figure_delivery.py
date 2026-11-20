"""Bounded legend and display corrections; no scientific calculation is run.

PDF data objects are reused. SVG edits are separate, editable layout objects.
Each changed region is recorded for pixel-level outside-region verification.
"""
from pathlib import Path
import csv, json, hashlib, re, shutil, zipfile, copy
from collections import Counter
import fitz
import numpy as np
from lxml import etree as ET

P = Path('/Users/taki/Desktop/泛血管稿件_批注修订_20261116')
OLD = Path('/Users/taki/Desktop/投稿包_最新版_20260919')
R = P/'图件编辑资料/全图核对与交付_20261118'
O = P/'Figure'
SVG = O/'Editable_SVG'
for x in [O,SVG,R/'previews_after',R/'backup',R/'checks']:
    x.mkdir(exist_ok=True, parents=True)
REG = '/System/Library/Fonts/Supplemental/Arial.ttf'
BOLD = '/System/Library/Fonts/Supplemental/Arial Bold.ttf'
NS='http://www.w3.org/2000/svg'; Q=lambda t:'{'+NS+'}'+t
INK='#22262B'
CLEAN_META={k:'' for k in ['title','author','subject','keywords','creator','producer','creationDate','modDate','trapped']}
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rows(path):
    with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def rgb(h):return tuple(int(h[i:i+2],16)/255 for i in [1,3,5])
source_hashes={str(f.relative_to(P)):sha(f) for f in (P/'Source_Data_Figures').rglob('*') if f.is_file()}

def svg_text(root,x,y,text,size=6,bold=False,rotate=None):
    e=ET.SubElement(root,Q('text'),x=str(x),y=str(y),fill=INK)
    e.set('style',f'font-family:Arial; font-size:{size}px;'+('font-weight:700;' if bold else ''))
    if rotate is not None:e.set('transform',f'rotate({rotate} {x} {y})')
    e.text=text
    return e
def font(page,bold=False):
    name='DeliveryArialB' if bold else 'DeliveryArial'
    page.insert_font(fontname=name,fontfile=BOLD if bold else REG,set_simple=False)
    return name
def text(page,root,x,y,t,size=6,bold=False):
    page.insert_text((x,y),t,fontsize=size,fontname=font(page,bold),color=rgb(INK))
    svg_text(root,x,y,t,size,bold)
def circle(page,root,x,y,r,color,filled=True):
    page.draw_circle((x,y),r,color=rgb(color),fill=rgb(color) if filled else (1,1,1),width=.65)
    ET.SubElement(root,Q('circle'),cx=str(x),cy=str(y),r=str(r),fill=color if filled else '#FFFFFF',stroke=color,**{'stroke-width':'.65'})
def line(page,root,x,y,w,dashed=False):
    page.draw_line((x,y),(x+w,y),color=rgb('#4B5157'),width=.75,dashes='[2 1.5] 0' if dashed else None)
    e=ET.SubElement(root,Q('line'),x1=str(x),x2=str(x+w),y1=str(y),y2=str(y),stroke='#4B5157',**{'stroke-width':'.75'})
    if dashed:e.set('stroke-dasharray','2 1.5')
def set_pdf_stream(doc,page,before,after):
    count=0
    for x in page.get_contents():
        s=doc.xref_stream(x);n=s.count(before)
        if n:doc.update_stream(x,s.replace(before,after));count+=n
    assert count==1,(before,count)
def clean_svg(root):
    for e in root.xpath('//*[local-name()="metadata"]'):e.getparent().remove(e)
def replace_base_panel_font(doc,page):
    old=[f for f in page.get_fonts(full=True) if not doc.extract_font(f[0])[3]]
    for f in old:
        assert f[3]=='Helvetica-Bold',f
        affected=[s for b in page.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans'] if s['font']=='Helvetica-Bold']
        assert all(s['text'] in list('abcde') for s in affected),affected
        new=page.insert_font(fontname='EmbeddedPanelArialBold',fontfile=BOLD,set_simple=True)
        doc.xref_copy(new,f[0])
    return [(s['bbox']) for s in affected] if old else []

names=[f'Figure{i}' for i in range(1,5)]+[f'ExtendedData_Figure{i}' for i in range(1,11)]
changes={};allowed={};manifest=[];numeric_checks=[]
data3=rows(P/'Source_Data_Figures/Figure3/Fig3a_L1.tsv')
assert len(data3)==15
for name in names:
    ismain=name.startswith('Figure')
    src=P/('Figures' if ismain else 'Extended_Data')/(name+'.pdf')
    d=fitz.open(src);pg=d[0]
    assert len(d)==1
    source_svg=P/('Figures' if ismain else 'Extended_Data')/(name+'.svg')
    if not source_svg.exists():
        n=int(name.split('Figure')[-1]);oldn=n if n<6 else n+1 if n<10 else 6
        source_svg=OLD/'图件编辑资料/editable_svg'/f'ExtendedData_Figure{oldn}.svg'
        if n==10:source_svg=None
    root=ET.parse(str(source_svg)).getroot() if source_svg else ET.fromstring(pg.get_svg_image(text_as_path=False).encode())
    groups=ET.SubElement(root,Q('g'),id='publication_legend_corrections_20261118')
    changes[name]=[];allowed[name]=[]
    if name=='Figure2':
        # Empty space above the lower factor ellipses; no paths/loadings moved.
        line(pg,groups,334,276,13)
        text(pg,groups,351,278.1,'P < 0.05',6)
        line(pg,groups,405,276,13,True)
        text(pg,groups,422,278.1,'P ≥ 0.05',6)
        allowed[name]+=[(330,267,465,283)]
        changes[name].append('2c: added solid/dashed path key using original nominal P < 0.05 threshold')
    if name=='Figure3':
        # The existing scatter, all 15 coordinates, size and filled status remain.
        cols=[]
        for i,row in enumerate(data3,1):
            use=root.xpath(f'//*[@id="identity_a_PathCollection_{i}"]//*[local-name()="use"]')[0]
            st=use.get('style','');color=re.search(r'stroke:\s*(#[0-9a-fA-F]+)',st).group(1)
            xx=float(use.get('x'));yy=float(use.get('y'))
            xexpect=104.951027+float(row['BETA_STD'])*(221.102362-104.951027)/.04
            yexpect=77.163803-float(row['primary_r0.1__assoc_mcz'])*(77.163803-31.921182)/6
            assert abs(xx-xexpect)<2e-5 and abs(yy-yexpect)<2e-5
            filled='fill: #ffffff' not in st
            expect=float(row['FDR_BH'])<.05 and float(row['primary_r0.1__assoc_mcp_BH'])<.05
            assert filled==expect
            cols.append((row['group'],color))
        for i,(label,color) in enumerate(cols):
            y=22+i*6.85
            circle(pg,groups,240.5,y-2.0,1.65,color)
            text(pg,groups,245.5,y,label,5.8)
        circle(pg,groups,118,132.3,2.2,'#747C86',False)
        text(pg,groups,124.5,134.4,'Not both FDR < 0.05',6)
        allowed[name]+=[(236,12,321,124),(113,126,225,138)]
        changes[name].append('3a/b: complete L1 colour key; open-symbol key added; all 15 plotted coordinates and joint FDR states verified')
        numeric_checks.append({'check':'Figure3a 15 coordinates and filled/open statuses','result':'PASS'})
    if name=='ExtendedData_Figure4':
        # One shared text matrix drives three vertical lineage names.
        set_pdf_stream(d,pg,b'0 -6 -6 0 185.804688 150.90625 Tm',b'0 -6 -6 0 169.804688 150.90625 Tm')
        moved=0
        for e in root.xpath('//*[local-name()="text"]'):
            tr=e.get('transform','')
            if 'translate(185.803062 ' in tr:
                e.set('transform',tr.replace('translate(185.803062 ','translate(169.803062 '));moved+=1
        assert moved==3,moved
        allowed[name]+=[(161,104,189,160),(161,300,189,396)]
        changes[name].append('4b/d: shifted only vertical lineage group labels 16 pt left; removed overlap with complete state names')
    if name=='ExtendedData_Figure5':
        # The old a lies under a later white map panel. Remove it and restore atop.
        set_pdf_stream(d,pg,b'8 0 0 -8 11.339844 14.230469 Tm\n(a)Tj',b'8 0 0 -8 11.339844 14.230469 Tm\n()Tj')
        es=[e for e in root.xpath('//*[local-name()="text"]') if ''.join(e.itertext()).strip()=='a']
        assert len(es)==1
        e=es[0];e.getparent().remove(e)
        text(pg,groups,11.339844,14.230469,'a',8,True)
        allowed[name]+=[(8,5,20,18)]
        changes[name].append('5a: restored obscured panel letter as a separate editable overlay')
    if name=='ExtendedData_Figure8':
        text(pg,groups,344.0,461.2,'Module matching similarity (cosine)',6)
        allowed[name]+=[(340,453,460,464)]
        changes[name].append('8d: labelled colour bar with module matching cosine similarity; original 18 values retained')
    if name in ['Figure4','ExtendedData_Figure9']:
        boxes=replace_base_panel_font(d,pg)
        allowed[name]+=[tuple(fitz.Rect(b)+(-1,-1,1,1)) for b in boxes]
        changes[name].append('Embedded Arial Bold for existing panel letters; numerical panels unchanged')
    # Strip personal/edit metadata in deliverables, never data.
    d.set_metadata(CLEAN_META);d.del_xml_metadata()
    d.subset_fonts()
    target=O/(name+'.pdf')
    if target.exists() and not (R/'backup'/target.name).exists():shutil.copy2(target,R/'backup'/target.name)
    d.save(target,garbage=4,deflate=True)
    clean_svg(root)
    if len(groups)==0:root.remove(groups)
    ET.ElementTree(root).write(str(SVG/(name+'.svg')),encoding='utf-8',xml_declaration=True)
    out=fitz.open(target);p2=out[0]
    mat=fitz.Matrix(2,2)
    a=fitz.open(src)[0].get_pixmap(matrix=mat,alpha=False)
    b=p2.get_pixmap(matrix=mat,alpha=False)
    assert (a.width,a.height)==(b.width,b.height)
    aa=np.frombuffer(a.samples,np.uint8).reshape(a.height,a.width,a.n)
    bb=np.frombuffer(b.samples,np.uint8).reshape(b.height,b.width,b.n)
    delta=np.any(aa!=bb,axis=2)
    for rect in allowed[name]:
        r=fitz.Rect(rect);x0=max(0,int(r.x0*2)-2);y0=max(0,int(r.y0*2)-2);x1=min(a.width,int(r.x1*2)+3);y1=min(a.height,int(r.y1*2)+3)
        delta[y0:y1,x0:x1]=False
    outside=int(delta.sum())
    assert outside==0,(name,'Unexpected pixels outside correction regions',outside)
    embedded=all(bool(out.extract_font(f[0])[3]) for f in p2.get_fonts(full=True))
    assert embedded,(name,'unembedded font')
    p2.get_pixmap(matrix=fitz.Matrix(1800/p2.rect.width,1800/p2.rect.width)).save(R/'previews_after'/(name+'.png'))
    manifest.append({'figure':name,'source_pdf':str(src),'source_sha256':sha(src),'delivered_pdf':str(target),'pdf_sha256':sha(target),'editable_svg':str(SVG/(name+'.svg')),'width_mm':round(p2.rect.width*25.4/72,3),'height_mm':round(p2.rect.height*25.4/72,3),'fonts_embedded':embedded,'outside_correction_pixels_changed':outside,'changes':'; '.join(changes[name]) or 'Original content reused; metadata removed'})

# Existing summary statistics are checked, not recomputed from expression or models.
coverage=rows(P/'Source_Data_Figures/Figure1_current/Figure1e_L2_coverage_plot_data.tsv')
assert len(coverage)==116 and len(set(r['state'] for r in coverage))==116
marker=rows(P/'Source_Data_Figures/Current_Display_20260920/marker_audit_20260921/original_all392_ordered.tsv')
assert len(marker)==392
est=rows(P/'Source_Data_Figures/Figure4_current/Figure5b_60_gene_unit_estimates.tsv')
assert Counter(r['gene'] for r in est)=={'FN1':34,'SCARB1':26}
assert len(set(r['biological_unit_id'] for r in est))==34
secs=rows(P/'Source_Data_Figures/Figure4_current/section_input_audit.tsv')
assert len(secs)==49 and len(set(r['biological_unit_id'] for r in secs))==35
all29=rows(P/'Source_Data_Figures/Figure4_current/Fig5c_all29.tsv');assert len(all29)==29
states=rows(P/'Source_Data_Figures/ExtendedData3_4_5/ED345_all_named_state_index.tsv')
counts=Counter(r['level1_annotation'] for r in states)
assert sum(counts[x] for x in ['Fibroblast','Fibromyocyte','Pericyte','Smooth muscle cell'])==41
assert counts['Endothelial cell']==9
assert counts['T cell']+counts['Natural killer cell']==16
assert sum(counts[x] for x in ['Dendritic cell','Macrophage','Monocyte','Neutrophil'])==38
numeric_checks.extend([{'check':k,'result':'PASS'} for k in [
 'Figure1g:116 distinct named states','Figure1h:392 stored cells in marker matrix (14 x 28)',
 'Figure4c:34 FN1 + 26 SCARB1 estimates; 34 distinct recorded IDs',
 'Figure4a:49 sections and 35 original analysis units','Figure4d:all29 candidates retained',
 'ED3/4/5:41 mural,9 EC,16 T/NK,38 myeloid named states']])

# Assemble native pages with no scaling or clipping.
merged=fitz.open()
for name in names:merged.insert_pdf(fitz.open(O/(name+'.pdf')))
merged.set_metadata(CLEAN_META);merged.save(O/'Figures_1-4_and_Extended_Data.pdf',garbage=4,deflate=True)
m=fitz.open(O/'Figures_1-4_and_Extended_Data.pdf')
merge_checks=[]
for i,name in enumerate(names):
    d=fitz.open(O/(name+'.pdf'))
    assert m[i].rect==d[0].rect
    assert m[i].get_text()==d[0].get_text()
    a=m[i].get_pixmap();b=d[0].get_pixmap()
    aa=np.frombuffer(a.samples,np.uint8);bb=np.frombuffer(b.samples,np.uint8)
    max_delta=int(np.abs(aa.astype(int)-bb.astype(int)).max())
    # PDF coordinate serialization can cause 1–2/255 antialias differences.
    # Native page rectangles and all text/figure content must remain unchanged.
    assert max_delta<=2,(name,max_delta)
    merge_checks.append({'figure':name,'native_page_size_unchanged':True,'extracted_text_identical':True,'max_render_channel_delta_255':max_delta})
assert source_hashes=={str(f.relative_to(P)):sha(f) for f in (P/'Source_Data_Figures').rglob('*') if f.is_file()}
(R/'checks/scientific_source_hashes.json').write_text(json.dumps(source_hashes,ensure_ascii=False,indent=2))
(R/'checks/numeric_checks.json').write_text(json.dumps(numeric_checks,ensure_ascii=False,indent=2))
(R/'checks/correction_regions.json').write_text(json.dumps(allowed,indent=2))
(R/'checks/merged_page_checks.json').write_text(json.dumps(merge_checks,indent=2))
(R/'delivery_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
with (O/'图件版本与尺寸.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(manifest[0]),delimiter='\t');w.writeheader();w.writerows(manifest)
print(json.dumps({'figures':len(names),'all_fonts_embedded':True,'all_outside_region_pixel_checks_passed':True,'merged_pages':len(m),'data_source_files_unchanged':len(source_hashes),'output':str(O)},ensure_ascii=False))
