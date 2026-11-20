"""Install a display-only revision of Figure 4a in the current package.

The current delivered PDF is the base. Only the spatial schematic is replaced;
the CoVarNet branch, separate panel letters and all data panels are preserved.
"""
from pathlib import Path
from lxml import etree as E
import copy, hashlib, json, shutil
import fitz
import numpy as np

P=Path('/Users/taki/Desktop/泛血管稿件_批注修订_20261116')
W=Path(__file__).resolve().parents[1]
B=W/'backup_before_plain_spatial_labels'
R=P/'图件编辑资料/全图核对与交付_20261118'
MM=72/25.4
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
original_hashes=json.loads((B/'input_hashes.json').read_text())
clean={k:'' for k in ['title','author','subject','keywords','creator','producer','creationDate','modDate','trapped']}

# Strip export metadata, retaining vector text and graphics.
standalone=W/'panels/Figure4a.pdf'
doc=fitz.open(standalone);doc.set_metadata(clean);doc.del_xml_metadata()
tmp=standalone.with_suffix('.clean.pdf');doc.save(tmp,garbage=4,deflate=True);doc.close();tmp.replace(standalone)
sr=E.parse(str(W/'panels/Figure4a.svg')).getroot()
for e in sr.xpath('//*[local-name()="metadata"]'):e.getparent().remove(e)
E.ElementTree(sr).write(str(W/'panels/Figure4a.svg'),encoding='utf-8',xml_declaration=True)
assert not [e for e in sr.xpath('//*[local-name()="text"]') if ''.join(e.itertext()).strip()=='a']

src=B/'Figure/Figure4.pdf'
old=fitz.open(src);d=fitz.open(src);pg=d[0]
clip=fitz.Rect(92*MM,0,183*MM,25*MM)
pg.add_redact_annot(clip,fill=(1,1,1))
pg.apply_redactions(images=2,graphics=1,text=0)
pa=fitz.open(standalone)
pg.show_pdf_page(clip,pa,0,clip=clip,keep_proportion=True)
d.set_metadata(clean);d.del_xml_metadata()
out=P/'Figure/Figure4.pdf'
d.save(out,garbage=4,deflate=True);d.close()
new=fitz.open(out)
assert new[0].rect==old[0].rect
text=new[0].get_text()
assert 'Adjusted associations' not in text and 'Spatial associations' not in text
assert 'Gene expression in tissue' in text and 'Compare within sections' in text
mat=fitz.Matrix(3,3)
o=old[0].get_pixmap(matrix=mat,alpha=False);n=new[0].get_pixmap(matrix=mat,alpha=False)
aa=np.frombuffer(o.samples,np.uint8).reshape(o.height,o.width,o.n)
bb=np.frombuffer(n.samples,np.uint8).reshape(n.height,n.width,n.n)
changed=np.any(aa!=bb,axis=2)
mask=np.ones(changed.shape,bool)
mask[:int(np.ceil(clip.y1*3))+2,int(np.floor(clip.x0*3))-2:]=False
outside=int(changed[mask].sum());assert outside==0,outside

# Replace only the SVG panel-a subtree; leave b-d and their text anchors intact.
dest=P/'Figure/Editable_SVG/Figure4.svg'
tree=E.parse(str(B/'Figure/Editable_SVG/Figure4.svg'));root=tree.getroot()
group=root.xpath('//*[@id="panel_a"]')[0]
protected_svg={i:E.tostring(root.xpath(f'//*[@id="panel_{i}"]')[0]) for i in 'bcd'}
for child in list(group):group.remove(child)
ids={el.get('id'):'a_'+el.get('id') for el in sr.iter() if el.get('id')}
for el in sr.iter():
    if el.get('id'):el.set('id',ids[el.get('id')])
    for key,value in list(el.attrib.items()):
        for key_old,key_new in ids.items():
            value=value.replace(f'url(#{key_old})',f'url(#{key_new})')
            if value==f'#{key_old}':value='#'+key_new
        el.set(key,value)
for child in sr:group.append(copy.deepcopy(child))
for i,v in protected_svg.items():assert E.tostring(root.xpath(f'//*[@id="panel_{i}"]')[0])==v
tree.write(str(dest),encoding='utf-8',xml_declaration=True)

# Synchronise the existing working copies and previews, not other figures.
for ext,src in [('pdf',out),('svg',dest)]:
    for target in [P/'Figures'/f'Figure4.{ext}',W/'review'/f'Figure4.{ext}']:
        shutil.copy2(src,target)
for target,dpi in [(W/'review/Figure4.png',180),(W/'review/Figure4_actual_size_96dpi.png',96),(R/'previews_after/Figure4.png',180)]:
    new[0].get_pixmap(matrix=fitz.Matrix(dpi/72,dpi/72),alpha=False).save(target)
new[0].get_pixmap(matrix=fitz.Matrix(3,3),clip=fitz.Rect(0,0,183*MM,25*MM),alpha=False).save(W/'review/Figure4a_preview.png')

merged_checks=[]
for path in [P/'Figure/Figures_1-4_and_Extended_Data.pdf',P/'Figures_and_ExtendedData_审阅.pdf']:
    if not path.exists():continue
    source=fitz.open(path);assert len(source)==14
    target=fitz.open()
    for i in range(14):target.insert_pdf(new if i==3 else source,from_page=0 if i==3 else i,to_page=0 if i==3 else i)
    target.set_metadata(clean);target.del_xml_metadata()
    tmp=path.with_suffix('.tmp.pdf');target.save(tmp,garbage=4,deflate=True);target.close()
    verify=fitz.open(tmp)
    for i in range(14):
        base=new[0] if i==3 else source[i]
        assert base.rect==verify[i].rect
        assert base.get_text()==verify[i].get_text()
        assert base.get_pixmap(alpha=False).samples==verify[i].get_pixmap(alpha=False).samples
    source.close();verify.close();tmp.replace(path)
    merged_checks.append(str(path))

allowed={str(out),str(dest),str(P/'Figures/Figure4.pdf'),str(P/'Figures/Figure4.svg'),*merged_checks}
for path,h in original_hashes.items():
    if path not in allowed and not path.startswith(str(W)):
        assert sha(path)==h, f'Unrelated file changed: {path}'
for letter in 'bcd':
    for ext in ['svg','pdf','png']:
        assert (W/'panels'/f'Figure4{letter}.{ext}').read_bytes()==(P/'图件编辑资料/Figure4/panels'/f'Figure4{letter}.{ext}').read_bytes()

manifest=json.loads((R/'delivery_manifest.json').read_text())
for item in manifest:
    if item['figure']=='Figure4':
        item['pdf_sha256']=sha(out);item['source_sha256']=sha(P/'Figures/Figure4.pdf')
        item['changes']='Spatial schematic uses concrete expression/cell-abundance comparison and pooling steps; CoVarNet and numerical panels unchanged'
        item['outside_correction_pixels_changed']=outside
(R/'delivery_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
report={'changed_panel':'Figure4a spatial branch','canvas_mm':[183,170],'panel_mm':[183,25],
        'outside_spatial_branch_changed_pixels':outside,'protected_bcd_byte_identical':True,
        'manuscript_tables_and_other_figures_unchanged':True,'scientific_analyses_run':False,
        'label_is_separate':True,'combined_previews_updated':merged_checks,
        'pdf_sha256':sha(out),'biorender_status':'PENDING_SYNC',
        'labels':{'Spatial associations':'Gene expression in tissue','Adjusted associations':'Compare within sections','Meta-analysis':'Combine results'}}
(W/'review/plain_spatial_labels_checks.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
(W/'README.md').write_text('''# Current Figure 4a

The spatial branch now reads: Gene expression in tissue → compare gene expression with estimated cell abundance within sections → combine results across 35 analysis units. Counts remain 49 sections and 29 candidates. Concrete RNA, cell and pooling symbols replace the abstract “Adjusted associations” label. These are schematic workflow motifs, not measured spatial patterns. The original conditional model, adjustment terms and pooling method remain in the existing Methods and Figure 4 legend.

Current source: `code/redraw_figure4a.py`. Current artwork: `panels/Figure4a.svg` and `.pdf`. Current full figure: `review/Figure4.pdf`. Native panel size: 183 × 25 mm; full figure: 183 × 170 mm. The CoVarNet branch and panels b–d are unchanged. The letter a remains separate in the assembly and BioRender.

BioRender: https://app.biorender.com/illustrations/6af982741be2d4f7c66d0f97
''')
print(json.dumps(report,ensure_ascii=False,indent=2))
