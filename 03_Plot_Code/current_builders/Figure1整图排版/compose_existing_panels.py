"""Recompose the current panel files without redrawing unaffected panels."""
from pathlib import Path
import json,hashlib
from lxml import etree as E
import fitz

O=Path(__file__).parent; MM=72/25.4
NS='http://www.w3.org/2000/svg';XL='http://www.w3.org/1999/xlink';Q=lambda x:'{'+NS+'}'+x
manifest=json.loads((O/'panel_placement.json').read_text())
protected=json.loads((O/'Figure1b_protected_panel_hashes.json').read_text())
assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in protected.items())
root=E.Element(Q('svg'),nsmap={None:NS,'xlink':XL},width='183mm',height='170mm',viewBox=f'0 0 {183*MM} {170*MM}')
doc=fitz.open();page=doc.new_page(width=183*MM,height=170*MM)
for name in ['a','b','c','f','e','d','g']:
    m=manifest[name];xx,yy=m['origin_mm'];w,h=m['width_mm'],m['height_mm']
    r=E.parse(str(O/'panels'/f'Figure1{name}.svg')).getroot()
    ids={e.get('id'):f'{name}_{e.get("id")}' for e in r.iter() if e.get('id')}
    for e in r.iter():
        for k,v in list(e.attrib.items()):
            if k=='id': e.set(k,ids[v])
            else:
                for old,newid in ids.items():
                    v=v.replace('url(#'+old+')','url(#'+newid+')')
                    v=('#'+newid) if v=='#'+old else v
                e.set(k,v)
    g=E.SubElement(root,Q('g'),id='panel_'+name,transform=f'translate({xx*MM} {yy*MM})')
    for e in list(r):g.append(e)
    with fitz.open(O/'panels'/f'Figure1{name}.pdf') as src:
        assert abs(src[0].rect.width-w*MM)<.1 and abs(src[0].rect.height-h*MM)<.1
        page.show_pdf_page(fitz.Rect(xx*MM,yy*MM,(xx+w)*MM,(yy+h)*MM),src,0)
page.insert_font(fontname='PanelArialBold',fontfile='/System/Library/Fonts/Supplemental/Arial Bold.ttf')
for name,m in manifest.items():
    xx,yy=m['origin_mm'];px=(xx+1)*MM;py=(yy+2.2)*MM+2.8
    page.insert_text((px,py),name,fontname='PanelArialBold',fontsize=8,color=(.133,.157,.192))
    t=E.SubElement(root,Q('text'),id='panel_letter_'+name,x=str(px),y=str(py),fill='#222831',
        **{'font-family':'Arial','font-weight':'bold','font-size':'8'});t.text=name
doc.save(O/'Figure1.pdf',garbage=4,deflate=True)
page.get_pixmap(dpi=200).save(O/'review/Figure1.png');doc.close()
E.ElementTree(root).write(str(O/'Figure1.svg'),xml_declaration=True,encoding='utf-8')
stats=json.loads((O/'preserved_data_checks.json').read_text())
stats['b']=json.loads((O/'review/Figure1b_layout_check.json').read_text())
(O/'preserved_data_checks.json').write_text(json.dumps(stats,indent=2))
manifest['b'].update({'axes':[], 'layout_source':'figure1b_layout.py',
 'primary_data_box_mm':[18,7,46,38.312824751456525],
 'tissue_data_box_mm':[80.5,7,25.5,21.238631112220465],
 'tissue_legend_box_mm':[80.5,30.7,25.5,18.7]})
(O/'panel_placement.json').write_text(json.dumps(manifest,indent=2))
print('Current preview recomposed. Unaffected panel assets unchanged:',len(protected))
