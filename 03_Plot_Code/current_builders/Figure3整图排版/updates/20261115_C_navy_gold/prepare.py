from pathlib import Path
from zipfile import ZipFile
from lxml import etree as E
import copy, json, hashlib, shutil
import fitz

P=Path('/Users/taki/Desktop/投稿包_最新版_20260919')
A=P/'图件编辑资料/BioRender_当前素材/Figure3整图排版'
V=P/'图件编辑资料/Figure3ef_优化预览_20261115'
R=Path(__file__).resolve().parent
S=R/'stage'
for d in ['panels','review','documents','source_data','upload_ready']:(S/d).mkdir(parents=True,exist_ok=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
meta=json.loads((V/'audit.json').read_text())
for v in meta['source_hashes'].values():
    assert sha(v['path'])==v['sha256'],v['path']
    shutil.copy2(v['path'],S/'source_data'/Path(v['path']).name)
for p,h in meta['current_formal_files_untouched'].items():assert sha(p)==h,('formal changed',p)
for k,src in [('e',V/'配色调整/Figure3e_C_navy_gold'),('f',V/'panels/Figure3f_regional_GWAS')]:
    for ext in ['svg','pdf','png']:shutil.copy2(src.with_suffix('.'+ext),S/'panels'/f'Figure3{k}.{ext}')
    shutil.copy2(S/'panels'/f'Figure3{k}.svg',S/'upload_ready'/f'Figure3{k}_C_20261115.svg')
MM=72/25.4
old=fitz.open(P/'Figures/Figure3.pdf');new=fitz.open();pg=new.new_page(width=old[0].rect.width,height=old[0].rect.height)
for clip in [fitz.Rect(0,0,183*MM,81*MM),fitz.Rect(0,81*MM,117*MM,190*MM)]:pg.show_pdf_page(clip,old,0,clip=clip)
for k,x,y,w,h in [('e',117,81,65,105),('f',2,190,181,55)]:pg.show_pdf_page(fitz.Rect(x*MM,y*MM,(x+w)*MM,(y+h)*MM),fitz.open(S/'panels'/f'Figure3{k}.pdf'),0)
new.save(S/'review/Figure3.pdf');pg.get_pixmap(dpi=240).save(S/'review/Figure3.png')
base=E.parse(str(P/'图件编辑资料/editable_svg/Figure3.svg'));rt=copy.deepcopy(base)
for k in 'ef':
    target=rt.xpath('//*[@id=$id]',id='adopted_panel_Figure3'+k)[0]
    for child in list(target):target.remove(child)
    part=E.parse(str(S/'panels'/f'Figure3{k}.svg')).getroot();target.set('viewBox',part.get('viewBox'))
    ids={x.get('id'):f'approved_C_{k}_'+x.get('id') for x in part.iter() if x.get('id')}
    for el in part.iter():
        for attr,v in list(el.attrib.items()):
            if attr=='id':el.set(attr,ids[v]);continue
            for before,after in ids.items():v=v.replace(f'url(#{before})',f'url(#{after})');v='#'+after if v=='#'+before else v
            el.set(attr,v)
    for el in part:target.append(copy.deepcopy(el))
for name in meta['a_to_d_svg_subtrees_retained']:
    assert E.tostring(rt.xpath('//*[@id=$id]',id=name)[0])==E.tostring(base.xpath('//*[@id=$id]',id=name)[0])
rt.write(str(S/'review/Figure3.svg'),encoding='utf-8',xml_declaration=True)

# Shared row review for the unchanged d and new e.
pair=fitz.open();page=pair.new_page(width=183*MM,height=105*MM)
for k,x,w,path in [('d',2,114,A/'panels/Figure3d.pdf'),('e',117,65,S/'panels/Figure3e.pdf')]:page.show_pdf_page(fitz.Rect(x*MM,0,(x+w)*MM,105*MM),fitz.open(path),0)
pair.save(S/'review/Figure3de_shared_rows.pdf');page.get_pixmap(dpi=240).save(S/'review/Figure3de_shared_rows.png')
NS='http://www.w3.org/2000/svg';pair_rt=E.Element('{'+NS+'}svg',nsmap={None:NS},width=f'{183*MM}pt',height=f'{105*MM}pt',viewBox=f'0 0 {183*MM} {105*MM}')
for k in 'de':
    el=copy.deepcopy(rt.xpath('//*[@id=$id]',id='adopted_panel_Figure3'+k)[0]);el.set('y','0');pair_rt.append(el)
E.ElementTree(pair_rt).write(str(S/'review/Figure3de_shared_rows.svg'),encoding='utf-8',xml_declaration=True)

# Preserve the full a-d caption and every other document paragraph.
caption=(V/'proposed_caption.md').read_text().strip().replace('Teal indicates','Blue indicates').replace('orange with a downward triangle','pale gold with a downward triangle').replace('\n\n',' ')
caption+=' Fib, fibroblast; EC, endothelial cell; Mac, macrophage; PB, pseudobulk; ABF, approximate Bayes factor; FDR, false discovery rate.'
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'};W='{'+ns['w']+'}'
checks=[]
for rel in ['Manuscript.docx','Nature_Submission/Figure_Legends.docx']:
    p=P/rel
    with ZipFile(p) as z:infos=z.infolist();entries={x.filename:z.read(x.filename) for x in infos}
    root=E.fromstring(entries['word/document.xml']);paras=root.findall('.//w:body/w:p',ns)
    texts=lambda q:''.join(q.xpath('.//w:t/text()',namespaces=ns))
    targets=[q for q in paras if texts(q).startswith('Fig. 3 ')]
    assert len(targets)==1;q=targets[0];oldtext=texts(q);assert ' e, ' in oldtext
    prefix=oldtext.split(' e, ',1)[0]+' '
    before=[E.tostring(v) for v in paras if v is not q]
    nodes=q.findall('.//w:t',ns);normal=None
    for t in nodes:
        rp=t.getparent().find('w:rPr',ns)
        if len(t.text or '')>80 and (rp is None or rp.find('w:b',ns) is None):normal=copy.deepcopy(rp);break
    n=0
    for t in nodes:
        text=t.text or '';keep=max(0,min(len(text),len(prefix)-n));t.text=text[:keep];n+=len(text)
    run=E.SubElement(q,W+'r')
    if normal is not None:run.append(normal)
    node=E.SubElement(run,W+'t');node.text=caption;node.set('{http://www.w3.org/XML/1998/namespace}space','preserve')
    assert texts(q)==prefix+caption
    assert [E.tostring(v) for v in paras if v is not q]==before
    entries['word/document.xml']=E.tostring(root,encoding='UTF-8',xml_declaration=True,standalone=True)
    dst=S/'documents'/p.name
    with ZipFile(dst,'w') as z:
        for info in infos:z.writestr(info,entries[info.filename])
    (S/'review'/f'{p.stem}_Figure3_caption.md').write_text(prefix+caption+'\n')
    checks.append({'path':str(p),'a_to_d_caption_preserved':True,'other_paragraphs_preserved':True,'baseline_sha256':sha(p)})
protected={str(p):sha(p) for b in ['独立小图_最新版/Figure3','图件编辑资料/可视化库重绘_当前/Figure3'] for p in (P/b).glob('*') if p.is_file() and p.stem in ['Figure3a','Figure3b','Figure3c','Figure3d']}
report={'status':'C_APPROVED_STAGED','palette':'C_navy_gold','source_hashes':meta['source_hashes'],'baseline_formal':meta['current_formal_files_untouched'],'protected_panel_hashes':protected,'a_to_d_svg_subtrees_preserved':meta['a_to_d_svg_subtrees_retained'],'caption_checks':checks,'e_displayed_cells':87,'e_constant_columns_in_caption':[29,29],'f_records':10,'scientific_analysis_rerun':False,'biorender_updated':False}
(R/'adoption_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
print(json.dumps({'stage':str(S),'palette':'C','a_to_d_retained':True,'captions_staged':len(checks)},ensure_ascii=False))
