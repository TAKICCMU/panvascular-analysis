"""Install the approved Figure 3 presentation, with scoped backups and checks."""
from pathlib import Path
from zipfile import ZipFile
from lxml import etree as E
import hashlib, json, shutil, csv, io, copy
import fitz

P=Path('/Users/taki/Desktop/投稿包_最新版_20260919')
R=Path(__file__).resolve().parent
B=R/'backup_before_adoption'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
changed=[]
def backup(p):
    q=B/p.relative_to(P)
    if p.exists() and not q.exists():
        q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
def install(src,dst):
    backup(dst);dst.parent.mkdir(parents=True,exist_ok=True)
    tmp=dst.with_name(dst.name+'.figure3_tmp');shutil.copy2(src,tmp);tmp.replace(dst);changed.append(dst)
def write(p,content):
    backup(p);tmp=p.with_name(p.name+'.figure3_tmp');tmp.write_text(content);tmp.replace(p);changed.append(p)
def table(p,fn):
    rows=list(csv.DictReader(io.StringIO(p.read_text()),delimiter='\t'));keys=list(rows[0])
    for row in rows:fn(row)
    out=io.StringIO();w=csv.DictWriter(out,fieldnames=keys,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows);write(p,out.getvalue())

meta=json.loads((R/'panel_placement.json').read_text())
for item in meta['source_hashes'].values():assert sha(item['path'])==item['sha256']
protected=list((P/'Figures').glob('*.pdf'))+list((P/'图件编辑资料/editable_svg').glob('*.svg'))
protected=[p for p in protected if p.stem!='Figure3']
for base in ['独立小图_最新版/Figure3','图件编辑资料/可视化库重绘_当前/Figure3']:
    protected += [p for p in (P/base).glob('Figure3*.*') if p.stem in ['Figure3a','Figure3b','Figure3c']]
protected_hashes={str(p):sha(p) for p in protected}

for src,dst in [('review/Figure3.pdf','Figures/Figure3.pdf'),('review/Figure3.svg','图件编辑资料/editable_svg/Figure3.svg'),('review/Figure3.svg','Public_Release/Figure_Code/frozen_svg/Figure3.svg'),('review/Figure3.pdf','图件编辑资料/可视化库重绘_当前/review/Figure3_整图审阅.pdf'),('review/Figure3.png','图件编辑资料/可视化库重绘_当前/review/Figure3_整图审阅.png')]:install(R/src,P/dst)
for k in 'def':
    for base in ['独立小图_最新版/Figure3','图件编辑资料/可视化库重绘_当前/Figure3']:
        for ext in ['pdf','svg','png']:install(R/'panels'/f'Figure3{k}.{ext}',P/base/f'Figure3{k}.{ext}')

tail=('d, Pseudobulk effects and 95% confidence intervals for all 29 Tier 1 gene–cell-state pairs. '
      'e, Evidence aligned across the same rows for primary MAGMA gene-level Bonferroni significance, primary pseudobulk analysis, MAGMA gene-level FDR significance excluding CAD, regional ABF in the primary scenario, and primary conditional spatial abundance association using q05 across all sections. Circles indicate support, downward triangles inverse association, and grey dashes no support under the original criterion. FDR<0.05 excluding CAD is stricter than the nominal-P eligibility criterion; grey does not remove a candidate. '
      'f, Five genes with strong ABF support in at least one scenario, using the context with the highest posterior among contexts meeting the strong-support criteria and showing both scenarios in that context. Paired stacked bars show PP.H4, PP.H3 and the sum of PP.H0–PP.H2; PP.H4 labels are rounded to four decimal places. Adjacent classifications retain the original posterior and regional signal criteria. Shared SNP counts are overlapping GWAS–QTL variants, not sample sizes. High PP.H4 alone does not establish strong support or mediation by a cell state. Leave-CAD means omitting the CAD phenotype, not excluding participants. '
      'Fib, fibroblast; EC, endothelial cell; Mac, macrophage; PB, pseudobulk; ABF, approximate Bayes factor; FDR, false discovery rate.')
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
W='{'+ns['w']+'}'
caption_checks=[]
for rel in ['Manuscript.docx','Nature_Submission/Figure_Legends.docx']:
    p=P/rel;backup(p)
    with ZipFile(p) as z:
        infos=z.infolist();entries={x.filename:z.read(x.filename) for x in infos}
    root=E.fromstring(entries['word/document.xml']);paras=root.findall('.//w:body/w:p',ns)
    texts=lambda q:''.join(q.xpath('.//w:t/text()',namespaces=ns))
    targets=[q for q in paras if texts(q).startswith('Fig. 3 ')]
    assert len(targets)==1
    q=targets[0];old=texts(q);prefix=old.split(' d, ',1)[0]+' '
    other_before=[E.tostring(v) for v in paras if v is not q]
    nodes=q.findall('.//w:t',ns);normal=None
    for t in nodes:
        rp=t.getparent().find('w:rPr',ns)
        if len(t.text or '')>80 and (rp is None or rp.find('w:b',ns) is None):normal=copy.deepcopy(rp);break
    n=0
    for t in nodes:
        text=t.text or '';keep=max(0,min(len(text),len(prefix)-n));t.text=text[:keep];n+=len(text)
    run=E.SubElement(q,W+'r')
    if normal is not None:run.append(normal)
    node=E.SubElement(run,W+'t');node.text=tail;node.set('{http://www.w3.org/XML/1998/namespace}space','preserve')
    assert texts(q)==prefix+tail
    assert [E.tostring(v) for v in paras if v is not q]==other_before
    entries['word/document.xml']=E.tostring(root,encoding='UTF-8',xml_declaration=True,standalone=True)
    tmp=p.with_name(p.name+'.figure3_tmp')
    with ZipFile(tmp,'w') as z:
        for info in infos:z.writestr(info,entries[info.filename])
    tmp.replace(p);changed.append(p)
    caption_checks.append({'path':str(p),'abc_caption_retained':True,'other_paragraphs_unchanged':True})
    (R/'review'/('Figure3_caption_manuscript.md' if rel.startswith('Manuscript') else 'Figure3_caption_legends.md')).write_text(prefix+tail+'\n')

pdf=P/'Figures/Figure3.pdf';svg=P/'图件编辑资料/editable_svg/Figure3.svg'
def figure_row(row):
    if row['figure']!='Figure3':return
    row.update(source_pdf=str(pdf),source_svg=str(svg),version='Approved Figure 3d–f redesign 20261115; formal a–c retained',status='APPROVED_SCHEME_ADOPTED_BIORENDER_SAVED',width_mm='183',height_mm='247',sha256=sha(pdf),svg_sha256=sha(svg),package_pdf=str(pdf),package_svg=str(svg),package_pdf_sha256=sha(pdf),package_svg_sha256=sha(svg),provenance_source_pdf=str(R/'review/Figure3.pdf'),provenance_source_svg=str(R/'review/Figure3.svg'))
table(P/'图件采用版本.tsv',figure_row)
changes={'d':'Approved compact forest; all 29 effects and stored 95% CI endpoints retained; no redundant title/footer',
         'e':'Approved 145-cell symbol matrix sharing d row labels and spacing; original categories retained; no redundant title/footer',
         'f':'Approved paired posterior compositions for the same five contexts and two scenarios; original H0–H4, SNP counts and classifications retained; no redundant title/footer'}
dim={r['panel']:(r['width_mm'],r['height_mm']) for r in meta['panels']}
def panels(row):
    if row.get('figure')!='Figure3' or row['panel'] not in 'def':return
    k=row['panel'];w,h=dim['Figure3'+k];base=P/'独立小图_最新版/Figure3'
    pg=fitz.open(base/f'Figure3{k}.pdf')[0]
    sizes=[s['size'] for b in pg.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans']]
    row.update(method='approved_redraw_existing_data',width_mm=str(w),height_mm=str(h),font_p10_pt=str(round(sorted(sizes)[len(sizes)//10],2)),source_files=meta['source_hashes'][k]['path'],source_objects='29 rows' if k=='d' else '145 frozen cells' if k=='e' else '10 original gene-context-scenario records',change=changes[k],pdf_sha256=sha(base/f'Figure3{k}.pdf'),svg_sha256=sha(base/f'Figure3{k}.svg'),data_integrity='source hashes checked; stored values and original classifications retained; no scientific analysis rerun',layout_legibility='d/e aligned shared rows; essential legends retained',visual_review='approved scheme adopted; local render and BioRender reload checked')
table(P/'独立小图_最新版/逐panel来源与修改清单.tsv',panels)
def current(row):
    if row['panel'] not in ['Figure3d','Figure3e','Figure3f']:return
    k=row['panel'][-1];w,h=dim[row['panel']]
    row.update(operation='approved_redraw_existing_data',source_files=meta['source_hashes'][k]['path'],change=changes[k],width_mm=str(w),height_mm=str(h))
table(P/'图件编辑资料/可视化库重绘_当前/逐小图_来源与修改.tsv',current)

for p,h in protected_hashes.items():assert sha(p)==h,p
for item in meta['source_hashes'].values():assert sha(item['path'])==item['sha256']
audit={'changed_files':[str(p) for p in changed],'protected_files_unchanged':protected_hashes,'source_data_hashes_unchanged':meta['source_hashes'],'caption_checks':caption_checks,'new_figure3_sha256':sha(pdf),'new_figure3_svg_sha256':sha(svg),'abc_original_svg_subtrees_preserved':meta['abc_svg_subtrees_unchanged'],'scientific_analysis_rerun':False}
(R/'review/adoption_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2))
print(json.dumps({'changed':len(changed),'protected_unchanged':len(protected_hashes),'caption_scope':'Figure 3 d/e/f only','audit':str(R/'review/adoption_audit.json')},ensure_ascii=False))
