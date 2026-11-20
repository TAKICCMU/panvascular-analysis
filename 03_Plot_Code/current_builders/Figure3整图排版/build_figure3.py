"""Adopt approved 3d/e/f and retain the original formal 3a/c objects and PDF content."""
from pathlib import Path
import copy, hashlib, json, re, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from lxml import etree as ET
import fitz

P=Path('/Users/taki/Desktop/投稿包_最新版_20260919')
R=Path(__file__).resolve().parent
S=P/'Source_Data_Figures'
for sub in ['panels','upload_ready','review','source_data','backup_before_adoption']:
    (R/sub).mkdir(parents=True,exist_ok=True)
MM=72/25.4; NS='http://www.w3.org/2000/svg'; Q=lambda s:'{'+NS+'}'+s
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
bp=R/'backup_before_adoption'
to_backup=[P/'Figures/Figure3.pdf',P/'图件编辑资料/editable_svg/Figure3.svg',
    P/'图件采用版本.tsv',P/'文件来源与校验.tsv',P/'独立小图_最新版/逐panel来源与修改清单.tsv',
    P/'图件编辑资料/可视化库重绘_当前/逐小图_来源与修改.tsv',
    P/'Public_Release/Figure_Code/frozen_svg/Figure3.svg']
for fld in ['独立小图_最新版/Figure3','图件编辑资料/可视化库重绘_当前/Figure3']:
    to_backup += list((P/fld).glob('Figure3*.*'))
for p in to_backup:
    dest=bp/p.relative_to(P)
    if p.exists() and not dest.exists():dest.parent.mkdir(exist_ok=True,parents=True);shutil.copy2(p,dest)
original_svg=bp/'图件编辑资料/editable_svg/Figure3.svg'
original_pdf=bp/'Figures/Figure3.pdf'
source_paths={'d':S/'Figure3/Fig3d_e_all29.tsv','e':P/'图件编辑资料/可视化库重绘_当前/data_used/Figure3e_145_frozen_cells.tsv',
    'f':S/'Current_Display_20260920/source_data/Figure3f_same5_contexts_both_scenarios.tsv'}
read=lambda p:pd.read_csv(p,sep='\t',float_precision='round_trip')
D,E,F=[read(source_paths[k]) for k in ['d','e','f']]
assert len(D)==29 and len(E)==145 and len(F)==10
handoff=P/'图件编辑资料/Figure3def_data_handoff_20261115.md'
assert all(sha(v) in handoff.read_text() for v in source_paths.values())
for k,p in source_paths.items():shutil.copy2(p,R/'source_data'/p.name)
shutil.copy2(handoff,R/'source_data/owning_chat_handoff.md')
C=json.loads((P/'图件编辑资料/semantic_colors_current.json').read_text())
INK='#26333D';MUTED='#65727B';NAVY='#3B6F91';TEAL=C['evidence']['support'];ORANGE=C['evidence']['inverse'];GUIDE='#E8EDF0'
plt.rcParams.update({'font.family':'Arial','font.size':6.3,'axes.labelsize':6.3,'xtick.labelsize':6,
    'ytick.labelsize':6.3,'svg.fonttype':'none','pdf.fonttype':42,'text.color':INK,'axes.labelcolor':INK,
    'axes.edgecolor':MUTED,'xtick.color':MUTED,'ytick.color':INK,'axes.linewidth':.5,'axes.unicode_minus':False})
def fig(w,h):return plt.figure(figsize=(w/25.4,h/25.4))
def ax(f,x,y,w,h):
    W,H=f.get_size_inches()*25.4;return f.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def tx(f,x,y,t,**kw):
    W,H=f.get_size_inches()*25.4;return f.text(x/W,1-y/H,t,**kw)
def clean(a):
    a.spines[['top','right','left']].set_visible(False);a.tick_params(axis='y',length=0);a.tick_params(axis='x',length=2,pad=2,width=.5)
def fonts(root):
    for e in root.iter():
        st=e.get('style','');m=re.search(r'font:\s*([^;]+)',st)
        if not m:continue
        v=re.match(r'(?:(italic|oblique|normal)\s+)?(?:(\d{3}|bold|normal)\s+)?([\d.]+)px\s+(.+)',m.group(1))
        if not v:continue
        fs,fw,size,fam=v.groups();s=f'font-size: {size}px; font-family: {fam}'
        if fs:s+=f'; font-style: {fs}'
        if fw:s+=f'; font-weight: {fw}'
        e.set('style',st[:m.start()]+s+st[m.end():])
def save(f,l):
    f.canvas.draw();rr=f.canvas.get_renderer();W,H=f.get_size_inches()*f.dpi
    for t in f.findobj(matplotlib.text.Text):
        if not t.get_visible() or not t.get_text():continue
        b=t.get_window_extent(rr);assert b.x0>=-1 and b.y0>=-1 and b.x1<=W+1 and b.y1<=H+1,(l,t.get_text(),b)
    for ext in ['pdf','svg']:f.savefig(R/'panels'/f'Figure3{l}.{ext}')
    rt=ET.parse(str(R/'panels'/f'Figure3{l}.svg'));fonts(rt.getroot());rt.write(str(R/'panels'/f'Figure3{l}.svg'),encoding='utf-8',xml_declaration=True)
    f.savefig(R/'panels'/f'Figure3{l}.png',dpi=400);plt.close(f)
    assert not fitz.open(R/'panels'/f'Figure3{l}.pdf')[0].get_images()

top=12;hh=85;pitch=hh/29
fd=fig(114,105);a=ax(fd,38,top,73,hh);a.set(xlim=(-3,2.4),ylim=(28.5,-.5));a.set_xticks([-2,-1,0,1,2]);a.set_yticks([])
for z in [-2,-1,0,1,2]:a.axvline(z,color='#A9B4BC' if z==0 else GUIDE,lw=.65 if z==0 else .4,zorder=0)
for i,r in D.iterrows():
    if i%2==0:a.axhspan(i-.5,i+.5,color='#FAFBFC',zorder=-2)
    y=top+(i+.5)*pitch
    tx(fd,1,y,r.gene_symbol,fontstyle='italic',va='center',fontsize=6.4)
    tx(fd,29.8,y,{'Fibroblast':'Fib','Macrophage':'Mac','Endothelial cell':'EC'}[r.cell_state],ha='center',va='center',color=C['L1'][r.cell_state],fontsize=6.1)
    ln,=a.plot([r.meta_CI_low,r.meta_CI_high],[i,i],c=NAVY,lw=.8);ln.set_gid(f'CI_{i}_{r.gene_symbol}')
    pt=a.scatter(r.meta_logFC,i,s=11,c=NAVY,ec='white',lw=.3,zorder=3);pt.set_gid(f'effect_{i}_{r.gene_symbol}')
    assert np.array_equal(ln.get_xdata(),[r.meta_CI_low,r.meta_CI_high]) and float(pt.get_offsets()[0,0])==r.meta_logFC
clean(a);a.set_xlabel('Meta log2 fold change (95% CI)',labelpad=3)
tx(fd,.5,3.4,'d',weight='bold',fontsize=8);tx(fd,1,9,'Gene',weight='bold');tx(fd,29.8,9,'Cell',ha='center',weight='bold')
save(fd,'d')

fe=fig(65,105);a=ax(fe,1,top,63,hh);a.set(xlim=(-.5,4.5),ylim=(28.5,-.5));a.set_xticks([]);a.set_yticks([]);a.spines[:].set_visible(False)
for i in range(29):
    if i%2==0:a.axhspan(i-.5,i+.5,color='#FAFBFC',zorder=-2)
for _,v in E.iterrows():
    if v.source_value==1:artist,=a.plot(v.col,v.row,marker='_',ms=3,c='#A6B1B9',mew=.8,ls='')
    else:artist=a.scatter(v.col,v.row,s=13 if v.source_value==3 else 16,marker='o' if v.source_value==3 else 'v',color=TEAL if v.source_value==3 else ORANGE,ec='none')
    artist.set_gid(f'evidence_r{v.row}_c{v.col}_code{v.source_value}')
for j,lab in enumerate(['Primary\nMAGMA','Primary\nPB','Leave-CAD\nMAGMA','Regional\nABF','Spatial']):
    tx(fe,1+(j+.5)*63/5,10,lab,ha='center',va='bottom',fontsize=5.9,linespacing=1.1)
tx(fe,.5,3.4,'e',weight='bold',fontsize=8)
fe.legend(handles=[Line2D([],[],marker='o',c=TEAL,ls='',ms=3.4,label='Support'),Line2D([],[],marker='v',c=ORANGE,ls='',ms=3.4,label='Inverse'),Line2D([],[],marker='_',c='#A6B1B9',ls='',ms=3.4,label='No support')],
    ncol=3,frameon=False,loc='lower center',bbox_to_anchor=(.5,.005),fontsize=5.7,columnspacing=.8,handletextpad=.35,handlelength=.7)
save(fe,'e')

ff=fig(181,55);topf=8;hf=35;a=ax(ff,51,topf,72,hf);a.set(xlim=(0,1),ylim=(9.5,-.5));a.set_yticks([]);a.set_xticks([0,.25,.5,.75,1]);clean(a);a.set_xlabel('Posterior probability',labelpad=2)
pairs=F[['gene','context']].drop_duplicates().reset_index(drop=True)
contexts={'Heart_Left_Ventricle':'Heart left ventricle','Whole_Blood':'Whole blood','Heart_Atrial_Appendage':'Heart atrial appendage','Cells_Cultured_fibroblasts':'Cultured fibroblasts'}
status={'strong_colocalisation':'Strong support','suggestive_colocalisation':'Suggestive support','no_regional_GWAS_signal':'GWAS criterion not met'}
for i,r in pairs.iterrows():
    rows=F[F.gene.eq(r.gene)&F.context.eq(r.context)].set_index('scenario');assert rows.n_overlap.nunique()==1
    y=topf+(2*i+1)*hf/10
    tx(ff,1,y-1.0,r.gene,fontsize=6.8,fontstyle='italic',weight='bold')
    tx(ff,1,y+1.7,contexts[r.context],fontsize=5.9,color=MUTED)
    tx(ff,172,y,str(int(rows.iloc[0].n_overlap)),ha='center',va='center',fontsize=6.4)
    for j,sc in enumerate(['primary','leave_one_CAD']):
        q=rows.loc[sc];k=2*i+j;vals=[q['PP.H4'],q['PP.H3'],sum(q[f'PP.H{n}'] for n in range(3))];left=0
        for tag,v,co in zip(['H4','H3','H0_H1_H2'],vals,[NAVY,'#97BDB9','#E0E5E8']):
            rect=a.barh(k,v,left=left,height=.70,color=co,ec='white',lw=.15)[0];rect.set_gid(f'{r.gene}_{sc}_{tag}');assert np.isclose(rect.get_width(),v,rtol=0,atol=1e-15);left+=v
        assert abs(left-1)<1e-12
        a.text(vals[0]/2,k,f'{vals[0]:.4f}',ha='center',va='center',fontsize=5.8,color='white')
        yy=topf+(k+.5)*hf/10
        tx(ff,49,yy,'Primary' if sc=='primary' else 'Leave-CAD',ha='right',va='center',fontsize=5.7,color=MUTED)
        tx(ff,128,yy,status[q.classification],va='center',fontsize=5.9,color=TEAL if q.classification=='strong_colocalisation' else MUTED,
           weight='bold' if q.classification=='strong_colocalisation' else 'normal')
tx(ff,.5,3.4,'f',weight='bold',fontsize=8);tx(ff,128,5.5,'Classification',fontsize=6.1,weight='bold');tx(ff,172,6,'Shared SNPs',ha='center',fontsize=6.1,weight='bold')
ff.legend(handles=[Patch(fc=NAVY,label='H4: shared variant'),Patch(fc='#97BDB9',label='H3: distinct variants'),Patch(fc='#E0E5E8',label='H0–H2: other hypotheses')],ncol=3,frameon=False,
    loc='lower center',bbox_to_anchor=(.52,-.015),fontsize=5.8,columnspacing=1.5,handletextpad=.4,handlelength=1)
save(ff,'f')

# Formal a–c are exported from the actual PDF, at exactly original physical scale.
old=fitz.open(original_pdf);clips={'a':fitz.Rect(0,0,258.5,138.7),'b':fitz.Rect(258.5,0,old[0].rect.width,138.7),'c':fitz.Rect(0,138.7,old[0].rect.width,223.2)}
placements=[]
for l,clip in clips.items():
    doc=fitz.open();pg=doc.new_page(width=clip.width,height=clip.height);pg.show_pdf_page(pg.rect,old,0,clip=clip)
    doc.save(R/'panels'/f'Figure3{l}.pdf');(R/'panels'/f'Figure3{l}.svg').write_text(pg.get_svg_image(text_as_path=True))
    pg.get_pixmap(matrix=fitz.Matrix(3,3)).save(R/'panels'/f'Figure3{l}.png')
    placements.append(dict(panel=f'Figure3{l}',x_mm=clip.x0/MM,y_mm=clip.y0/MM,width_mm=clip.width/MM,height_mm=clip.height/MM,source='original formal PDF, unchanged scale/content'))
for l,x,y,w,h in [('d',2,81,114,105),('e',117,81,65,105),('f',2,190,181,55)]:placements.append(dict(panel=f'Figure3{l}',x_mm=x,y_mm=y,width_mm=w,height_mm=h,source='approved real-data redraw, redundant titles/notes removed'))
for l in 'abcdef':shutil.copy2(R/'panels'/f'Figure3{l}.svg',R/'upload_ready'/f'Figure3{l}.svg')

# Assemble vector PDF: top region comes directly from original, without redraw.
doc=fitz.open();page=doc.new_page(width=183*MM,height=247*MM);clip=fitz.Rect(0,0,old[0].rect.width,223.2)
page.show_pdf_page(clip,old,0,clip=clip)
for r in placements[3:]:
    src=fitz.open(R/'panels'/f'{r["panel"]}.pdf');rect=fitz.Rect(r['x_mm']*MM,r['y_mm']*MM,(r['x_mm']+r['width_mm'])*MM,(r['y_mm']+r['height_mm'])*MM)
    page.show_pdf_page(rect,src,0)
doc.save(R/'review/Figure3.pdf');page.get_pixmap(matrix=fitz.Matrix(3,3)).save(R/'review/Figure3.png')

# Assemble SVG retaining exact original a–c object subtrees and coordinate transforms.
rt=ET.parse(str(original_svg)).getroot();orig=copy.deepcopy(rt)
abc_ids=['identity_a_figure_1','figure3_b_content','legend_2','c_magma','c_pseudobulk','text_82','text_83','text_103']
subhash=lambda root,name:hashlib.sha256(ET.tostring(root.xpath('//*[@id=$v]',v=name)[0])).hexdigest()
abc_before={k:subhash(orig,k) for k in abc_ids}
frozen=copy.deepcopy(orig.xpath('//*[@id="FigureYa_vector_matrix_Figure3e"]')[0])
for name in ['axes_6','figure3d_compact_forest','figure3e_three_symbol_key','figure3f_three_evidence_columns','text_106','text_107','text_108']:
    for el in rt.xpath('//*[@id=$v]',v=name):el.getparent().remove(el)
rt.set('width',f'{183*MM}pt');rt.set('height',f'{247*MM}pt');rt.set('viewBox',f'0 0 {183*MM} {247*MM}')
# Keep frozen scientific evidence attributes in a non-rendered defs record for provenance and existing readers.
defs=ET.SubElement(rt,Q('defs'),id='preserved_frozen_evidence_data');defs.append(frozen)
for r in placements[3:]:
    rr=ET.parse(str(R/'panels'/f'{r["panel"]}.svg')).getroot();ids={e.get('id'):'adopted_'+r['panel']+'_'+e.get('id') for e in rr.iter() if e.get('id')}
    for e in rr.iter():
        for k,v in list(e.attrib.items()):
            if k=='id':e.set(k,ids[v]);continue
            for before,after in ids.items():v=v.replace(f'url(#{before})',f'url(#{after})');v=('#'+after) if v=='#'+before else v
            e.set(k,v)
    g=ET.SubElement(rt,Q('svg'),id='adopted_panel_'+r['panel'],x=str(r['x_mm']*MM),y=str(r['y_mm']*MM),width=str(r['width_mm']*MM),height=str(r['height_mm']*MM),viewBox=rr.get('viewBox'),overflow='hidden')
    for el in rr:g.append(copy.deepcopy(el))
assert all(subhash(rt,k)==abc_before[k] for k in abc_ids)
ET.ElementTree(rt).write(str(R/'review/Figure3.svg'),encoding='utf-8',xml_declaration=True)
spans=lambda pg:[(s['text'],s['bbox'],s['size']) for bl in pg.get_text('dict')['blocks'] for ln in bl.get('lines',[]) for s in ln['spans'] if s['bbox'][3]<=223.2]
before_spans,after_spans=spans(old[0]),spans(page)
assert [s[0] for s in before_spans]==[s[0] for s in after_spans]
max_text_delta=max(abs(u-v) for s1,s2 in zip(before_spans,after_spans) for u,v in zip(s1[1],s2[1]))
assert max_text_delta<.0001
source_hashes={k:{'path':str(p),'sha256':sha(p)} for k,p in source_paths.items()}
metadata={'page_width_mm':183,'page_height_mm':247,'panels':placements,'abc_svg_subtrees_unchanged':abc_before,
    'abc_pdf_original_vector_content_reused':True,'abc_text_position_max_error_pt':max_text_delta,'source_records':{'d_effects_and_CI':29,'e_frozen_cells':145,'f_records':10},'source_hashes':source_hashes,
    'scientific_analysis_rerun':False,'redundant_titles_and_notes_removed':True,'status':'PREPARED_FOR_BIORENDER_IMPORT'}
(R/'panel_placement.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2))
print(json.dumps({'root':str(R),'panels':len(placements),'abc_original_vectors_retained':True,'size_mm':[183,247]},ensure_ascii=False))
