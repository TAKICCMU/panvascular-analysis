"""Figure2 display preparation only; preserves adopted sources and scientific results.

Exact quantitative plots are never redrawn by generative image tools. BioRender
MCP draft provides a composition reference; the compact workflow is an exact
editable SVG. Standalone assets deliberately contain NO panel letters.
"""
from pathlib import Path
import os, copy, io, json, hashlib, shutil, base64
os.environ['DYLD_FALLBACK_LIBRARY_PATH']='/opt/homebrew/lib'
from lxml import etree as E
import numpy as np
import pandas as pd
from PIL import Image
import cairosvg
import pymupdf as fitz
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap

P=Path('/Users/taki/Desktop/投稿包_最新版_20260919')
R=Path(__file__).resolve().parent
O=R/'panels'; O.mkdir(exist_ok=True)
S=P/'Source_Data_Figures/Figure2'
NS='http://www.w3.org/2000/svg'; Q=lambda x:'{'+NS+'}'+x
MM=72/25.4
INK='#20252A'; BLUE='#345A7D'; TEAL='#168E8A'; GOLD='#B79150'
records=[]

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def export(root,name,source,change):
    p=O/name
    E.ElementTree(root).write(str(p)+'.svg',encoding='UTF-8',xml_declaration=True)
    cairosvg.svg2pdf(bytestring=E.tostring(root),write_to=str(p)+'.pdf')
    doc=fitz.open(str(p)+'.pdf'); page=doc[0]
    page.get_pixmap(matrix=fitz.Matrix(4,4)).save(str(p)+'.png')
    records.append(dict(panel=name,source=[str(x) for x in source],source_sha256={str(x):sha(x) for x in source if Path(str(x)).is_file()},change=change,width_mm=page.rect.width/MM,height_mm=page.rect.height/MM,panel_letter_baked=False,svg_sha256=sha(str(p)+'.svg')))

def native_panel(l):
    source=P/'独立小图_最新版/Figure2'/f'Figure2{l}.svg'
    root=E.parse(str(source)).getroot()
    # Remove only the single outside panel letter, never in-chart labels.
    found=[x for x in root.iter(Q('text')) if ''.join(x.itertext()).strip()==l]
    assert len(found)==1,(l,len(found))
    found[0].getparent().remove(found[0])
    if l=='b':
        # Cairo does not honour SVG pixelated for this 8x8 image. Preserve its
        # literal RGBA bytes as 64 vector rectangles instead of smoothing cells.
        for im in list(root.iter(Q('image'))):
            if im.get('width')=='8' and im.get('height')=='8':
                arr=np.asarray(Image.open(io.BytesIO(base64.b64decode(im.get('{http://www.w3.org/1999/xlink}href').split(',')[1]))).convert('RGBA'))
                assert arr.shape==(8,8,4)
                g=E.Element(Q('g'),transform=im.get('transform'))
                for i in range(8):
                    for j in range(8):
                        r,gc,b,a=map(int,arr[i,j]);E.SubElement(g,Q('rect'),x=str(j),y=str(i),width='1',height='1',fill=f'#{r:02x}{gc:02x}{b:02x}',**{'fill-opacity':str(a/255),'stroke':'none'})
                im.getparent().replace(im,g)
    b=cairosvg.svg2png(bytestring=E.tostring(root),scale=2,background_color='white')
    a=np.asarray(Image.open(io.BytesIO(b)).convert('RGB')); yy,xx=np.where(a.min(axis=2)<248)
    v=list(map(float,root.get('viewBox').split()))
    x0=v[0]+max(0,xx.min()-5)/a.shape[1]*v[2];y0=v[1]+max(0,yy.min()-5)/a.shape[0]*v[3]
    x1=v[0]+min(a.shape[1],xx.max()+6)/a.shape[1]*v[2];y1=v[1]+min(a.shape[0],yy.max()+6)/a.shape[0]*v[3]
    root.set('viewBox',f'{x0} {y0} {x1-x0} {y1-y0}');root.set('width',f'{x1-x0}pt');root.set('height',f'{y1-y0}pt')
    export(root,'Figure2'+l,[source],'Outside letter removed; visible content cropped with margin; all native data geometry retained')

def compact_workflow():
    root=E.Element(Q('svg'),nsmap={None:NS},width='180mm',height='34mm',viewBox='0 0 510.236220 96.377953')
    def rect(x,y,w,h,fill,stroke='#91A5B5',rx=1.6):
        return E.SubElement(root,Q('rect'),x=str(x*MM),y=str(y*MM),width=str(w*MM),height=str(h*MM),rx=str(rx*MM),fill=fill,stroke=stroke,**{'stroke-width':'.55'})
    def text(x,y,s,size=6.2,bold=False,anchor='start',fill=INK):
        el=E.SubElement(root,Q('text'),x=str(x*MM),y=str(y*MM),fill=fill,**{'font-family':'Arial','font-size':str(size),'font-weight':'700' if bold else '400','text-anchor':anchor});el.text=s
    def arrow(x1,y1,x2,y2,color='#56616B'):
        E.SubElement(root,Q('path'),d=f'M{x1*MM},{y1*MM} L{x2*MM},{y2*MM}',fill='none',stroke=color,**{'stroke-width':'.65'})
        if y1==y2:pts=[(x2,y2),(x2-1.1,y2-.65),(x2-1.1,y2+.65)]
        else:pts=[(x2,y2),(x2-.65,y2-1.1),(x2+.65,y2-1.1)]
        E.SubElement(root,Q('polygon'),points=' '.join(f'{x*MM},{y*MM}' for x,y in pts),fill=color)
    for x,w,col in [(0,48,'#EEF3F7'),(51,32,'#EEF3F7'),(86,34,'#F4F6F8'),(123,57,'#EAF4F3')]:rect(x,.25,w,33.5,col)
    for x1,x2 in [(48,51),(83,86),(120,123)]:arrow(x1,17,x2,17)
    text(24,4,'8 vascular GWAS',7.5,True,'middle')
    text(4,7.3,'Phenotype',5.8);text(44,7.3,'N: cases + controls',5.8,anchor='end')
    names=['AD','AAA','TAA','CAD','Stroke','PAD','CA','CAVS'];ns=[464256,443400,19646,547261,1296908,449548,48434,956682]
    for i,(name,n) in enumerate(zip(names,ns)):
        y=10.5+i*2.95;text(4,y,name);text(44,y,f'{n:,}',anchor='end')
    text(67,4,'LDSC +',7.2,True,'middle');text(67,7.3,'Genomic SEM',7.2,True,'middle')
    rect(57,10.4,20,5.3,'#C9DAE8');text(67,14,'Shared G',6.5,True,'middle')
    rect(57,17.5,20,5.3,'#E0E8F1');text(67,21.1,'Aortic',6.3,anchor='middle')
    rect(53.5,24.6,27,5.3,'#E0E8F1');text(67,28.2,'Atherosclerotic',6.3,anchor='middle')
    text(67,32.3,'Orthogonal bifactor',5.7,anchor='middle')
    text(103,4,'Shared factor GWAS',7.1,True,'middle')
    text(103,9,'Regional Q-SNP',6.2,anchor='middle');text(103,12,'exclusion guided by LD',6.2,anchor='middle')
    arrow(103,13.2,103,16.2)
    rect(89,17,28,5.5,'#E5ECF3');text(103,20.8,'FUMA + MAGMA',6.6,True,'middle')
    text(103,27.1,'155 primary loci',6.6,True,'middle');text(103,31.1,'67 without CAD',6.2,anchor='middle')
    text(151.5,4,'Evidence integration',7.3,True,'middle')
    rect(125,6.5,22,7.5,'#DBE7F0');text(136,9.6,'Genetic',6.2,anchor='middle');text(136,12.3,'evidence',6.2,anchor='middle')
    rect(150,6.5,28,7.5,'#F5E9D2',stroke='#B59F73');text(164,9.6,'Pseudobulk',6.2,anchor='middle');text(164,12.3,'meta-analysis',6.2,anchor='middle')
    # Undirected convergence branches to the evaluable input set; not a causal diagram.
    E.SubElement(root,Q('path'),d=f'M{136*MM},{14*MM} L{136*MM},{16*MM} L{151.5*MM},{16*MM} L{164*MM},{16*MM} L{164*MM},{14*MM}',fill='none',stroke='#637581',**{'stroke-width':'.6'})
    arrow(151.5,16,151.5,18)
    text(151.5,20.4,'1,730 genes · 15,639 evaluable pairs',6.1,True,'middle')
    text(151.5,24.4,'Five joint criteria',6.2,anchor='middle')
    arrow(151.5,25,151.5,27)
    rect(129,27.5,45,5.2,'#C5DDE0',stroke='#8AAFB4');text(151.5,31.2,'29 prioritised pairs',7.3,True,'middle')
    source=P/'图件编辑资料/BioRender_当前素材/GWAS_8_phenotypes_source.tsv'
    export(root,'Figure2a',[source,P/'Manuscript.docx'],'MCP-informed compact one-row workflow; exact editable SVG text, 34 mm tall; thresholds retained in existing legend, no new scientific calculation')

def robustness():
    src=S/'Figure2_Fig3e_leave_one_axis_correlations.tsv';d=pd.read_csv(src,sep='\t').set_index('analysis')
    order=['full','leave_one_AA','leave_one_AAA','leave_one_TAAD','leave_one_CAD','leave_one_Stroke','leave_one_PAD','leave_one_CA','leave_one_CAVS']
    d=d.loc[order].copy();d['display_label']=['Full','−AD','−AAA','−TAA','−CAD','−Stroke','−PAD','−CA','−CAVS']
    d['CI_low']=d.STD_Genotype-1.96*d.STD_Genotype_SE;d['CI_high']=d.STD_Genotype+1.96*d.STD_Genotype_SE
    d.to_csv(R/'source_data/Figure2g_existing_effects_and_CI.tsv',sep='\t')
    plt.rcParams.update({'font.family':'Arial','font.size':6.2,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.55,'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK})
    f=plt.figure(figsize=(65/25.4,46/25.4));a=f.add_axes([13/65,10/46,49/65,32/46])
    vals=d.STD_Genotype.to_numpy();se=d.STD_Genotype_SE.to_numpy();y=np.arange(9)
    a.axvline(vals[0],lw=.55,ls=':',c='#A7AFB7',zorder=0)
    for i in range(9):
        c=BLUE if i in [0,4] else '#747C86'
        a.errorbar(vals[i],i,xerr=1.96*se[i],fmt='D' if i==0 else 'o',ms=3 if i in [0,4] else 2.6,lw=.8,color=c,capsize=1.5)
    a.set(ylim=(8.55,-.55),xlim=(0,.51));a.set_yticks(y,d.display_label);a.set_xticks([0,.2,.4]);a.tick_params(axis='y',length=0,pad=3);a.tick_params(axis='x',length=2,pad=2)
    a.set_xlabel('Factor correlation (95% CI)',labelpad=3)
    a.spines['left'].set_visible(False)
    f.savefig(O/'Figure2g.svg');plt.close(f)
    root=E.parse(str(O/'Figure2g.svg')).getroot()
    production=P/'图件编辑资料/全流程代码审核_20260922/Panvascular_code_review/analysis_code/retained_genetics/Genetics/full_pipeline/03_correlated_two_factor_leave_one_CFA.R'
    export(root,'Figure2g',[src,production],'Proposed main-panel reuse of ED7c: same nine estimates, original Wald 95% CI definition, no refitting; full plus all eight omissions')

if __name__=='__main__':
    for l in 'bcdef': native_panel(l)
    compact_workflow();robustness()
    (R/'panel_manifest.json').write_text(json.dumps(records,ensure_ascii=False,indent=2))
