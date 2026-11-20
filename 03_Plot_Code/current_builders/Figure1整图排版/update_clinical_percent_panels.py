"""Current review panels: f clinical percentages, g L2 study coverage.

No source statistics or current BioRender canvas are changed.  Lowercase
panel letters remain separate layout objects.  Uses Table 02 donor records
already extracted by extract_figure1e_clinical.py; no clinical imputation.
"""
from pathlib import Path
import copy,csv,hashlib,json,shutil
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from lxml import etree as E
import cairosvg
import pymupdf

O=Path(__file__).parent;P=O.parents[2];R=O/'review';MM=72/25.4
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
source=P/'Supplementary_Tables.xlsx'
source_hash=sha(source)
meta=json.loads((R/'Figure1e_clinical_metadata_check.json').read_text())
assert meta['source_sha256']==source_hash
d=pd.read_csv(R/'Figure1e_clinical_donor_records.tsv',sep='\t')
assert len(d)==172 and d.donor_id.nunique()==172
protected=list((O/'panels').glob('*'))+[O/'Figure1.pdf',O/'Figure1.svg',P/'Manuscript.docx',source]
protected_hash={str(p):sha(p) for p in protected if p.is_file()}
backup=P/'图件编辑资料/轻量备份/Figure1临床百分比与编号调整前_20261115'
for f in [R/'Figure1e_L2_and_clinical_prototype.pdf',R/'Figure1e_L2_and_clinical_prototype.svg',
          R/'Figure1e_L2_and_clinical_prototype.png',R/'Figure1e_CURRENT_PROTOTYPE.md',
          O/'panel_placement.json']:
    target=backup/f.relative_to(O)
    if f.exists() and not target.exists():target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,target)

denominator=len(d)
age_bins=np.array([0,20,40,60,80,100])
age_counts=np.r_[np.histogram(d.age_years.dropna(),bins=age_bins)[0],d.age_years.isna().sum()]
sex_counts=d.sex.fillna('Not reported').value_counts().reindex(['male','female','Not reported']).to_numpy()
assert age_counts.sum()==sex_counts.sum()==denominator
age_percent=age_counts/denominator*100
sex_percent=sex_counts/denominator*100
assert np.isclose(age_percent.sum(),100) and np.isclose(sex_percent.sum(),100)
age_labels=['<20','20–39','40–59','60–79','≥80','Not\nreported']
sex_labels=['Male','Female','Not\nreported']
age_colors=['#718997']*5+['#D9DDE1']
sex_colors=['#647F92','#AD8D75','#D9DDE1']
mpl.rcParams.update({'font.family':'Arial','font.size':7.5,'text.color':'#222831',
 'axes.labelcolor':'#222831','axes.edgecolor':'#87939D','axes.linewidth':.5,
 'xtick.color':'#48515A','ytick.color':'#48515A','pdf.fonttype':42,'svg.fonttype':'none',
 'savefig.facecolor':'white'})
W,H=122.,39.
f=plt.figure(figsize=(W/25.4,H/25.4))
def axes(box):
 x,y,w,h=box
 return f.add_axes([x/W,(H-y-h)/H,w/W,h/H])
ax_age=axes([12.,3.,62.,23.])
ax_sex=axes([87.,3.,33.,23.])
for ax,values,labels,colors in [(ax_age,age_percent,age_labels,age_colors),(ax_sex,sex_percent,sex_labels,sex_colors)]:
    x=np.arange(len(values))
    ax.bar(x,values,width=.72,color=colors,edgecolor='none',zorder=3)
    ax.set_ylim(0,62);ax.set_xlim(-.65,len(values)-.35)
    ax.set_yticks([0,20,40,60]);ax.yaxis.set_major_formatter(PercentFormatter(100,decimals=0))
    ax.set_xticks(x,labels)
    ax.tick_params(axis='x',length=0,pad=3,labelsize=7.5)
    ax.tick_params(axis='y',length=2,pad=2,labelsize=7)
    ax.grid(axis='y',color='#E5E9EC',lw=.45,zorder=0)
    ax.spines[['top','right']].set_visible(False)
    for xx,value in zip(x,values):
        ax.text(xx,float(value)+1.8,f'{value:.1f}%',ha='center',va='bottom',fontsize=7.7)
ax_age.set_ylabel('Donors (%)',fontsize=7.5,labelpad=3)
ax_sex.set_yticklabels([]);ax_sex.tick_params(axis='y',length=0);ax_sex.spines['left'].set_visible(False)
overlay=axes([0,0,W,H]);overlay.set(xlim=(0,W),ylim=(H,0));overlay.set_axis_off()
overlay.text(43.,36.3,'Age (years)',ha='center',va='center',fontsize=7.8)
overlay.text(103.5,36.3,'Sex',ha='center',va='center',fontsize=7.8)
f.canvas.draw();renderer=f.canvas.get_renderer();outside=[]
for ax in f.axes:
    texts=list(ax.texts)
    if ax.axison:texts+=[ax.xaxis.label,ax.yaxis.label]+ax.get_xticklabels()+ax.get_yticklabels()
    for t in texts:
        if not t.get_visible() or not t.get_text():continue
        bb=t.get_window_extent(renderer).transformed(f.transFigure.inverted())
        if min(bb.x0,bb.y0)<0 or max(bb.x1,bb.y1)>1:outside.append(t.get_text())
assert not outside,outside
f_prefix=R/'Figure1f_clinical_percent'
for ext in ['pdf','svg']:f.savefig(f_prefix.with_suffix('.'+ext))
f.savefig(f_prefix.with_suffix('.png'),dpi=300);plt.close(f)

# Separate the approved left SVG without retaining hidden clinical objects.
svg=E.parse(str(R/'Figure1e_L2_and_clinical_prototype.svg'))
ns={'s':'http://www.w3.org/2000/svg'}
root=svg.getroot();figure=root.find('.//s:g[@id="figure_1"]',ns)
for gid in ['axes_2','axes_3']:
    node=figure.find(f's:g[@id="{gid}"]',ns);assert node is not None;figure.remove(node)
ov=figure.find('s:g[@id="axes_4"]',ns)
remove_text={'172 recorded donor IDs','Age (years)','Age reported: 111; missing: 61','Sex (donor IDs)'}
for child in list(ov):
    texts=child.xpath('.//s:text/text()',namespaces=ns)
    if texts and set(texts)<=remove_text:ov.remove(child)
assert not any(t in remove_text for t in root.xpath('.//s:text/text()',namespaces=ns))
root.set('width',f'{132*MM}pt');root.set('height',f'{55*MM}pt');root.set('viewBox',f'0 0 {132*MM} {55*MM}')
background=figure.find('s:g[@id="patch_1"]/s:path',ns)
background.set('d',f'M 0 {55*MM} L {132*MM} {55*MM} L {132*MM} 0 L 0 0 z')
g_prefix=R/'Figure1g_L2_coverage'
svg.write(str(g_prefix.with_suffix('.svg')),xml_declaration=True,encoding='utf-8')
cairosvg.svg2pdf(url=str(g_prefix.with_suffix('.svg')),write_to=str(g_prefix.with_suffix('.pdf')))
cairosvg.svg2png(url=str(g_prefix.with_suffix('.svg')),write_to=str(g_prefix.with_suffix('.png')),dpi=300)

summary=[]
for field,labels,counts,percents in [('age',age_labels,age_counts,age_percent),('sex',sex_labels,sex_counts,sex_percent)]:
    for label,n,pct in zip(labels,counts,percents):summary.append({'variable':field,'category':label.replace('\n',' '),
        'stored_count':int(n),'denominator_recorded_donor_ids':denominator,'percentage_unrounded':float(pct),'display_percentage':f'{pct:.1f}%'})
pd.DataFrame(summary).to_csv(R/'Figure1f_clinical_percent_data.tsv',sep='\t',index=False)

mapping=[('a','a','Atlas construction',O/'panels/Figure1a.pdf'),
 ('b','b','L1 UMAP and Tissue',O/'panels/Figure1b.pdf'),
 ('c','c','Cells per lineage',O/'panels/Figure1c.pdf'),
 ('f','d','Studies per lineage',O/'panels/Figure1f.pdf'),
 ('d','e','Tissue composition',O/'panels/Figure1d.pdf'),
 ('new clinical prototype','f','Recorded donor age and sex (%)',f_prefix.with_suffix('.pdf')),
 ('new L2 coverage prototype','g','Studies per L2 state',g_prefix.with_suffix('.pdf')),
 ('g','h','L1 marker dotplot',O/'panels/Figure1g.pdf')]
pd.DataFrame([{'old_panel':old,'new_panel':new,'content':name,'source_pdf':str(path),
 'embedded_panel_letter':False,'status':'CURRENT_REVIEW_MAPPING_NOT_YET_APPLIED_TO_BIORENDER'} for old,new,name,path in mapping]).to_csv(
 R/'Figure1_current_panel_mapping.tsv',sep='\t',index=False)
for old,new,name,source_pdf in mapping:
    if (old,new) in [('f','d'),('d','e'),('g','h')]:
        for ext in ['pdf','svg','png']:
            shutil.copy2(source_pdf.with_suffix('.'+ext),R/f'Figure1{new}_current.{ext}')

notes='''# 本轮 Figure1 审阅材料

临床图现为Figure1f_clinical_percent；L2覆盖图现为Figure1g_L2_coverage。
原f（Studies per lineage）→d；原d（组织组成）→e；原g（marker dotplot）→h；a–c不变。旧UMAP放大图不再进入拟采用主图。
所有小图不内嵌panel字母。字母在BioRender中作为独立可编辑文字添加。本轮没有导出BioRender，也没有修改正式正文/整图。

临床图去掉总标题、172 IDs标题及重复的缺失说明。Age (years)、Sex和Donors (%)是必要坐标标签。
年龄和性别统一使用172个既有供者ID为分母，未报告单列灰色；只把显示数量转换为百分比。
实际尺寸122×39 mm。年龄有效图区62×23 mm，原37×14.5 mm；普通文字7–7.8 pt。
L2覆盖图为132×55 mm；保留全部116个状态、原研究数、原类别色及排布，没有转成临床百分比。

待采用图注：f, Recorded age and sex, expressed as percentages of 172 donor identifiers represented by 247 sample identifiers in Supplementary Table 2. Each donor ID contributes once. Not reported values remain in the denominator and are shown in grey. Percentages are rounded to one decimal place and may not sum to exactly 100% after rounding. g, Study coverage of 116 named L2 states; each dot represents one state, vertical position gives the number of contributing studies, and numbers above the categories indicate the number of states.

原表临床字段来源限制继续保留在Figure1e_clinical_metadata_check.json；本轮不独立核验年龄/性别，也不改动临床值或供者身份。
'''
(R/'Figure1_current_panel_notes.md').write_text(notes)
(R/'Figure1e_CURRENT_PROTOTYPE.md').write_text('当前审阅入口已拆分为 Figure1f_clinical_percent.pdf 和 Figure1g_L2_coverage.pdf。\n完整编号对应见 Figure1_current_panel_mapping.tsv。\n本轮没有替换正式整图或BioRender。\n')

pdf_checks={}
for key,prefix in [('f',f_prefix),('g',g_prefix)]:
    with pymupdf.open(prefix.with_suffix('.pdf')) as doc:
        p=doc[0];text=p.get_text()
        p.get_pixmap(dpi=96).save(prefix.with_name(prefix.name+'_native_96dpi').with_suffix('.png'))
        spans=[s for b in p.get_text('dict')['blocks'] if 'lines' in b for line in b['lines'] for s in line['spans'] if s['text'].strip()]
        pdf_checks[key]={'dimensions_mm':[round(p.rect.width/MM,2),round(p.rect.height/MM,2)],
            'font_range_pt':[round(min(s['size'] for s in spans),2),round(max(s['size'] for s in spans),2)]}
        if key=='f':assert '172' not in text and '111' not in text and '94' not in text
        if key=='g':assert not any(v in text for v in ['Age (years)','172 recorded','Female','Sex (donor'])
assert all(sha(Path(k))==v for k,v in protected_hash.items())
check={'status':'REVIEW_PANELS_UPDATED_BIORENDER_AND_FORMAL_ASSEMBLY_NOT_CHANGED',
 'denominator':172,'age_percentage_sum':float(age_percent.sum()),'sex_percentage_sum':float(sex_percent.sum()),
 'clinical_percentages':summary,'pdfs':pdf_checks,'source_sha256':source_hash,'source_and_formal_files_unchanged':True,
 'text_outside_page':outside,'no_embedded_panel_letters':True,
 'mapping':{old:new for old,new,_,_ in mapping},
 'l2_svg_preservation':'Original axes_1 and the L2-only elements of axes_4 retained without transform; clinical axes and text removed, not hidden.'}
(R/'Figure1_clinical_percent_and_panel_check.json').write_text(json.dumps(check,indent=2))
print(json.dumps({k:v for k,v in check.items() if k not in ['clinical_percentages','source_sha256']},indent=2))
