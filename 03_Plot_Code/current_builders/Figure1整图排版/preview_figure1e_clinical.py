"""Review prototype: unchanged L2 study coverage + recorded age and sex.

First run extract_figure1e_clinical.py (openpyxl). This drawing uses only its
stored Table 02 extraction and the existing state summaries. No models run.
"""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pymupdf
from PIL import Image

O=Path(__file__).parent;P=O.parents[2];R=O/'review';MM=72/25.4
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
sources=[P/'Source_Data_Figures/ExtendedData3_4_5/ED345_all_named_state_index.tsv',
 P/'Source_Data_Figures/Figure1/L2_cross_study_recurrence.tsv',O/'L1_display_palette.tsv',
 P/'Supplementary_Tables.xlsx']
before={str(p):sha(p) for p in sources}
protected=list((O/'panels').glob('Figure1*.*'))+[O/'Figure1.pdf',O/'Figure1.svg',O/'panel_placement.json']
protected_hash={str(p):sha(p) for p in protected if p.is_file()}
audit=json.loads((R/'Figure1e_clinical_metadata_check.json').read_text())
assert audit['source_sha256']==sha(sources[3])
donors=pd.read_csv(R/'Figure1e_clinical_donor_records.tsv',sep='\t',dtype={'donor_id':str})
assert len(donors)==172 and donors.donor_id.nunique()==172
assert donors.age_years.notna().sum()==audit['age_reported']

named=pd.read_csv(sources[0],sep='\t');coverage=pd.read_csv(sources[1],sep='\t')
coverage[['level1_annotation','level2_annotation']]=coverage.state.str.split(' | ',regex=False,expand=True)
keys=['level1_annotation','level2_annotation']
for table in (named,coverage):
 for key in keys:table[key]=table[key].str.strip()
z=named.merge(coverage,on=keys,how='left',validate='one_to_one',suffixes=('_named',''))
assert len(z)==116 and z.studies.notna().all() and z.n_cells.eq(z.n_cells_named).all()
palette=pd.read_csv(sources[2],sep='\t').sort_values('order')
order=[k for k in palette.L1 if k in set(z.level1_annotation)]
colors=palette.set_index('L1').hex.to_dict()
display={'Endothelial cell':'Endothelial','Fibroblast':'Fibroblast','Smooth muscle cell':'Smooth muscle',
 'Fibromyocyte':'Fibromyocyte','Pericyte':'Pericyte','Macrophage':'Macrophage','Monocyte':'Monocyte',
 'Neutrophil':'Neutrophil','Dendritic cell':'Dendritic','T cell':'T cell','Natural killer cell':'Natural killer',
 'B cell':'B cell','Plasma cell':'Plasma cell','Mast cell':'Mast cell'}
mpl.rcParams.update({'font.family':'Arial','font.size':7,'text.color':'#222831',
 'axes.labelcolor':'#222831','axes.edgecolor':'#87939D','axes.linewidth':.5,
 'xtick.color':'#48515A','ytick.color':'#48515A','pdf.fonttype':42,'svg.fonttype':'none',
 'savefig.facecolor':'white'})
W,H=183.,55.;f=plt.figure(figsize=(W/25.4,H/25.4))
def axes(box):
 x,y,w,h=box
 return f.add_axes([x/W,(H-y-h)/H,w/W,h/H])
left_box=[10.5,8.5,117.,25.5]
a=axes(left_box)
a.set_xlim(-.5,13.5);a.set_ylim(5,27.5);a.set_yticks([5,10,15,20,26])
a.set_ylabel('Studies per L2 state',fontsize=7,labelpad=3)
a.set_xticks(np.arange(14),[display[v] for v in order],rotation=55,ha='right',rotation_mode='anchor')
a.tick_params(axis='x',length=0,pad=3.5,labelsize=7)
a.tick_params(axis='y',length=2,pad=2,labelsize=6.7)
a.spines[['top','right']].set_visible(False);a.spines['bottom'].set_visible(False)
a.grid(axis='y',color='#E5E9EC',lw=.45,zorder=0)
dot_diameter=1.05;dot_step=1.12;category_mm=left_box[2]/14;records=[]
for i,parent in enumerate(order):
 if i%2==0:a.axvspan(i-.5,i+.5,fc='#F7F8F9',ec='none',zorder=-1)
 g=z[z.level1_annotation.eq(parent)]
 for count,tied in g.groupby('studies',sort=True):
  tied=tied.sort_values('level2_annotation')
  xx=i+(np.arange(len(tied))-(len(tied)-1)/2)*dot_step/category_mm
  a.scatter(xx,np.full(len(tied),count),s=(dot_diameter*MM)**2,c=colors[parent],edgecolors='none',linewidths=0,zorder=3)
  for (_,row),x in zip(tied.iterrows(),xx):
   r=row.to_dict();r.update(display_label=display[parent],plot_x_category=float(x),plot_y_studies=int(count),display_hex=colors[parent]);records.append(r)
 a.plot([i-.34,i+.34],[5.15,5.15],color=colors[parent],lw=1.35,solid_capstyle='butt',clip_on=True)

# Two small clinical panels. One observation per recorded donor ID in both.
clinical_colors={'age':'#718997','male':'#647F92','female':'#AD8D75','Not reported':'#D9DDE1'}
age=donors.age_years.dropna().to_numpy()
bins=[0,20,40,60,80,100]
age_counts,_=np.histogram(age,bins=bins)
assert age_counts.sum()==len(age)==111
ab=axes([144.,10.0,37.,14.5])
ab.bar(np.arange(5),age_counts,width=.74,color=clinical_colors['age'],edgecolor='none',zorder=2)
ab.set_xlim(-.6,4.6);ab.set_ylim(0,60)
ab.set_xticks(np.arange(5),['<20','20–39','40–59','60–79','≥80'])
ab.set_yticks([0,20,40])
ab.tick_params(axis='both',length=2,pad=1.5,labelsize=6.2)
ab.spines[['top','right']].set_visible(False)
ab.set_ylabel('Donor IDs',fontsize=6.3,labelpad=2)
ab.grid(axis='y',color='#E5E9EC',lw=.4,zorder=0)
for x,n in enumerate(age_counts):ab.text(x,float(n)+1.5,str(int(n)),ha='center',va='bottom',fontsize=6.5)

sex_names=['male','female','Not reported']
sex_counts=donors.sex.fillna('Not reported').value_counts().reindex(sex_names,fill_value=0).astype(int)
assert sex_counts.sum()==172
sb=axes([155.2,37.,24.,13.2])
sb.barh(np.arange(3),sex_counts,height=.5,color=[clinical_colors[k] for k in sex_names],edgecolor='none')
sb.set_xlim(0,100);sb.set_ylim(2.65,-.65)
sb.set_yticks(np.arange(3),['Male','Female','Not reported'])
sb.set_xticks([]);sb.tick_params(axis='y',length=0,pad=2,labelsize=6.5)
sb.spines[['top','bottom','left','right']].set_visible(False)
for y,n in enumerate(sex_counts):sb.text(float(n)+2,y,str(int(n)),va='center',ha='left',fontsize=6.7)

o=axes([0,0,W,H]);o.set(xlim=(0,W),ylim=(H,0));o.set_axis_off()
o.text(1.,3.,'L2 states',fontsize=6.7,va='center')
for i,parent in enumerate(order):
 x=left_box[0]+(i+.5)*category_mm
 o.text(x,3.,str(int(z.level1_annotation.eq(parent).sum())),fontsize=7.4,fontweight='bold',ha='center',va='center')
o.scatter([74.],[51.6],s=(dot_diameter*MM)**2,c='#61717C',edgecolors='none')
o.text(76.,51.6,'One dot = one named L2 state',ha='left',va='center',fontsize=6.7)
o.text(158.5,2.9,'172 recorded donor IDs',fontsize=7,ha='center',va='center',fontweight='bold')
o.text(162.5,7.4,'Age (years)',fontsize=6.7,ha='center',va='center')
o.text(158.5,30.2,'Age reported: 111; missing: 61',fontsize=6.3,ha='center',va='center')
o.text(158.5,34.1,'Sex (donor IDs)',fontsize=6.7,ha='center',va='center')

f.canvas.draw();renderer=f.canvas.get_renderer();outside=[]
for axis in f.axes:
 texts=list(axis.texts)
 if axis.axison:texts += [axis.xaxis.label,axis.yaxis.label]+axis.get_xticklabels()+axis.get_yticklabels()
 for t in texts:
  if not t.get_visible() or not t.get_text():continue
  bb=t.get_window_extent(renderer).transformed(f.transFigure.inverted())
  if min(bb.x0,bb.y0)<0 or max(bb.x1,bb.y1)>1:outside.append(t.get_text())
assert not outside,outside
prefix=R/'Figure1e_L2_and_clinical_prototype'
for ext in ['pdf','svg']:f.savefig(prefix.with_suffix('.'+ext))
f.savefig(prefix.with_suffix('.png'),dpi=300);plt.close(f)

# Exact left-side preservation compared with the preceding accepted prototype.
prev=Image.open(R/'Figure1e_compact_coverage_and_size.png').convert('RGB')
curr=Image.open(prefix.with_suffix('.png')).convert('RGB')
assert prev.size==curr.size
crop=(0,0,int(132/25.4*300),prev.height)
unchanged=np.array_equal(np.asarray(prev.crop(crop)),np.asarray(curr.crop(crop)))
assert unchanged,'Approved left graphic changed'
summary=[]
for label,lo,hi,n in zip(['<20','20–39','40–59','60–79','≥80'],bins[:-1],bins[1:],age_counts):
 summary.append({'variable':'Age','category':label,'count':int(n),'denominator':172,'reported_denominator':111,'rule':f'{lo} <= age < {hi}'})
summary.append({'variable':'Age','category':'Not reported','count':61,'denominator':172,'reported_denominator':111,'rule':'Missing stored age'})
for label,n in sex_counts.items():summary.append({'variable':'Sex','category':label,'count':int(n),'denominator':172,'reported_denominator':134,'rule':'Stored category; missing shown separately'})
pd.DataFrame(summary).to_csv(R/'Figure1e_clinical_plot_summary.tsv',sep='\t',index=False)
pd.DataFrame(records).to_csv(R/'Figure1e_L2_clinical_left_data.tsv',sep='\t',index=False)
with pymupdf.open(prefix.with_suffix('.pdf')) as d:
 assert len(d)==1
 page=d[0];page.get_pixmap(dpi=96).save(R/'Figure1e_L2_and_clinical_native_96dpi.png')
 spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans'] if s['text'].strip()]
 assert 'Cells per L2 state' not in page.get_text()
 font_range=[min(s['size'] for s in spans),max(s['size'] for s in spans)]
assert all(sha(Path(k))==v for k,v in before.items())
assert all(sha(Path(k))==v for k,v in protected_hash.items())
check={'status':'REVIEW_ONLY_NOT_APPLIED_TO_BIORENDER_OR_FORMAL_FIGURE',
 'size_mm':[W,H],'left_132mm_pixel_identical':unchanged,'named_states':116,
 'recorded_donor_ids':172,'age_reported':111,'age_missing':61,'sex_counts':sex_counts.to_dict(),
 'minimum_font_pt':min(font_range),'maximum_font_pt':max(font_range),
 'no_embedded_panel_letter':True,'text_outside_page':outside,'source_sha256':before,
 'source_and_adopted_figures_unchanged':True,'clinical_colors':clinical_colors}
(R/'Figure1e_L2_clinical_layout_check.json').write_text(json.dumps(check,indent=2))
note='''# Figure1e：L2覆盖与临床信息样板

左侧沿用已认可的116个L2状态研究覆盖图；132 mm左侧渲染与上一样板逐像素一致。
右侧以最新版Supplementary_Tables.xlsx的Table 02为准，247个样本标识按172个既有Donor ID去重。
同一Donor ID的已报告年龄、性别没有冲突。不按细胞数加权，不把多份样本当作多名供者。

年龄111个ID有记录，61个缺失；性别男性94、女性40、未报告38。年龄按20岁区间展示，所有已报告值均保留。
这些是原表临床信息的描述性汇总，不是新增独立核验。GSE131778的7个用户确认独立供者ID仍沿用Table 02已有年龄/性别值，本轮没有重新确认其临床来源。

未绘制疾病/健康供者比例：此前指定的“健康分析组”包含部分患病供者的参照组织，不能直接改称健康供者。
本轮不改任何来源字段、注释、科学结果、正文或正式图。图件为183×55 mm审阅样板，没有内嵌panel字母。

待采用英文图注：Study coverage of the 116 named L2 states, grouped by parent lineage. Each dot denotes one state and vertical position indicates the recorded number of contributing studies. Numbers above the groups indicate the number of named states. Right, recorded age and sex from Supplementary Table 2, summarised once per donor identifier across 247 sample identifiers (172 recorded donor IDs). Age is available for 111 IDs and missing for 61; sex is recorded as male for 94 IDs, female for 40 and not reported for 38. Repeated sample records with the same donor ID are counted once. The clinical panels describe the stored metadata rather than a new verification of demographic information or donor independence.
'''
(R/'Figure1e_clinical_notes.md').write_text(note)
print(json.dumps({k:v for k,v in check.items() if k!='source_sha256'},indent=2))
