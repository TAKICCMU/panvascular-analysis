"""Display-only Figure 2 alternatives from existing results; no model fitting."""
from pathlib import Path
import csv, hashlib, json, subprocess, textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.font_manager import FontProperties
import fitz

P = Path('/Users/taki/Desktop/投稿包_最新版_20260919')
R = Path(__file__).resolve().parent
OLD = Path('/Users/taki/Desktop/Codex任务1/Figure2_genetics_redesign_v3/source_data')
for d in ['panels','source_data','review']: (R/d).mkdir(parents=True,exist_ok=True)
SOURCES = {
 'loadings':P/'Source_Data_Figures/Figure2/Figure2_Fig3d_bifactor_loadings.tsv',
 'ecv':OLD/'server_export/CFA_02__bifactor_ECV.tsv',
 'workbook':P/'Supplementary_Tables.xlsx',
 'manuscript':P/'Manuscript.docx',
 'credible_sets':OLD/'local_existing/credible_sets_with_acceptance.tsv',
 'current_regions':P/'Source_Data_Figures/Figure2/all_region_summary_with_acceptance.tsv',
 'retained_regions':OLD/'local_existing/all_region_summary_with_acceptance.tsv',
 'share_formula':P/'图件编辑资料/全流程代码审核_20260922/Panvascular_code_review/analysis_code/retained_genetics/Genetics/full_pipeline/02_three_model_CFA_comparison.R',
}
hashes={k:hashlib.sha256(f.read_bytes()).hexdigest() for k,f in SOURCES.items()}
assert SOURCES['current_regions'].read_bytes()==SOURCES['retained_regions'].read_bytes()

# Read the final workbook directly, keeping all five data records and original order.
extract="""import json,openpyxl,sys
w=openpyxl.load_workbook(sys.argv[1],read_only=True,data_only=True)
s=w['Table 10']; h=[c.value for c in s[3]]
rows=[dict(zip(h,r)) for r in s.iter_rows(min_row=4,values_only=True) if isinstance(r[0],str) and r[0].startswith('rs')]
print(json.dumps(rows,ensure_ascii=False))
"""
PY='/Users/taki/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3'
novel=pd.DataFrame(json.loads(subprocess.run([PY,'-c',extract,str(SOURCES['workbook'])],check=True,capture_output=True,text=True).stdout))
assert novel['Lead variant'].tolist()==['rs7525267','rs75135846','rs78051210','rs333929','rs7236320']
novel.to_csv(R/'source_data/option_B_final_Table10_rows.tsv',sep='\t',index=False)

load=pd.read_csv(SOURCES['loadings'],sep='\t')
order=['AA','AAA','TAAD','CAD','Stroke','PAD','CA','CAVS']
names={'AA':'AD','TAAD':'TAA'}
general=load[load.lhs.eq('G')].set_index('rhs')['STD_Genotype'].loc[order]
specific=load[~load.lhs.eq('G')].set_index('rhs')['STD_Genotype'].loc[order]
share=general**2/(general**2+specific**2)
ecv=float(pd.read_csv(SOURCES['ecv'],sep='\t').ECV_general.iloc[0])
assert abs(ecv-(general**2).sum()/((general**2).sum()+(specific**2).sum()))<1e-12
pd.DataFrame({'trait':order,'display_label':[names.get(x,x) for x in order],
 'general_loading':general.values,'specific_loading':specific.values,
 'general_common_variance_share':share.values}).to_csv(R/'source_data/option_A_display_shares.tsv',sep='\t',index=False)

cs=pd.read_csv(SOURCES['credible_sets'],sep='\t')
cs=cs[cs.scenario.eq('primary') & cs.method.eq('finite502_EB_z_only') & cs.credible_set_analysis_ready.eq(True)].copy()
assert len(cs)==156 and cs.region_id.nunique()==136
assert not cs.duplicated(['scenario','region_id','method','cs_id']).any()
cs.to_csv(R/'source_data/option_C_all_accepted_primary_sets.tsv',sep='\t',index=False)

INK='#222B35'; MUTED='#697480'; BLUE='#355A80'; TEAL='#278B88'; GREY='#DCE2E8'; GRID='#E8EDF1'
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':7,'xtick.labelsize':6.5,
 'ytick.labelsize':7,'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,
 'axes.edgecolor':'#68737D','axes.linewidth':.55,'svg.fonttype':'none','pdf.fonttype':42,
 'savefig.facecolor':'white','svg.hashsalt':'Figure2_existing_results_options'})
mm=1/25.4
W,H=90,72
def fig(): return plt.figure(figsize=(W*mm,H*mm),facecolor='white')
def axes(f,x,y,w,h): return f.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def text(f,x,y,s,**kw): return f.text(x/W,1-y/H,s,va=kw.pop('va','top'),**kw)
def clean(a,left=False):
 for sp in ['top','right']: a.spines[sp].set_visible(False)
 if not left: a.spines['left'].set_visible(False)
 a.tick_params(axis='y',length=0)
 a.tick_params(axis='x',length=2.5,width=.5,pad=2)
def save(f,name):
 for ext in ['svg','pdf','png']:
  f.savefig(R/'panels'/f'{name}.{ext}',dpi=300)
 plt.close(f)

# A: the same loading-square ratio used by the original CFA display code.
# It is NOT total heritability, individual disease risk or a causal fraction.
f=fig()
text(f,1,1,'Shared factor',fontsize=7.3,fontweight='bold')
text(f,89,1.5,f'{100*ecv:.1f}%',fontsize=13,fontweight='bold',ha='right',color=BLUE)
text(f,1,7.5,'of modelled common variance across all eight traits',fontsize=6.5,color=MUTED)
a=axes(f,16,20,63,42)
y=np.arange(8)
a.barh(y,100*share.values,color=BLUE,height=.62,zorder=3)
a.barh(y,100*(1-share.values),left=100*share.values,color=GREY,height=.62,zorder=2)
a.set_yticks(y,[names.get(t,t) for t in order]);a.set_ylim(7.7,-.7)
a.set_xlim(0,100);a.set_xticks([0,25,50,75,100]);a.set_xlabel('Share of common variance (%)',labelpad=3)
a.grid(axis='x',color='white',lw=.6,zorder=4);clean(a)
for yi,v in zip(y,100*share.values):
 a.text(102,yi,f'{v:.0f}%',va='center',fontsize=6.2,color=BLUE,clip_on=False)
for x,c,s in [(16,BLUE,'Shared G'),(47,GREY,'Domain factor')]:
 f.patches.append(Rectangle((x/W,1-16.4/H),2.6/W,2.6/H,transform=f.transFigure,fc=c,ec='none'))
 text(f,x+3.6,14,s,fontsize=6.7)
save(f,'Option_A_shared_variance')

# B: exact final Table10 P and PIP; all five shortlisted loci, genomic order.
f=fig()
text(f,1,1,'Lead variant / regional genes',fontsize=6.5,color=MUTED)
text(f,51,1,'Shared GWAS',fontsize=6.8,fontweight='bold',ha='center')
text(f,77.5,1,'Variant PIP',fontsize=6.8,fontweight='bold',ha='center')
genes=['C1orf21 / COLGALT2','LINC01105 / SOX11','EPB41L2','MAP2K6 region','RAB27B']
axp=axes(f,41,10,20,48); axi=axes(f,68,10,20,48)
yp=np.arange(5)
vals=-np.log10(novel['P'].astype(float).to_numpy());pip=novel['Primary posterior inclusion probability (PIP)'].astype(float).to_numpy()
for a in [axp,axi]:
 a.set_ylim(4.6,-.6);a.set_yticks([]);clean(a)
 for j in range(5):
  if j%2==0:a.axhspan(j-.5,j+.5,color='#F6F8FA',zorder=0,lw=0)
axp.set_xlim(0,10);axp.set_xticks([0,5,10]);axp.set_xlabel('−log10(P)',labelpad=3)
axp.axvline(-np.log10(5e-8),ls=(0,(2,2)),lw=.6,color=MUTED,zorder=2)
axp.barh(yp,vals,color=BLUE,height=.30,zorder=3)
axi.set_xlim(-.02,1.04);axi.set_xticks([0,.5,1]);axi.set_xlabel('PIP',labelpad=3)
axi.hlines(yp,0,pip,color='#AED5D3',lw=1.5,zorder=2)
axi.scatter(pip,yp,s=18,c=TEAL,edgecolors='none',zorder=3)
for j,(rs,g) in enumerate(zip(novel['Lead variant'],genes)):
 yy=10+(j+.6)/5.2*48
 text(f,1,yy-2.5,rs,fontsize=7,fontweight='bold',va='center')
 text(f,1,yy+.8,g,fontsize=6.2,color=MUTED,va='center')
 text(f,64.7,yy,f'{pip[j]:.3f}',fontsize=6.3,ha='center',va='center',color=TEAL)
# Dashed line is a threshold key, not a new panel title or result.
f.add_artist(Line2D([.455,.49],[.035,.035],transform=f.transFigure,ls=(0,(2,2)),color=MUTED,lw=.6))
text(f,45,69.1,r'$P = 5\times10^{-8}$',fontsize=6.2,va='center')
save(f,'Option_B_five_candidate_loci')

# C: each point is an already accepted credible set; exact coordinates, no jitter.
f=fig()
text(f,1,1,'156 accepted credible sets',fontsize=8,fontweight='bold')
text(f,1,7,'across 136 of 155 primary regions',fontsize=6.7,color=MUTED)
a=axes(f,13,18,74,42)
a.scatter(cs.cs_size,cs.lead_PIP,s=14,color=BLUE,alpha=.65,edgecolors='none')
a.set_xscale('log');a.set_xlim(.8,300);a.set_ylim(-.035,1.04)
a.set_xticks([1,5,10,50,100,250],[1,5,10,50,100,250]);a.minorticks_off()
a.set_yticks([0,.25,.5,.75,1]);a.set_xlabel('Variants per credible set (log scale)',labelpad=3)
a.set_ylabel('Lead variant PIP',labelpad=3)
a.grid(color=GRID,lw=.45);a.set_axisbelow(True);clean(a,left=True)
text(f,13,69.2,'One point = one accepted credible set',fontsize=6.3,color=MUTED,va='center')
save(f,'Option_C_finemapping_resolution')

# Review contact sheet only. Separate panel files have no baked panel letters.
cn=FontProperties(fname='/Library/Fonts/Arial Unicode.ttf')
CW,CH=294,117
f=plt.figure(figsize=(CW*mm,CH*mm),facecolor='white')
def ct(x,y,s,**kw): return f.text(x/CW,1-y/CH,s,fontproperties=cn,va='top',**kw)
ct(5,3,'Figure 2：三个可替代敏感性图的内容方案',fontsize=12,color=INK)
ct(5,10.5,'均使用已有真实结果；仅生成审阅样板，当前 BioRender Figure 2 保持不变。',fontsize=7,color=MUTED)
options=[
 ('A','共享遗传成分有多大？','Option_A_shared_variance',
  '突出正文40.9%的结果。比例分母是模型共同方差，\n不是总体疾病风险或全遗传力；仅显示点估计。'),
 ('B','发现了哪些候选位点？  优先推荐','Option_B_five_candidate_loci',
  '关联强度与原PIP分开显示；“新关联”限于原检索范围。\n附近基因不是已确认效应基因，低PIP位点同样保留。'),
 ('C','精细定位的分辨率如何？','Option_C_finemapping_resolution',
  '保留全部156个获接受可信集合，不挑高PIP结果。\n散点展示定位分辨率，不将PIP称为因果证明。'),
]
for i,(key,title,name,note) in enumerate(options):
 x=5+i*98
 ct(x,18,f'方案{key}  {title}',fontsize=8.2,color=INK)
 with fitz.open(R/'panels'/f'{name}.pdf') as d:
  pix=d[0].get_pixmap(matrix=fitz.Matrix(3,3),alpha=False)
  im=np.frombuffer(pix.samples,dtype=np.uint8).reshape(pix.height,pix.width,3)
 ax=f.add_axes([x/CW,1-((26+72)/CH),90/CW,72/CH]);ax.imshow(im);ax.axis('off')
 ct(x,102,note,fontsize=6.6,color=MUTED,linespacing=1.65)
for ext in ['png','pdf']:
 f.savefig(R/'review'/f'Figure2_options_comparison.{ext}',dpi=220)
plt.close(f)

manifest={'scope':'display-only alternatives; no model fit or new significance test',
 'current_BioRender_unchanged':True,'source_hashes':{k:{'path':str(SOURCES[k]),'sha256':v} for k,v in hashes.items()},
 'A':{'formula':'lambda_G^2/(lambda_G^2+lambda_domain^2); original CFA code lines 95-101',
      'global_ECV_stored':ecv,'display_transform_only':True,'not_total_heritability':True},
 'B':{'worksheet':'Table 10','rows':[4,5,6,7,8],'all_five_preserved':True,'regional_display_gene_names':'shortened labels follow the current manuscript; full original Table10 fields retained in source_data','novelty_scope':'historical component-disease search, not refreshed online or global novelty'},
 'C':{'method':'finite502_EB_z_only','scenario':'primary','filter':'credible_set_analysis_ready == True',
      'sets':156,'regions':136,'region_summary_byte_identical_to_current_package':True},
 'panel_size_mm':[90,72],'letters_baked_into_panel_files':False,
 'protected_source_hashes_unchanged':all(hashlib.sha256(f.read_bytes()).hexdigest()==hashes[k] for k,f in SOURCES.items())}
(R/'source_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
(R/'README.md').write_text('''# Figure 2 补充内容样板

本目录仅供选择展示方案，尚未替换 BioRender 中当前的 g，也未修改正式图、正文或源表。

- A：共享/领域成分占模型共同方差的比例。40.9%读取既有ECV文件；逐表型比例按照原CFA代码的平方载荷比例做展示转换。不是重新拟合，也不代表总体疾病风险或全部遗传力。未计算新的CI/P值。
- B（优先推荐）：当前Supplementary Table10的全部5个候选位点，展示原P和原PIP。这5个位点与29个gene–state候选不是同一集合。保留最低PIP的rs7236320；不因视觉效果筛选。候选新关联仅在原文限定的历史检索范围内，未做新的在线新颖性审计；附近基因不等于已确认效应基因。显示基因名按正文简写，完整区域基因字段保留在source_data。
- C：既有主分析、z-only方法、通过原接受标志的156个可信集合；逐个绘出cs_size和lead_PIP。沿用原接受规则，不是把所有155区域都称为有可信集合。

独立PDF/SVG位于panels；三方案比较位于review；真实作图数据位于source_data。每个单图为90×72mm，不包含正式panel字母。B需要足够横向宽度，不宜直接压入原61.5mm的g小槽。

本轮未运行Genomic SEM、SuSiE、GWAS、LD计算或新显著性检验。source_manifest.json记录原文件身份和显示转换。
''')
print(json.dumps({'output':str(R),'ECV':ecv,'novel_variants':novel['Lead variant'].tolist(),'accepted_sets':len(cs),'sources_unchanged':manifest['protected_source_hashes_unchanged']},ensure_ascii=False,indent=2))
