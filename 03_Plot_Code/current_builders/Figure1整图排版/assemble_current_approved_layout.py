"""Approved in-place Figure1 layout. Plot stored summaries only; no inference.

Lowercase panel letters are assembly-level objects, never in panel assets.
The local mirror is prepared for the existing BioRender canvas, not an export.
"""
from pathlib import Path
import copy, hashlib, json, shutil
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
from matplotlib.ticker import PercentFormatter
from lxml import etree as E
import pymupdf as fitz

O=Path(__file__).parent; P=O.parents[2]; R=O/'review'; A=R/'current_layout_assets'; A.mkdir(exist_ok=True)
MM=72/25.4
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
backup=P/'图件编辑资料/轻量备份/Figure1临床百分比与编号调整前_20261115/BioRender合入前'
for file in list((O/'panels').glob('Figure1*'))+[O/'Figure1.pdf',O/'Figure1.svg',O/'panel_placement.json']:
    dst=backup/file.relative_to(O)
    if not dst.exists():dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(file,dst)
protected=[P/'Manuscript.docx',P/'Supplementary_Tables.xlsx']
protected_hash={str(p):sha(p) for p in protected}
palette=pd.read_csv(O/'L1_display_palette.tsv',sep='\t').sort_values('order')
pal=palette.set_index('L1').hex.to_dict()
short={k:k.replace(' cell','') for k in pal}
mpl.rcParams.update({'font.family':'Arial','font.size':7,'text.color':'#222831','axes.labelcolor':'#222831',
    'axes.edgecolor':'#87939D','axes.linewidth':.5,'xtick.color':'#48515A','ytick.color':'#48515A',
    'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
def new(w,h):
    f=plt.figure(figsize=(w/25.4,h/25.4)); f._mm=(w,h); return f
def axes(f,box):
    W,H=f._mm; x,y,w,h=box
    return f.add_axes([x/W,(H-y-h)/H,w/W,h/H])
def tx(f,x,y,text,size=7,**kw):
    w,h=f._mm;return f.text(x/w,1-y/h,text,fontsize=size,va='center',**kw)
checks={}
def save(f,name):
    f.canvas.draw(); rr=f.canvas.get_renderer(); out=[]
    for ax in f.axes:
        texts=list(ax.texts)
        if ax.axison:texts += [ax.xaxis.label,ax.yaxis.label]+ax.get_xticklabels()+ax.get_yticklabels()
        for t in texts:
            if not t.get_visible() or not t.get_text():continue
            bb=t.get_window_extent(rr).transformed(f.transFigure.inverted())
            if min(bb.x0,bb.y0)<-.001 or max(bb.x1,bb.y1)>1.001:out.append(t.get_text())
    assert not out,(name,out)
    for ext in ['pdf','svg']:f.savefig(A/f'Figure1{name}.{ext}')
    f.savefig(A/f'Figure1{name}.png',dpi=900)
    checks[name]={'size_mm':list(f._mm),'text_outside':out,'no_panel_letter':True}
    plt.close(f)

# Copy untouched panels from the lightweight pre-relabel backup.
for old,newlabel in [('a','a'),('b','b'),('c','c'),('f','d'),('d','e')]:
    for ext in ['pdf','svg','png']:
        src=backup/f'panels/Figure1{old}.{ext}'
        shutil.copy2(src,A/f'Figure1{newlabel}.{ext}')

# f: use the already computed percentages, excluding only the missing display rows.
clinical=pd.read_csv(R/'Figure1f_clinical_percent_data.tsv',sep='\t')
assert clinical.denominator_recorded_donor_ids.eq(172).all()
f=new(58,49)
age=clinical[(clinical.variable=='age')&(clinical.category!='Not reported')]
sex=clinical[(clinical.variable=='sex')&(clinical.category!='Not reported')]
assert len(age)==5 and len(sex)==2
a=axes(f,[10,4,46,20]); x=np.arange(5)
a.bar(x,age.percentage_unrounded,width=.72,color='#718997',edgecolor='none',zorder=3)
a.set(xlim=(-.6,4.6),ylim=(0,33));a.set_yticks([0,10,20,30]);a.yaxis.set_major_formatter(PercentFormatter(100,decimals=0))
a.set_xticks(x,age.category,fontsize=7)
a.tick_params(length=2,pad=2,labelsize=7);a.tick_params(axis='x',length=0)
a.set_ylabel('Donors (%)',fontsize=7,labelpad=2)
a.spines[['top','right']].set_visible(False);a.grid(axis='y',color='#E5E9EC',lw=.45,zorder=0)
for xx,row in zip(x,age.itertuples()):a.text(xx,row.percentage_unrounded+1,row.display_percentage,ha='center',va='bottom',fontsize=7.3)
tx(f,33,29,'Age (years)',7.4,ha='center')
s=axes(f,[15,33,35,10]); s.barh([0,1],sex.percentage_unrounded,height=.55,color=['#647F92','#AD8D75'],edgecolor='none')
s.set(xlim=(0,60),ylim=(1.6,-.6));s.set_yticks([0,1],sex.category,fontsize=7.3);s.set_xticks([]);s.tick_params(length=0,pad=3)
for sp in s.spines.values():sp.set_visible(False)
for yy,row in enumerate(sex.itertuples()):s.text(row.percentage_unrounded+1,yy,row.display_percentage,va='center',fontsize=7.3)
tx(f,32,47,'Sex · donors (%)',7.4,ha='center')
save(f,'f')
assert 'Not reported' not in fitz.open(A/'Figure1f.pdf')[0].get_text()

# g: every stored L2 study count retained; categorical tie offsets are display only.
z=pd.read_csv(R/'Figure1e_L2_clinical_left_data.tsv',sep='\t')
assert len(z)==116
order=[k for k in palette.L1 if k in set(z.level1_annotation)]
display=z.drop_duplicates('level1_annotation').set_index('level1_annotation').display_label.to_dict()
f=new(122,40);box=[10.5,6.2,108,17.2];a=axes(f,box)
a.set(xlim=(-.5,13.5),ylim=(5,27.5));a.set_yticks([5,10,15,20,26]);a.set_ylabel('Studies per L2 state',fontsize=7,labelpad=3)
a.set_xticks(range(14),[display[k] for k in order],rotation=50,ha='right',rotation_mode='anchor',fontsize=7)
a.tick_params(axis='x',length=0,pad=3);a.tick_params(axis='y',length=2,pad=2,labelsize=6.7)
a.spines[['top','right','bottom']].set_visible(False);a.grid(axis='y',color='#E5E9EC',lw=.45,zorder=0)
for i,parent in enumerate(order):
    if i%2==0:a.axvspan(i-.5,i+.5,fc='#F7F8F9',ec='none',zorder=-1)
    group=z[z.level1_annotation==parent]
    for study,tied in group.groupby('studies',sort=True):
        xs=i+(np.arange(len(tied))-(len(tied)-1)/2)*1.03/(box[2]/14)
        a.scatter(xs,np.full(len(tied),study),s=(.88*MM)**2,c=pal[parent],edgecolors='none',lw=0,zorder=3)
    a.plot([i-.34,i+.34],[5.15,5.15],color=pal[parent],lw=1.35,solid_capstyle='butt')
    tx(f,box[0]+(i+.5)*box[2]/14,2.8,str(len(group)),7.4,ha='center',weight='bold')
# The count-row definition is in the caption; reserve the corner for native g.
save(f,'g')

# h: stored B scaled mean and detection fractions. Layout-only re-export.
src=P/'Source_Data_Figures/Current_Display_20260920/marker_audit_20260921'
V=pd.read_csv(src/'per_gene_scaled_mean.tsv',sep='\t',index_col=0)
F=pd.read_csv(src/'original_detection_fraction.tsv',sep='\t',index_col=0)
assert V.shape==(14,28) and V.index.equals(F.index) and V.columns.equals(F.columns)
cm=LinearSegmentedColormap.from_list('blue_teal_sequential',['#EFF4F5','#BAD4D8','#76A8B4','#367789','#173F5C'])
f=new(122,44);a=axes(f,[22,3,83,28]);x,y=np.meshgrid(range(28),range(14))
smax=(.70*min(83/28,28/14)*MM)**2;valid=V.notna().to_numpy()&F.notna().to_numpy()&(F.to_numpy()>0)
a.scatter(x[valid],y[valid],s=smax*F.to_numpy()[valid],c=V.to_numpy()[valid],cmap=cm,norm=Normalize(0,1),edgecolors='none',lw=0)
a.set(xlim=(-.5,27.5),ylim=(13.5,-.5));a.set_xticks(range(28),V.columns,rotation=90,fontstyle='italic',fontsize=6.2)
a.set_yticks(range(14),[short[k] for k in V.index],fontsize=6.2);a.tick_params(length=0,pad=2);a.tick_params(axis='y',pad=5)
for sp in a.spines.values():sp.set_visible(False)
for j,state in enumerate(V.index):a.add_patch(Rectangle((-.86,j-.37),.12,.74,fc=pal[state],ec='none',clip_on=False))
tx(f,113.4,5,'Per-gene scaled\nmean expression',6,ha='center',linespacing=1.05)
cb=f.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap=cm),cax=axes(f,[109,10,2,12]))
cb.set_ticks([0,.5,1]);cb.ax.tick_params(labelsize=6,length=2,pad=2);cb.outline.set_linewidth(.35)
tx(f,113.5,26,'Detected\ncells (%)',6.1,ha='center',linespacing=1)
la=axes(f,[107,30,14,12]);la.set(xlim=(0,14),ylim=(12,0));la.axis('off')
for j,q in enumerate([.25,.5,.75,1]):la.scatter(2,1+j*3,s=smax*q,c='#617E89',edgecolors='none');la.text(5.6,1+j*3,str(int(q*100)),fontsize=6.1,va='center')
save(f,'h')

positions={'a':[0,0,183,29],'b':[0,31,107,51],'c':[109,31,60,51],'d':[170,31,13,51],
           'e':[0,84,58,35],'f':[0,121,58,49],'g':[61,84,122,40],'h':[61,126,122,44]}
NS='http://www.w3.org/2000/svg';XL='http://www.w3.org/1999/xlink';Q=lambda v:'{'+NS+'}'+v
root=E.Element(Q('svg'),nsmap={None:NS,'xlink':XL},width='183mm',height='170mm',viewBox=f'0 0 {183*MM} {170*MM}')
doc=fitz.open();page=doc.new_page(width=183*MM,height=170*MM)
for label,(x,y,w,h) in positions.items():
    svg=E.parse(str(A/f'Figure1{label}.svg')).getroot(); ids={e.get('id'):f'{label}_{e.get("id")}' for e in svg.iter() if e.get('id')}
    for e in svg.iter():
        for k,v in list(e.attrib.items()):
            if k=='id':e.set(k,ids[v])
            else:
                for old,newid in ids.items():v=v.replace('url(#'+old+')','url(#'+newid+')');v='#'+newid if v=='#'+old else v
                e.set(k,v)
    g=E.SubElement(root,Q('g'),id='panel_'+label,transform=f'translate({x*MM} {y*MM})')
    for child in list(svg):g.append(child)
    with fitz.open(A/f'Figure1{label}.pdf') as d:
        assert np.allclose([d[0].rect.width/MM,d[0].rect.height/MM],[w,h],atol=.02)
        page.show_pdf_page(fitz.Rect(x*MM,y*MM,(x+w)*MM,(y+h)*MM),d,0)
page.insert_font(fontname='PanelArialBold',fontfile='/System/Library/Fonts/Supplemental/Arial Bold.ttf')
for label,(x,y,w,h) in positions.items():
    xx,yy=x+1,y+2.5
    page.insert_text((xx*MM,yy*MM),label,fontname='PanelArialBold',fontsize=8,color=(.05,.05,.05))
    t=E.SubElement(root,Q('text'),x=str(xx*MM),y=str(yy*MM),fill='#111111',**{'font-family':'Arial','font-size':'8','font-weight':'bold','data-role':'separate-panel-letter'});t.text=label
doc.save(R/'Figure1_current_preview.pdf',garbage=4,deflate=True)
page.get_pixmap(dpi=220).save(R/'Figure1_current_preview.png')
page.get_pixmap(dpi=96).save(R/'Figure1_current_native_preview.png')
E.ElementTree(root).write(str(R/'Figure1_current_preview.svg'),xml_declaration=True,encoding='utf-8')
checks.update({'status':'LOCAL_MIRROR_PREPARED_AWAITING_BIORENDER_UPDATE','positions_mm':positions,
    'clinical_denominator':172,'hidden_missing_counts':{'age':61,'sex':38},'L2_states':116,
    'marker_shape':[14,28],'marker_area_rule':'area proportional to original detection fraction; zero area at zero',
    'unchanged_sources':{str(p):sha(p) for p in [src/'per_gene_scaled_mean.tsv',src/'original_detection_fraction.tsv',R/'Figure1e_L2_clinical_left_data.tsv',R/'Figure1f_clinical_percent_data.tsv']}})
assert all(sha(Path(k))==v for k,v in protected_hash.items())
(R/'Figure1_current_layout_check.json').write_text(json.dumps(checks,indent=2))
print(json.dumps({'sizes':positions,'clinical_denominator':172,'not_reported_drawn':False,'named_states':len(z),'marker_cells':V.size},indent=2))
