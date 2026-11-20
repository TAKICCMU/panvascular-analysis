"""Display-only Figure 1b: frozen cell coordinates, labels and colours.

No panel letter is baked into the image; BioRender owns that text object.
The six tissue classes reuse the reviewed sample mapping without inference.
"""
from pathlib import Path
import json, shutil
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

MM = 72 / 25.4

def draw_panel(D, pal, TP, LIM, source, out):
    plt.rcParams.update({'font.family':'Arial','font.size':6.2,
        'svg.fonttype':'none','pdf.fonttype':42,'text.color':'#222831',
        'savefig.facecolor':'white'})
    f=plt.figure(figsize=(107/25.4,51/25.4))
    def axis(x,y,w,h): return f.add_axes([x/107,1-(y+h)/51,w/107,h/51])
    Z=D[D.level1_annotation.ne('Unknown')]
    # Larger primary map; the inset retains the same limits and equal aspect.
    for x,y,w,h,key in [(18,7,46,38.312824751456525,'level1_annotation'),(80.5,7,25.5,21.238631112220465,'Tissue')]:
        a=axis(x,y,w,h)
        a.scatter(Z.UMAP1,Z.UMAP2,c=Z[key].map(pal if key=='level1_annotation' else TP),
            s=.035,lw=0,alpha=.85,rasterized=True)
        a.set(xlim=LIM[:2],ylim=LIM[2:]);a.set_aspect('equal');a.axis('off')
        if key=='level1_annotation': main=a
    overlay=axis(0,0,107,51);overlay.set(xlim=(0,107),ylim=(51,0));overlay.axis('off')
    overlay.text(93.25,3.8,'Tissue',fontsize=6.5,ha='center',va='center')
    # Short, fixed label callouts: full biological identities are unchanged.
    positions={
        'Schwann cell':(17.1,12,'right'), 'Plasma cell':(17.1,17.4,'right'),
        'Pericyte':(17.1,22.8,'right'), 'Fibroblast':(17.1,28.2,'right'),
        'Fibromyocyte':(17.1,34.0,'right'), 'Endothelial cell':(17.1,41.5,'right'),
        'B cell':(32.5,13.5,'center'), 'Macrophage':(47,5.7,'center'),
        'Neutrophil':(56,12,'left'), 'Monocyte':(56,18,'left'),
        'Dendritic cell':(56,24,'left'), 'Natural killer cell':(64.3,34,'left'),
        'T cell':(57,43,'center'), 'Smooth muscle cell':(33,44.8,'center'),
        'Mast cell':(44.0,35.6,'center')}
    names={'Natural killer cell':'Natural killer\ncell','Smooth muscle cell':'Smooth muscle\ncell',
           'Endothelial cell':'Endothelial\ncell'}
    A=pd.read_csv(source/'Current_Display_20260920/Figure1b_display15_direct_labels_20260921.tsv',sep='\t')
    assert set(A.L1)==set(positions)
    f.canvas.draw()
    xy=main.transData.transform(A[['UMAP1','UMAP2']].to_numpy())
    xy=f.transFigure.inverted().transform(xy);xy=np.c_[xy[:,0]*107,(1-xy[:,1])*51]
    texts=[]; routes=[]
    for i,r in A.iterrows():
        x,y,ha=positions[r.L1]
        t=overlay.text(x,y,names.get(r.L1,r.L1),fontsize=6.3,ha=ha,va='center',
            linespacing=1.05,zorder=4,bbox=dict(fc='white',ec='none',pad=.18))
        texts.append(t)
        A.loc[i,'label_x_mm']=x;A.loc[i,'label_y_mm']=y
        A.loc[i,'display_label']=names.get(r.L1,r.L1);A.loc[i,'display_hex']=pal[r.L1]
    f.canvas.draw();renderer=f.canvas.get_renderer()
    boxes=[]
    for i,(r,t) in enumerate(zip(A.itertuples(),texts)):
        bb=t.get_window_extent(renderer).transformed(overlay.transData.inverted())
        x0,x1=sorted([bb.x0,bb.x1]);y0,y1=sorted([bb.y0,bb.y1]);boxes.append([x0,y0,x1,y1])
        axx,ayy=xy[i];ha=positions[r.L1][2]
        if ha=='right':end=np.array([x1+.9,positions[r.L1][1]]);elbow=end+np.array([1.8,0])
        elif ha=='left':end=np.array([x0-.9,positions[r.L1][1]]);elbow=end-np.array([1.8,0])
        else:
            end=np.array([np.clip(axx,x0+.5,x1-.5),y1+.75 if ayy>y1 else y0-.75])
            elbow=end+np.array([0,1.15 if ayy>y1 else -1.15])
        route=np.vstack([xy[i],elbow,end]);routes.append(route)
        overlay.plot(route[:,0],route[:,1],color='#919CA5',lw=.28,zorder=2,solid_capstyle='round')
        # A short coloured terminal helps match the full name to its category.
        v=end-elbow;v=v/(np.linalg.norm(v) or 1)
        p=end-v*.65
        overlay.plot([p[0],end[0]],[p[1],end[1]],color=pal[r.L1],lw=.9,zorder=3,solid_capstyle='round')
    # The tissue legend uses only the six approved labels and exact existing colours.
    tissue_labels={k:k for k in TP};tissue_labels['Other systemic arteries']='Other systemic\narteries'
    ty=[31.3,34.4,37.5,40.6,43.7,47.3]
    for (k,color),y in zip(TP.items(),ty):
        overlay.add_patch(Rectangle((80.5,y-.55),1.25,1.1,fc=color,ec='none'))
        overlay.text(82.7,y,tissue_labels[k],fontsize=6.0,ha='left',va='center',linespacing=1.0)
    def cross(a,b,c,d):
        def cr(x,y,z):u,v=y-x,z-x;return u[0]*v[1]-u[1]*v[0]
        return cr(a,b,c)*cr(a,b,d)<-1e-8 and cr(c,d,a)*cr(c,d,b)<-1e-8
    crossings=[]
    for i in range(len(routes)):
        for j in range(i):
            if any(cross(a,b,c,d) for a,b in zip(routes[i][:-1],routes[i][1:]) for c,d in zip(routes[j][:-1],routes[j][1:])):
                crossings.append([A.L1.iloc[i],A.L1.iloc[j]])
    collisions=[]
    for i,b in enumerate(boxes):
        for j,c in enumerate(boxes[:i]):
            if min(b[2],c[2])>max(b[0],c[0]) and min(b[3],c[3])>max(b[1],c[1]):collisions.append([A.L1.iloc[i],A.L1.iloc[j]])
    f.canvas.draw();renderer=f.canvas.get_renderer()
    text_outside=[]
    for t in overlay.texts:
        bb=t.get_window_extent(renderer).transformed(overlay.transData.inverted())
        if min(bb.x0,bb.x1)<0 or max(bb.x0,bb.x1)>107 or min(bb.y0,bb.y1)<0 or max(bb.y0,bb.y1)>51:
            text_outside.append(t.get_text())
    assert not crossings and not collisions and not text_outside, (crossings,collisions,text_outside)
    out=Path(out); (out/'panels').mkdir(exist_ok=True)
    for ext in ['pdf','svg','png']:f.savefig(out/'panels'/f'Figure1b.{ext}',dpi=900)
    f.savefig(out/'review/Figure1b.png',dpi=300)
    plt.close(f)
    A.to_csv(out/'Figure1b_label_anchors.tsv',sep='\t',index=False)
    pd.DataFrame([{'Tissue':k,'display_label':tissue_labels[k],'hex':v} for k,v in TP.items()]).to_csv(out/'Figure1b_Tissue_legend.tsv',sep='\t',index=False)
    check={'cells_displayed':len(Z),'Unknown_not_displayed':len(D)-len(Z),
       'full_coordinate_limits':LIM,'line_crossings':len(crossings),'crossing_names':crossings,
       'label_collisions':collisions,'text_outside_panel':text_outside,'tissue_classes':list(TP),'palette_unchanged':True,
       'panel_letter_embedded':False,'label_font_pt':6.3,'tissue_legend_font_pt':6.0}
    (out/'review/Figure1b_layout_check.json').write_text(json.dumps(check,indent=2))
    shutil.copy2(out/'panels/Figure1b.png',out/'BioRender_upload/Figure1b.png')
    return check

if __name__=='__main__':
    out=Path(__file__).parent; P=out.parents[2]; source=P/'Source_Data_Figures'
    palette=pd.read_csv(out/'L1_display_palette.tsv',sep='\t').set_index('L1').hex.to_dict()
    tissue=pd.read_csv(out/'Tissue_display_palette.tsv',sep='\t').set_index('Tissue').hex.to_dict()
    D=pd.read_csv(source/'Figure1/full_cells_frozen_umap.tsv.gz',sep='\t')
    T=pd.read_csv(source/'Current_Display_20260920/Figure1b_tissue_sample_mapping.tsv',sep='\t')
    assert T.sample_uid.is_unique and len(D)==1175115
    D['Tissue']=D.sample_uid.map(T.set_index('sample_uid').Tissue)
    assert D.Tissue.notna().all() and set(D.Tissue)==set(tissue)
    limits=json.loads((out/'preserved_data_checks.json').read_text())['b']['full_coordinate_limits']
    print(json.dumps(draw_panel(D,palette,tissue,limits,source,out),indent=2))
