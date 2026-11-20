from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907')
O=Path(__file__).resolve().parents[1];E=R/'04_annotation/level1_marker_evidence'
t=pd.read_csv(O/'labels/primary_supervision_labels.tsv.gz',sep='\t')
if t.columns[0]!='cell_id':
    t=t.rename(columns={t.columns[0]:'cell_id'})
    t.to_csv(O/'labels/primary_supervision_labels.tsv.gz',sep='\t',index=False)
cell_ids=pd.read_csv(E/'cell_ids.tsv.gz',sep='\t').cell_id.to_numpy()
assert np.array_equal(cell_ids[t.source_row_position.to_numpy()],t.cell_id.to_numpy())
m=json.loads((E/'level1_marker_evidence_manifest.json').read_text());genes=m['marker_gene_order']
x=np.load(E/'level1_marker_log1p_cpm_expression.npy',mmap_mode='r')
c=pd.read_csv(E/'level1_marker_program_contract.tsv',sep='\t')
means={};dets={}
for label,z in t.groupby('reliable_L1',sort=True):
    v=np.asarray(x[z.source_row_position.to_numpy()]);means[label]=v.mean(axis=0);dets[label]=(v>0).mean(axis=0)
mean=pd.DataFrame(means,index=genes).T;det=pd.DataFrame(dets,index=genes).T
mean.to_csv(O/'labels/primary_marker_mean.tsv',sep='\t');det.to_csv(O/'labels/primary_marker_detection.tsv',sep='\t')
for k in range(0,len(c),5):
    panel=c.iloc[k:k+5];gg=[g for s in panel.present_markers for g in s.split(';')];v=mean[gg].to_numpy();fraction=det[gg].to_numpy();norm=v/np.maximum(v.max(axis=0),1e-9)
    fig,ax=plt.subplots(figsize=(max(12,len(gg)*.31),7));xx,yy=np.meshgrid(np.arange(len(gg)),np.arange(len(mean)))
    im=ax.scatter(xx.ravel(),yy.ravel(),s=fraction.ravel()*100,c=norm.ravel(),cmap='Reds',vmin=0,vmax=1)
    ax.set_xticks(range(len(gg)),gg,rotation=90);ax.set_yticks(range(len(mean)),[f'{l} (n={int((t.reliable_L1==l).sum()):,})' for l in mean.index]);ax.invert_yaxis();ax.set_ylim(len(mean)-.5,-2)
    start=0
    for row in panel.itertuples():
        n=len(row.present_markers.split(';'));ax.text(start+(n-1)/2,-1.15,row.program,ha='center',fontsize=8);ax.axvline(start-.5,c='grey',lw=.4);start+=n
    fig.colorbar(im,ax=ax,label='Mean expression / gene maximum across labels')
    ax.set_title('Primary reliable L1 review | dot area = detection fraction\nUnknown retained; labels are marker-supported candidates',pad=25)
    for f in [.25,.5,1.]:ax.scatter([],[],s=100*f,c='grey',label=f'{f:.0%}')
    ax.legend(title='Raw > 0',loc='upper left',bbox_to_anchor=(1.2,1))
    fig.tight_layout()
    for ext in ['png','pdf']:fig.savefig(O/f'labels/marker_support_panel{k//5+1}.{ext}',dpi=180,bbox_inches='tight')
    plt.close(fig)
print('MARKER_SUPPORT_COMPLETE')
