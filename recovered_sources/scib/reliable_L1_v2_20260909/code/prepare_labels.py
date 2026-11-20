from pathlib import Path
import json, hashlib, shutil
import numpy as np
import pandas as pd
import anndata as ad

R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907')
O=Path(__file__).resolve().parents[1]
W=R/'02_scib_benchmark/studywise_sampleonly_v1'
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def nh(x):return hashlib.sha256(b''.join(str(v).encode()+b'\0' for v in x)).hexdigest()
for d in ['inputs','labels','embeddings','models','logs','metrics','audit']: (O/d).mkdir(exist_ok=True)
assert not (O/'inputs/primary.h5ad').exists(), 'Do not overwrite completed inputs'
p=W/'inputs/studywise_marker_preannotation_all_cells.tsv.gz'
t=pd.read_csv(p,sep='\t',index_col='cell_id')
assert t.index.is_unique
trusted=(t.marker_top_program_n_markers_detected>=3)&(t.marker_top_program_detection_fraction>=.4)&(t.marker_top_to_second_score_ratio>=2)&(t.studywise_marker_level1!='Unknown')
t['reliable_L1']=np.where(trusted,t.studywise_marker_level1,'Unknown')
t['supervision_label']=np.where(trusted,t.studywise_marker_level1,'__unlabeled__')
t['trusted']=trusted
t['decision_reason']=np.select([t.marker_top_program_n_markers_detected<3,t.marker_top_program_detection_fraction<.4,t.marker_top_to_second_score_ratio<2],['insufficient_marker_count','insufficient_marker_fraction','competing_programs'],default='trusted')
t.to_csv(O/'labels/all_cells_L1_evidence.tsv.gz',sep='\t')
contract={'minimum_detected_markers':3,'minimum_program_detection_fraction':.4,'minimum_top_second_score_ratio':2,'unknown_is_unlabeled':True,'thresholds_fixed_before_new_scores':True,'L2_generated':False,'source':str(p),'source_sha256':sha(p),'marker_vocabulary':'unchanged 14 primary canonical programs; rare programs remain audit-only','limitations':'marker-supported provisional labels, not author-curated ground truth; overlapping Fib/SMC/FMC programs may be unresolved'}
(O/'labels/LABEL_CONTRACT.json').write_text(json.dumps(contract,indent=2))
a=ad.read_h5ad(W/'inputs/scib_primary_sampleonly_cap1500_hvg2000.h5ad')
assert nh(a.obs_names)=='e2ac0ad6c47a2349b3b26a129287868dd211b6b43a72b2644a9a5a1ed08f88d4'
q=t.loc[a.obs_names]
for c in ['reliable_L1','supervision_label','decision_reason']:a.obs[c]=pd.Categorical(q[c].to_numpy())
a.obs['trusted']=q.trusted.to_numpy()
q.to_csv(O/'labels/primary_supervision_labels.tsv.gz',sep='\t',index_label='cell_id')
a.write_h5ad(O/'inputs/primary.h5ad')
for scope,z in [('all',t),('primary',q)]:
    s=z.groupby('reliable_L1').agg(n_cells=('trusted','size'),median_markers=('marker_top_program_n_markers_detected','median'),median_detection=('marker_top_program_detection_fraction','median'),median_ratio=('marker_top_to_second_score_ratio','median'))
    s['fraction']=s.n_cells/len(z);s.to_csv(O/f'labels/{scope}_L1_summary.tsv',sep='\t')
    z.groupby(['study_file','reliable_L1']).size().rename('n_cells').to_csv(O/f'labels/{scope}_study_L1_counts.tsv',sep='\t')
records=[]
for method,config in [('pca','baseline'),('harmony','theta1'),('scvi','latent20_layers2'),('scanorama','dim30'),('liger','k20_lambda5'),('scalex','iter30000_bs512')]:
    src=W/f'primary/embeddings/{method}__{config}.npy';meta=json.loads(src.with_suffix('.json').read_text());assert meta['obs_names_sha256']==nh(a.obs_names)
    arr=np.load(src,mmap_mode='r');assert arr.shape[0]==len(a) and np.isfinite(arr).all()
    dest=O/f'embeddings/{method}.npy';dest.symlink_to(src)
    records.append({'method':method,'source':str(src),'sha256':sha(src),'metadata':meta,'retrained':False})
(O/'audit/REUSED_EMBEDDINGS.json').write_text(json.dumps(records,indent=2))
(O/'inputs/PREPARE_COMPLETE.json').write_text(json.dumps({'n_cells':len(a),'n_trusted':int(q.trusted.sum()),'n_unknown':int((~q.trusted).sum()),'cell_ids_sha256':nh(a.obs_names),'input_sha256':sha(O/'inputs/primary.h5ad')},indent=2))
print(q.reliable_L1.value_counts().to_string(),flush=True)
