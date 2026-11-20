"""Display-only paired UMAPs from frozen official trusted-cell 15-NN graphs.
No integration, no clustering, no effect on scores. CPU only for both UMAPs.
"""
from pathlib import Path
import os,json,sys,time
os.environ['CUDA_VISIBLE_DEVICES']='';os.environ['JAX_PLATFORMS']='cpu'
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='8'
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
from sklearn.utils import check_random_state
from umap.umap_ import fuzzy_simplicial_set
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];O=R.parent/'reliable_L1_v2_20260909';method=sys.argv[1]
assert method in ['harmony_cpu','harmony_gpu']
source=ad.read_h5ad(R/'inputs/primary.h5ad',backed='r');obs=source.obs.copy();source.file.close();mask=obs.trusted.to_numpy(bool);obs=obs.loc[mask].copy();z=np.asarray(np.load(R/f'embeddings/{method}.npy',mmap_mode='r')[mask],dtype=np.float32)
a=ad.AnnData(obs=obs);a.obsm['X_harmony']=z
coords=R/f'figures/{method}_umap.npy'
if coords.exists():xy=np.load(coords)
else:
 graphpath=R/'metrics/harmony_cpu.neighbors.npz' if method=='harmony_cpu' else O/'metrics/harmony.neighbors.npz'
 with np.load(graphpath) as g:idx=g['indices'][:,:15].copy();dist=g['distances'][:,:15].copy()
 assert len(idx)==len(obs)
 graph=fuzzy_simplicial_set(z,15,check_random_state(0),'euclidean',knn_indices=idx,knn_dists=dist)[0].tocsr()
 a.obsp['connectivities']=graph
 a.uns['neighbors']={'connectivities_key':'connectivities','params':{'n_neighbors':15,'method':'umap','metric':'euclidean','use_rep':'X_harmony','random_state':0}}
 sc.tl.umap(a,min_dist=.3,spread=1.,maxiter=500,init_pos='spectral',n_components=2,alpha=1.,negative_sample_rate=5,random_state=0)
 xy=a.obsm['X_umap'];assert np.isfinite(xy).all();np.save(coords,xy);sparse.save_npz(R/f'metrics/{method}_visual_connectivities_n15.npz',graph)
fields=['reliable_L1','study_file' if 'study_file' in obs else 'study','sample_uid']
for field in fields:
 vals=obs[field].astype(str);cats=sorted(vals.unique());colors=plt.colormaps['tab20'](np.linspace(0,1,len(cats))) if len(cats)<=20 else plt.colormaps['turbo'](np.linspace(0,1,len(cats)));pal=dict(zip(cats,colors));order=np.random.default_rng(0).permutation(len(obs))
 fig,ax=plt.subplots(figsize=(10,8));ax.scatter(xy[order,0],xy[order,1],c=np.array([pal[v] for v in vals])[order],s=.2,linewidths=0,rasterized=True);ax.set(title=f'{method} | {field} | n={len(obs):,}',xlabel='UMAP 1',ylabel='UMAP 2');ax.spines[['top','right']].set_visible(False)
 if len(cats)<=30:
  handles=[plt.Line2D([],[],marker='o',linestyle='',color=pal[c],label=c,markersize=4) for c in cats];ax.legend(handles=handles,bbox_to_anchor=(1.01,1),loc='upper left',fontsize=6,frameon=False)
 fig.savefig(R/f'figures/{method}_{field}.png',dpi=240,bbox_inches='tight');fig.savefig(R/f'figures/{method}_{field}.pdf',bbox_inches='tight');plt.close(fig)
 pd.DataFrame(colors,index=cats,columns=['red','green','blue','alpha']).to_csv(R/f'tables/{field}_palette.tsv',sep='\t') if method=='harmony_cpu' else None
pd.DataFrame(xy,index=obs.index,columns=['UMAP1','UMAP2']).to_csv(R/f'tables/{method}_umap_coordinates.tsv.gz',sep='\t',index_label='cell_id')
(R/f'audit/{method}_visualization.json').write_text(json.dumps({'CPU_only':True,'n_cells':len(obs),'subset':'identical frozen trusted evaluation cells','neighbors':15,'source_graph':'official frozen GPU graph reused; newly computed CPU graph' ,'min_dist':.3,'maxiter':500,'seed':0,'Leiden':'not run: original benchmark uses KMeans NMI/ARI','scope':'display only, no score change'},indent=2))
