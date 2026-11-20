from pathlib import Path
import argparse, json, time, fcntl, os
import numpy as np
import anndata as ad
from benchmark_common import set_seed, environment_record, obs_names_sha256

R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907')
O=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('method',choices=['scanvi','scpoli']);args=p.parse_args()
assert not (O/f'embeddings/{args.method}.npy').exists()
set_seed(0)
lock=open(R/'06_reanalysis_takeover_20260908/GPU_HEAVY_STAGE.lock','a');fcntl.flock(lock,fcntl.LOCK_EX)
import torch
assert torch.cuda.is_available()
torch.set_num_threads(8);torch.set_float32_matmul_precision('high')
free,total=torch.cuda.mem_get_info();assert free>12*1024**3, 'Insufficient free GPU memory; no silent CPU fallback'
torch.cuda.set_per_process_memory_fraction(min(.60,(free-4*1024**3)/total))
a=ad.read_h5ad(O/'inputs/primary.h5ad');assert obs_names_sha256(a.obs_names)=='e2ac0ad6c47a2349b3b26a129287868dd211b6b43a72b2644a9a5a1ed08f88d4'
assert np.array_equal(a.obs.supervision_label.astype(str).eq('__unlabeled__').to_numpy(),~a.obs.trusted.to_numpy())
started=time.time()
if args.method=='scanvi':
    import scvi
    scvi.settings.seed=0
    src=R/'02_scib_benchmark/studywise_sampleonly_v1/primary/models/latent20_layers2/scvi__latent20_layers2'
    base=scvi.model.SCVI.load(src,adata=a,accelerator='gpu',device='auto')
    existing=np.load(O/'embeddings/scvi.npy',mmap_mode='r')
    ii=np.linspace(0,len(a)-1,1000,dtype=int)
    actual=base.get_latent_representation(indices=ii)
    np.testing.assert_allclose(actual,existing[ii],rtol=1e-4,atol=1e-4)
    model=scvi.model.SCANVI.from_scvi_model(base,adata=a,labels_key='supervision_label',unlabeled_category='__unlabeled__')
    model.train(max_epochs=50,accelerator='gpu',devices=1,early_stopping=True,check_val_every_n_epoch=5,batch_size=1024,datasplitter_kwargs={'num_workers':8})
    emb=model.get_latent_representation()
    model.save(O/'models/scanvi',save_anndata=False,overwrite=False)
    params={'max_epochs':50,'batch_size':1024,'early_stopping':True,'check_val_every_n_epoch':5,'scvi_reused':str(src),'scvi_latent_reproduction_checked_cells':1000}
else:
    from scarches.models.scpoli import scPoli
    a.X=a.layers['counts'].copy()
    model=scPoli(a,condition_keys='sample_uid',cell_type_keys='supervision_label',unknown_ct_names=['__unlabeled__'],labeled_indices=np.flatnonzero(a.obs.trusted.to_numpy()).tolist(),latent_dim=10,embedding_dims=10,recon_loss='nb')
    assert '__unlabeled__' not in model.cell_types_
    assert len(model.labeled_indices_)==int(a.obs.trusted.sum())
    model.train(n_epochs=50,pretraining_epochs=40,batch_size=512,n_workers=8,eta=5.,prototype_training=True,unlabeled_prototype_training=False)
    emb=model.get_latent(a,mean=True)
    model.save(str(O/'models/scpoli'),overwrite=False)
    params={'epochs':50,'pretraining_epochs':40,'batch_size':512,'n_workers':8,'eta':5.,'latent_dim':10,'embedding_dims':10,'recon_loss':'nb','unlabeled_prototype_training':False}
emb=np.asarray(emb,dtype=np.float32);assert emb.shape[0]==len(a) and np.isfinite(emb).all()
np.save(O/f'embeddings/{args.method}.npy',emb)
(O/f'embeddings/{args.method}.json').write_text(json.dumps({'status':'complete','method':args.method,'seed':0,'input_cell_ids_sha256':obs_names_sha256(a.obs_names),'parameters':params,'n_unlabeled':int((~a.obs.trusted).sum()),'training_label_key':'supervision_label','evaluation_label_key':'reliable_L1','batch_key':'sample_uid','seconds':time.time()-started,'environment':environment_record(['torch','scvi-tools','scarches']),'GPU':torch.cuda.get_device_name()},indent=2))
print('COMPLETE',args.method,emb.shape,flush=True)
