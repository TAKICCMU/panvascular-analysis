"""Prepare CellHint primary harmonization once, before meta-neighbor selection.

Cache identical confusion-matrix requests within the 15-way native threshold
scan. Copies are returned, so native downstream mutation cannot corrupt cache.
All native parameter grids, dataset ordering, filtering and alignment rules stay
unchanged. Original pilot processes remain untouched.
"""
import os
os.environ.update(OMP_NUM_THREADS='4', MKL_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4',
                  NUMEXPR_NUM_THREADS='4', NUMBA_NUM_THREADS='4', LOKY_MAX_CPU_COUNT='4',
                  PYTHONHASHSEED='0')
if hasattr(os, 'sched_getaffinity'):
    os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:4])
import argparse, hashlib, json, random, time
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import cellhint
from cellhint.distance import Distance

R = Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907')
W = R / '02_scib_benchmark/studywise_sampleonly_v1'
OUT = R / '06_reanalysis_takeover_20260908/scib'


def fh(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        while b:=f.read(8*1024*1024):h.update(b)
    return h.hexdigest()


def nh(v):
    h=hashlib.sha256()
    for x in v:h.update(str(x).encode());h.update(b'\0')
    return h.hexdigest()


class ConfusionCache:
    def __init__(self):
        self.original=Distance.to_multi_confusion
        self.key=None; self.value=None; self.hits=0; self.misses=0
    def install(self):
        owner=self
        def cached(instance, relation, D, check=True):
            key=(id(instance),str(D),bool(check),tuple(relation.columns),tuple(relation.index),
                 tuple(map(tuple,relation.to_numpy(dtype=object))))
            if key!=owner.key:
                owner.value=owner.original(instance,relation,D,check)
                owner.key=key;owner.misses+=1
            else:owner.hits+=1
            return tuple(x.copy(deep=True) for x in owner.value)
        Distance.to_multi_confusion=cached
        return self
    def restore(self):Distance.to_multi_confusion=self.original


def smoke(out):
    rng=np.random.default_rng(0)
    n=360; batch=np.repeat(['s1','s2','s3','s4'],90)
    label=np.tile(np.repeat(['A','B','C'],30),4)
    x=rng.normal(0,.3,(n,8))+np.array([{'A':0,'B':3,'C':6}[z] for z in label])[:,None]
    a=ad.AnnData(x.astype(np.float32),obs=pd.DataFrame({'sample_uid':batch,'studywise_marker_level1':label},index=[f'c{i}' for i in range(n)]))
    a.obsm['X_pca_omicverse']=x.astype(np.float32)
    kw=dict(dataset='sample_uid',cell_type='studywise_marker_level1',use_rep='X_pca_omicverse',use_pct=False,filter_cells=True)
    original=cellhint.harmonize(a,**kw)
    cache=ConfusionCache().install()
    try:optimized=cellhint.harmonize(a,**kw)
    finally:cache.restore()
    pd.testing.assert_frame_equal(original.relation,optimized.relation)
    pd.testing.assert_frame_equal(original.reannotation,optimized.reannotation)
    np.testing.assert_array_equal(original.minimum_unique_percents,optimized.minimum_unique_percents)
    np.testing.assert_array_equal(original.minimum_divide_percents,optimized.minimum_divide_percents)
    result={'status':'passed','n_cells':n,'n_samples':4,'identical_relation':True,'identical_reannotation':True,'identical_selected_thresholds':True,'cache_hits':cache.hits,'cache_misses':cache.misses}
    (out/'CACHE_EQUIVALENCE_SMOKE.json').write_text(json.dumps(result,indent=2));print(result,flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--smoke',action='store_true');args=parser.parse_args()
    out=OUT/'primary/shared_harmonization';out.mkdir(parents=True,exist_ok=True)
    if args.smoke:smoke(out);return
    if not (out/'CACHE_EQUIVALENCE_SMOKE.json').exists():raise RuntimeError('Equivalence smoke required')
    if (out/'SHARED_HARMONIZATION_COMPLETE.json').exists():print('Already complete');return
    random.seed(0);np.random.seed(0)
    p=W/'inputs/scib_primary_sampleonly_cap1500_hvg2000.h5ad'
    pp=W/'omicverse_pca/primary/omicverse_pca.npy';pm=pp.with_name('omicverse_pca.metadata.json')
    a=ad.read_h5ad(p,backed='r');obs=a.obs.copy();names=a.obs_names.copy();a.file.close()
    obs['sample_uid']=obs['sample_uid'].astype(str).astype('category')
    obs['studywise_marker_level1']=obs['studywise_marker_level1'].astype(str).astype('category')
    meta=json.loads(pm.read_text());pca=np.load(pp,mmap_mode='r')
    assert meta['obs_names_sha256']==nh(names) and meta['artifacts']['pca_sha256']==fh(pp)
    assert len(names)==315949 and obs.sample_uid.nunique()==247
    # Harmonize explicitly consumes use_rep, so expression storage is unnecessary.
    data=ad.AnnData(X=np.zeros((len(obs),1),dtype=np.float32),obs=obs)
    data.obsm['X_pca_omicverse']=np.asarray(pca,dtype=np.float32)
    contract={'input':str(p),'input_sha256':fh(p),'pca':str(pp),'pca_sha256':fh(pp),'obs_names_sha256':nh(names),
              'n_cells':len(names),'batch_key':'sample_uid','label_key':'studywise_marker_level1','n_samples':247,
              'parameters':{'use_pct':False,'filter_cells':True,'use_rep':'X_pca_omicverse'},
              'optimization':'single-entry confusion matrix memoization; deep copies; all native harmonization thresholds/order unchanged',
              'reusable_for':['meta3','meta5'],'cpu_affinity':sorted(os.sched_getaffinity(0))}
    (out/'SHARED_HARMONIZATION_RUNNING.json').write_text(json.dumps(contract,indent=2))
    started=time.time();cache=ConfusionCache().install()
    try:
        alignment=cellhint.harmonize(data,dataset='sample_uid',cell_type='studywise_marker_level1',use_rep='X_pca_omicverse',use_pct=False,filter_cells=True)
    finally:cache.restore()
    alignment.reannotation.to_csv(out/'cellhint_primary_reannotation.tsv.gz',sep='\t',compression='gzip',index_label='cell_id')
    alignment.relation.to_csv(out/'cellhint_primary_relation.tsv',sep='\t',index=False)
    np.savez_compressed(out/'selected_harmonization_thresholds.npz',minimum_unique_percents=alignment.minimum_unique_percents,minimum_divide_percents=alignment.minimum_divide_percents)
    result=contract|{'status':'completed','elapsed_seconds':time.time()-started,'cache_hits':cache.hits,'cache_misses':cache.misses,'reannotated_cells':len(alignment.reannotation),
                     'reannotation_sha256':fh(out/'cellhint_primary_reannotation.tsv.gz')}
    (out/'SHARED_HARMONIZATION_COMPLETE.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
