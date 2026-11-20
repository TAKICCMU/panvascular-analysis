"""Read-only extraction from existing source-count artifacts. No model or hypothesis test."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import anndata as ad,numpy as np,pandas as pd,hashlib,sys
root=Path('/data/anzhen_vhd/scRNAseq/codex/panvascular_downstream_storage_v1/expanded_parallel_20260916_v1/overnight_manuscript_integration_20260917')
joint=root/'spatial_target_excluded_joint_v1';source=root/'spatial_target_excluded_v1/outputs'
roster=pd.read_csv(joint/'target_model_roster.tsv',sep='\t');pairs=roster[['gene','cell_context']].drop_duplicates();genes=pairs.gene.drop_duplicates().tolist();assert len(pairs)==29
spots=pd.read_csv(joint/'spot_joint_source.tsv.gz',sep='\t');rows=[]
for gse in roster.GSE.unique():
 path=source/(gse+'_full_atlas_gene_UMI.h5ad');a=ad.read_h5ad(path);source_sha=hashlib.sha256(path.read_bytes()).hexdigest()
 assert a.obs_names.is_unique
 present=[g for g in genes if g in a.var_names];X=a[:,present].X;X=X.toarray() if hasattr(X,'toarray') else np.asarray(X)
 assert np.all(np.isfinite(X)) and np.all(X>=0) and np.allclose(X,np.rint(X))
 vals={g:X[:,i] for i,g in enumerate(present)}
 for sid in a.obs.section_id.unique():
  expected=spots.loc[spots.section_id.eq(sid)].set_index('spot_id');idx=a.obs_names.get_indexer(expected.index);assert (idx>=0).all();assert len(idx)==int(a.obs.section_id.eq(sid).sum())
  library=a.obs.iloc[idx].total_counts_original.to_numpy();assert np.array_equal(library,expected.library_size_original.to_numpy())
  for g in genes:
   old=expected['expr__'+g].to_numpy();ok=g in vals and np.isfinite(old).any();v=vals[g][idx] if ok else None;norm=np.log1p(v*10000/library) if ok else None;err=float(np.nanmax(np.abs(norm-old))) if ok else None
   if ok:assert err<1e-10,(g,sid,err)
   else:
    assert np.isnan(old).all()
    if g in vals:assert np.all(vals[g][idx]==0),(g,sid,'missing source but nonzero union matrix')
   rows.append(dict(gene=g,section_id=sid,study=gse,n_spots=len(idx),gene_available=ok,n_positive=int((v>0).sum()) if ok else None,positive_fraction=float((v>0).mean()) if ok else None,target_total_UMI=int(v.sum()) if ok else None,raw_count_value_status='原始counts已核对' if ok else '未提供_原流程记为该切片基因缺失_不把union填零当测量零',normalised_max_abs_difference=err,source_raw_h5ad=str(path),source_raw_h5ad_sha256=source_sha,spot_identity_match='EXACT',library_denominator_match='EXACT',gene_id_in_raw_artifact='SYMBOL_ONLY',source_artifact_n_genes=a.n_vars))
 del a,X,vals
pd.DataFrame(rows).to_csv(sys.stdout,sep='\t',index=False,na_rep='NA')
