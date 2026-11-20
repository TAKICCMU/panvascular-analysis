"""Small synthetic tests; no server dataset or OmicVerse execution."""
import importlib.util
from pathlib import Path
import tempfile
from types import SimpleNamespace
import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py

path=Path(__file__).with_name('l1_cpu_local_refinement.py')
spec=importlib.util.spec_from_file_location('l1refine_test',path)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
assert len(m.APPROVED)==15 and len(m.LOCAL_TARGETS)==13
assert set(m.APPROVED).isdisjoint(m.LOCAL_TARGETS)
assert set(m.EXPECTED_COUNTS)==set(m.APPROVED)|set(m.LOCAL_TARGETS)|{'28'}
assert sum(m.EXPECTED_COUNTS.values())==1175115
assert sum(m.EXPECTED_COUNTS[k] for k in m.APPROVED)==824927
assert sum(m.EXPECTED_COUNTS[k] for k in m.LOCAL_TARGETS)==350118
assert sum(m.EXPECTED_L1.values())==1175115
for name,n in m.EXPECTED_L1.items():
    actual=sum(v for k,v in m.EXPECTED_COUNTS.items() if m.APPROVED.get(k,'Unknown')==name)
    assert actual==n,(name,actual,n)
X=sp.csr_matrix(np.array([[4,0,2],[0,2,0],[1,3,0],[0,0,5]],dtype=np.int32))
rawhash=m.sparse_sha256(X)
ids=pd.Index(['a','b','c','d']); genes=pd.Index(['G1','G2','G3'])
# all-gene per-cell normalization on a copy
lib=np.asarray(X.sum(axis=1)).ravel()
log=X.astype(np.float64).multiply((10000/lib)[:,None]).tocsr();log.data=np.log1p(log.data)
A=SimpleNamespace(X=X,layers={'log1p':log},obs_names=ids,var_names=genes,n_obs=4)
mean,frac,raw,n=m.aggregate_markers(A,['C2_r0p5_0','C2_r0p5_0','C2_r0p5_1','C2_r0p5_1'],genes.tolist())
np.testing.assert_allclose(mean.to_numpy(),np.vstack([log.toarray()[:2].mean(0),log.toarray()[2:].mean(0)]))
np.testing.assert_allclose(frac.to_numpy(),np.array([[.5,.5,.5],[.5,.5,.5]]))
np.testing.assert_allclose(raw.to_numpy(),np.array([[2,1,1],[.5,1.5,2.5]]))
assert m.sparse_sha256(X)==rawhash
m.lock_sparse(X)
assert all(not getattr(X,k).flags.writeable for k in ('data','indices','indptr'))
assert m.sparse_sha256(X)==rawhash
with tempfile.TemporaryDirectory() as td:
    p=Path(td)/'test.h5ad'
    with h5py.File(p,'w') as f:
        g=f.create_group('X');g.attrs['encoding-type']='csr_matrix';g.attrs['shape']=X.shape
        for key in ('data','indices','indptr'):g.create_dataset(key,data=getattr(X,key))
    assert m.h5_sparse_sha256(p)==rawhash
    m.readonly_file(p)
    assert p.stat().st_mode & 0o222==0
print('PASS: reviewed counts/mapping, disjoint membership, raw-hash/HDF5 equivalence, exact all-cell means/detection, and in-memory/file write guards.')
