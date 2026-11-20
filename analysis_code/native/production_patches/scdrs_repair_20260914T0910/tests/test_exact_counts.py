import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
import numpy as np,pandas as pd
from scipy import sparse
from scdrs_adapter import prepare_scdrs_counts
x=sparse.csr_matrix([[1,3],[2,4]],dtype=np.float32);g=pd.DataFrame({'feature_id':['A','B'],'scdrs_gene_symbol':['A','B'],'ensembl_id':['ENSG1',None]})
y,n,a=prepare_scdrs_counts(x,g,['A','B']);assert np.shares_memory(x.data,y.data);assert np.array_equal(x.toarray(),y.toarray());assert a['source_symbol_only_features']==1
h=g.copy();h.scdrs_gene_symbol=['B','A'];z,_,_=prepare_scdrs_counts(x,h,['A','B']);assert np.array_equal(z.toarray(),x.toarray()[:,::-1]);assert not np.shares_memory(x.data,z.data)
assert np.array_equal(x.sum(axis=1),y.sum(axis=1))
print('PASS identity no-copy, nonidentity exact ordering, source-only preservation and row-total checks')
