
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,hashlib,sys
import pandas as pd,anndata as ad
s=Path(__file__).resolve().parents[1]
for dtype in ['category','object','string']:
 a=pd.Series(['A',None,'','NA'],dtype=dtype)
 assert a.astype('string').fillna('').astype(str).tolist()==['A','','','NA']
try:pd.Series(['A',None],dtype='category').fillna('')
except TypeError:pass
else:raise AssertionError('Old failure not reproduced')
n=0
for f in sorted((s/'outputs/full_umi').glob('*.h5ad')):
 a=ad.read_h5ad(f,backed='r')
 maps=list((s/'evidence').glob('*__'+f.stem+'__gene_map.tsv.gz'));assert len(maps)==1
 gm=pd.read_csv(maps[0],sep='\t',keep_default_na=False)
 for col in ['source_id','source_symbol','atlas_gene']:
  assert a.var[col].astype('string').fillna('').astype(str).tolist()==gm[col].astype(str).tolist()
 a.file.close();n+=1
assert n==49,n
(s/'validation_patch_v2/TEST_PASS.json').write_text(json.dumps({'status':'PASS','categorical_bug_reproduced':True,'string_object_category_fixture':True,'real_sections_gene_map_checked':n,'no_model_or_input_modification':True},indent=2))
print('PASS',n)
