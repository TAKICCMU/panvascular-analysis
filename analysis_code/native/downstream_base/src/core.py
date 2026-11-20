"""Read-only release contracts and bounded sparse access. No integration operations."""
from __future__ import annotations
import hashlib, json, os, pathlib, tempfile, datetime
import numpy as np
import pandas as pd
import h5py
from scipy import sparse
from anndata.io import read_elem

class Blocked(RuntimeError): pass
MISSING={'','na','nr','nan','none','null','not_available','not available','unknown','<na>'}
def missing(v):
    return pd.isna(v) or str(v).strip().lower() in MISSING

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def digest_ids(values):
    return hashlib.sha256(('\n'.join(map(str,values))+'\n').encode()).hexdigest()

def atomic_json(path,obj):
    path=pathlib.Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.partial')
    temp.write_text(json.dumps(obj,indent=2,ensure_ascii=False,default=str));os.replace(temp,path)

def read_frame(path,slot='obs',columns=None):
    with h5py.File(path,'r') as f:
        g=f[slot]; key=g.attrs.get('_index','_index'); idx=np.asarray(read_elem(g[key])).astype(str)
        cols=columns if columns is not None else [c for c in g if c!=key]
        absent=set(cols)-set(g)
        if absent:raise Blocked(f'Missing {slot} fields: {sorted(absent)}')
        return pd.DataFrame({c:read_elem(g[c]) for c in cols},index=pd.Index(idx,name='cell_id' if slot=='obs' else 'gene_id'))

def sparse_block(f,slot,start,stop):
    g=f[slot]
    if not isinstance(g,h5py.Group):
        # Never dense-read a full atlas; a caller supplied bounded block only.
        return sparse.csr_matrix(g[start:stop])
    enc=g.attrs.get('encoding-type','')
    if enc!='csr_matrix':raise Blocked(f'{slot}: expected CSR row-block access, got {enc}')
    ptr=g['indptr'][start:stop+1];base=int(ptr[0]);end=int(ptr[-1])
    return sparse.csr_matrix((g['data'][base:end],g['indices'][base:end],ptr-base),shape=(stop-start,int(g.attrs['shape'][1])))

def expression_blocks(path,slot='X',block_rows=4096):
    with h5py.File(path,'r') as f:
        g=f[slot];n=int(g.attrs['shape'][0]) if isinstance(g,h5py.Group) else g.shape[0]
        for start in range(0,n,block_rows):
            stop=min(n,start+block_rows);yield start,stop,sparse_block(f,slot,start,stop)

def check_annotations(obs,ontology,l1='level1_annotation',l2='level2_annotation',allow_l2_alias=False):
    if 'L2_annotation' in obs:
        if l2 in obs and not np.array_equal(obs[l2].astype(str),obs['L2_annotation'].astype(str)):
            raise Blocked('Conflicting L2_annotation / level2_annotation')
        if l2 not in obs:
            if not allow_l2_alias:raise Blocked('Explicit L2_annotation migration approval required')
            obs=obs.copy();obs[l2]=obs['L2_annotation']
    if l1 not in obs or l2 not in obs:raise Blocked('Required formal annotation columns missing')
    bad=obs[l2].map(lambda x:missing(x) and str(x)!='Unknown') | obs[l2].astype(str).str.contains('pending|unresolved',case=False,regex=True)
    if bad.any():raise Blocked(f'Unapproved/missing L2 values: {int(bad.sum())}')
    if not (obs.loc[obs[l1]=='Unknown',l2]=='Unknown').all():raise Blocked('L1 Unknown must remain L2 Unknown')
    if not (obs.loc[obs[l1]=='Schwann cell',l2]=='Not_assessed').all():raise Blocked('Schwann must be Not_assessed')
    if ((obs[l2]=='Not_assessed') & (obs[l1]!='Schwann cell')).any():raise Blocked('Not_assessed outside Schwann')
    for c in [l1,l2,'approved']:
        if c not in ontology:raise Blocked(f'Ontology lacks {c}')
    named=ontology.loc[~ontology[l2].isin(['Unknown','Not_assessed'])].copy()
    if named[l2].duplicated().any():raise Blocked('Duplicate named L2 in ontology')
    if not named['approved'].astype(str).str.lower().isin(['true','1']).all():raise Blocked('Ontology not approved')
    expected=dict(zip(named[l2],named[l1])); pairs=obs[[l1,l2]].drop_duplicates()
    for parent,state in pairs.itertuples(index=False,name=None):
        if state not in ['Unknown','Not_assessed'] and expected.get(state)!=parent:raise Blocked(f'Invalid parent for {state}')
    used=set(obs[l2])-{'Unknown','Not_assessed'}
    if used!=set(expected):raise Blocked('Ontology named-state universe differs from observed labels')
    if 'n_cells' in named:
        actual=obs[l2].value_counts()
        if any(int(row.n_cells)!=int(actual[row[l2]]) for _,row in named.iterrows()):raise Blocked('Ontology cell counts differ')
    return obs

def preflight(config,output,verify_hash=True):
    """Production release verification; caller must gate before expression checks."""
    a=config['atlas'];p=pathlib.Path(a['path']);mp=pathlib.Path(a['manifest'])
    if not p.exists() or not mp.exists():raise Blocked('WAITING_FOR_FINAL_L2_ATLAS')
    m=json.loads(mp.read_text())
    if m.get('release_status')!='FINAL_APPROVED' or m.get('ready_for_downstream') is not True or m.get('level1_frozen') is not True or m.get('level2_frozen') is not True:
        raise Blocked('WAITING_FOR_FINAL_L2_ATLAS: final approval absent')
    if pathlib.Path(m.get('h5ad_path','')).resolve()!=p.resolve():raise Blocked('Manifest names a different atlas')
    if m.get('release_version')!=a['expected_release_version']:raise Blocked('Wrong release version')
    expected=m.get('h5ad_sha256');pinned=a.get('input_sha256')
    if not expected or not pinned or expected!=pinned:raise Blocked('Pin actual reviewed release SHA in config')
    if verify_hash and sha(p)!=expected:raise Blocked('Atlas checksum mismatch')
    checksum_file=pathlib.Path(a['checksum']);
    if not checksum_file.exists() or checksum_file.read_text().split()[0]!=expected:raise Blocked('Release SHA sidecar missing/inconsistent')
    contract=json.loads(pathlib.Path(a['schema_contract']).read_text())
    support=pathlib.Path(config['outputs']['release_support'])
    for name in contract['release_sidecars']:
        side=support/name;dig=m.get('release_sidecar_sha256',{}).get(name)
        if not side.exists() or not dig or sha(side)!=dig:raise Blocked(f'Missing/unverified release sidecar {name}')
    compat = bool(a.get('release_compatibility_required') or m.get('support_revision'))
    compat_checks = None
    if compat:
        from .release_compat import validate_release_sources
        compat_checks = validate_release_sources(m,config)
    with h5py.File(p,'r') as f:
        if tuple(f['X'].attrs['shape'])!=tuple(a['expected_shape']):raise Blocked('Shape mismatch')
        if 'log1p' not in f['layers'] or tuple(f['layers/log1p'].attrs['shape'])!=tuple(a['expected_shape']):raise Blocked('Missing/wrong all-gene log1p')
        uns=read_elem(f['uns'])
        for k,v in {'release_status':'FINAL_APPROVED','level1_frozen':True,'level2_frozen':True,'atlas_release_version':a['expected_release_version']}.items():
            if uns.get(k)!=v:raise Blocked(f'In-file release declaration inconsistent: {k}')
        if uns.get('annotation_keys')!={'level1':a['level1_key'],'level2':a['level2_key']}:raise Blocked('Wrong annotation key mapping')
        keys=m.get('global_key_map',{})
        for semantic in ['pca','harmony','umap']:
            key=keys.get(semantic)
            if not key or key not in f['obsm'] or f['obsm'][key].shape[0]!=a['expected_shape'][0]:raise Blocked(f'Missing registered global {semantic}')
        graph_uns=keys.get('neighbors_uns')
        if not graph_uns or graph_uns not in uns:raise Blocked('Registered global neighbors absent')
        for typ in ['connectivities_key','distances_key']:
            k=uns[graph_uns].get(typ)
            if not k or k not in f['obsp']:raise Blocked(f'Registered {typ} absent')
    fields=list(dict.fromkeys([a['level1_key'],a['level2_key']]+contract['data_layout']['metadata_columns_expected_from_v4']))
    with h5py.File(p,'r') as f:
        if 'L2_annotation' in f['obs']:fields.append('L2_annotation')
    obs=read_frame(p,columns=fields);var=read_frame(p,'var')
    if not obs.index.is_unique or not var.index.is_unique:raise Blocked('Duplicate cell/gene IDs')
    if obs[a['sample_key']].map(missing).any() or obs[a['sample_key']].nunique()!=contract['expected_identity']['n_sample_uid']:raise Blocked('Sample identity mismatch')
    if obs[a['level1_key']].value_counts().to_dict()!=contract['expected_identity']['expected_level1_counts']:raise Blocked('Frozen L1 changed')
    ontology=pd.read_csv(support/'level2_ontology.tsv',sep='\t',keep_default_na=False)
    obs=check_annotations(obs,ontology,a['level1_key'],a['level2_key'])
    for name,values in [('cell_id_order_sha256',obs.index),('gene_id_order_sha256',var.index)]:
        from .release_compat import release_digest
        if m.get(name)!=release_digest(values,m,name,compat_required=compat):raise Blocked(f'{name} mismatch under declared release encoding')
    if compat:
        for name,key in [('annotation_sha256',a['level2_key']),('level1_annotation_sha256',a['level1_key'])]:
            if m.get(name)!=release_digest(obs[key].astype(str),m,name,compat_required=True):raise Blocked(f'{name} mismatch')
    # Semantic spot check on spread-out rows, not a claim of full numerical equivalence.
    with h5py.File(p,'r') as f:
        for i in np.unique(np.linspace(0,len(obs)-1,min(32,len(obs)),dtype=int)):
            x=sparse_block(f,'X',int(i),int(i)+1).astype('float64');z=sparse_block(f,'layers/log1p',int(i),int(i)+1)
            if not np.isfinite(x.data).all() or (x.data<0).any() or not np.equal(x.data,np.floor(x.data)).all():raise Blocked('X not finite integer counts')
            if x.sum()>0:x.data*=10000/float(x.sum())
            x.data=np.log1p(x.data);d=(z-x).tocsr()
            if len(d.data) and np.max(np.abs(d.data))>1e-5:raise Blocked('log1p is not CP10k natural log1p on checked rows')
    # Exact release-sidecar metadata and annotation identity, no fuzzy matching.
    ann=pd.read_csv(support/'cell_annotations.tsv.gz',sep='\t',index_col='cell_id',keep_default_na=False)
    if not ann.index.equals(obs.index) or not np.array_equal(ann[[a['level1_key'],a['level2_key']]].astype(str),obs[[a['level1_key'],a['level2_key']]].astype(str)):raise Blocked('Annotation sidecar differs')
    sm=pd.read_csv(support/'sample_metadata.csv',keep_default_na=False)
    if sm[a['sample_key']].duplicated().any() or set(sm[a['sample_key']])!=set(obs[a['sample_key']]):raise Blocked('Sample sidecar keys differ')
    for c in set(sm)-{a['sample_key']}:
        if c in obs:
            vals=obs[[a['sample_key'],c]].copy();vals[c]=vals[c].map(lambda v:'' if missing(v) else str(v))
            if vals.groupby(a['sample_key'],observed=True)[c].nunique().max()>1:raise Blocked(f'Clinical field varies within sample: {c}')
            source=sm.set_index(a['sample_key'])[c].map(lambda v:'' if missing(v) else str(v))
            actual=vals.drop_duplicates(a['sample_key']).set_index(a['sample_key'])[c]
            # numeric age formatting is compared numerically
            equal=np.allclose(pd.to_numeric(actual.reindex(source.index),errors='coerce'),pd.to_numeric(source,errors='coerce'),equal_nan=True) if c=='age' else np.array_equal(actual.reindex(source.index),source)
            if not equal:raise Blocked(f'Sample metadata sidecar mismatch: {c}')
    result={'status':'PASS_RELEASE_PREFLIGHT','sha256':expected,'manifest_sha256':sha(mp),'shape':a['expected_shape'],'n_samples':obs[a['sample_key']].nunique(),'log1p_check':'32 spread-out rows; entire file authenticated by SHA','global_key_map':m['global_key_map']}
    if compat:result['compatibility_checks']=compat_checks
    atomic_json(output,result);return obs,ontology,result
