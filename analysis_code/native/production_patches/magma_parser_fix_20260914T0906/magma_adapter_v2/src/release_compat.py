"""Bound external release adapters only; never rewrite atlas data or generic ID hashes."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import h5py
from anndata.io import read_elem
from .core import Blocked,sha,digest_ids,read_frame

LENGTH_ENCODING='length_prefixed_utf8_uint64_le_v1'
NEWLINE_ENCODING='newline_terminated_utf8_v1'

def encoded_digest(values,encoding):
    if encoding==NEWLINE_ENCODING:return digest_ids(values)
    if encoding!=LENGTH_ENCODING:raise Blocked('Unknown declared release hash encoding: '+str(encoding))
    h=hashlib.sha256()
    for v in values:
        b=str(v).encode('utf-8');h.update(len(b).to_bytes(8,'little'));h.update(b)
    return h.hexdigest()

def bound_file(record,cache=None):
    cache={} if cache is None else cache
    if not isinstance(record,dict) or not record.get('path') or not record.get('sha256'):raise Blocked('Missing file binding')
    p=Path(record['path'])
    if not p.is_file():raise Blocked('Bound source file absent: '+str(p))
    # Invocation-local cache only: no persistent stat-only authentication.
    key=str(p.resolve())
    if key not in cache:cache[key]=sha(p)
    if cache[key]!=record['sha256']:raise Blocked('Bound source SHA mismatch: '+str(p))
    return p

def validate_manifest_binding(m,config):
    if m.get('support_revision')!='acceptance_r1':raise Blocked('Unknown/missing compatibility support revision')
    if m.get('downstream_execution_authorized') is not False:raise Blocked('Compatibility release cannot grant execution')
    original=bound_file(m.get('original_manifest',{}));old=json.loads(original.read_text())
    if m.get('original_manifest_sha256')!=sha(original):raise Blocked('Original manifest binding differs')
    for k in ['h5ad_path','h5ad_sha256','release_version','cell_id_order_sha256','gene_id_order_sha256','annotation_sha256','level1_annotation_sha256','global_key_map','level1_counts_verified','shape_verified']:
        if m.get(k)!=old.get(k):raise Blocked('Compatibility changed protected manifest field: '+k)
    publisher=bound_file(m.get('publisher_script',{}))
    expected_publisher=Path(old['finalization_run'])/'audit/finalize_level2.py'
    if publisher.resolve()!=expected_publisher.resolve():raise Blocked('Publisher script not from original committed run')
    encoding=m.get('hash_encoding',{})
    for k in ['cell_id_order_sha256','gene_id_order_sha256','annotation_sha256','level1_annotation_sha256']:
        if encoding.get(k)!=LENGTH_ENCODING:raise Blocked('Missing/unknown bound encoding for '+k)
    if m.get('release_support_dir')!=config['outputs']['release_support']:raise Blocked('Support path differs')
    return old

def release_digest(values,m,key,compat_required=False):
    if m.get('support_revision'):
        enc=m.get('hash_encoding',{}).get(key)
        if enc!=LENGTH_ENCODING:raise Blocked('Undeclared release encoding '+key)
        return encoded_digest(values,enc)
    if compat_required:raise Blocked('Compatibility manifest with explicit encoding required')
    # Unchanged legacy contract. Compatibility mode is never inferred by trying algorithms.
    return digest_ids(values)

def validate_metadata_sources(path):
    meta=json.loads(Path(path).read_text());cache={}
    if meta.get('schema')!='metadata_review_sources_acceptance_r1':raise Blocked('Metadata sources schema differs')
    required={'imported_metadata_csv','import_code','import_run','bundle_checksums','sample_metadata_snapshot','metadata_provenance','user_confirmation','audit_A','audit_B','merge_record','public_GSA_evidence'}
    roles={x.get('role') for x in meta.get('files',[])}
    if not required<=roles:raise Blocked('Missing actual metadata evidence roles: '+repr(sorted(required-roles)))
    for record in meta['files']:
        if not record.get('evidence_type') or not record.get('scope') or not record.get('version'):raise Blocked('Metadata evidence lacks scope/type/version')
        bound_file(record,cache)
    if meta.get('import_link_validation',{}).get('status')!='PASS':raise Blocked('Metadata import linkage unverified')
    # Exact recorded imported CSV is compared to release sample snapshot, independently of claims.
    def single(role):
        rows=[x for x in meta['files'] if x['role']==role]
        if len(rows)!=1:raise Blocked('Nonunique '+role)
        return rows[0]
    patch=single('imported_metadata_csv');sample=single('sample_metadata_snapshot');bundle=single('bundle_checksums')
    sums=json.loads(Path(bundle['path']).read_text())
    if sums.get(meta['import_link_validation']['bundle_member'])!=patch['sha256']:raise Blocked('Imported CSV not bound by recorded bundle checksum')
    a=pd.read_csv(patch['path'],dtype=str,keep_default_na=False).set_index('sample_uid');b=pd.read_csv(sample['path'],dtype=str,keep_default_na=False).set_index('sample_uid')
    if a.index.duplicated().any() or b.index.duplicated().any() or set(a.index)!=set(b.index):raise Blocked('Metadata source sample keys differ')
    for c in a:
        if c not in b:raise Blocked('Final metadata missing imported field '+c)
        x=a[c];y=b.reindex(a.index)[c]
        if c=='age':
            if not np.allclose(pd.to_numeric(x,errors='coerce'),pd.to_numeric(y,errors='coerce'),equal_nan=True):raise Blocked('Imported age differs')
        else:
            # CSV blanks and final empty/NA serialization only; identities never numeric/trimmed.
            y=y.mask(y.isin(['nan','NaN','None','<NA>']), '')
            if not np.array_equal(x,y):raise Blocked('Imported metadata differs from final sample snapshot: '+c)
    return {'verified_files':len(cache),'sample_count':len(a),'fields_checked':list(a),'external_web_fetched_this_run':False}

def validate_registry(path):
    reg=pd.read_csv(path,sep='\t',dtype=str,keep_default_na=False);cache={};records=[]
    required={'role','path','sha256','version','n_obs','n_vars','label_key','graph_key','embedding_key','mapping_path','mapping_sha256','version_source','version_source_sha256'}
    if not required<=set(reg):raise Blocked('Registry missing fields: '+repr(sorted(required-set(reg))))
    for _,r in reg.iterrows():
        if r.role!='lineage_analysis':raise Blocked('Unrecognized registry role: '+r.role)
        p=bound_file({'path':r.path,'sha256':r.sha256},cache);mp=bound_file({'path':r.mapping_path,'sha256':r.mapping_sha256},cache)
        vp=bound_file({'path':r.version_source,'sha256':r.version_source_sha256},cache)
        if not r.version or r.version!=str(vp.parent):raise Blocked('Version must be recorded originating run directory')
        with h5py.File(p,'r') as f:
            g=f['obs'];ids=np.asarray(read_elem(g[g.attrs['_index']])).astype(str)
            v=f['var'];ng=len(v[v.attrs['_index']])
            if (len(ids),ng)!=(int(r.n_obs),int(r.n_vars)):raise Blocked('Registry shape mismatch')
            if r.label_key not in g:raise Blocked('Invalid registry label key '+r.label_key)
            if r.graph_key not in f['uns'] or r.embedding_key not in f['obsm']:raise Blocked('Invalid registry graph/embedding key')
            graph=read_elem(f['uns'][r.graph_key])
            for key in ['connectivities_key','distances_key']:
                if graph.get(key) not in f['obsp']:raise Blocked('Registered graph slot absent')
                if tuple(f['obsp'][graph[key]].attrs['shape'])!=(len(ids),len(ids)):raise Blocked('Registry graph shape mismatch')
            if f['obsm'][r.embedding_key].shape[0]!=len(ids):raise Blocked('Registry embedding shape mismatch')
            labels=np.asarray(read_elem(g[r.label_key])).astype(str)
        mapping=pd.read_csv(mp,sep='\t',dtype=str,keep_default_na=False)
        if 'cell_id' not in mapping or r.label_key not in mapping or not np.array_equal(mapping.cell_id,ids) or not np.array_equal(mapping[r.label_key],labels):raise Blocked('Registry mapping does not match exact source cells/labels')
        if 'label_usage_path' in reg:
            up=bound_file({'path':r.label_usage_path,'sha256':r.label_usage_sha256},cache)
            usage=pd.read_csv(up,sep='\t',dtype=str,keep_default_na=False)
            if not {'branch','label_key','selection_rule','original_cluster_uid'}<=set(usage):raise Blocked('Registry usage schema missing')
            selected=usage[(usage.branch==r.branch)&(usage.label_key==r.label_key)]
            if any(v and v not in set(labels) for v in selected.original_cluster_uid):raise Blocked('Registry reviewed cluster UID not in source partition')
        records.append({'branch':r.get('branch'),'label_key':r.label_key,'n_obs':len(ids),'n_vars':ng})
    return {'rows':len(reg),'unique_files_authenticated':len(cache),'records':records}

def validate_release_sources(m,config):
    validate_manifest_binding(m,config)
    support=Path(config['outputs']['release_support'])
    meta=validate_metadata_sources(support/'metadata_review_sources.json')
    reg=validate_registry(support/'source_registry.tsv')
    return {'metadata':meta,'registry':reg}
