# -*- coding: utf-8 -*-
"""Approved working L2 writeback + evidence-preserving metadata + limited CPU review.

This module never re-integrates/re-clusters, never changes frozen level1_annotation,
never fabricates clinical values/donor identities, and never writes the input h5ad.
Run from the existing notebook with run_from_notebook(); all scientific work runs
in fresh CPU-only processes using the same Python executable. Existing notebook
objects, including Macrophage_L2, are not mutated. See README_CN.md for scope.
"""
from __future__ import annotations
import argparse, base64, csv, gc, hashlib, html, itertools, json, os, re
import shutil, stat, subprocess, sys, time, traceback, zipfile
from pathlib import Path
from datetime import datetime

CODE_VERSION = 'L2_targeted_followup_metadata_v4_1.0'
ROOT = Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908')
BUNDLE = Path(__file__).resolve().parent
L1 = 'level1_annotation'
L2 = 'L2_annotation'
EXPECTED_L1 = {'Endothelial cell':104415,'Fibroblast':251371,'Fibromyocyte':29347,
 'Smooth muscle cell':243346,'Pericyte':4513,'Macrophage':210967,'Monocyte':31046,
 'Neutrophil':18509,'T cell':156349,'Natural killer cell':25719,'B cell':10496,
 'Plasma cell':4802,'Mast cell':12929,'Dendritic cell':10556,'Schwann cell':2832,'Unknown':57918}
EXTRA_MARKERS = {
 'Macrophage':['KYNU','ATP13A3','CD83','RGS1','SGK1','CYP27A1','PSD3','PCSK6','NCKAP5','SELENOP','IGF1','SLC40A1','PHLDA1','CCL20','FCN1','VCAN','PDE4C','LRRC23','TTN','TNNT2','MYH6','PRX','SOX10'],
 'Fibroblast':['ITGA8','LTBP2','RSAD2','CDH19','SOX10','ERBB3','S100B','GFAP','SCARA5','PI16','CLDN1','THBS4','FAP','TNC','PTGDS'],
 'SMC':['PRKG1','PDE4D','CACNA1C','DMD','MYOCD','CARMN','DLX5','SOST','RUNX2','SP7','CDH19','SOX10','CNTNAP5','DPP10','EFEMP1'],
 'Endothelial':['ACKR1','ZNF385D','OMD','NRG1','BMPER','SMOC1','PGF','APLN','SFTPC','SFTPB','PROX1','SEMA3G','AKR1C1'],
 'Fibromyocyte':['ITGA10','CD9','PRKG1','MYOCD','LAMA2','LGALS1','C7'],
 'Pericyte':['CD36','STEAP4','PLXDC1','COL4A1','PRKG1','COX4I2','RGS5','PDGFRA','PI15'],
 'T_cell':['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2','TRDC','TRGC1','TRGC2','GIMAP7','BACH2','THEMIS','SLC4A10','KLRB1','PDE4C','LRRC23'],
 'NK_cell':['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2','TRDC','TRGC1','TRGC2','FGFBP2','RAP1GAP2','LDB2','CXCR4','IL12RB2'],
 'Dendritic':['CD1C','FCER1A','CLEC10A','CLEC9A','XCR1','AXL','SIGLEC6','LILRA4','CLEC4C','FCN1','C1QC','MSR1','KYNU','LAMP3','CCR7'],
 'Monocyte':['FCGR3A','CX3CR1','LYPD2','MERTK','C1QC','FCN1','DOCK4','SERPINB2','CXCL10','S100A12'],
 'B_cell':['AICDA','BCL6','RGS13','SERPINA9','TCL1A','IGHD','BANK1','NR4A1','S100A4','CD27','PDE4C','LRRC23'],
 'Plasma':['MZB1','JCHAIN','XBP1','DERL3','FKBP11','ELL2','CREB3L2','SEL1L3','ALAS2','AHSP','SLC4A1','HBB','HBA1','HBA2','MKI67','TOP2A','CD3D','TRAC','CDKN1A','RGS1','IGHA1','IGHG1','IGHM'],
 'Mast':['TPSAB1','TPSB2','KIT','HDC','HPGDS','P2RY14','PRG2','CTSG','CMA1','CPA3','NFKB1','NFKBIA','ABCC4','PPARD'],
 'Neutrophil':['FCGR3B','CSF3R','CXCR2','CXCR4','S100A12','CXCL8','C1QC','CSF1R','DOCK4','IRAK2','SORL1','MPO','ELANE'],
}

def require(ok, message):
    if not bool(ok): raise RuntimeError(message)
def now(): return datetime.now().strftime('%Y%m%d_%H%M%S_%f')
def read_json(p): return json.loads(Path(p).read_text(encoding='utf-8'))
def write_json(p, x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_name(p.name+'.tmp')
    q.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str),encoding='utf-8');q.replace(p)
def sha_file(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(16*1024*1024),b''):h.update(b)
    return h.hexdigest()
def io_funcs():
    try: from anndata.io import read_elem, write_elem, sparse_dataset
    except ImportError:
        from anndata._io.specs import read_elem, write_elem
        from anndata.experimental import sparse_dataset
    return read_elem,write_elem,sparse_dataset

def _value_json(v):
    import numpy as np
    if isinstance(v,bytes):return v.decode('utf-8')
    if isinstance(v,np.ndarray):return [_value_json(x) for x in v.tolist()]
    if isinstance(v,np.generic):return v.item()
    return v

def hash_h5_tree(node):
    """Logical full-data hash, chunked. No loading of a full expression matrix."""
    import h5py, numpy as np
    h=hashlib.sha256()
    def walk(n, rel):
        h.update(rel.encode());h.update(json.dumps({k:_value_json(v) for k,v in n.attrs.items()},sort_keys=True,default=str).encode())
        if isinstance(n,h5py.Dataset):
            require(not n.is_virtual,'Virtual dataset is not a self-contained snapshot.')
            h.update(str(n.dtype).encode());h.update(str(n.shape).encode())
            if n.shape==(): chunks=[n[()]]
            else:
                step=max(1,int(16*1024*1024/max(1,n.dtype.itemsize*int(np.prod(n.shape[1:])))))
                chunks=(n[i:i+step] for i in range(0,len(n),step))
            for x in chunks:
                a=np.asarray(x)
                if a.dtype.kind in 'OUS':
                    for v in a.ravel():
                        b=str(_value_json(v)).encode('utf8');h.update(len(b).to_bytes(8,'little'));h.update(b)
                else:h.update(memoryview(np.ascontiguousarray(a)).cast('B'))
        else:
            for k in sorted(n):
                require(isinstance(n.get(k,getlink=True),h5py.HardLink),'External/soft links require separate inspection.')
                walk(n[k],rel+'/'+k)
    walk(node,'');return h.hexdigest()

def protected_h5(p):
    import h5py
    with h5py.File(p,'r') as f:
        return {k:hash_h5_tree(f[k]) for k in ('X','layers','var','varm','varp','obsm','obsp','raw') if k in f}

def load_metadata(p):
    import pandas as pd
    m=pd.read_csv(p,encoding='utf-8-sig',dtype=str,keep_default_na=False)
    require('sample_uid' in m and m.sample_uid.is_unique,'Metadata needs unique sample_uid.')
    require(len(m)==247,'Metadata must have 247 sample rows.')
    a=pd.to_numeric(m.age.replace('',float('nan')),errors='raise')
    require(a.dropna().between(0,120).all(),'Invalid age value.')
    require(m.sex.isin(['','male','female']).all(),'Unexpected sex encoding.')
    require(not m.donor_uid.str.contains('PROVISIONAL|not_available|unknown',case=False,regex=True).any(),
            'Unresolved donor aliases must not be promoted into verified donor_uid.')
    return m

DROP_META_RE = re.compile(r'(^donor_(?:identity_resolved|resolution_status|evidence|evidence_class|status)$|'
 r'^(?:metadata_confidence|metadata_sources|verification_note|verification_reference|mapping_limitations|'
 r'replicate_type|replicate_evidence_status|technical_audit_note|independent_unit_for_inference|'
 r'primary_inference_constraint|downstream_repeat_policy|clinical_context_scope|analysis_capture_uid_scope)$|'
 r'^verified_|^donor_uid(?:_v2|_curated|_analysis)$|__status$|__values_counts_json$)',re.I)

def overlay_metadata(obs, meta):
    """Values from verified/provided samples only. NA stays NA; no mean imputation."""
    import pandas as pd,numpy as np
    require('sample_uid' in obs and obs.sample_uid.notna().all(),'Missing cell-to-sample linkage.')
    require(set(obs.sample_uid.astype(str))<=set(meta.sample_uid),'Unmatched sample in source.')
    cols=[c for c in obs if DROP_META_RE.search(str(c)) or (str(c).lower().startswith('donor') and c not in meta.columns)]
    out=obs.drop(columns=cols).copy(); mt=meta.set_index('sample_uid')
    old_ids=obs.index.copy();sample=obs.sample_uid.astype(str)
    for col in meta.columns:
        if col=='sample_uid':continue
        v=sample.map(mt[col])
        if col=='age':out[col]=pd.to_numeric(v.replace('',np.nan),errors='raise').to_numpy(dtype=float)
        else:out[col]=pd.Categorical(v.replace('',np.nan).to_numpy())
    require(out.index.equals(old_ids),'Metadata join altered cell order.')
    return out,cols

def resolve_paths(contracts, root, mac_run=None, other_run=None):
    paths={}
    for b,c in contracts.items():
        exact=Path(c['analysis_file'])
        if b=='Macrophage' and mac_run is not None:run=Path(mac_run)
        elif b!='Macrophage' and other_run is not None:run=Path(other_run)/'lineages'/b
        else:
            # Remap ONLY the known root prefix, never pick the most recent unrelated file.
            rel=exact.relative_to(ROOT);run=(Path(root)/rel).parent.parent
        analysis=run/'data'/exact.name
        mp=run/'tables/cell_cluster_mapping.tsv.gz'
        require(analysis.is_file() and mp.is_file(),f'Missing {b} source. Supply the exact Macrophage/Other run directory: {run}')
        require(sha_file(mp)==c['mapping_sha256'],f'{b} membership differs from reviewed input.')
        res=read_json(run/'RESULTS.json')
        require(res.get('analysis_sha256')==c['analysis_sha256'],f'{b} analysis manifest differs.')
        paths[b]={'run':str(run),'analysis':str(analysis),'mapping':str(mp)}
    return paths

def build_label_vector(obs, decisions, contracts, paths):
    import numpy as np,pandas as pd
    require(len(obs)==1175115 and obs.index.is_unique,'Unexpected mother cell IDs/count.')
    require(L1 in obs and obs[L1].notna().all(),'Frozen level1_annotation missing.')
    require(obs[L1].astype(str).value_counts().to_dict()==EXPECTED_L1,'Frozen L1 composition changed.')
    labels=pd.Series('Pending',index=obs.index,dtype=object)
    labels.loc[obs[L1].astype(str).eq('Unknown')]='Unknown'
    labels.loc[obs[L1].astype(str).eq('Schwann cell')]='Not_assessed'
    allseen=np.zeros(len(obs),dtype=bool);audit=[]
    for b,c in contracts.items():
        mp=pd.read_csv(paths[b]['mapping'],sep='\t',dtype=str,keep_default_na=False)
        require(len(mp)==c['n_cells'] and mp.cell_id.is_unique,f'{b}: wrong member count.')
        require(mp.cell_id.eq(mp.parent_cell_id).all(),f'{b}: parent cell linkage differs.')
        ix=obs.index.get_indexer(mp.cell_id);require((ix>=0).all(),f'{b}: cell absent from mother.')
        require(not allseen[ix].any(),'Cell appears in two L1 branches.');allseen[ix]=True
        require(obs[L1].iloc[ix].astype(str).eq(c['lineage']).all(),f'{b}: parent L1 mismatch.')
        require(np.array_equal(obs.sample_uid.iloc[ix].astype(str),mp.sample_uid.astype(str)),f'{b}: sample members changed.')
        dd=decisions.loc[decisions.branch.eq(b)]
        require(dd.cluster_uid.is_unique,f'{b}: duplicate decisions.')
        require(mp[c['selected_key']].value_counts().to_dict()==dict(zip(dd.cluster_uid,dd.n_cells.astype(int))),f'{b}: reviewed clusters/counts mismatch.')
        rowmap=dd.set_index('cluster_uid');v=mp[c['selected_key']].map(rowmap[L2])
        require(v.notna().all(),f'{b}: missing decision.')
        labels.iloc[ix]=v.to_numpy()
        audit.append(pd.DataFrame({'cell_id':mp.cell_id,'level1_annotation':c['lineage'],
            'source_branch':b,'cluster_key':c['selected_key'],'source_cluster_uid':mp[c['selected_key']],
            'working_L2_annotation':v,'decision_action':mp[c['selected_key']].map(rowmap.action)}))
    require(int(allseen.sum())==1114365,'Wrong total analysed cells.')
    if L2 in obs:
        current=obs[L2].astype(object)
        conflict=current.notna() & ~current.isin(['','Pending']) & current.ne(labels)
        require(not conflict.any(),'Latest source contains conflicting newer L2 labels; reconcile explicitly instead of overwriting.')
    return labels,pd.concat(audit,ignore_index=True)

def readonly(p):
    p=Path(p);p.chmod(stat.S_IMODE(p.stat().st_mode)&~0o222)

def save_updated_h5(source, dest, new_obs, run, expected_protected):
    import h5py,pandas as pd,numpy as np
    read_elem,write_elem,_=io_funcs()
    require(not dest.exists(),f'Refusing to overwrite {dest}')
    stage=dest.with_name(dest.stem+'.partial.h5ad')
    require(not stage.exists(),f'Incomplete file already exists: {stage}')
    require(shutil.disk_usage(dest.parent).free>source.stat().st_size+8*1024**3,'Insufficient disk space for full new h5ad.')
    shutil.copyfile(source,stage)
    with h5py.File(stage,'r+') as f:
        old_obs=read_elem(f['obs'])
        require(old_obs.index.equals(new_obs.index),'Cell order changed before write.')
        require(np.array_equal(old_obs[L1].astype(str),new_obs[L1].astype(str)),'L1 must remain identical.')
        del f['obs'];write_elem(f,'obs',new_obs,dataset_kwargs={'compression':'lzf'})
        u=f.require_group('uns')
        vals={'L2_started':True,'L2_frozen':False,'L2_annotation_key':L2,
          'L2_status':'PARTIAL_LABELS_APPLIED_TARGETED_REVIEW_PENDING',
          'level1_frozen':True,'L1_frozen':True,'active_annotation_key':L1,
          'metadata_version':'v4_evidence_preserving_clean','metadata_primary_donor_key':'donor_uid',
          'L2_followup_directory':str(run),'current_full_h5ad':str(dest),
          'metadata_evidence_archive':str(run/'audit'),'L2_approval_table':str(run/'inputs/approved_L2_mapping.tsv'),
          'L2_partial_semantics':{'Pending':'state not yet assigned','Unknown':'L1 Unknown or approved L2 identity exclusion',
                                 'Not_assessed':'Schwann excluded from L2; L1 remains Schwann cell'}}
        for k,v in vals.items():
            if k in u:del u[k]
            write_elem(u,k,v)
        f.flush()
    require(protected_h5(stage)==expected_protected,'New file expression/coordinates/graphs differ from source.')
    with h5py.File(stage,'r') as f:
        obs_saved=read_elem(f['obs'])
        pd.testing.assert_frame_equal(obs_saved,new_obs,check_dtype=False,check_categorical=False)
    stage.rename(dest);readonly(dest)
    return sha_file(dest)

def prepare(root=ROOT, source=None, mac_run=None, other_run=None, resume_run=None, threads=40):
    """Save one new complete metadata+partial-L2 checkpoint. Does not touch kernel adata."""
    import numpy as np,pandas as pd,h5py
    root=Path(root);src=Path(source) if source else root/'panvascular_CPU_theta2_level1_frozen_v1.h5ad'
    require(src.is_file(),f'Mother file not found: {src}')
    contracts=read_json(BUNDLE/'inputs/source_contracts.json')
    for rel,h in read_json(BUNDLE/'BUNDLE_CHECKSUMS.json').items():
        require(sha_file(BUNDLE/rel)==h,f'Bundle file changed: {rel}')
    if resume_run:
        run=Path(resume_run); st=read_json(run/'RUN.json')
        require(st['source']==str(src.resolve()),'Resume source differs.')
        dest=Path(st['output_h5ad'])
        require(dest.exists() and sha_file(dest)==st['output_sha256'],'Resume snapshot missing or changed.')
        return run,dest
    paths=resolve_paths(contracts,root,mac_run,other_run)
    source_sha=sha_file(src);source_stat=(src.stat().st_size,src.stat().st_mtime_ns)
    print('Hashing protected expression/coordinates once ...',flush=True)
    protected=protected_h5(src)
    original=root/'panvascular_CPU_theta2_level1_frozen_v1.h5ad'
    expected_source=next(iter(contracts.values()))['source_sha256']
    if source_sha!=expected_source:
        require(src.resolve()!=original.resolve(),'The original frozen mother SHA256 has changed.')
        require(original.exists() and sha_file(original)==expected_source,'Cannot verify the original frozen ancestor.')
        require(protected==protected_h5(original),'Selected newer mother is not a metadata-only descendant of the frozen L1.')
    read_elem,_,_=io_funcs()
    with h5py.File(src,'r') as f:
        obs=read_elem(f['obs']);var=read_elem(f['var']);uns=read_elem(f['uns'])
    require(len(var)==37336 and var.index.is_unique,'Wrong gene set.')
    require(bool(uns.get('level1_frozen',uns.get('L1_frozen',False))),'Source L1 is not frozen.')
    dec=pd.read_csv(BUNDLE/'inputs/working_L2_decisions.tsv',sep='\t',keep_default_na=False)
    labels,audit=build_label_vector(obs,dec,contracts,paths)
    meta=load_metadata(BUNDLE/'inputs/metadata_h5ad_patch_v4.csv')
    expected_counts=pd.read_csv(BUNDLE/'inputs/sample_counts_expected.csv',encoding='utf-8-sig').set_index('sample_uid').n_cells
    require(obs.sample_uid.astype(str).value_counts().to_dict()==expected_counts.to_dict(),'Metadata sample sizes/membership are not for this atlas.')
    new_obs,removed=overlay_metadata(obs,meta);new_obs[L2]=pd.Categorical(labels.to_numpy())
    require(new_obs[L1].astype(str).equals(obs[L1].astype(str)),'Metadata changed L1.')
    run=root/('L2_targeted_followup_metadata_v4_'+now());run.mkdir(parents=True)
    for sub in ('inputs','audit','tables','branches','reviews'): (run/sub).mkdir()
    shutil.copytree(BUNDLE/'inputs',run/'inputs',dirs_exist_ok=True)
    shutil.copytree(BUNDLE/'audit',run/'audit',dirs_exist_ok=True)
    shutil.copy2(__file__,run/'audit/executed_driver.py')
    # Full old observation metadata are external, not discarded silently.
    obs.to_csv(run/'audit/obs_before_metadata_update.tsv.gz',sep='\t',index_label='cell_id')
    write_json(run/'audit/removed_observation_columns.json',removed)
    write_json(run/'audit/protected_source_datasets.json',protected)
    write_json(run/'inputs/resolved_paths.json',paths)
    audit.to_csv(run/'tables/L2_cell_mapping_audit.tsv.gz',sep='\t',index=False)
    new_obs[['sample_uid',L1,L2]].to_csv(run/'tables/all_cells_L2_partial.tsv.gz',sep='\t',index_label='cell_id')
    composition=new_obs.groupby([L1,L2],observed=True).size().rename('n_cells').reset_index()
    composition.to_csv(run/'tables/L2_working_composition.tsv',sep='\t',index=False)
    dest=root/('panvascular_CPU_level1_L2_partial_metadata_v4_'+now()+'.h5ad')
    print('Saving full new h5ad (raw/log1p/coordinates unchanged):',dest,flush=True)
    output_sha=save_updated_h5(src,dest,new_obs,run,protected)
    require(source_stat==(src.stat().st_size,src.stat().st_mtime_ns),'Source changed during snapshot creation.')
    record={'code_version':CODE_VERSION,'source':str(src.resolve()),'source_sha256':source_sha,
      'output_h5ad':str(dest),'output_sha256':output_sha,'L1_frozen':True,'L2_frozen':False,
      'L1_changes':0,'no_uncertain_donor_ids_fabricated':True,'threads':threads,'created_at':datetime.now().isoformat(),
      'status':'PARTIAL_L2_AND_METADATA_SAVED_FOLLOWUP_PENDING','protected_datasets':protected,
      'clinical_missing_counts':read_json(BUNDLE/'inputs/expected_summary.json')}
    write_json(run/'RUN.json',record);write_json(dest.with_suffix('.manifest.json'),record)
    dest.with_suffix('.h5ad.sha256').write_text(output_sha+'  '+dest.name+'\n')
    (run/'PROJECT_STATE_L2_WORKING.md').write_text('# L2 working snapshot\n\nL1 frozen = TRUE; L2 frozen = FALSE.\n\n'
      f'Full h5ad: {dest}\nSHA256: {output_sha}\n\nOnly approved existing UIDs have working names/Unknown. '
      'Pending remains distinct. Follow-up produces evidence, not automatic new labels. '
      'Clinical missing values are not invented; evidence notes are outside obs.\n',encoding='utf8')
    print(composition.to_string(index=False),flush=True)
    del obs,new_obs,audit;gc.collect();return run,dest

# -------------------- Follow-up-only scientific helpers --------------------
def clean_sparse(x):
    import scipy.sparse as sp
    z=sp.csr_matrix(x,copy=True);z.sum_duplicates();z.eliminate_zeros();z.sort_indices();return z

def aggregated(x, labels):
    import numpy as np,pandas as pd,scipy.sparse as sp
    lab=np.asarray(labels,dtype=str);order=list(dict.fromkeys(lab));cat=pd.Categorical(lab,categories=order)
    ns=np.bincount(cat.codes,minlength=len(order));a=sp.csr_matrix((np.ones(len(lab)),(cat.codes,np.arange(len(lab)))),shape=(len(order),len(lab)))
    return (a@x).toarray()/ns[:,None],ns,order

def marker_values(raw, log, obs, genes, varnames, key):
    import numpy as np,pandas as pd
    ii=varnames.get_indexer(genes);require((ii>=0).all(),'Missing requested marker.')
    r=clean_sparse(raw[:,ii]);l=clean_sparse(log[:,ii]);bi=r.copy();bi.data=np.ones(bi.nnz)
    mean,n,order=aggregated(l,obs[key]);frac,_,_=aggregated(bi,obs[key]);rm,_,_=aggregated(r,obs[key])
    return pd.DataFrame(mean,index=order,columns=genes),pd.DataFrame(frac,index=order,columns=genes),pd.DataFrame(rm,index=order,columns=genes),n

def save_dotplot(mean,frac,name,out):
    import numpy as np,pandas as pd,anndata as ad,scanpy as sc,matplotlib.pyplot as plt
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    for start in range(0,mean.shape[1],48):
        m=mean.iloc[:,start:start+48];f=frac.reindex(index=m.index,columns=m.columns)
        mock=ad.AnnData(X=np.zeros(m.shape,dtype=np.float32),obs=pd.DataFrame({'group':pd.Categorical(m.index,categories=list(m.index))},index=m.index),var=pd.DataFrame(index=m.columns))
        for rel in (False,True):
            colors=m.copy()
            if rel:
                den=colors.max()-colors.min();colors=(colors-colors.min()).div(den.replace(0,1),axis=1)
            dp=sc.pl.DotPlot(mock,var_names=list(m.columns),groupby='group',dot_color_df=colors,dot_size_df=f,
                categories_order=list(m.index),standard_scale=None,figsize=(max(9,len(m.columns)*.21),max(4,len(m)*.27+2)))
            dp=dp.style(dot_max=1.0,dot_min=0)
            dp.legend(colorbar_title='Gene-wise relative mean' if rel else 'Mean log1p(CP10k)',size_title='Raw counts > 0 (%)')
            dp.make_figure();fig=dp.get_axes()['mainplot_ax'].figure
            stem=f'{name}_{"relative" if rel else "mean"}_{start//48+1:02d}'
            try:
                fig.savefig(out/(stem+'.png'),dpi=160,bbox_inches='tight')
                fig.savefig(out/(stem+'.pdf'),dpi=200,bbox_inches='tight')
            finally:plt.close(fig)


def existing_umap(xy,obs,key,out):
    import numpy as np,pandas as pd,anndata as ad,scanpy as sc,matplotlib.pyplot as plt
    order=np.random.default_rng(0).permutation(len(obs))
    a=ad.AnnData(obs=pd.DataFrame({'label':pd.Categorical(obs[key].astype(str).to_numpy()[order])},index=obs.index[order]),obsm={'X_plot':xy[order].copy()})
    fig,ax=plt.subplots(figsize=(10,8))
    try:
        sc.pl.embedding(a,basis='X_plot',color='label',show=False,save=False,ax=ax,frameon=False,size=1,alpha=1,title='Existing embedding; '+key)
        ax.set_aspect('equal',adjustable='box')
        for c in ax.collections:c.set_linewidths(0);c.set_edgecolor('none');c.set_rasterized(True)
        fig.savefig(out,dpi=180,bbox_inches='tight')
    finally:plt.close(fig)


def source_summary(obs,key,out):
    import pandas as pd,numpy as np
    out=Path(out);out.mkdir(parents=True,exist_ok=True);rows=[]
    for uid,g in obs.groupby(key,observed=True,sort=False):
        r={'cluster_uid':uid,'n_cells':len(g)}
        for f in ('sample_uid','study','donor_uid','tissue_detail','vascular_bed','disease_context','comparison_group','assay_type'):
            if f not in g:continue
            v=g[f].dropna().astype(str);v=v[v.ne('')];vc=v.value_counts()
            r['n_'+f]=len(vc);r['missing_'+f]=len(g)-len(v)
            r['top_'+f]=str(vc.index[0]) if len(vc) else '';r['top_'+f+'_fraction']=float(vc.iloc[0]/len(g)) if len(vc) else np.nan
        for f in ('L2_total_raw_counts','L2_detected_genes','L2_pct_counts_mt'):
            if f in g:r['median_'+f]=float(g[f].median())
        rows.append(r)
    pd.DataFrame(rows).to_csv(out/'cluster_sources_updated_metadata.tsv',sep='\t',index=False)
    for f in ('sample_uid','study','assay_type','vascular_bed','donor_uid','disease_context','comparison_group'):
        if f in obs:pd.crosstab(obs[key].astype(str),obs[f].astype('string').fillna('Missing')).to_csv(out/f'cluster_by_{f}.tsv',sep='\t')
    depth=obs[['sample_uid',key,'L2_total_raw_counts','L2_detected_genes']].copy()
    depth['within_sample_depth_quintile']=depth.groupby('sample_uid',observed=True).L2_total_raw_counts.rank(pct=True).mul(5).sub(1e-8).astype(int).clip(0,4)+1
    pd.crosstab(depth[key],depth.within_sample_depth_quintile).to_csv(out/'cluster_by_within_sample_depth_quintile.tsv',sep='\t')
    return depth


def choose_pairs(obs,z,c):
    import numpy as np,pandas as pd
    key=c['selected_key'];pairs={tuple(sorted(p)):'PREDEFINED_REVIEW_COMPARISON' for p in c['manual_pairs']}
    available=set(obs[key].astype(str));require(all(a in available and b in available for a,b in pairs),'Missing predefined cluster.')
    cent=pd.DataFrame(z).groupby(obs[key].astype(str).to_numpy()).mean()
    for uid in c['pending_uids']:
        # A peer selection aid ONLY, not a classifier. No new label follows this choice.
        candidates=[x for x in c['named_uids'] if x!=uid]
        if not candidates:continue
        d=((cent.loc[candidates]-cent.loc[uid])**2).sum(axis=1)
        peer=str(d.idxmin());pairs.setdefault(tuple(sorted((uid,peer))),'NEAREST_NAMED_CENTROID_FOR_REVIEW_ONLY')
    return [{'key':key,'a':a,'b':b,'selection':v} for (a,b),v in sorted(pairs.items())]+[
      {'key':'lineage_L2_r1p2','a':p[0],'b':p[1],'selection':'EXISTING_HIGHER_RESOLUTION_REVIEW_ONLY'} for p in c.get('highres_pairs',[])]


def paired_sample_markers(raw,log,obs,varnames,key,a,b,genes,out):
    """Within-sample descriptive differences, no pseudoreplicated p-values."""
    import numpy as np,pandas as pd
    sub=obs[key].astype(str).isin([a,b]).to_numpy();o=obs.loc[sub].copy()
    gi=varnames.get_indexer(genes);r=raw[:,gi][sub,:];l=log[:,gi][sub,:];o['_stratum']=o.sample_uid.astype(str)+'|||'+o[key].astype(str)
    m,f,_,n=marker_values(r,l,o,genes,pd.Index(genes),'_stratum')
    samples=[];summary=[]
    for sid in o.sample_uid.astype(str).unique():
        sa=sid+'|||'+a;sb=sid+'|||'+b
        if sa not in m.index or sb not in m.index:continue
        na=int((o['_stratum']==sa).sum());nb=int((o['_stratum']==sb).sum())
        if min(na,nb)<20:continue
        ga=o.loc[o['_stratum'].eq(sa)];gb=o.loc[o['_stratum'].eq(sb)]
        vals=(m.loc[sa]-m.loc[sb]).to_numpy()
        samples.extend({'sample_uid':sid,'gene':gene,'n_a':na,'n_b':nb,'mean_a':m.loc[sa,gene],
           'mean_b':m.loc[sb,gene],'mean_difference':diff,'fraction_a':f.loc[sa,gene],'fraction_b':f.loc[sb,gene],
           'median_raw_counts_a':ga.L2_total_raw_counts.median(),'median_raw_counts_b':gb.L2_total_raw_counts.median(),
           'donor_uid':str(ga.donor_uid.dropna().iloc[0]) if ga.donor_uid.notna().any() else '',
           'assay_type':str(ga.assay_type.dropna().iloc[0]) if ga.assay_type.notna().any() else ''}
           for gene,diff in zip(genes,vals))
    df=pd.DataFrame(samples)
    df.to_csv(out/'within_sample_marker_differences.tsv.gz',sep='\t',index=False)
    if len(df):
        for gene,g in df.groupby('gene',sort=False):
            summary.append({'gene':gene,'n_samples':len(g),'n_samples_a_higher':int((g.mean_difference>0).sum()),
               'median_difference':float(g.mean_difference.median()),'fraction_a_higher':float((g.mean_difference>0).mean()),
               'interpretation':'descriptive only; samples can share donors'})
    pd.DataFrame(summary).to_csv(out/'within_sample_direction_summary.tsv',sep='\t',index=False)


def pair_de(raw,log,obs,varnames,spec,review_genes,out):
    import numpy as np,pandas as pd,scanpy as sc,anndata as ad
    from threadpoolctl import threadpool_limits
    a,b,key=spec['a'],spec['b'],spec['key'];folder=out/(a+'__VS__'+b);folder.mkdir(exist_ok=True,parents=True)
    signature=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()
    if (folder/'DONE.json').exists():
        st=read_json(folder/'DONE.json');require(st['signature']==signature,'Pair cache mismatch.');return st
    mask=obs[key].astype(str).isin([a,b]).to_numpy(); o=obs.loc[mask].copy()
    n=o[key].astype(str).value_counts();require(set(n.index)=={a,b},'Invalid pair.')
    if min(n)<2:
        st={'signature':signature,'status':'SKIPPED_TOO_FEW_CELLS','a':a,'b':b};write_json(folder/'DONE.json',st);return st
    work=ad.AnnData(X=log[mask,:].copy(),obs=pd.DataFrame({'pair':pd.Categorical(o[key].astype(str))},index=o.index),var=pd.DataFrame(index=varnames))
    work.uns['log1p']={'base':None}
    gmask=np.asarray((work.X>0).sum(axis=0)).ravel()>=3
    require(gmask.any(),'No genes eligible in pair.')
    sc.tl.rank_genes_groups(work,groupby='pair',groups=[a],reference=b,use_raw=False,
        method='wilcoxon',tie_correct=True,mask_var=gmask,n_genes=int(gmask.sum()),
        pts=True,rankby_abs=False,corr_method='benjamini-hochberg')
    rank=sc.get.rank_genes_groups_df(work,group=a)
    ii=varnames.get_indexer(rank.names)
    aa=raw[obs[key].astype(str).eq(a).to_numpy(),:]
    rank['pct_a']=(np.asarray(aa.getnnz(axis=0)).ravel()/aa.shape[0])[ii];del aa
    bb=raw[obs[key].astype(str).eq(b).to_numpy(),:]
    rank['pct_b']=(np.asarray(bb.getnnz(axis=0)).ravel()/bb.shape[0])[ii];del bb
    rank['pct_difference']=rank.pct_a-rank.pct_b;rank.insert(0,'cluster_b',b);rank.insert(0,'cluster_a',a)
    rank.to_csv(folder/'FindMarkers_pair_all_genes.tsv.gz',sep='\t',index=False)
    topa=rank.sort_values('scores',ascending=False).head(100)
    reverse=rank.copy();reverse['cluster_a']=b;reverse['cluster_b']=a
    reverse['scores']=-rank['scores'];reverse['logfoldchanges']=-rank['logfoldchanges']
    reverse['pct_a']=rank['pct_b'];reverse['pct_b']=rank['pct_a'];reverse['pct_difference']=-rank['pct_difference']
    topb=reverse.sort_values('scores',ascending=False).head(100)
    topa.to_csv(folder/'top100_a_vs_b.tsv',sep='\t',index=False);topb.to_csv(folder/'top100_b_vs_a.tsv',sep='\t',index=False)
    # Add evidence-driven anchors; never use these choices to assign labels.
    choose=rank.loc[~rank.names.str.match(r'^(MT-|RPL|RPS)|^(MALAT1|NEAT1)$')]
    anchors=list(dict.fromkeys(list(choose.nlargest(6,'scores').names)+list(choose.nsmallest(6,'scores').names)))
    display=list(dict.fromkeys([g for g in review_genes if g in varnames]+anchors))
    gi=varnames.get_indexer(display)
    m,f,rm,ns=marker_values(raw[:,gi][mask,:],log[:,gi][mask,:],o,display,pd.Index(display),key)
    m.to_csv(folder/'marker_mean_log1p.tsv',sep='\t');f.to_csv(folder/'marker_detection_fraction.tsv',sep='\t')
    # Compact pair figure uses discovered anchors, not a 400-gene wide plot.
    save_dotplot(m[anchors],f[anchors],'pair_markers',folder/'figures')
    del work;gc.collect()
    paired_sample_markers(raw,log,obs,varnames,key,a,b,display,folder)
    st={'status':'DONE','signature':signature,'a':a,'b':b,'n_a':int(n[a]),'n_b':int(n[b]),
       'n_tested_genes':int(gmask.sum()),'interpretation':'exploratory cluster-marker contrast; not donor-level DE'}
    write_json(folder/'DONE.json',st)
    gc.collect();return st


def joint_programs(raw,obs,varnames,key,panels,out):
    import numpy as np,pandas as pd
    records={}
    for title,genes in panels.items():
        genes=[g for g in genes if g in varnames]
        if not genes:continue
        idx=varnames.get_indexer(genes);r=raw[:,idx]
        records[title]=np.asarray((r>0).sum(axis=1)).ravel().astype('int16')
    pc=pd.DataFrame(records,index=obs.index); pc.insert(0,'cell_id',obs.index);pc.insert(1,'cluster_uid',obs[key].astype(str).to_numpy())
    pc.to_csv(out/'per_cell_program_detection_counts.tsv.gz',sep='\t',index=False)
    groups=list(records);summ=[]
    for uid,ix in obs.groupby(key,observed=True).indices.items():
        for g in groups:
            for threshold in [1,2,3]:summ.append({'cluster_uid':str(uid),'program':g,'n_cells':len(ix),'min_genes':threshold,'fraction':float((records[g][ix]>=threshold).mean())})
    pd.DataFrame(summ).to_csv(out/'program_multigene_detection.tsv',sep='\t',index=False)
    pairs=[('CD3D','TRAC'),('CD3E','TRAC'),('MYH11','DCN'),('MYH11','VWF'),('C1QC','DCN'),('C1QC','CD1C'),('MZB1','ALAS2'),('MZB1','CD3D'),('RGS5','MYH11')]
    j=[]
    for ga,gb in pairs:
        if ga not in varnames or gb not in varnames:continue
        a=np.asarray(raw[:,varnames.get_loc(ga)].toarray()).ravel();b=np.asarray(raw[:,varnames.get_loc(gb)].toarray()).ravel()
        for uid,ix in obs.groupby(key,observed=True).indices.items():
            for cut in [1,2]:j.append({'cluster_uid':str(uid),'gene_a':ga,'gene_b':gb,'raw_cutoff':cut,'n_cells':len(ix),
              'frac_a':float((a[ix]>=cut).mean()),'frac_b':float((b[ix]>=cut).mean()),'frac_joint':float(((a[ix]>=cut)&(b[ix]>=cut)).mean())})
    pd.DataFrame(j).to_csv(out/'direct_same_cell_joint_detection.tsv',sep='\t',index=False)


def branch_worker(run,branch,threads):
    require(os.environ.get('CUDA_VISIBLE_DEVICES')=='','Worker GPU must be hidden before startup.')
    import numpy as np,pandas as pd,h5py,scanpy as sc,anndata as ad,scipy.sparse as sp
    import matplotlib;matplotlib.use('Agg')
    from threadpoolctl import threadpool_limits
    import numba
    require(sc.__version__=='1.11.5','This continuation was designed for the audited Scanpy 1.11.5 environment; do not silently change versions.')
    numba.set_num_threads(min(threads,numba.config.NUMBA_NUM_THREADS));sc.settings.n_jobs=threads;sc.settings.autoshow=False
    run=Path(run);c=read_json(run/'inputs/source_contracts.json')[branch];path=read_json(run/'inputs/resolved_paths.json')[branch]
    out=run/'branches'/branch;out.mkdir(exist_ok=True,parents=True)
    status=lambda s:write_json(out/'STATUS.json',{'status':s,'branch':branch,'time':datetime.now().isoformat()})
    status('VERIFYING_EXISTING_BRANCH');require(sha_file(path['analysis'])==c['analysis_sha256'],'Source branch h5ad changed.')
    read_elem,_,sparse_dataset=io_funcs()
    with h5py.File(path['analysis'],'r') as f:
        obs=read_elem(f['obs']);var=read_elem(f['var'])
        raw=read_elem(f['X']);log=read_elem(f['layers/log1p'])
        hkey='X_harmony_mac_L2' if branch=='Macrophage' else 'X_harmony_lineage_L2'
        ukey='X_umap_mac_L2' if branch=='Macrophage' else 'X_umap_lineage_L2'
        require(hkey in f['obsm'] and ukey in f['obsm'],'Expected stored embedding keys missing.')
        z=read_elem(f['obsm'][hkey]);xy=read_elem(f['obsm'][ukey])
    require(sp.issparse(raw) and sp.issparse(log),'Expected sparse raw/log1p.')
    require(len(obs)==c['n_cells'] and obs[L1].astype(str).eq(c['lineage']).all(),'Branch L1 or membership changed.')
    mp=pd.read_csv(path['mapping'],sep='\t',dtype=str,keep_default_na=False)
    require(np.array_equal(obs.index.astype(str),mp.cell_id),'Branch matrix membership differs.')
    for key in c['existing_keys']:require(np.array_equal(obs[key].astype(str),mp[key]),'Stored cluster changed.')
    meta=load_metadata(run/'inputs/metadata_h5ad_patch_v4.csv');obs,_=overlay_metadata(obs,meta)
    raw=clean_sparse(raw);log=clean_sparse(log)
    for x in (raw,log):
        require(np.isfinite(x.data).all() and (x.data>=0).all(),'Invalid expression values.')
        for a in (x.data,x.indices,x.indptr):a.flags.writeable=False
    totals=np.asarray(raw.sum(axis=1,dtype=np.float64)).ravel();obs['L2_total_raw_counts']=totals
    obs['L2_detected_genes']=raw.getnnz(axis=1)
    mt=var.index.str.upper().str.startswith('MT-');obs['L2_pct_counts_mt']=100*np.asarray(raw[:,mt].sum(axis=1)).ravel()/np.maximum(totals,1)
    dec=pd.read_csv(run/'inputs/working_L2_decisions.tsv',sep='\t',keep_default_na=False)
    d=dec.loc[dec.branch.eq(branch)].set_index('cluster_uid');obs[L2]=pd.Categorical(obs[c['selected_key']].astype(str).map(d[L2]))
    panels=read_json(run/f'inputs/{branch}_marker_panels.json')
    extras=list(EXTRA_MARKERS.get(branch,[]))
    extras+=list(itertools.chain.from_iterable(c['existing_top_anchors'].values()))
    extras+=list(itertools.chain.from_iterable(x.split(';') for x in d.marker_candidates))
    base=list(dict.fromkeys(itertools.chain.from_iterable(list(panels['state'].values())+list(panels['boundary'].values()))))
    genes=[g for g in dict.fromkeys(base+extras) if g and g in var.index]
    write_json(out/'marker_panel_used.json',{'state':panels['state'],'boundary':panels['boundary'],'additional_candidates':extras,'not_found':[g for g in base+extras if g and g not in var.index]})
    status('UPDATED_METADATA_AND_MARKER_REVIEW')
    with threadpool_limits(limits=threads):
        # Two available partitions reviewed separately; no merging of resolution assignments.
        keys=[c['selected_key']]
        if branch!='Macrophage':keys+=['lineage_L2_r1p2']
        for key in keys:
            t=out/key;t.mkdir(exist_ok=True,parents=True);m,f,rm,n=marker_values(raw,log,obs,genes,var.index,key)
            for df,fn in [(m,'marker_mean_log1p.tsv'),(f,'marker_detection_fraction.tsv'),(rm,'marker_mean_raw_counts.tsv')]:df.rename_axis('cluster_uid').to_csv(t/fn,sep='\t')
            source_summary(obs,key,t);existing_umap(xy,obs,key,t/'existing_UMAP.png')
            for name,panel in [('state',panels['state']),('boundary',panels['boundary'])]:
                gg=[g for g in dict.fromkeys(itertools.chain.from_iterable(panel.values())) if g in m]
                save_dotplot(m[gg],f[gg],name,t/'figures')
            gg=[g for g in dict.fromkeys(extras) if g and g in m]
            save_dotplot(m[gg],f[gg],'additional_markers',t/'figures')
            joint_programs(raw,obs,var.index,key,panels['boundary'],t)
        # Previous L1 expression-only references are NOT mixed into clustering or DE.
        ref_file=run/'tables/L1_reference_mean_log1p.tsv'
        if ref_file.exists():
            ref=pd.read_csv(ref_file,sep='\t',index_col=0);rf=pd.read_csv(run/'tables/L1_reference_detection_fraction.tsv',sep='\t',index_col=0)
            m,f,_,_=marker_values(raw,log,obs,genes,var.index,c['selected_key'])
            common=[g for g in genes if g in ref.columns]
            refs=ref.reindex([x for x in c['reference_lineages'] if x in ref.index])[common];refs.index='REF_'+refs.index.astype(str)
            fr=rf.reindex(c['reference_lineages'])[common].dropna(how='all');fr.index='REF_'+fr.index.astype(str)
            bm=pd.concat([m[common],refs]);bf=pd.concat([f[common],fr])
            bm.to_csv(out/'target_plus_L1_reference_mean.tsv',sep='\t');bf.to_csv(out/'target_plus_L1_reference_detection.tsv',sep='\t')
            save_dotplot(bm,bf,'target_plus_L1_reference',out/'reference_figures')
        pairs=choose_pairs(obs,z,c);pd.DataFrame(pairs).to_csv(out/'planned_pairwise_comparisons.tsv',sep='\t',index=False)
        status('TARGETED_PAIRWISE_MARKERS')
        pair_results=[]
        for i,p in enumerate(pairs):
            print(f'[{branch}] pair {i+1}/{len(pairs)}: {p["a"]} vs {p["b"]}',flush=True)
            # classic+explicit anchors for paired sample audit; direct DE still tests whole gene set.
            rg=list(dict.fromkeys(base+EXTRA_MARKERS.get(branch,[])+c['existing_top_anchors'].get(p['a'],[])+c['existing_top_anchors'].get(p['b'],[])))
            rg=[g for g in rg if g in var.index]
            pair_results.append(pair_de(raw,log,obs,var.index,p,rg,out/'pairwise'))
    obs[['sample_uid',L1,L2]+c['existing_keys']].to_csv(out/'unchanged_membership_and_working_names.tsv.gz',sep='\t',index_label='cell_id')
    make_review_html(out,branch)
    info={'status':'PAUSED_FOR_TARGETED_L2_REVIEW','branch':branch,'n_cells':len(obs),
      'n_pairwise_comparisons':len(pair_results),'hvg_pca_harmony_neighbors_leiden_umap_recomputed':False,
      'new_identity_assignments_from_followup':False,'metadata_updated_for_descriptive_outputs':True,
      'source_h5ad_unchanged':sha_file(path['analysis'])==c['analysis_sha256']}
    require(info['source_h5ad_unchanged'],'Source local h5ad changed.')
    write_json(out/'DONE.json',info);status('PAUSED_FOR_TARGETED_L2_REVIEW');pack_review(out,run/'reviews'/branch)
    return info


def make_review_html(out,branch):
    title=f'{branch}: targeted L2 evidence; no automatic relabelling'
    page=['<!doctype html><meta charset="utf-8"><style>body{font:14px sans-serif;margin:24px}img{max-width:100%}table{border-collapse:collapse}td,th{padding:5px;border:1px solid #ccc}</style>',f'<h1>{html.escape(title)}</h1>',
      '<p>Existing coordinates and partitions reused. Pairwise tests are exploratory cluster-marker evidence, not independent donor-level tests. Pending is not Unknown. Metadata provenance remains outside the h5ad observation columns.</p>']
    for p in sorted(out.rglob('*.png')):
        page.extend(['<h3>'+html.escape(str(p.relative_to(out)))+'</h3>','<img src="'+html.escape(str(p.relative_to(out)))+'">'])
    (out/'REVIEW_INDEX.html').write_text('\n'.join(page),encoding='utf8')


def pack_review(out,prefix,max_mb=35):
    """Independent <=35MiB target ZIPs; never include full matrices or full-genome DE."""
    out=Path(out);prefix=Path(prefix);prefix.parent.mkdir(exist_ok=True,parents=True)
    files=[p for p in sorted(out.rglob('*')) if p.is_file() and p.suffix in ('.png','.json','.tsv','.html') and p.name!='FindMarkers_pair_all_genes.tsv.gz']
    # Key compressed summaries are small; omit per-cell expression to reduce upload size.
    files += [p for p in sorted(out.rglob('*.tsv.gz')) if p.name in ('within_sample_marker_differences.tsv.gz','unchanged_membership_and_working_names.tsv.gz')]
    part=1;z=None;written=[]
    try:
        for p in files:
            if z is None or (z.fp.tell()+p.stat().st_size>max_mb*1024**2 and len(z.filelist)>0):
                if z:z.close()
                name=prefix.parent/(prefix.name+f'_review_{part:02d}.zip');part+=1
                z=zipfile.ZipFile(name,'w',zipfile.ZIP_DEFLATED,compresslevel=6);written.append(str(name))
            z.write(p,arcname=str(p.relative_to(out)))
    finally:
        if z:z.close()
    return written


def prepare_l1_references(run):
    """Marker-only L1 reference summaries, read-only and blockwise; no reference DE."""
    import numpy as np,pandas as pd,h5py,scipy.sparse as sp
    run=Path(run);p=run/'tables/L1_reference_mean_log1p.tsv'
    if p.exists():return
    record=read_json(run/'RUN.json');src=Path(record['output_h5ad'])
    contracts=read_json(run/'inputs/source_contracts.json')
    genes=[]
    for b in contracts:
        pan=read_json(run/f'inputs/{b}_marker_panels.json')
        genes+=list(itertools.chain.from_iterable(pan['state'].values()))+list(itertools.chain.from_iterable(pan['boundary'].values()))+EXTRA_MARKERS.get(b,[])
    read_elem,_,sparse_dataset=io_funcs()
    with h5py.File(src,'r') as f:
        obs=read_elem(f['obs']);var=read_elem(f['var']);genes=[g for g in dict.fromkeys(genes) if g in var.index]
        idx=var.index.get_indexer(genes);labels=obs[L1].astype(str);order=list(EXPECTED_L1);code=pd.Categorical(labels,categories=order).codes
        sums=np.zeros((len(order),len(genes)));det=sums.copy();n=np.bincount(code,minlength=len(order))
        raw=sparse_dataset(f['X']);log=sparse_dataset(f['layers/log1p'])
        for lo in range(0,len(obs),10000):
            hi=min(lo+10000,len(obs));r=clean_sparse(raw[lo:hi,:][:,idx]);l=clean_sparse(log[lo:hi,:][:,idx]);r.data=np.ones(r.nnz)
            membership=sp.csr_matrix((np.ones(hi-lo),(code[lo:hi],np.arange(hi-lo))),shape=(len(order),hi-lo))
            sums+=(membership@l).toarray();det+=(membership@r).toarray()
    pd.DataFrame(sums/n[:,None],index=order,columns=genes).to_csv(p,sep='\t')
    pd.DataFrame(det/n[:,None],index=order,columns=genes).to_csv(run/'tables/L1_reference_detection_fraction.tsv',sep='\t')


def cpu_environment(threads):
    env=os.environ.copy();env['CUDA_VISIBLE_DEVICES']='';env['MPLBACKEND']='Agg';env['PYTHONUNBUFFERED']='1'
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:env[key]=str(threads)
    return env


def run_followups(run,threads=40,max_workers=2,min_available_gb=32):
    import pandas as pd
    run=Path(run);contracts=read_json(run/'inputs/source_contracts.json')
    per=max(1,threads//max_workers);pending=list(contracts);running={};states={}
    # One scheduler per run. An explicit foreground wait avoids silently abandoned jobs.
    lock=run/'SCHEDULER.lock'
    fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.write(fd,str(os.getpid()).encode());os.close(fd)
    try:
        print('Preparing marker-only L1 comparison references ...',flush=True);prepare_l1_references(run)
        while pending or running:
            while pending and len(running)<max_workers:
                try:
                    import psutil
                    free=psutil.virtual_memory().available/(1024**3)
                except ImportError:free=float('inf')
                if free<min_available_gb:
                    if not running:raise RuntimeError(f'Available RAM {free:.1f}GiB < {min_available_gb}GiB; checkpoint saved; free memory then resume.')
                    break
                b=pending.pop(0);out=run/'branches'/b;out.mkdir(parents=True,exist_ok=True)
                if (out/'DONE.json').exists():
                    states[b]=read_json(out/'DONE.json')['status'];continue
                log=(out/('worker_'+now()+'.log')).open('w')
                proc=subprocess.Popen([sys.executable,str(BUNDLE/'run_L2_followup_and_update.py'),'--worker',str(run),b,'--threads',str(per)],env=cpu_environment(per),stdout=log,stderr=subprocess.STDOUT)
                running[b]=(proc,log);states[b]='RUNNING';print(f'Started {b}, PID={proc.pid}, threads={per}',flush=True)
            for b,(p,log) in list(running.items()):
                rc=p.poll()
                if rc is not None:
                    log.close();states[b]='PAUSED_FOR_TARGETED_L2_REVIEW' if rc==0 else 'FAILED_CHECK_LOG'
                    del running[b];print(b,states[b],flush=True)
            pd.DataFrame([{'branch':b,'status':states.get(b,'QUEUED')} for b in contracts]).to_csv(run/'batch_status.csv',index=False)
            if running:time.sleep(3)
        status='PAUSED_FOR_GLOBAL_L2_TARGETED_REVIEW' if all(v=='PAUSED_FOR_TARGETED_L2_REVIEW' for v in states.values()) else 'PARTIAL_FAILURE_REVIEW_AVAILABLE'
        write_json(run/'FOLLOWUP_STATUS.json',{'status':status,'branches':states,'new_annotations_from_followup_applied':False})
        links=['<!doctype html><meta charset="utf-8"><h1>Targeted L2 follow-up</h1>',f'<p>{status}</p>']
        for b in contracts:links.append(f'<p><a href="branches/{b}/REVIEW_INDEX.html">{b}</a> — {states.get(b)}</p>')
        (run/'GLOBAL_TARGETED_REVIEW.html').write_text('\n'.join(links),encoding='utf8')
    finally:
        # Interrupt stops spawned processes, not existing Macrophage or earlier jobs.
        for p,log in running.values():
            if p.poll() is None:p.terminate()
            log.close()
        lock.unlink(missing_ok=True)
    return states


def run_from_notebook(*,root=ROOT,source=None,macrophage_run=None,other_run=None,
                      total_threads=40,max_workers=2,resume_run=None,run_targeted=True):
    """Single-cell launcher. No scientific imports/changes in caller notebook."""
    root=Path(root)
    request=root/('L2_followup_request_'+now()+'.json')
    write_json(request,{'root':str(root),'source':str(source) if source else None,
       'mac_run':str(macrophage_run) if macrophage_run else None,'other_run':str(other_run) if other_run else None,
       'resume_run':str(resume_run) if resume_run else None,'threads':total_threads,'workers':max_workers,'run_targeted':run_targeted})
    command=[sys.executable,str(BUNDLE/'run_L2_followup_and_update.py'),'--driver',str(request)]
    log=request.with_suffix('.log')
    with log.open('w') as handle:
        p=subprocess.Popen(command,env=cpu_environment(total_threads),stdout=handle,stderr=subprocess.STDOUT)
        print(f'Running CPU driver PID {p.pid}; log: {log}',flush=True)
        pos=0
        try:
            while p.poll() is None:
                time.sleep(3)
                with log.open(errors='replace') as f:f.seek(pos);txt=f.read();pos=f.tell()
                if txt:print(txt,end='',flush=True)
            with log.open(errors='replace') as f:f.seek(pos);txt=f.read()
            if txt:print(txt,end='',flush=True)
        except KeyboardInterrupt:
            print('Driver is a separate process; check its PID/log before relaunching.',flush=True);raise
    require(p.returncode==0,f'Driver failed; inspect {log}. Completed snapshot/branch outputs are retained.')
    record=read_json(request.with_suffix('.result.json'))
    return Path(record['run']),Path(record['h5ad'])


def cli():
    ap=argparse.ArgumentParser();ap.add_argument('--driver');ap.add_argument('--worker',nargs=2);ap.add_argument('--threads',type=int,default=20);a=ap.parse_args()
    try:
        if a.worker:branch_worker(a.worker[0],a.worker[1],a.threads)
        elif a.driver:
            cfg=read_json(a.driver)
            run,final=prepare(root=cfg['root'],source=cfg['source'],mac_run=cfg['mac_run'],other_run=cfg['other_run'],resume_run=cfg['resume_run'],threads=cfg['threads'])
            # Persist recovery route BEFORE lengthy follow-up computations.
            write_json(Path(a.driver).with_suffix('.result.json'),{'run':str(run),'h5ad':str(final)})
            print('RECOVERY_RUN =',run,flush=True)
            if cfg['run_targeted']:run_followups(run,cfg['threads'],cfg['workers'])
            print('Full new h5ad:',final,flush=True);print('Follow-up reviews:',run/'reviews',flush=True)
        else:ap.print_help()
    except BaseException:
        traceback.print_exc();raise
if __name__=='__main__':cli()
