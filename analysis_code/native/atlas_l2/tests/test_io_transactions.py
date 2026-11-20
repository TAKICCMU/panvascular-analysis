"""Small local checks: no AnnData/Scanpy, no server files or science models."""
from pathlib import Path
import importlib.util, tempfile, json, os, shutil
from unittest import mock
import numpy as np
import pandas as pd
import h5py
P=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('freeze_test_module',P/'finalize_level2.py')
M=importlib.util.module_from_spec(s);s.loader.exec_module(M)
results=[]
def passed(name): results.append({'test':name,'passed':True})
def reject(name, f, kinds=(RuntimeError,)):
    try:f()
    except kinds: passed(name);return
    raise AssertionError('Expected rejection: '+name)
def h5(p):
    with h5py.File(p,'w') as f:
        for key in ('X','layers/log1p'):
            g=f.create_group(key);g.attrs['encoding-type']='csr_matrix';g.attrs['shape']=[3,2]
            g.create_dataset('data',data=np.array([1,2,3],dtype=np.int32 if key=='X' else np.float32))
            g.create_dataset('indices',data=[0,1,0]);g.create_dataset('indptr',data=[0,1,2,3])
        f.create_group('var').create_dataset('names',data=np.array(['A','B'],dtype=h5py.string_dtype()))
        f.create_group('obsm').create_dataset('X_umap_cpu',data=np.array([[0.,1.],[2.,3.],[4.,5.]]))
        for slot in ('varm','varp','obsp','obs','uns'): f.create_group(slot)
with tempfile.TemporaryDirectory(prefix='L2_final_io_test_') as td:
    root=Path(td);p=root/'raw.h5ad';h5(p)
    old=M.protected_h5(p)
    M.verify_raw_storage(p,(3,2));passed('HDF5 CSR raw/log1p shape and integer/nonnegative storage check')
    with h5py.File(p,'r+') as f:f['obs'].create_dataset('level2',data=np.array(['A','B','C'],dtype=h5py.string_dtype()))
    assert old==M.protected_h5(p);passed('Metadata-only HDF5 edit leaves all protected logical hashes identical')
    with h5py.File(p,'r+') as f:f['X/data'][0]=4
    assert old['X']!=M.protected_h5(p)['X'];passed('Changed raw value is detected by protected hash')
    with h5py.File(p,'r+') as f:f['layers/log1p/data'][0]=-2.
    reject('Invalid log1p negative value rejected',lambda:M.verify_raw_storage(p,(3,2)))
    with h5py.File(p,'r+') as f:f['obsm']['foreign']=h5py.SoftLink('/X')
    reject('HDF5 soft/external links rejected in protected slots',lambda:M.protected_h5(p))
    # Ensure helper convention agrees with original already-run program.
    h5(root/'same.h5ad')
    original=Path('/mnt/data/_final_L2_work/L2_annotation_update_next_review/review_utils.py')
    if original.exists():
        spec=importlib.util.spec_from_file_location('old_hash_helpers',original);v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
        assert M.protected_h5(root/'same.h5ad')==v.protected_h5(root/'same.h5ad')
        passed('All logical dataset hashes match the preceding update program convention')
    obs=pd.DataFrame({'sample_uid':['s1','s1','s2'],'donor_uid':['d1','d1','d2'],'age':[None,None,60.],
                      'sex':[None,None,'male'],'library_size_raw_cpu':[10,100,200], 'level1_annotation':['T cell']*3})
    meta=M.write_sample_metadata(obs,root/'samples.csv')
    assert len(meta)==2 and meta.age.isna().sum()==1 and 'library_size_raw_cpu' not in meta
    passed('Clinical missing values preserved; per-cell QC excluded from sample conflict checks')
    bad=obs.copy();bad['age']=[40.,50.,60.]
    reject('True within-sample clinical conflict rejected without first-row imputation',lambda:M.write_sample_metadata(bad,root/'bad.csv'))
    # Normal full-file replacement with a shared-inode canonical alias.
    stage=root/'candidate.h5ad';stage.write_bytes(b'new-file');stage.chmod(0o444)
    target=root/'work.h5ad';target.write_bytes(b'original-file');archive=root/'backup';archive.mkdir()
    can=root/'canonical.h5ad';sig=M.file_stat(target)
    back=M.transactional_publish(stage,target,can,archive,sig)
    assert target.read_bytes()==b'new-file' and back.read_bytes()==b'original-file'
    passed('Original full file retained by backup before explicit working-path replacement')
    assert os.path.samefile(target,can) and not can.is_symlink()
    passed('Canonical and working names share one fully embedded hard-linked file, not two copies')
    assert not(target.stat().st_mode&0o222)
    passed('Published file retains read-only permissions')
    stage2=root/'s2';stage2.write_bytes(b'next');sig=M.file_stat(target)
    target.chmod(0o644);target.write_bytes(b'changed')
    reject('Concurrent target change rejected before overwrite',lambda:M.transactional_publish(stage2,target,root/'c2',archive,sig))
    reject('Existing canonical release cannot be silently overwritten',lambda:M.transactional_publish(stage2,target,can,archive,M.file_stat(target)))
    # Rollback after the target replacement but before canonical link install.
    r=root/'rollback';r.mkdir();a=r/'backup';a.mkdir();t=r/'work';t.write_bytes(b'old');s=r/'stage';s.write_bytes(b'new');c=r/'canon'
    original_replace=os.replace
    def fail_second(src,dst):
        if Path(dst)==c: raise OSError('injected canonical installation failure')
        return original_replace(src,dst)
    with mock.patch.object(M.os,'replace',side_effect=fail_second):
        reject('Canonical install failure rolls back working file',lambda:M.transactional_publish(s,t,c,a,M.file_stat(t)),(OSError,))
    assert t.read_bytes()==b'old' and (a/'work').read_bytes()==b'old' and not c.exists()
    passed('Rollback preserves readable original and backup')
    # Post-h5ad publication recovery, SHA checked. No AnnData required here.
    r=root/'recovery';r.mkdir();run=r/'run';run.mkdir();(run/'backup').mkdir();(run/'audit').mkdir();(run/'release_staging').mkdir()
    t=r/'work.h5ad';t.write_bytes(b'new-data');c=r/'final.h5ad';os.link(t,c)
    backup=run/'backup/work.h5ad';backup.write_bytes(b'old-data')
    release=r/'atlas_release/L1L2_final_v1';support=run/'release_staging'
    M.dump_json(support/'validation_report.json',{'status':'PASS'})
    (support/'level2_ontology.tsv').write_text('level2_state_id\tlevel2_annotation\nL2_01\tCD68 macrophage\n')
    manifest={'working_path_alias':str(t),'h5ad_path':str(c),'release_support_dir':str(release),'h5ad_sha256':M.sha_file(t),
              'backup_h5ad':str(backup),'source_h5ad_sha256':M.sha_file(backup),
              'release_sidecar_sha256':{p.name:M.sha_file(p) for p in support.iterdir()},
              'shape_verified':[100,4],'level2_special_counts_verified':{'Unknown':10,'Not_assessed':2},'named_level2_states_verified':2}
    M.dump_json(run/'FREEZE_STATUS.json',{'status':'H5AD_COMMITTED_SIDECARS_PENDING'})
    M.finish_publication(manifest,run)
    assert release.exists() and t.with_suffix('.manifest.json').exists() and c.with_suffix('.h5ad.sha256').exists()
    assert M.load_json(run/'FREEZE_STATUS.json')['status']=='COMMITTED_FINAL_APPROVED'
    passed('Interrupted sidecar publication can be completed from the verified staged manifest')
    M.finish_publication(manifest,run)
    passed('Finishing a published release again is idempotent and does not rewrite expression')
    (release/'level2_ontology.tsv').write_text('changed')
    reject('Changed final sidecar rejected on release recheck',lambda:M.finish_publication(manifest,run))
    # Check flock conventions used by existing focus scheduler.
    import fcntl
    rev=root/'active_review';rev.mkdir();lock=rev/'focus_scheduler.lock'
    with lock.open('w') as handle:
        fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        reject('Running previous focus scheduler rejected via its actual flock',lambda:M.require_idle_review(rev))
        fcntl.flock(handle,fcntl.LOCK_UN)
    M.dump_json(rev/'FOCUS_STATUS.json',{'status':'PARTIAL_FAILURE_LABELS_SAVED'})
    reject('Incomplete prior focused analysis rejected',lambda:M.require_idle_review(rev))
    assert M.marker_anchor('PDE4C-high macrophage')=='PDE4C'
    assert M.marker_anchor('FCGR3A-expressing monocyte')=='FCGR3A'
    assert M.marker_anchor('HLA-DRA antigen-presenting macrophage')=='HLA-DRA'
    passed('Ontology marker parsing preserves real hyphenated genes and strips only defined display suffixes')
    reject('Conflicting formal L2 alias columns rejected',lambda:M.labels_from_obs(pd.DataFrame({'level2_annotation':['A'],'L2_annotation':['B']})))

out={'n_tests':len(results),'tests':results,'native_anndata_io_tested_locally':False,'server_h5ad_accessed':False}
(P/'tests/io_transaction_checks.json').write_text(json.dumps(out,indent=2,ensure_ascii=False))
print('PASS:',len(results),'I/O / transaction checks')
