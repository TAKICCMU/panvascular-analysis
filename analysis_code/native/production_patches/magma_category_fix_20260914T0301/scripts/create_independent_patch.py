from pathlib import Path
import shutil,difflib,ast
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908');B=R/'downstream_CPU_L1L2_v2_1_release_compat_20260914';W=R/'production_launches/magma_category_fix_20260914T0301';C=W/'magma_adapter_v1'
for d in ['scripts','logs','evidence','configs','approvals','outputs']:(W/d).mkdir(parents=True,exist_ok=True)
for d in ['src','scripts','external_legacy']:shutil.copytree(B/d,C/d,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
f=C/'src/external_modules.py';old=f.read_text();new=old.replace('def specificity_from_counts(','''def identity_codes(values, mapping):
    """Map exact identity values to numeric codes on a disposable object Series.
    Missing/unmapped identities become -1; no category mutation or string coercion.
    Leading zeros and distinct string identities remain distinct.
    """
    return pd.Series(values, copy=False).astype(object).map(mapping).fillna(-1).to_numpy(int)

def specificity_from_counts(''').replace("gi=labels.map(dict(zip(groups,range(len(groups))))).fillna(-1).to_numpy(int)","gi=identity_codes(labels,dict(zip(groups,range(len(groups)))))").replace("si=obs[study_key].map(dict(zip(studies,range(len(studies))))).fillna(-1).to_numpy(int)","si=identity_codes(obs[study_key],dict(zip(studies,range(len(studies)))))")
assert new!=old;assert new.count('identity_codes(')==3;ast.parse(new);f.write_text(new);(W/'evidence/EXACT_SOURCE_CHANGE.diff').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile=str(B/'src/external_modules.py'),tofile=str(f))))
print(str(W))
