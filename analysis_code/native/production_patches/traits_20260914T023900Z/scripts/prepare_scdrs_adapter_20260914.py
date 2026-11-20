from pathlib import Path
import shutil,difflib,json,hashlib,ast
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908');C=R/'downstream_CPU_L1L2_v2_1_release_compat_20260914';W=R/'production_launches/traits_20260914T023900Z';N=W/'scdrs_adapter_v1';N.mkdir(exist_ok=False)
for folder in ['src','scripts','external_legacy']:shutil.copytree(C/folder,N/folder,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
f=N/'scripts/external_run_module.py';old=f.read_text();new=old.replace("hp=out/'scdrs_raw_symbols.h5ad';b.write_h5ad(hp)",'''# Preserve the exact frozen global graph for native scDRS heterogeneity.
            # No reclustering or neighbor recomputation; b retains exact original cell order.
            ng=a.uns['neighbors_cpu_n30']
            ck=ng['connectivities_key']; dk=ng['distances_key']
            if not b.obs_names.equals(a.obs_names): raise ValueError('Frozen graph cell order differs')
            if a.obsp[ck].shape != (b.n_obs,b.n_obs): raise ValueError('Frozen graph shape differs')
            b.obsp['connectivities']=a.obsp[ck].copy()
            b.obsp['distances']=a.obsp[dk].copy()
            b.uns['neighbors']={'connectivities_key':'connectivities','distances_key':'distances','params':dict(ng['params'])}
            hp=out/'scdrs_raw_symbols.h5ad';b.write_h5ad(hp,compression='lzf')''').replace("for key in group_keys:\n                    call([cfg['executable']", "for key in group_keys:\n                    (so/key).mkdir(parents=True,exist_ok=False)\n                    call([cfg['executable']")
assert new!=old;assert "b.write_h5ad(hp)" not in new;ast.parse(new);f.write_text(new)
(W/'evidence/SCDRS_ADAPTER.diff').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile=str(C/'scripts/external_run_module.py'),tofile=str(f))))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
(W/'evidence/SCDRS_ADAPTER_STATUS.json').write_text(json.dumps({'status':'NOT_READY_FOR_PRODUCTION','code_dir':str(N),'syntax_pass':True,'baseline_source_unchanged_sha256':sha(C/'scripts/external_run_module.py'),'new_source_sha256':sha(f),'fixes':['lzf compression only on new raw-symbol object','exact frozen n30 graph copied into expected native keys; no recomputation','native group output directory explicitly created'],'native_tiny_compute':'PASS_nctrl1000','native_tiny_group_with_graph_and_existing_directory':'PASS','full_adapter_tested':False,'blockers':['735 unresolved canonical symbols remain; no fallback guessing','Exact-cell raw-derived covariates not generated yet','Need full derived-object compressed disk estimate incl graph and three full score files','Original old scDRS environment compatibility with full-derived obs schema not yet verified'],'no_production_receipt':True},indent=2))
print(str(N))
