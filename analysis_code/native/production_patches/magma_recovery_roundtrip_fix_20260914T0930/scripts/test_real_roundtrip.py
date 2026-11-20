import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='2'
import sys,json,hashlib,shutil,ast,time
from pathlib import Path
from types import SimpleNamespace
import h5py,numpy as np,pandas as pd
from anndata.io import read_elem,sparse_dataset
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908');W=R/'production_launches/magma_recovery_roundtrip_fix_20260914T0930';C=W/'magma_adapter_v3';S=W/'smoke';S.mkdir();inputs=S/'input';inputs.mkdir();out=S/'actual_recovery_entry';out.mkdir();sys.path.insert(0,str(C/'src'))
from magma_recovery import run_recovery,verify_serialized_frame,_read_and_validate
from external_modules import specificity_from_counts,sha256
old=R/'production_launches/magma_category_fix_20260914T0301/outputs/runs/20260913T190814_808850_7ccf71427e40/magma';mapping=pd.read_csv(old/'L1_variables.tsv',sep='\t');gmfile=R/'production_launches/traits_20260914T023900Z/inputs/gene_map_MAGMA_READY.tsv';gm=pd.read_csv(gmfile,sep='\t')
t=time.time()
with h5py.File(R/'panvascular_CPU_theta2_L1L2_final_v1.h5ad','r') as f:
 obs=pd.DataFrame({'study':read_elem(f['obs/study']),'level1_annotation':read_elem(f['obs/level1_annotation'])},index=read_elem(f['obs/_index']))
 chosen=[]
 for group in mapping.group:
  ii=np.flatnonzero(obs.level1_annotation.eq(group).to_numpy());counts=obs.iloc[ii].study.value_counts();studies=counts[counts>=20].index[:3];assert len(studies)==3
  selected=[]
  for study in studies:selected.extend(ii[obs.iloc[ii].study.eq(study).to_numpy()][:20].tolist())
  rest=ii[~np.isin(ii,selected)];selected.extend(rest[:500-len(selected)].tolist());assert len(selected)==500;chosen.extend(selected)
 chosen=np.array(sorted(chosen));assert len(np.unique(chosen))==7500;X=sparse_dataset(f['X'])[chosen,:].tocsr();smallobs=obs.iloc[chosen].copy()
smallobs['_state']=smallobs.level1_annotation.astype(object) # Plumbing-only fixture: L2 deliberately mirrors L1; never used in scientific output.
pd.DataFrame({'cell_id':smallobs.index,'L1':smallobs.level1_annotation,'study':smallobs.study}).to_csv(S/'selected_real_cells.tsv',sep='\t',index=False)
print('Real count subset',X.shape,X.nnz,'seconds',time.time()-t,flush=True)
cell,study,mp,audit=specificity_from_counts(X,smallobs,gm.ensembl_id,smallobs.level1_annotation)
assert len(mp)==15;assert cell.columns.name=='variable';cell.attrs['gene_mapping_audit'].to_csv(inputs/'L1_gene_mapping_audit.tsv',sep='\t',index=False);cell.to_csv(inputs/'L1_cell_weighted_specificity.tsv',sep='\t',index=False);mp.to_csv(inputs/'L1_variables.tsv',sep='\t',index=False);audit.to_csv(inputs/'L1_eligibility.tsv',sep='\t',index=False)
serialized=pd.read_csv(inputs/'L1_cell_weighted_specificity.tsv',sep='\t');assert serialized.columns.name is None
verify_serialized_frame(serialized,cell);passed=['actual_7500_real_cells_37336_genes_specificity_TSV_readback_same_reuse_validator_pass']
renamed=cell.copy().rename_axis(index='rows_not_stored_by_CSV',columns='variable_not_stored_by_CSV');verify_serialized_frame(serialized,renamed);passed.append('axis_names_not_persisted_by_CSV_only_normalized')
def reject(name,bad):
 try:verify_serialized_frame(bad,cell)
 except AssertionError:passed.append(name);return
 raise AssertionError(name+' unexpectedly accepted')
bad=serialized.copy();bad.iloc[0,1]+=1e-3;reject('numerical_content_change_rejected',bad)
bad=serialized.copy();bad=bad.rename(columns={bad.columns[2]:'WRONG_STATE'});reject('column_identity_change_rejected',bad)
reject('column_order_change_rejected',serialized.iloc[:,[1,0]+list(range(2,len(serialized.columns)))])
reject('row_order_change_rejected',serialized.iloc[::-1].reset_index(drop=True))
bad=serialized.copy();bad['Average']=bad.Average.astype(object);reject('unexpected_dtype_change_rejected',bad)
for name in ['primary_r0.1_L1_cell_weighted.gsa.out','primary_r0.1_L1_cell_weighted.log']:shutil.copyfile(old/name,inputs/name)
target,conditions,models=_read_and_validate(inputs/'primary_r0.1_L1_cell_weighted.gsa.out',mp,'test_actual_gsa_parser');assert len(target)==len(conditions)==len(models)==15;passed.append('real_native_gsa_full_parser_with_matching_real_L1_mapping_pass')
dummy=inputs/'SMOKE_ONLY.raw';dummy.write_text('No GWAS analysis occurs in this smoke fixture.\n')
cert={'atlas_sha256':'SMOKE_ONLY_REAL_COUNT_SUBSET','annotation_and_study_keys':{'level1_key':'level1_annotation','level2_key':'level2_annotation','study_key':'study'},'specificity_source_sha256':sha256(C/'src/external_modules.py'),'files':{p.name:{'path':str(p),'sha256':sha256(p)} for p in inputs.iterdir() if p.is_file()}};manifest=inputs/'reuse_manifest.json';manifest.write_text(json.dumps(cert,indent=2));a=SimpleNamespace(X=X,obs=smallobs,var_names=pd.Index(gm.feature_id));cfg={'recovery_manifest':str(manifest),'recovery_manifest_sha256':sha256(manifest),'executable':'SMOKE_NATIVE_FIXTURE_NOT_EXECUTED','scenarios':[{'name':name,'gene_results':str(dummy),'sha256':sha256(dummy)} for name in ['primary_r0.1','primary_r0.6','leave_one_CAD_r0.1']]};ac={'input_sha256':'SMOKE_ONLY_REAL_COUNT_SUBSET','level1_key':'level1_annotation','level2_key':'level2_annotation','study_key':'study'};calls=[]
def native_fixture(cmd):
 assert cmd[0]=='SMOKE_NATIVE_FIXTURE_NOT_EXECUTED';prefix=Path(cmd[cmd.index('--out')+1]);assert prefix.parent==out;shutil.copyfile(old/'primary_r0.1_L1_cell_weighted.gsa.out',str(prefix)+'.gsa.out');Path(str(prefix)+'.log').write_text('SMOKE ONLY: stored native output fixture, executable was not run.\n');calls.append(cmd)
run_recovery(a,ac,cfg,{'gene_map':str(gmfile)},out,native_fixture)
assert len(calls)==11;assert len(list(out.glob('*_with_BH.tsv')))==12;assert len(list(out.glob('*_condition_rows.tsv')))==12;assert json.loads((out/'RECOVERY_PROGRESS.json').read_text())['stage']=='ALL_12_ANALYSES_COMPLETE';passed.append('actual_run_recovery_entry_complete_all_12_output_paths_with_real_raw_subset_and_stored_native_fixtures')
# Resource bootstrap branch: only parent magma worker is reduced; scientific external loader remains configured 12.
text=(C/'scripts/run_downstream.py').read_text();assign=next(n for n in ast.parse(text).body if isinstance(n,ast.Assign) and any(isinstance(x,ast.Name) and x.id=='threads' for x in n.targets));expr=compile(ast.Expression(assign.value),'<bootstrap-test>','eval')
for worker,expected in [('magma','1'),('cohort','12'),(None,'12')]:assert eval(expr,{'args':SimpleNamespace(worker=worker),'c':{'execution':{'per_worker_threads':12}}})==expected
passed.append('parent_magma_worker_one_thread_other_worker_budget_unchanged')
sys.path.insert(0,str(C));from src.dispatcher import code_digest
report={'status':'PASS','passed':len(passed),'tests':passed,'real_subset_cells':7500,'real_subset_genes':37336,'real_subset_nnz':int(X.nnz),'original_minimums_retained':{'min_cells':500,'min_studies':3,'min_study_cells':20},'actual_native_gsa_parser':True,'full_recovery_entry_exercised':True,'new_native_BINARY_calls_in_this_test':0,'stored_native_output_fixture_calls':11,'test_only_L2_mirrors_L1':True,'elapsed_seconds':time.time()-t,'code_digest':code_digest(C),'test_scientific_results_not_valid_for_inference':True};(W/'evidence/REAL_ROUNDTRIP_TEST_RESULTS.json').write_text(json.dumps(report,indent=2));(S/'SMOKE_ONLY.json').write_text(json.dumps({'not_scientific_results':True,'purpose':'Real raw-count roundtrip and recovery/parser plumbing test'}));print(json.dumps(report),flush=True)
