"""Build new small external sidecars from existing records; never open atlas writable."""
from pathlib import Path
import sys,json,hashlib,shutil,datetime,os
sys.dont_write_bytecode=True
B=Path(__file__).resolve().parents[1];R=B.parent;sys.path.insert(0,str(B))
import pandas as pd,numpy as np,h5py
from anndata.io import read_elem
from src.core import sha,atomic_json
from src.release_compat import validate_metadata_sources,LENGTH_ENCODING
S=R/'atlas_release/L1L2_final_v1_acceptance_r1';OLD=R/'atlas_release/L1L2_final_v1';M=R/'panvascular_CPU_theta2_L1L2_final_v1.acceptance_r1.manifest.json';OM=R/'panvascular_CPU_theta2_L1L2_final_v1.manifest.json'
assert not S.exists() and not M.exists(),'No overwrite: choose separately reviewed suffix'
m=json.loads(OM.read_text());assert m['h5ad_sha256']=='00fc74c1ddaa5e4b32d8664f134e0ec5a0c762f5e9a9fb45b96acbd40804a310'
S.mkdir();D=R/'L2_targeted_followup_metadata_v4_20260913_154850_913989';C=R/'L2_followup_metadata_v4_code';A=R/'metadata_provenance_impact_audit_20260913_013529'
protected={str(p):sha(p) for root in [OLD,R/'downstream_CPU_L1L2_v2_1_revision_20260913'] for p in root.rglob('*') if p.is_file() and '__pycache__' not in str(p)};protected[str(OM)]=sha(OM)
atomic_json(B/'evidence/PROTECTED_SMALL_FILES_BEFORE.json',protected)
for name,h in m['release_sidecar_sha256'].items():
 p=OLD/name;assert sha(p)==h,name;shutil.copy2(p,S/name)
files=[]
def rec(path,role,scope,kind='saved_provenance_record',version='20260913',uri=None):
 p=Path(path);assert p.is_file(),p
 x={'role':role,'path':str(p),'sha256':sha(p),'version':version,'scope':scope,'evidence_type':kind}
 if uri:x['evidence_uri']=uri
 files.append(x);return x
patch=rec(D/'inputs/metadata_h5ad_patch_v4.csv','imported_metadata_csv','all247samples;30 imported source columns',version='v4_evidence_preserving_clean_user_donor_confirmed')
assert patch['sha256']==sha(C/'inputs/metadata_h5ad_patch_v4.csv')
rec(D/'audit/executed_driver.py','import_code','load_metadata/overlay_metadata and saved RUN source-output chain')
run=rec(D/'RUN.json','import_run','actual successful metadata import run; includes donor confirmation summary')
rr=json.loads(Path(run['path']).read_text());assert rr['output_h5ad']==m['source_h5ad_path'] and rr['clinical_missing_counts']['donor_sample_rows_updated_by_user_confirmation']==8
bundle=rec(C/'BUNDLE_CHECKSUMS.json','bundle_checksums','current user-confirmed bundle; exact copied run inputs match')
assert json.loads(Path(bundle['path']).read_text())['inputs/metadata_h5ad_patch_v4.csv']==patch['sha256']
rec(OLD/'sample_metadata.csv','sample_metadata_snapshot','actual release metadata from final obs',version='L1L2_final_v1')
rec(D/'audit/metadata_provenance_v4.csv','metadata_provenance','all247samples;field/sample source classes and limits')
rec(D/'audit/donor_update_user_confirmation_20260913.csv','user_confirmation','GSE166676/GSE184073;8samples;24donor/biosample field updates','user_confirmation_not_independent_public_verification')
rec(D/'audit/protected_source_datasets.json','import_invariance','raw/log1p/coordinates hashes recorded before import')
rec(D/'audit/removed_observation_columns.json','import_column_inventory','removed prior audit/provisional columns archived externally')
rec(D/'audit/obs_before_metadata_update.tsv.gz','before_import_obs','historical observation snapshot; not reconstructed from final')
rec(D/'audit/web_verification_sources.json','historical_web_locators','prior recorded GSE237230 checks; URL locators not new page hashes')
rec(A/'METADATA_AUDIT_REPORT_CN.md','audit_A','metadata provenance/impact audit; proposals not themselves applied')
rec(A/'GSE213128_sample_donor_evidence.csv','audit_A','5aortic valve samples; GSA Individual explicit mappings')
rec(A/'GSE213128_mapping_review.md','audit_A','GEO/GSA name mismatch and limitations retained')
for p in sorted((A/'evidence/public').glob('*')):
 if p.stem.startswith(('HRS43993','HRI27492')) or p.name=='GSE213128_family.soft.txt':
  rec(p,'public_GSA_evidence','GSE213128 saved original sample/Individual/GEO evidence','saved_public_record','20260913')
assert sha(R/'source_B_metadata_curated_v2.csv')==sha(B/'evidence/audit_B/sample_study_clinical_metadata_curated_v2.csv')
rec(R/'source_B_metadata_curated_v2.csv','audit_B','inherited247sample curation; donor/tissue/clinical fields; not independently revalidated now')
rec(B/'evidence/audit_B/SAMPLE_METADATA_COMPLETENESS_AND_DUPLICATION_REPORT.md','audit_B','original auditB report copied byte-identically from local audit deliverable')
rec(B/'evidence/audit_B/SHA256SUMS.txt','audit_B','original auditB checksum list')
for n in ['metadata_changes_v2_to_v3.csv','metadata_import_contract.json','sample_study_clinical_metadata_curated_v3.csv']:
 rec(R/n,'merge_record','historical merged A/B correction and predecessor-v3 policy; NOT claimed final imported CSV')
rec(Path(m['finalization_run'])/'FREEZE_STATUS.json','finalization_record','unchanged metadata and protected datasets through final L2 commit','release_completion','20260914')
rec(OLD/'validation_report.json','final_invariance','final L1 and metadata unchanged from imported working object','release_validation','20260914')
# Bind and verify exact CSV->final sample table relation, then independent native preflight verifies sample->obs.
meta={'schema':'metadata_review_sources_acceptance_r1','files':files,'import_link_validation':{'status':'PASS','bundle_member':'inputs/metadata_h5ad_patch_v4.csv','import_run':str(D),'executed_code_lines':'250-263 load exact bundle CSV, copy inputs, archive old obs; RUN.json records output and user confirmation','csv_sha256':patch['sha256'],'chain':'bundle CSV == copied run CSV -> recorded successful imported working H5AD -> finalization metadata_unchanged -> actual final sample_metadata -> final obs (native preflight)','evidence_limitation':'Historical import RUN records final paths and invariance, not a separate timestamped CSV digest; bundle checksum plus copied run CSV and exact semantic comparison are the retained evidence. No new independent web verification.'},'user_confirmations_are_public_verification':False,'web_policy':'Only existing saved files receive SHA; online URLs are locators only.'}
atomic_json(S/'metadata_review_sources.json',meta);proof=validate_metadata_sources(S/'metadata_review_sources.json');atomic_json(B/'evidence/METADATA_SOURCE_LINK_VALIDATION.json',proof)
# External precise usage records supplement each actual local h5ad registry row.
reg=pd.read_csv(OLD/'source_registry.tsv',sep='\t',dtype=str,keep_default_na=False)
contracts=json.loads((D/'inputs/source_contracts.json').read_text());dec=pd.read_csv(OLD/'approved_level2_mapping.tsv',sep='\t',dtype=str,keep_default_na=False)
rows=[];usage=[]
for _,row in reg.iterrows():
 p=Path(row.path);mp=pd.read_csv(row.mapping_path,sep='\t',dtype=str,keep_default_na=False);ct=contracts[row.branch]
 assert sha(row.mapping_path)==row.mapping_sha256
 with h5py.File(p,'r') as f:
  g=f['obs'];ids=np.asarray(read_elem(g[g.attrs['_index']])).astype(str);vg=f['var'];ng=len(vg[vg.attrs['_index']]);assert np.array_equal(ids,mp.cell_id)
  assert len(ids)==ct['n_cells'] and row.sha256==ct['analysis_sha256']
  # Unique registered neighbor object, not a copied global key.
  gs=[]
  for k in f['uns']:
   if isinstance(f['uns'][k],h5py.Group) and 'connectivities_key' in f['uns'][k] and 'distances_key' in f['uns'][k]:gs.append(k)
  assert len(gs)==1,(row.branch,gs)
  ek='X_umap_mac_L2' if row.branch=='Macrophage' else 'X_umap_lineage_L2';assert ek in f['obsm']
  used={ct['selected_key']};branchdec=dec[dec.branch==row.branch]
  for _,d in branchdec.iterrows():
   keys=[k for k in ct['existing_keys'] if d.original_cluster_uid and d.original_cluster_uid in set(mp[k])]
   if d.cluster_key:keys.append(d.cluster_key)
   assert keys,(row.branch,d.review_unit)
   for k in sorted(set(keys)):
    assert k in g and k in mp
    used.add(k);usage.append({'branch':row.branch,'review_unit':d.review_unit,'original_cluster_uid':d.original_cluster_uid,'label_key':k,'selection_rule':d.selection_rule,'level2_annotation':d.level2_annotation,'n_cells_reviewed':d.n_cells,'decision_basis':d.decision_basis})
  task=p.parent.parent/'task.json';assert task.is_file()
  for k in sorted(used):
   assert np.array_equal(np.asarray(read_elem(g[k])).astype(str),mp[k])
   out=row.to_dict();out.update(version=str(task.parent),version_source=str(task),version_source_sha256=sha(task),n_obs=len(ids),n_vars=ng,label_key=k,graph_key=gs[0],embedding_key=ek,role_schema='local_lineage_with_graph_v1',label_usage_path=str(S/'registry_label_usage.tsv'),label_usage_policy='selected partition plus exact registered final decision keys; residual overrides retain selection rules')
   rows.append(out)
pd.DataFrame(usage).to_csv(S/'registry_label_usage.tsv',sep='\t',index=False)
for row in rows:row['label_usage_sha256']=sha(S/'registry_label_usage.tsv')
pd.DataFrame(rows).to_csv(S/'source_registry.tsv',sep='\t',index=False)
# Declare all publisher annotation/axis hashes, retain every original digest.
m.update(support_revision='acceptance_r1',original_manifest={'path':str(OM),'sha256':sha(OM)},original_manifest_sha256=sha(OM),publisher_script={'path':str(Path(m['finalization_run'])/'audit/finalize_level2.py'),'sha256':sha(Path(m['finalization_run'])/'audit/finalize_level2.py')},hash_encoding={k:LENGTH_ENCODING for k in ['cell_id_order_sha256','gene_id_order_sha256','annotation_sha256','level1_annotation_sha256']},release_support_dir=str(S),validation_report_path=str(S/'validation_report.json'),downstream_execution_authorized=False,compatibility_only=True)
m['release_sidecar_sha256']={p.name:sha(p) for p in S.iterdir() if p.is_file()};atomic_json(M,m)
cfg=json.loads((B/'evidence/phase1_config_PROPOSED.yaml').read_text());cfg['atlas'].update(manifest=str(M),release_compatibility_required=True);cfg['outputs'].update(release_support=str(S),root=str(B/'phase1_outputs_AFTER_APPROVAL'));cfg['execution'].update(run_full_analysis=False,approval_receipt=None)
cfg['statistical_units']['comparison_field']='comparison_group'
# Explicit semantic crosswalk, no change to old design or final metadata.
atomic_json(B/'configs/phase1_config_PROPOSED.yaml',cfg)
atomic_json(B/'evidence/BUILD_COMPAT_SUPPORT.json',{'status':'BUILT_NOT_YET_NATIVE_ACCEPTED','metadata_proof':proof,'registry_rows':len(rows),'registry_usage_rows':len(usage),'manifest':str(M),'manifest_sha256':sha(M),'h5ad_not_copied':True})
print(json.dumps({'metadata_source_files':len(files),'registry_rows':len(rows),'manifest':str(M)}),flush=True)
