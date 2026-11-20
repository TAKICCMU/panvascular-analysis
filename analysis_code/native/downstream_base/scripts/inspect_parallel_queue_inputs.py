"""Read matrix headers/features and executable versions only. No spatial model."""
from pathlib import Path
import sys,json,os,gzip,csv,subprocess,datetime
sys.dont_write_bytecode=True
B=Path(__file__).resolve().parents[1];R=B.parent;sys.path.insert(0,str(B))
from src.core import read_frame,sha,atomic_json
import pandas as pd
cfg=json.loads((R/'downstream_CPU_L1L2_v2_1_revision_20260913/configs/analysis_config.yaml').read_text());p=Path(cfg['external']['spatial']['manifest']);atlas=Path(cfg['atlas']['path']);var=read_frame(atlas,'var');ids=set(var.index.astype(str));symbols=set(var['gene_symbol'].dropna().astype(str)) if 'gene_symbol' in var else ids
rows=[]
def op(path):return gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)
for r in csv.DictReader(p.open(),delimiter='\t'):
 row={k:r[k] for k in ['GSE','sample_id','section_id','biological_unit_id','include_mapping','include_group_contrast']};missing=[k for k in ['matrix','barcodes','features','positions'] if not Path(r[k]).is_file()];row['missing_files']=';'.join(missing)
 if not missing:
  with op(r['matrix']) as f:
   header=next(line.strip() for line in f if not line.startswith('%') and line.strip());dims=list(map(int,header.split()))
  with op(r['features']) as f:features=[line.rstrip('\n').split('\t') for line in f]
  with op(r['barcodes']) as f:nbar=sum(1 for line in f if line.strip())
  gid=[x[0] for x in features];sym=[x[1] if len(x)>1 else x[0] for x in features]
  row.update(n_features=len(features),n_barcodes=nbar,matrix_shape_nnz=json.dumps(dims),dimensions_match=dims[:2]==[len(features),nbar],exact_ID_overlap=len(ids&set(gid)),exact_symbol_overlap=len(symbols&set(sym)),duplicate_feature_symbols=len(sym)-len(set(sym)),features_sha256=sha(r['features']),barcodes_sha256=sha(r['barcodes']),matrix_bytes=Path(r['matrix']).stat().st_size,matrix_whole_sha='NOT_RECOMPUTED_HEADER_ONLY')
 rows.append(row)
pd.DataFrame(rows).to_csv(B/'evidence/SPATIAL_INPUT_HEADER_AND_GENE_CHECK.tsv',sep='\t',index=False)
resources=[]
for module,key in [('magma','gene_results'),('scdrs','geneset')]:
 for x in cfg['external'][module]['scenarios']:
  q=Path(x[key]);h=sha(q) if q.is_file() else None;resources.append({'module':module,'trait':x['name'],'path':str(q),'exists':q.is_file(),'sha256':h,'matches_declared':h==x.get('sha256')})
for module in ['magma','scdrs']:
 q=Path(cfg['external'][module]['executable']);resources.append({'module':module,'executable':str(q),'exists':q.is_file(),'executable_bit':os.access(q,os.X_OK)})
atomic_json(B/'evidence/CPU_TRAIT_INPUT_CHECK.json',resources)
q={'timestamp':datetime.datetime.now(datetime.timezone.utc).isoformat(),'spatial_manifest':str(p),'spatial_manifest_sha256':sha(p),'sections':len(rows),'headers_valid':sum(x.get('dimensions_match',False) for x in rows),'missing_rows':sum(bool(x['missing_files']) for x in rows),'exact_symbol_overlap_range':[min(x.get('exact_symbol_overlap',0) for x in rows),max(x.get('exact_symbol_overlap',0) for x in rows)],'target_exclusion_rule_approved':False,'target_exclusion_list':cfg['external']['spatial'].get('excluded_targets',[]),'forced_genes_reviewed':cfg['external']['spatial'].get('forced_gene_policy_reviewed',False),'whole_spatial_matrices_authenticated':False,'training_started':False,'GPU_visibility_instruction':'Separate GPU process environment, do not inherit CPU CUDA_VISIBLE_DEVICES empty mask','prior_spatial_parameters_PROPOSED_NOT_APPROVED':{k:cfg['external']['spatial'].get(k) for k in ['reference_max_epochs','mapping_max_epochs','N_cells_per_location','detection_alpha','max_cells_per_state','min_cells_per_state','posterior_samples','reference_mode']}}
atomic_json(B/'evidence/PARALLEL_QUEUE_INPUT_CHECK.json',q);print(json.dumps(q),flush=True)
