from pathlib import Path
import csv,gzip,json,hashlib,shutil,sys
from exclusion import token,row_tokens
D=Path('/data/anzhen_vhd/scRNAseq/codex/panvascular_downstream_storage_v1');W=D/'expanded_parallel_20260916_v1';O=W/'overnight_manuscript_integration_20260917';J=O/'spatial_target_excluded_v1';old=D/'spatial_gpu_data_v1/shared_gpu_r3';resume=W/'spatial_gpu_resume_v2'
for x in ['code','configs','evidence','outputs','smoke']:(J/x).mkdir(parents=True,exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def rows(p):
 with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:return list(csv.DictReader(f,delimiter='\t'))
target=O/'candidate_prioritization/Tier1_spatial_interface.tsv';tr=rows(target);ids=sorted({token(r['gene_id']) for r in tr});assert len(ids)==29
excluded={token(r[k]) for r in tr for k in ['gene_id','gene_symbol','SYMBOL'] if r[k]}
universe=old/'outputs/reference_gene_universe.tsv';ur=rows(universe)
maps=sorted((resume/'evidence').glob('*__gene_map.tsv.gz'));assert maps
records=[(str(universe),r.get('',next(iter(r.values()))),r) for r in ur]
for p in maps:records.extend((str(p),str(r['source_row']),r) for r in rows(p))
# Only curated identity columns, never fuzzy symbols or numerical metadata.
changed=True
while changed:
 before=len(excluded)
 for p,i,r in records:
  ts=row_tokens(i,r)
  # source_row numeric index is not a gene identity
  ts={t for t in ts if not t.isdigit()}
  if ts&excluded:excluded.update(ts)
 changed=len(excluded)>before
with open(J/'evidence/target_alias_resolution.tsv','w') as f:
 f.write('source\tfeature\tidentity_tokens\n')
 for p,i,r in records:
  ts={t for t in row_tokens(i,r) if not t.isdigit()}
  if ts&excluded:f.write(p+'\t'+i+'\t'+';'.join(sorted(ts))+'\n')
with open(J/'evidence/target_29_gene_audit.tsv','w') as f:
 f.write('gene_id\tgene_symbol\treference_matching_features\n')
 for ident in ids:
  rs=[r for r in tr if token(r['gene_id'])==ident];sy=sorted({token(r['gene_symbol']) for r in rs});matches=[r.get('','') for r in ur if ident in row_tokens(r.get('',''),r)]
  f.write(ident+'\t'+';'.join(sy)+'\t'+';'.join(matches)+'\n')
c=json.loads((resume/'configs/config.json').read_text());c.update(reference_mode='target_excluded_both_stages',excluded_targets=sorted({token(r['gene_symbol']) for r in tr}),target_ids=ids,excluded_feature_tokens=sorted(excluded),reference_counts=str(old/'outputs/reference_counts.h5ad'),storage_output_root=str(J),target_list=str(target),forced_genes=[])
(J/'configs/config.json').write_text(json.dumps(c,indent=2))
shutil.copy2(resume/'evidence/SPATIAL_CONTENT_VALIDATION.tsv',J/'evidence/SPATIAL_CONTENT_VALIDATION.tsv')
for p in maps:shutil.copy2(p,J/'evidence'/p.name)
manifest=Path(c['spatial_manifest']);mr=[r for r in rows(manifest) if r['include_mapping'].lower() in ('true','1','yes')]
files=[target,universe,old/'outputs/reference_counts.h5ad',old/'outputs/reference_factors.tsv',manifest,J/'configs/config.json',J/'evidence/SPATIAL_CONTENT_VALIDATION.tsv']+list((J/'code').glob('*.py'))+list((J/'evidence').glob('*gene_map.tsv.gz'))
files+=sorted({Path(r[k]) for r in mr for k in ['matrix','barcodes','features','positions']})
receipt={'scope':['target_excluded_reference','target_excluded_mapping'],'cpu_threads':4,'bound_files':{str(p):sha(p) for p in files},'target_id_count':29,'alias_scope':'exact version-normalized Ensembl and symbol closure over reference var metadata and all verified spatial gene maps; no external alias dictionary claimed','reference_shared_parameters_require_retrain':True,'status':'PREPARED_NOT_LAUNCHED'}
(J/'configs/PRODUCTION_RECEIPT.json').write_text(json.dumps(receipt,indent=2))
print(json.dumps({'job':str(J),'targets':len(ids),'excluded_tokens':len(excluded),'bound_files':len(files),'status':'PREPARED_NOT_LAUNCHED'}))
