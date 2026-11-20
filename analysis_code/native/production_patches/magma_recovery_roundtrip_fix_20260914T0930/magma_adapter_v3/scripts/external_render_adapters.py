#!/usr/bin/env python3
"""Generate auditable complete inherited code, without executing any analysis."""
import argparse,json,re,sys,difflib
from pathlib import Path
import yaml
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from external_modules import render_spatial_adapter,sha256
p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);args=p.parse_args()
c=yaml.safe_load(Path(args.config).read_text());o=Path(args.out);o.mkdir(parents=True,exist_ok=False)
root=Path(__file__).resolve().parents[1];ext=c['external'];status=[]
sp=ext.get('spatial',{})
if sp.get('manifest') and sp.get('forced_gene_policy_reviewed'):
    source=root/'external_legacy/spatial/02_spatial_v4_cell2location_GPU.py';before=source.read_text()
    ac=c['atlas'];after=render_spatial_adapter(before,ac['path'],o/'spatial_run',sp['manifest'],
       ac.get('level1_key','level1_annotation'),ac.get('level2_key','level2_annotation'),ac.get('study_key','study'),
       sp.get('forced_genes',[]),ac.get('input_sha256') or 'UNPUBLISHED_BLOCKED',Path(args.config).resolve(),training_config=sp)
    (o/'spatial_CPU_dynamic.py').write_text(after)
    (o/'spatial_source.diff').write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True))))
    status.append({'module':'spatial','status':'CODE_RENDERED_NOT_EXECUTED','source_sha256':sha256(source)})
else:status.append({'module':'spatial','status':'BLOCKED_RESOURCE_OR_FORCED_GENE_POLICY_REVIEW'})
# Exact literal mapping for all inherited genomics R/Python resources; retain every algorithm line.
# Missing old path resolutions block rendering, rather than retaining unusable /path/to strings.
for job in ext.get('inherited_jobs',[]):
    source=root/'external_legacy'/job['source'];before=source.read_text();after=before
    for old,new in job['literal_path_map'].items():
        if old not in after:raise ValueError(f'Unmatched audited path mapping: {old}')
        after=after.replace(old,new)
    unresolved=sorted(set(re.findall(r'/path/to/[^\"\n]+',after)))
    if unresolved:
        status.append({'module':job['name'],'status':'BLOCKED_UNRESOLVED_RESOURCES','unresolved':unresolved});continue
    dest=o/(job['name']+source.suffix)
    # The generated source is provided for review, never called automatically.
    dest.write_text(after)
    (o/(job['name']+'.diff')).write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True))))
    status.append({'module':job['name'],'status':'CODE_RENDERED_NOT_EXECUTED','source_sha256':sha256(source),'rendered_sha256':sha256(dest)})
(o/'adapter_render_status.json').write_text(json.dumps(status,indent=2))
