from pathlib import Path
import subprocess,json,hashlib,pandas as pd
O=Path('/Users/taki/Desktop/Panvascular_Spatial_Evidence_Assessment_20260920');S=O/'code/native_sources';S.mkdir(exist_ok=True)
remote={'sample_manifest.tsv':'/home/anzhen_vhd/scRNAseq/codex/spatial_validation_v3/00_inputs/spatial_sample_manifest_curated.tsv','pathology_availability.tsv':'/home/anzhen_vhd/scRNAseq/codex/spatial_validation_v3/00_inputs/public_pathology_annotations/public_pathology_label_availability.tsv','pathology_proxy.tsv':'/home/anzhen_vhd/scRNAseq/codex/spatial_validation_v3/00_inputs/public_pathology_annotations/GSE274572_HPcluster_pathology_proxy_lookup.tsv','historical_residual_moran.tsv.gz':'/home/anzhen_vhd/scRNAseq/codex/spatial_validation_v5_locked71_spatial_statistics_audit_v1/01_target_NB_Poisson_residual_Moran_by_section.tsv.gz','historical_manifest.json':'/home/anzhen_vhd/scRNAseq/codex/spatial_validation_v5_locked71_spatial_statistics_audit_v1/run_manifest.json'}
for name,path in remote.items():
 data=subprocess.check_output(['ssh','-o','BatchMode=yes','codex-shpc-gpu','cat',path]);(S/name).write_bytes(data)
(S/'source_paths.json').write_text(json.dumps(remote,indent=2))
