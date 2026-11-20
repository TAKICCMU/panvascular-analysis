import json,sys
from pathlib import Path
import numpy as np,pandas as pd,h5py
from anndata._io.specs import read_elem
J=Path(__file__).resolve().parents[1];c=json.loads((J/'configs/config.json').read_text())
with h5py.File(c['atlas'],'r') as f:
 obs=read_elem(f['obs']);nnz=np.diff(f['X/indptr'][:])
l1=obs.level1_annotation.astype(str);l2=obs.level2_annotation.astype(str);lab=(l1+' | '+l2).where(~l2.isin(['Unknown','Not_assessed']),'Other | '+l1);upper=0;nc=0;nf=0
for g,ix in lab.groupby(lab,observed=True).indices.items():
 if len(ix)<40:continue
 k=min(1500,len(ix));upper+=int(np.sort(nnz[ix])[-k:].sum());nc+=k;nf+=1
sp=pd.read_csv(J/'evidence/SPATIAL_CONTENT_VALIDATION.tsv',sep='\t');sparse_ref=upper*8+(nc+1)*8
ref_2=2*sparse_ref;sp_2=int(sp.raw_matrix_bytes.sum())*2;summary=int(sp.qc_retained.sum())*nf*4*5
# Additional 2 GiB covers checkpoints, annotation columns, signatures, misc compressed output.
total=ref_2+sp_2+summary+2*1024**3
r={'reference_maximum_selected_cells':nc,'reference_factors_max':nf,'exact_upper_bound_selected_nnz_from_largest_cells_per_factor':upper,'reference_two_uncompressed_sparse_objects_bytes':ref_2,'mapping_two_original_sparse_objects_bytes':sp_2,'posterior_summary_bytes_5_summaries':summary,'checkpoint_metadata_margin_GiB':2,'estimated_conservative_total_GiB':total/1024**3,'reserved_GiB':16,'reservation_sufficient':total<=16*1024**3,'raw_atlas_counts_read':False,'method':'Read X indptr and group maximum 1500 cell nnz upper bounds; no compression benefit assumed; no persisted posterior draws'}
(J/'evidence/DISK_RESERVATION.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2));assert r['reservation_sufficient']
