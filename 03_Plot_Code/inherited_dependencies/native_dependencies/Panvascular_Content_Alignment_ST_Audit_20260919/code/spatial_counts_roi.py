from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from matplotlib.path import Path as PolygonPath
O=Path('/Users/taki/Desktop/Panvascular_Content_Alignment_ST_Audit_20260919');R=Path('/Users/taki/Desktop/Codex任务1/spatial_visium_geo_local/extracted')
roi_path=O/'audit/spatial_roi_definitions.json';config=json.loads(roi_path.read_text());summaries=[];allrows=[]
def describe(a,prefix):
 a=pd.to_numeric(a,errors='coerce').dropna();return {prefix+'_n_available':len(a),prefix+'_sum':float(a.sum()) if len(a) else None,prefix+'_min':float(a.min()) if len(a) else None,prefix+'_median':float(a.median()) if len(a) else None,prefix+'_max':float(a.max()) if len(a) else None,prefix+'_mean':float(a.mean()) if len(a) else None}
for roi in config['rois']:
 gene=roi['gene'];gse='GSE248608' if gene=='LIPA' else 'GSE314851';z=pd.read_csv(O/f'audit/spatial_counts_{gene}_all_barcode_universe.tsv',sep='\t');z=z[z.in_positions].copy();vertices=np.asarray(roi['polygon'],dtype=float);points=z[['source_hires_pixel_x','source_hires_pixel_y']].to_numpy();a=vertices;b=np.roll(vertices,-1,axis=0);ab=b-a
 t=np.sum((points[:,None,:]-a[None,:,:])*ab[None,:,:],axis=2)/np.sum(ab*ab,axis=1)[None,:];closest=a[None,:,:]+np.clip(t,0,1)[:,:,None]*ab[None,:,:];dist=np.linalg.norm(points[:,None,:]-closest,axis=2).min(axis=1)
 inside=PolygonPath(np.vstack([vertices,vertices[0]])).contains_points(points)|(dist<=1e-9)
 sf=json.loads((R/gse/(roi['section']+'_scalefactors_json.json')).read_text());radius=sf['spot_diameter_fullres']*sf['tissue_hires_scalef']/2
 z.insert(0,'roi_id',roi['id']);z.insert(1,'roi_role',roi['role']);z['center_in_roi']=inside;z['distance_to_polygon_boundary_hires_pixels']=dist;z['signed_distance_hires_pixels_positive_inside']=np.where(inside,dist,-dist);z['physical_spot_radius_hires_pixels']=radius;z['entire_physical_spot_disk_inside_roi']=inside&(dist>=radius-1e-9);z['physical_disk_intersects_roi']=inside|(dist<=radius+1e-9);z['raw_gene_CP10k']=z.raw_gene_count/z.matrix_total_UMI_all_features*1e4;z['raw_gene_log1p_CP10k']=np.log1p(z.raw_gene_CP10k)
 z.to_csv(O/f'audit/spatial_counts_roi_{roi["id"]}_all_positions.tsv',sep='\t',index=False)
 q=z[z.center_in_roi].copy();q.to_csv(O/f'audit/spatial_counts_roi_{roi["id"]}_center_included.tsv',sep='\t',index=False);allrows.append(q)
 for subset,mask in [('center_included',q.center_in_roi),('whole_disk_included',q.entire_physical_spot_disk_inside_roi)]:
  v=q[mask];r={'roi_id':roi['id'],'gene':gene,'section':roi['section'],'role':roi['role'],'subset':subset,'physical_spot_radius_hires_pixels':radius,'n_positions':len(v),'n_in_matrix':int(v.in_matrix.sum()),'n_missing_matrix':int((~v.in_matrix).sum()),'n_plotted':int(v.in_plot.sum()),'n_not_plotted':int((~v.in_plot).sum()),'n_in_tissue_1':int(v.in_tissue.eq(1).sum()),'n_in_tissue_0':int(v.in_tissue.eq(0).sum()),'n_raw_gene_positive_available':int(v.raw_gene_count.gt(0).sum()),'n_raw_gene_zero_available':int(v.raw_gene_count.eq(0).sum()),'n_whole_disk_inside':int(v.entire_physical_spot_disk_inside_roi.sum())}
  for col,pref in [('raw_gene_count','raw_gene_UMI'),('matrix_total_UMI_all_features','total_UMI'),('raw_gene_CP10k','CP10k'),('raw_gene_log1p_CP10k','log1p_CP10k'),('distance_to_polygon_boundary_hires_pixels','boundary_distance_pixels')]:r.update(describe(v[col],pref))
  summaries.append(r)
pd.concat(allrows,ignore_index=True).to_csv(O/'audit/spatial_counts_roi_all_center_included.tsv',sep='\t',index=False);pd.DataFrame(summaries).to_csv(O/'audit/spatial_counts_roi_summary.tsv',sep='\t',index=False)
(O/'audit/spatial_counts_roi_summary.json').write_text(json.dumps({'roi_file_sha256':hashlib.sha256(roi_path.read_bytes()).hexdigest(),'selection_basis':config['selection_basis'],'inclusion_rule':'Boundary included; center tested against polygon. Entire disk iff center inside and minimum Euclidean distance to polygon edges >= physical spot radius. No marker display radius used.','missing_values':'No matrix means unavailable counts/library/CP10k; never imputed as zero.','summaries':summaries},indent=2,ensure_ascii=False))
print(pd.DataFrame(summaries)[['roi_id','subset','n_positions','n_plotted','n_missing_matrix','n_in_tissue_1','n_raw_gene_positive_available','raw_gene_UMI_sum','raw_gene_UMI_max','total_UMI_median','log1p_CP10k_median','physical_spot_radius_hires_pixels']].to_string(index=False))
