from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy.io import mmread
from PIL import Image
O=Path('/Users/taki/Desktop/Panvascular_Content_Alignment_ST_Audit_20260919')
R=Path('/Users/taki/Desktop/Codex任务1/spatial_visium_geo_local/extracted')
D=Path('/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Nature_Reconstruction_20260918/SourceData/spatial')
params=json.loads((D/'Figure5a_display_parameters.json').read_text());results=[];hashes={}
for ind,(gse,stem,gene,expected) in enumerate([('GSE248608','GSM7917895_CA2','LIPA',1254),('GSE314851','GSM9414344_LCA','STAT3',2072)]):
 base=R/gse;paths={k:base/(stem+'_'+k) for k in ['matrix.mtx','features.tsv','barcodes.tsv','tissue_positions_list.csv','scalefactors_json.json','tissue_hires_image.png']};paths['plot_spots']=D/f'Figure5a_{gene}_spots.tsv'
 for key,path in paths.items():hashes[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
 features=pd.read_csv(paths['features.tsv'],sep='\t',header=None,names=['feature_id','symbol','feature_type']);features['matrix_row_1based']=np.arange(1,len(features)+1)
 features[features.duplicated('symbol',keep=False)|features.duplicated('feature_id',keep=False)].to_csv(O/f'audit/spatial_counts_{gene}_duplicate_features.tsv',sep='\t',index=False)
 target=features[features.symbol.eq(gene)];target.to_csv(O/f'audit/spatial_counts_{gene}_target_features.tsv',sep='\t',index=False)
 assert len(target)==1,('Ambiguous symbol',gene,target.to_dict('records'))
 assert features.feature_id.eq(target.iloc[0].feature_id).sum()==1
 bc=pd.read_csv(paths['barcodes.tsv'],sep='\t',header=None,names=['barcode']);assert not bc.barcode.duplicated().any()
 mat=mmread(paths['matrix.mtx']).tocsr();assert mat.shape==(len(features),len(bc));assert np.all(mat.data>=0) and np.all(mat.data==np.floor(mat.data))
 raw=np.asarray(mat[target.index[0],:].todense()).ravel().astype('int64');library=np.asarray(mat.sum(axis=0)).ravel().astype('int64')
 bc['matrix_column_1based']=np.arange(1,len(bc)+1);bc['in_matrix']=True;bc['raw_gene_count']=raw;bc['matrix_total_UMI_all_features']=library;bc['matrix_detected_features']=np.asarray((mat>0).sum(axis=0)).ravel();bc['audit_log1p_1e4_raw_over_library']=np.log1p(np.divide(raw*1e4,library,out=np.zeros(len(raw),dtype=float),where=library>0))
 pos=pd.read_csv(paths['tissue_positions_list.csv'],header=None,names=['barcode','in_tissue','array_row','array_col','source_pixel_row_y','source_pixel_col_x']);assert not pos.barcode.duplicated().any();pos['in_positions']=True
 p=pd.read_csv(paths['plot_spots'],sep='\t');p['barcode']=p.spot_id.str.rsplit('::',n=1).str[-1];assert not p.barcode.duplicated().any();p['in_plot']=True;p['plot_order_1based']=np.arange(1,len(p)+1);p=p.rename(columns={'x':'plot_source_x','y':'plot_source_y','library_size_original':'plotted_original_library_size','expr__'+gene:'plotted_normalised_value'})
 z=pos.merge(bc,on='barcode',how='outer',validate='one_to_one').merge(p,on='barcode',how='outer',validate='one_to_one')
 for c in ['in_matrix','in_positions','in_plot']:z[c]=z[c].fillna(False).astype(bool)
 sf=json.loads(paths['scalefactors_json.json'].read_text());par=next(x for x in params if x['gene']==gene);scale=sf['tissue_hires_scalef'];im=Image.open(paths['tissue_hires_image.png']);assert scale==par['scale']
 z['gene']=gene;z['feature_id']=target.iloc[0].feature_id;z['hires_scale']=scale;z['source_hires_pixel_x']=z.source_pixel_col_x*scale;z['source_hires_pixel_y']=z.source_pixel_row_y*scale;z['plot_hires_pixel_x']=z.plot_source_x*scale;z['plot_hires_pixel_y']=z.plot_source_y*scale
 z['raw_library_minus_plotted_library']=z.matrix_total_UMI_all_features-z.plotted_original_library_size;z['audit_normalised_minus_plotted']=z.audit_log1p_1e4_raw_over_library-z.plotted_normalised_value;z['fullres_x_difference']=z.plot_source_x-z.source_pixel_col_x;z['fullres_y_difference']=z.plot_source_y-z.source_pixel_row_y
 x0,x1=par['xlim'];y1,y0=par['ylim'];z['plot_left_fraction']=(z.plot_hires_pixel_x-x0)/(x1-x0);z['plot_top_fraction']=(z.plot_hires_pixel_y-y0)/(y1-y0)
 # Existing 39x39 mm axes with aspect='equal': derive archived-page point centres, no rendering changes.
 ratio=(x1-x0)/(y1-y0);aw=min(39.,39.*ratio);ah=min(39.,39./ratio);left=5+47*ind+(39-aw)/2;top=10+(39-ah)/2
 z['archived_page_x_mm']=left+aw*z.plot_left_fraction;z['archived_page_y_mm']=top+ah*z.plot_top_fraction
 z['not_plotted_status']=np.select([z.in_plot,~z.in_matrix,~z.in_positions,z.in_tissue.eq(0)],['plotted','not_in_matrix','not_in_positions','outside_tissue'],default='in_tissue_matrix_not_plotted')
 z.to_csv(O/f'audit/spatial_counts_{gene}_all_barcode_universe.tsv',sep='\t',index=False)
 plotted=z[z.in_plot].sort_values('plot_order_1based');plotted.to_csv(O/f'audit/spatial_counts_{gene}_plotted_spots.tsv',sep='\t',index=False);z[~z.in_plot].to_csv(O/f'audit/spatial_counts_{gene}_not_plotted.tsv',sep='\t',index=False);z[(~z.in_matrix)|(~z.in_positions)].to_csv(O/f'audit/spatial_counts_{gene}_missing_matrix_or_positions.tsv',sep='\t',index=False)
 r={'gene':gene,'section':stem,'matrix_shape':list(mat.shape),'matrix_nnz':int(mat.nnz),'feature_types':features.feature_type.value_counts().to_dict(),'duplicate_gene_symbol_rows':int(features.symbol.duplicated(keep=False).sum()),'target_feature':target.to_dict('records'),'n_position_spots':len(pos),'n_matrix_spots':len(bc),'n_plotted_spots':len(plotted),'expected_plotted':expected,'n_positions_in_tissue':int(pos.in_tissue.eq(1).sum()),'plotted_in_tissue_counts':plotted.in_tissue.value_counts(dropna=False).to_dict(),'n_plotted_not_in_matrix':int((~plotted.in_matrix).sum()),'n_plotted_not_in_positions':int((~plotted.in_positions).sum()),'n_in_tissue_matrix_not_plotted':int((z.in_tissue.eq(1)&z.in_matrix&~z.in_plot).sum()),'n_matrix_not_plotted':int((z.in_matrix&~z.in_plot).sum()),'n_positions_not_in_matrix':int((z.in_positions&~z.in_matrix).sum()),'n_raw_positive_plotted':int(plotted.raw_gene_count.gt(0).sum()),'n_normalised_positive_plotted':int(plotted.plotted_normalised_value.gt(0).sum()),'raw_count_sum_plotted':int(plotted.raw_gene_count.sum()),'raw_count_max_plotted':int(plotted.raw_gene_count.max()),'max_abs_library_difference':float(plotted.raw_library_minus_plotted_library.abs().max()),'max_abs_normalised_difference':float(plotted.audit_normalised_minus_plotted.abs().max()),'normalised_matches_atol_1e_12':bool(np.allclose(plotted.audit_log1p_1e4_raw_over_library,plotted.plotted_normalised_value,atol=1e-12,rtol=0)),'max_abs_fullres_x_difference':float(plotted.fullres_x_difference.abs().max()),'max_abs_fullres_y_difference':float(plotted.fullres_y_difference.abs().max()),'hires_scale':scale,'image_size_pixels':list(im.size),'xlim':par['xlim'],'ylim':par['ylim'],'n_outside_display_limits':int(((plotted.plot_left_fraction<0)|(plotted.plot_left_fraction>1)|(plotted.plot_top_fraction<0)|(plotted.plot_top_fraction>1)).sum()),'coordinate_convention':'x=source pixel column; y=source pixel row; multiplied by positive hires scale; imshow upper origin and descending y limits. No swap, mirror, rotation or shift in audited plot transform.','histology_alignment':'NOT established by numerical coordinate audit; parent performing independent visual/ROI audit'}
 results.append(r);print(json.dumps(r,ensure_ascii=False,indent=2),flush=True)
(O/'audit/spatial_counts_summary.json').write_text(json.dumps(results,ensure_ascii=False,indent=2));(O/'audit/spatial_counts_input_sha256.json').write_text(json.dumps(hashes,ensure_ascii=False,indent=2))
