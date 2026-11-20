from pathlib import Path
import os,json,hashlib
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import numpy as np,pandas as pd,matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image
from scipy.io import mmread
O=Path('/Users/taki/Desktop/Panvascular_Spatial_Evidence_Assessment_20260920');A=Path('/Users/taki/Desktop/Panvascular_Content_Alignment_ST_Audit_20260919/audit');B=Path('/Users/taki/Desktop/Panvascular_Marker_Spatial_Review_20260919');RAW=Path('/Users/taki/Desktop/Codex任务1/spatial_visium_geo_local/extracted');JOINT=Path('/Users/taki/Desktop/Codex任务1/overnight_manuscript_integration_20260917/spatial_target_excluded_joint_v1/spot_joint_source.tsv.gz')
manifest=json.loads((B/'audit/Figure5a_display.json').read_text())['data_regions'];cols=['spot_id','section_id','library_size_original','x','y','expr__LIPA','expr__STAT3','abundance__Fibroblast','program__fibrous_cap_core'];joint=pd.concat([c[c.section_id.isin(['GSE248608::GSM7917895_CA2','GSE314851::GSM9414344_LCA'])] for c in pd.read_csv(JOINT,sep='\t',usecols=cols,chunksize=50000)],ignore_index=True);assert joint.spot_id.is_unique
mpl.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,'pdf.fonttype':42,'svg.fonttype':'none'})
cm=mpl.colors.LinearSegmentedColormap.from_list('frozen_target_sequence',['#482878','#3B528B','#2196C4','#35B779','#FDE725','#F68B32','#B2182B']);audit=[];feature_rows=[]
with PdfPages(O/'review/same_section_diagnostic_2pages.pdf') as pp:
 for spec in manifest:
  target=spec['gene'];study='GSE248608' if target=='LIPA' else 'GSE314851';stem=spec['section'];d=RAW/study;z=pd.read_csv(A/f'spatial_counts_{target}_plotted_spots.tsv',sep='\t',float_precision='round_trip');n=1254 if target=='LIPA' else 2072;assert len(z)==n and z.barcode.is_unique
  bc=pd.read_csv(d/f'{stem}_barcodes.tsv',sep='\t',header=None)[0];ft=pd.read_csv(d/f'{stem}_features.tsv',sep='\t',header=None,names=['id','symbol','type']);mat=mmread(d/f'{stem}_matrix.mtx').tocsr();ix=z.matrix_column_1based.to_numpy(dtype=int)-1;assert np.array_equal(bc.iloc[ix].to_numpy(),z.barcode.to_numpy());den=z.matrix_total_UMI_all_features.to_numpy();assert np.array_equal(np.asarray(mat[:,ix].sum(0)).ravel(),den)
  assert np.array_equal(np.asarray((mat[:,ix]>0).sum(0)).ravel(),z.matrix_detected_features.to_numpy())
  for g in [target,'DCN','COL1A1','LST1','ACTA2']:
   rows=np.where(ft.symbol.eq(g)&ft.type.eq('Gene Expression'))[0]
   if len(rows)>1:raise ValueError(f'Duplicate requested marker {g}: unresolved')
   vals=np.asarray(mat[rows[0],ix].toarray()).ravel() if len(rows) else np.full(n,np.nan)
   z['raw__'+g]=vals;z['log1p_CP10k__'+g]=np.log1p(1e4*vals/den);feature_rows.append(dict(section=stem,gene=g,found=bool(len(rows)),feature_id=ft.iloc[rows[0]].id if len(rows) else None,matrix_row_1based=int(rows[0]+1) if len(rows) else None))
  assert np.allclose(z['log1p_CP10k__'+target],z.plotted_normalised_value,atol=1e-12)
  z=z.merge(joint.drop(columns='section_id'),on='spot_id',how='left',validate='one_to_one',indicator=True);assert z['_merge'].eq('both').all();assert np.allclose(z.library_size_original,den);assert np.allclose(z.x,z.source_pixel_col_x);assert np.allclose(z.y,z.source_pixel_row_y);assert np.allclose(z['expr__'+target],z['log1p_CP10k__'+target])
  z['histology_uncertain_original_flag']=z.barcode.isin(spec['flagged_barcodes']);assert z.histology_uncertain_original_flag.sum()==(7 if target=='LIPA' else 0)
  z.to_csv(O/'tables'/f'diagnostic_{target}_same_section_source.tsv',sep='\t',index=False,na_rep='NA')
  im=Image.open(d/f'{stem}_tissue_hires_image.png');W,H=im.size;x=z.source_hires_pixel_x.to_numpy();y=z.source_hires_pixel_y.to_numpy();flag=z.histology_uncertain_original_flag.to_numpy()
  f,axs=plt.subplots(2,5,figsize=(350/25.4,180/25.4));f.subplots_adjust(left=.015,right=.985,top=.90,bottom=.14,wspace=.14,hspace=.35)
  f.suptitle(f'{target} | {study} / {stem} | all {n:,} frozen plotted spots | full source image field',fontsize=12,y=.98)
  panels=[('H&E (unmodified full image)',None),('Total UMI','matrix_total_UMI_all_features'),('Detected genes','matrix_detected_features'),(target,'log1p_CP10k__'+target),('DCN (Fib reference)','log1p_CP10k__DCN'),('COL1A1 (Fib reference)','log1p_CP10k__COL1A1'),('LST1 (myeloid reference)','log1p_CP10k__LST1'),('ACTA2 (SMC reference)','log1p_CP10k__ACTA2'),('Fib abundance (q05; model-inferred)','abundance__Fibroblast'),('Fibrous-cap expression program','program__fibrous_cap_core')]
  records=[]
  for a,(title,key) in zip(axs.ravel(),panels):
   a.set_title(title,pad=5);a.set_xlim(0,W);a.set_ylim(H,0);a.set_aspect('equal');a.set_xticks([]);a.set_yticks([])
   if key is None:a.imshow(im,origin='upper',extent=[0,W,H,0],interpolation='nearest');continue
   a.set_facecolor('#050708');v=z[key].to_numpy(float);finite=np.isfinite(v);expression=key.startswith('log1p');zero=finite&(v==0)
   if key=='log1p_CP10k__'+target:lo,hi=spec['vmin'],spec['vmax'];palette=cm;label='log1p(CP10k); frozen limits'
   elif expression:lo,hi=0,float(np.nanmax(v)) if finite.any() else 1;palette=cm;label='log1p(CP10k)'
   else:lo=0 if key.startswith('matrix_') else float(np.nanmin(v));hi=float(np.nanmax(v));palette='viridis';label='UMIs' if 'total_UMI' in key else 'Genes' if 'detected_features' in key else 'q05 abundance' if key=='abundance__Fibroblast' else 'Existing program score'
   if hi<=lo:hi=lo+1
   norm=mpl.colors.Normalize(lo,hi,clip=False);cmap=mpl.colormaps[palette] if isinstance(palette,str) else palette
   normal=finite&~flag;colours=cmap(norm(v));
   if expression:colours[zero]=mpl.colors.to_rgba('#8B8F95')
   a.scatter(x[normal],y[normal],s=1.5,c=colours[normal],lw=0,rasterized=False)
   a.scatter(x[flag&finite],y[flag&finite],s=9,facecolors='none',edgecolors=colours[flag&finite],lw=.6)
   a.scatter(x[~finite],y[~finite],s=5,c='#f0f0f0',marker='x',lw=.4)
   ext='both' if finite.any() and v[finite].min()<lo and v[finite].max()>hi else 'max' if finite.any() and v[finite].max()>hi else 'neither'
   cb=f.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmap),ax=a,orientation='horizontal',fraction=.045,pad=.025,extend=ext);cb.ax.tick_params(labelsize=7,length=2);cb.set_label(label,fontsize=7,labelpad=2)
   records.append(dict(panel=title,key=key,n=n,n_finite=int(finite.sum()),n_zero=int(zero.sum()),n_missing=int((~finite).sum()),vmin=lo,vmax=hi,colourbar_extend=ext))
  ref=axs.ravel()[1].get_position(); hp=axs.ravel()[0].get_position(); axs.ravel()[0].set_position([hp.x0+hp.width/2-ref.width/2,ref.y0,ref.width,ref.height])
  # Equal physical image-to-page scale across H&E and all spot/QC maps.
  assert all(np.allclose(a.get_position().bounds[2:],ref.bounds[2:],atol=1e-8) for a in axs.ravel())
  f.text(.02,.075,'H&E, QC and expression are separate views of the same full field; identical orientation and coordinate limits. No smoothing, imputation or spot exclusion.',fontsize=8)
  f.text(.02,.05,'Grey expression dots = observed zero; x = missing. Hollow spots retain the 7 original LIPA histology-uncertain flags; all STAT3 boundary spots remain.',fontsize=8)
  f.text(.02,.025,'Reference markers are descriptive. q05 abundance is model-inferred, not staining; fibrous-cap program is expression-based, not a histological compartment.',fontsize=8)
  f.savefig(O/'review'/f'diagnostic_{target}_full_field.pdf');f.savefig(O/'review'/f'diagnostic_{target}_full_field.png',dpi=180);pp.savefig(f);plt.close(f)
  audit.append(dict(gene=target,section=stem,n_plotted=n,n_original_flags=int(flag.sum()),full_image_size=[W,H],xlim=[0,W],ylim=[H,0],source_matrix=str(d/f'{stem}_matrix.mtx'),source_audit=str(A/f'spatial_counts_{target}_plotted_spots.tsv'),join_source=str(JOINT),join_key='spot_id',join_complete=True,normalization='log1p(10000 * raw / original all-feature library)',target_and_qc_exact_source_check='PASS',panels=records))
  print(target,n,flush=True)
(O/'tables/diagnostic_gene_feature_dictionary.tsv').write_text(pd.DataFrame(feature_rows).to_csv(sep='\t',index=False,na_rep='NA'));(O/'notes/diagnostic_render_manifest.json').write_text(json.dumps(audit,indent=2))
