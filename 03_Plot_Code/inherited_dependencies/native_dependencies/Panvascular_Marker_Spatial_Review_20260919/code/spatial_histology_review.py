from pathlib import Path
import json
import numpy as np,pandas as pd
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle,Polygon,Rectangle
from matplotlib.colors import LinearSegmentedColormap,Normalize
B=Path('/Users/taki/Desktop/Panvascular_Content_Alignment_ST_Audit_20260919');O=Path('/Users/taki/Desktop/Panvascular_Marker_Spatial_Review_20260919');R=Path('/Users/taki/Desktop/Codex任务1/spatial_visium_geo_local/extracted');D=Path('/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Nature_Reconstruction_20260918/SourceData/spatial');params=json.loads((D/'Figure5a_display_parameters.json').read_text());rois=json.loads((B/'audit/spatial_roi_definitions.json').read_text())['rois']
plt.rcParams.update({'font.family':'Arial','font.size':8,'pdf.fonttype':42,'savefig.facecolor':'white'})
cm=LinearSegmentedColormap.from_list('original_expression',['#F2B447','#D46A22','#A32C36','#541B43']);manifest=[]
def source(gene):
 stem='GSM7917895_CA2' if gene=='LIPA' else 'GSM9414344_LCA';gse='GSE248608' if gene=='LIPA' else 'GSE314851';return Image.open(R/gse/(stem+'_tissue_hires_image.png')),next(p for p in params if p['gene']==gene)
def footer(fig,text):fig.text(.04,.045,text,fontsize=7,va='bottom',linespacing=1.4)
def emit(pdf,fig,name):
 pdf.savefig(fig);fig.savefig(O/f'review/{name}.png',dpi=180);plt.close(fig)
with PdfPages(O/'review/Histology_ROI_three_view_review.pdf') as pdf:
 # Original whole tissue context; only fixed ROI outlines, no expression data.
 f=plt.figure(figsize=(10,7));a=f.add_axes([.08,.13,.58,.79]);im,par=source('LIPA');a.imshow(im);a.set_axis_off()
 for i,r in enumerate([r for r in rois if r['gene']=='LIPA' and r['role']=='suspected_acellular_space'],1):
  v=np.asarray(r['polygon']);a.add_patch(Polygon(v,closed=True,fill=False,ec='#BD2850' if i==1 else '#267E8E',lw=1));a.annotate(str(i),v.mean(axis=0),xytext=(v[:,0].max()+55,v[:,1].min()-25),color='#222222',fontsize=10,arrowprops={'arrowstyle':'-','lw':.7})
 f.text(.69,.77,'LIPA · GSM7917895_CA2',fontsize=12,weight='bold');f.text(.69,.68,'Original whole-tissue H&E\n\n1  Fixed central clear-space ROI\n\nROI coordinates and selection unchanged.',fontsize=9,linespacing=1.7)
 footer(f,'Review regions only; not pathologist-confirmed lumen labels.\nNo expression overlay on this localization page.');emit(pdf,f,'histology_00_LIPA_localization')
 for j,r in enumerate([r for r in rois if r['role']=='suspected_acellular_space'],1):
  g=r['gene'];im,par=source(g);v=np.asarray(r['polygon'],float);p=pd.read_csv(B/f'audit/spatial_counts_roi_{r["id"]}_all_positions.tsv',sep='\t');q=p[p.center_in_roi].copy().sort_values(['source_hires_pixel_y','source_hires_pixel_x']).reset_index(drop=True);q['review_spot_number']=np.arange(1,len(q)+1);q.to_csv(O/f'review/{r["id"]}_number_key.tsv',sep='\t',index=False)
  radius=float(p.physical_spot_radius_hires_pixels.iloc[0]);pad=max(12,radius*1.6);x0,x1=v[:,0].min()-pad,v[:,0].max()+pad;y0,y1=v[:,1].min()-pad,v[:,1].max()+pad
  f=plt.figure(figsize=(12,7.7));f.suptitle(f'{j}. {r["id"]}',x=.035,ha='left',fontsize=12,weight='bold',y=.965)
  f.text(.035,.91,'Fixed ROI on original H&E; comparison role: '+r['role'].replace('_',' '),fontsize=8)
  axes=[]
  for k,title in enumerate(['Unmodified H&E','H&E + estimated capture disks','Original expression values (overlay)']):
   a=f.add_axes([.045+k*.323,.37,.285,.48]);a.imshow(im,alpha=1 if k<2 else .68);a.set(xlim=(x0,x1),ylim=(y1,y0),aspect='equal');a.set_title(title,fontsize=9,pad=9);a.set_xlabel('Hires pixel x');a.set_ylabel('');a.tick_params(labelsize=7);a.tick_params(axis='y',left=False,labelleft=False);axes.append(a)
   if k: a.add_patch(Polygon(v,closed=True,fill=False,ec='#197C90',lw=.7,ls='--'))
  for row in q.itertuples():
   x=row.source_hires_pixel_x;y=row.source_hires_pixel_y;num=str(row.review_spot_number)
   axes[1].add_patch(Circle((x,y),radius,fill=False,ec='#BC2448',lw=.9,ls='-' if row.in_plot else '--'))
  # Identical outside-frame numbering in all views; fine leaders do not cover with labels.
  for side in [0,1]:
   group=q[(q.source_hires_pixel_x<(x0+x1)/2) if side==0 else (q.source_hires_pixel_x>=(x0+x1)/2)].sort_values('source_hires_pixel_y')
   last=-1.
   for row in group.itertuples():
    frac=(row.source_hires_pixel_y-y0)/(y1-y0);frac=max(frac,last+.075);frac=min(frac,.97);last=frac
    for a in axes:a.annotate(str(row.review_spot_number),xy=(row.source_hires_pixel_x,row.source_hires_pixel_y),xytext=(-.04 if side==0 else 1.04,1-frac),textcoords='axes fraction',ha='right' if side==0 else 'left',va='center',fontsize=7,color='#101010',arrowprops={'arrowstyle':'-','lw':.32,'color':'#444444'},annotation_clip=False,zorder=7)
  # Recreate original expression point geometry in hires units; this is explicitly not capture geometry.
  X=par['xlim'][1]-par['xlim'][0];Y=par['ylim'][0]-par['ylim'][1];width_mm=min(39.,39.*X/Y);px_per_pt=X/(width_mm*72/25.4);norm=Normalize(par['positive_cap02'],par['positive_cap98'],clip=True)
  visible=p[p.in_plot & p.source_hires_pixel_x.between(x0-20,x1+20)&p.source_hires_pixel_y.between(y0-20,y1+20)]
  for row in visible.itertuples():
   val=row.plotted_normalised_value;rr=np.sqrt(2.65 if val>0 else 1.45)/2*px_per_pt;color=cm(norm(val)) if val>0 else '#CED3D6';axes[2].add_patch(Circle((row.source_hires_pixel_x,row.source_hires_pixel_y),rr,fc=color,ec='#542037' if val>0 else 'none',lw=.05,zorder=2))
  # Move number objects over expression markers.
  for t in axes[2].texts:t.set_zorder(5)
  cbax=f.add_axes([.728,.305,.17,.013]);cb=f.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cm),cax=cbax,orientation='horizontal');cb.set_ticks([par['positive_cap02'],par['positive_cap98']]);cb.set_ticklabels([f'≤{par["positive_cap02"]:.3g}',f'≥{par["positive_cap98"]:.3g}']);cb.ax.tick_params(labelsize=7,length=2);cb.set_label('Original log1p(CP10k); positive-only scale',fontsize=7,labelpad=2)
  f.text(.04,.265,f'Estimated physical capture radius: {radius:.3f} hires px. Red circles use original spot_diameter_fullres × hires_scale / 2.\nSolid / dashed circle: plotted / not plotted. Teal dashed polygon: fixed review ROI. Display dots (right) are not physical capture disks.',fontsize=7,va='top')
  # Compact numbered original-values key; no summary statistics recomputed.
  lines=['No.   Barcode                Raw UMI    Library     In plot    Full disk in ROI']
  for row in q.itertuples():lines.append(f'{row.review_spot_number:>2}    {row.barcode:<20} {row.raw_gene_count:>5g}       {row.matrix_total_UMI_all_features:>6g}       {str(row.in_plot):<5}      {str(row.entire_physical_spot_disk_inside_roi)}')
  cols=[lines[:8],(['No.   Barcode                Raw UMI    Library     In plot    Full disk in ROI']+lines[8:]) if len(lines)>8 else []]
  for k,ls in enumerate(cols):
   if ls:f.text(.04+k*.49,.19,'\n'.join(ls),fontsize=6.1,fontfamily='DejaVu Sans Mono',va='top',linespacing=1.35)
  footer(f,'H&E carries only outside-frame IDs and fine leaders; identical IDs map to barcodes in the companion TSV.\nOriginal-value overlay is a same-crop reconstruction, not a baseline screenshot. Display point size is not capture size; status remains unresolved.')
  emit(pdf,f,f'histology_{j:02}_{r["id"]}');manifest.append({'page':j+1,'roi':r['id'],'crop_hires_xyxy':[x0,y0,x1,y1],'physical_radius_px':radius,'original_expression_display':True,'value_source':str(B/f'audit/spatial_counts_roi_{r["id"]}_all_positions.tsv')})
(O/'notes/histology_review_manifest.json').write_text(json.dumps(manifest,indent=2))
print('Written 5-page PDF and 5 PNGs')
