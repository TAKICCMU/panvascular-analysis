"""Resolution-only rerender of Figure1b's frozen coordinate table; no embedding computation."""
from redraw_data_panels import *
import io,base64
from PIL import Image
from export_native_panels import Q
from reflow_native_maps import reflow
p=S/'Figure1/full_cells_frozen_umap.tsv.gz';d=read(p);assert len(d)==1175115
mfile=S/'Current_Display_20260920/Figure1b_tissue_sample_mapping.tsv';m=read(mfile);assert m.sample_uid.is_unique
assert d.sample_uid.isin(m.sample_uid).all();d['Tissue']=d.sample_uid.map(m.set_index('sample_uid').Tissue)
shown=d[d.level1_annotation.ne('Unknown')];assert len(shown)==1117197
xmin,xmax=d.UMAP1.min(),d.UMAP1.max();ymin,ymax=d.UMAP2.min(),d.UMAP2.max();dx=xmax-xmin;dy=ymax-ymin
r=E.parse(str(P/'图件编辑资料/editable_svg/Figure1.svg')).getroot();audit=[]
for id,x,y,w,h,field,pal in [('axes_2',34,27,49,33,'level1_annotation',C['L1']),('axes_3',125,23,51,43,'Tissue',C['tissue'])]:
 cols=shown[field].map(pal);assert cols.notna().all(),shown[field].unique()
 f=fig(183,170);a=axmm(f,x,y,w,h);a.set(xlim=(xmin-.04*dx,xmax+.04*dx),ylim=(ymin-.04*dy,ymax+.04*dy));a.set_aspect('equal');a.axis('off');a.scatter(shown.UMAP1,shown.UMAP2,c=cols,s=.055,linewidths=0,alpha=.8,rasterized=True)
 buf=io.BytesIO();f.savefig(buf,format='svg',dpi=600);plt.close(f);new=E.fromstring(buf.getvalue()).find('.//'+Q('image'));old=r.xpath('//*[@id=$v]',v=id)[0].find(Q('image'));href='{http://www.w3.org/1999/xlink}href'
 image=Image.open(io.BytesIO(base64.b64decode(new.get(href).split(',')[1])))
 # Raster crop bounds from the identical axes/data transform are taken from this600dpi rendering.
 for k in ['x','y','width','height','transform',href]:old.set(k,new.get(k))
 audit.append({'object':id,'n_cells':len(shown),'omitted_Unknown':len(d)-len(shown),'pixel_dimensions':image.size,'dpi':600,'point_area_pt2':.055,'alpha':.8,'ordered_coordinate_source':str(p),'counts_labels_embedding_unchanged':True});print(id,image.size,flush=True)
out=O/'review/Figure1_highres_native.svg';E.ElementTree(r).write(str(out),encoding='utf-8',xml_declaration=True)
rec=reflow('Figure1b',[['axes_2','axes_4'],['axes_3','text_148']],2,[255,195],source_override=out);rec['mode']='frozen_coordinates_high_resolution_rerender';rec['source_files']=[str(p),str(mfile)];rec['change']='600dpi point raster from identical coordinates, order, mask, area, alpha and palette; existing label anchors retained';rec['details']=audit
(O/'review/highres_manifest.json').write_text(json.dumps([rec],ensure_ascii=False,indent=2));(O/'review/highres_source_hashes.json').write_text(json.dumps(sources,ensure_ascii=False,indent=2))
