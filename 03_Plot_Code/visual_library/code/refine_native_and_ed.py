"""Layout only: existing vector objects, saved marker statistics, original H&E.
No counts, model, integration, region, or expression computation is performed.
"""
from pathlib import Path
import json, re, copy, hashlib, base64
import numpy as np
import pandas as pd
from lxml import etree as E
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
import redraw_main as m
import export_native_panels as native
import pymupdf as fitz
import cairosvg

R=m.ROOT;P=m.P;S=m.S;Q=native.Q;MM=72/25.4

def emit(root,path):
 path.parent.mkdir(exist_ok=True)
 E.ElementTree(root).write(str(path)+'.svg',encoding='utf-8',xml_declaration=True)
 cairosvg.svg2pdf(bytestring=E.tostring(root),write_to=str(path)+'.pdf')
 with fitz.open(str(path)+'.pdf') as d:d[0].get_pixmap(matrix=fitz.Matrix(2,2)).save(str(path)+'.png')

def affine(s):
 M=np.eye(3)
 for name,args in re.findall(r'(matrix|translate|scale|rotate)\(([^)]+)\)',s or ''):
  v=[float(x) for x in re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?',args)];T=np.eye(3)
  if name=='translate':T[0,2]=v[0];T[1,2]=v[1] if len(v)>1 else 0
  elif name=='scale':T[0,0]=v[0];T[1,1]=v[-1]
  elif name=='matrix':T=np.array([[v[0],v[2],v[4]],[v[1],v[3],v[5]],[0,0,1]])
  elif name=='rotate':
   rad=np.deg2rad(v[0]);T[:2,:2]=[[np.cos(rad),-np.sin(rad)],[np.sin(rad),np.cos(rad)]]
   if len(v)==3:T[:2,2]=np.array(v[1:])-T[:2,:2]@np.array(v[1:])
  M=M@T
 return M

def clipbox(root,ax):
 el=next(e for e in ax.iter() if e.get('clip-path'))
 ref=el.get('clip-path')[5:-1];rr=root.xpath('//*[@id=$id]',id=ref)[0].find(Q('rect'))
 x,y,w,h=[float(rr.get(k)) for k in ['x','y','width','height']]
 M=np.eye(3)
 for a in reversed([ax]+list(ax.iterancestors())):M=M@affine(a.get('transform'))
 points=M@np.array([[x,x+w],[y,y+h],[1,1]])
 return [points[0,0],points[1,0],points[0,1],points[1,1]],(x,y,w,h)

def spatial():
 src=R/'review/本轮修改前/Figure5/Figure5a.svg';root=E.parse(str(src)).getroot()
 table=S/'Current_Display_20260920/Figure5a_current_sections.tsv';meta=pd.read_csv(table,sep='\t')
 # Replace only the tiny embedded display thumbnails with byte-identical
 # original H&E files. The inherited complete field rectangle is retained.
 proofs=[]
 for n,(_,row) in zip([1,6],meta.iterrows()):
  ax=root.xpath('//*[@id=$id]',id=f'r2_axes_{n}')[0];bb,(x,y,w,h)=clipbox(root,ax)
  file=Path(row.image);raw=file.read_bytes();assert hashlib.sha256(raw).hexdigest()==row.image_sha256
  im=next(ax.iter(Q('image')));im.attrib.pop('transform',None)
  for k,v in zip(['x','y','width','height'],[x,y,w,h]):im.set(k,str(v))
  im.set('{http://www.w3.org/1999/xlink}href','data:image/png;base64,'+base64.b64encode(raw).decode())
  proofs.append({'section':row.full_section_id,'image_sha256':row.image_sha256,'full_field_rectangle':list((x,y,w,h)),'spots':int(row.n_existing_spots),'image_enhancement':False})
 out=E.Element(Q('svg'),nsmap=root.nsmap,width=f'{153*MM}pt',height=f'{105*MM}pt',viewBox=f'0 0 {153*MM} {105*MM}')
 defs=E.SubElement(out,Q('defs'));seen=set()
 def add(ids,box,x,y,scale):
  part=native.extract(root,ids)
  for dd in part.findall(Q('defs')):
   for ch in list(dd):
    key=ch.get('id') or E.tostring(ch)
    if key not in seen:defs.append(ch);seen.add(key)
  g=E.SubElement(out,Q('g'),transform=f'translate({x} {y}) scale({scale}) translate({-box[0]} {-box[1]})')
  for ch in list(part):
   if ch.tag!=Q('defs'):g.append(ch)
 def text(x,y,s,size=8,anchor='middle'):
  e=E.SubElement(out,Q('text'),x=str(x*MM),y=str(y*MM),style=f'font:{size}px Arial;fill:#20252A;text-anchor:{anchor}');e.text=s
 text(1.5,4,'a',10,'start')
 for row,(ids,barids,top) in enumerate([([1,2,4],[3,5],13),([6,7,9],[8,10],59)]):
  boxes=[clipbox(root,root.xpath('//*[@id=$id]',id=f'r2_axes_{n}')[0])[0] for n in ids]
  scale=36*MM/(boxes[0][3]-boxes[0][1]);assert all(abs((b[3]-b[1])*scale-36*MM)<1e-6 for b in boxes)
  for j,(n,b,x) in enumerate(zip(ids,boxes,[3,49,105])):
   add([f'r2_axes_{n}'],b,x*MM,top*MM,scale)
   centre=x+(b[2]-b[0])*scale/MM/2
   title=f'H&E · {meta.iloc[row].study} / '+('CA1' if row==0 else 'T1077') if j==0 else ['','Inferred Fib abundance','FN1'][j]
   text(centre,top-3,title)
   if row==1 and j>0:text(centre,99,['','Σ marginal q05','log1p(CP10k)'][j])
  for n,x in zip(barids,[87,143]):
   sub=native.extract(root,[f'r2_axes_{n}']);b=native.content_bbox(sub)
   # Keep bar and its original tick labels together; no value/scale changes.
   add([f'r2_axes_{n}'],b,x*MM,(top+3)*MM,1.0)
 text(20,101,'Zero',8)
 E.SubElement(out,Q('circle'),cx=str(11*MM),cy=str(100*MM),r='1.7',fill='#8b8f95')
 before_uses=[tuple(e.get(k) for k in ['x','y','style','{http://www.w3.org/1999/xlink}href']) for e in root.iter(Q('use'))]
 after_uses=[tuple(e.get(k) for k in ['x','y','style','{http://www.w3.org/1999/xlink}href']) for e in out.iter(Q('use'))]
 # Exclude only the old legend glyph; every original spot glyph is retained.
 from collections import Counter
 for key,count in Counter(after_uses).items():assert Counter(before_uses)[key]>=count
 for n in [2,4,7,9]:
  a=root.xpath('//*[@id=$id]',id=f'r2_axes_{n}')[0];b=out.xpath('//*[@id=$id]',id=f'r2_axes_{n}')[0]
  assert E.tostring(a)==E.tostring(b)
 emit(out,R/'Figure5/Figure5a')
 (R/'review/Figure5a_layout_image_check.json').write_text(json.dumps({'source':str(src),'source_table':str(table),'rows':proofs,'all_four_spot_axes_byte_identical':True,'mapping':'same full frame for all 3 views in each row; isotropic placement','color_limits_changed':False},indent=2))
 return {'figure':'Figure5','panel':'a','mode':'reflowed_native_geometry','source_files':[str(src),str(table)]+meta.image.tolist(),'change':'Larger 36mm views; compact columns; original high resolution H&E; unchanged full spot axes and colour bars','details':{'H_E_original_bytes_verified':True,'spot_axes_unchanged':True}}

def marker_ed():
 src=P/'Public_Release/Figure_Code/native_dependencies/Panvascular_Marker_Communication_Enrichment_20260920/source_data/display_marker_values.tsv';D=m.read(src)
 idx=m.read(S/'ExtendedData3_4_5/ED345_all_named_state_index.tsv').set_index('level2_annotation')
 mappings=[]
 for n,letter,ai in [(2,'b',1),(3,'b',1),(4,'b',1),(4,'d',4),(5,'b',1)]:
  name=f'ExtendedData_Figure{n}{letter}';genes,labels=m.h.axis_labels(f'ExtendedData_Figure{n}',f'newED{n}_axes_{ai}')
  z=D[D.display_label.isin(labels)&D.gene_id.isin(genes)].copy();assert len(z)==len(labels)*len(genes)
  states=[z.loc[z.display_label.eq(l),'state_id'].iloc[0] for l in labels]
  mean=z.pivot(index='state_id',columns='gene_id',values='mean_log1p').reindex(index=states,columns=genes)
  fraction=z.pivot(index='state_id',columns='gene_id',values='detection_fraction').reindex(index=states,columns=genes)
  caps=z[['vmin','vmax']].drop_duplicates();assert len(caps)==1
  lo,hi=caps.iloc[0].tolist();nr,nc=mean.shape;rh=3.7;cw=3.4
  fp=FontProperties(family='Arial',size=7.5)
  labelwidth=max(TextPath((0,0),s,prop=fp).get_extents().width/MM for s in labels)
  grouped=nr>20;left=labelwidth+(13 if grouped else 5);right=18;top=4;bottom=24
  W=left+cw*nc+right;H=top+rh*nr+bottom;f=m.fig(W,H);a=m.axmm(f,left,top,cw*nc,rh*nr)
  xx,yy=np.meshgrid(np.arange(nc),np.arange(nr));valid=mean.notna().to_numpy()&fraction.notna().to_numpy();smax=(.72*min(cw,rh)*MM)**2
  # Stored statistics and bounds only. 0 detection is exactly 0 area.
  artist=a.scatter(xx[valid],yy[valid],s=(smax*fraction.to_numpy())[valid],c=mean.to_numpy()[valid],cmap=m.h.BL,vmin=lo,vmax=hi,ec='none',lw=0)
  if (~valid).any():a.scatter(xx[~valid],yy[~valid],s=9,marker='x',c=m.GREY,lw=.6)
  a.set(xlim=(-.5,nc-.5),ylim=(nr-.5,-.5));a.set_yticks(range(nr),labels,fontsize=7.5);a.set_xticks(range(nc),genes,fontsize=7.5,rotation=90,fontstyle='italic');a.tick_params(length=0,pad=3);a.tick_params(axis='y',pad=8);a.spines[:].set_visible(False)
  parents=[]
  for i,st in enumerate(states):
   parent=st if st in m.L1 else idx.loc[st,'level1_annotation'];parents.append(parent)
   co=m.L1[st] if st in m.L1 else m.C['L2'][m.C['L2_name_to_id'][st]]
   a.add_patch(m.Rectangle((-.84,i-.43),.13,.86,fc=co,ec='none',clip_on=False))
  if grouped:
   for pa in dict.fromkeys(parents):
    ids=np.where(np.array(parents)==pa)[0];t=top+rh*ids.min();hh=rh*len(ids)
    f.add_artist(m.Rectangle((6/W,1-(t+hh-.3)/H),.75/W,(hh-.6)/H,transform=f.transFigure,fc=m.L1[pa],ec='none'))
    f.text(3/W,1-(t+hh/2)/H,pa.replace(' cell',''),rotation=90,ha='center',va='center',fontsize=7.5)
   for i in range(1,nr):
    if parents[i]!=parents[i-1]:a.axhline(i-.5,c='#E2E6E8',lw=.4,zorder=0)
  cb=f.colorbar(artist,cax=m.axmm(f,left+nc*cw+5,top+3,2.3,min(28,nr*rh-3)));cb.set_label('Mean log1p(CP10k)',fontsize=7.5,labelpad=4);cb.ax.tick_params(labelsize=7.5);cb.outline.set_visible(False)
  la=m.axmm(f,(W-72)/2,top+nr*rh+17,72,6);la.set(xlim=(0,100),ylim=(0,1));la.axis('off');la.text(0,.75,'Detected (%)',va='center',fontsize=7.5)
  for x,v in zip([44,61,78,95],[.25,.5,.75,1]):la.scatter(x,.75,s=smax*v,c=m.GREY,ec='none');la.text(x,-.05,str(round(v*100)),ha='center',va='top',fontsize=7.5)
  details={'rows':nr,'columns':nc,'cells':nr*nc,'vmin':float(lo),'vmax':float(hi),'area_map':'smax * original detection fraction; zero area at zero','area_max_pt2':smax,'row_pitch_mm':rh,'column_pitch_mm':cw,'missing_cells':int((~valid).sum()),'zero_cells':int((fraction==0).sum().sum()),'statistic':'saved mean_log1p; no rescaling','palette':'unchanged adopted Blues(.15,.95)','row_order':states,'gene_order':genes}
  m.save(f,name,[src,P/f'图件编辑资料/editable_svg/ExtendedData_Figure{n}.svg'],'Compact grouped row labels; larger 3.7mm rows / 3.4mm columns; same full marker matrix and original colour bounds',details)
  z.to_csv(R/'data_used'/f'{name}_unchanged_values.tsv',sep='\t',index=False)
  mappings.extend({'panel':name,'display_label':label,'state_id':st,'parent':pa} for label,st,pa in zip(labels,states,parents))
 pd.DataFrame(mappings).to_csv(R/'data_used/ED_marker_identity_map.tsv',sep='\t',index=False)

def main():
 rec=json.loads((R/'review/panel_manifest.json').read_text());extra=spatial();rec=[r for r in rec if not(r['figure']=='Figure5' and r['panel']=='a')]+[extra]
 marker_ed();rec+=m.RECORDS
 # Copy unchanged companion panels for self-contained ED2–5 review sheets.
 for n,letters in [(2,'a'),(3,'a'),(4,'ac'),(5,'a')]:
  for letter in letters:
   name=f'ExtendedData_Figure{n}{letter}';old=P/f'独立小图_最新版/ExtendedData_Figure{n}'
   dst=R/f'ExtendedData_Figure{n}';dst.mkdir(exist_ok=True)
   import shutil
   for ext in ['pdf','svg','png']:shutil.copy2(old/f'{name}.{ext}',dst/f'{name}.{ext}')
   rec.append({'figure':f'ExtendedData_Figure{n}','panel':letter,'mode':'preserved_native_geometry','source_files':[str(old/f'{name}.svg')],'change':'Existing companion map/heatmap preserved with exact data and colours','details':{}})
 (R/'review/panel_manifest.json').write_text(json.dumps(rec,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
