"""Restore ordinary marker dots from frozen SVG paired tile/path objects.
No analysis, sorting, coordinate edits, scaling or panel centering.
"""
from pathlib import Path
import os,sys,re,json,copy,io,hashlib
os.environ.setdefault('DYLD_FALLBACK_LIBRARY_PATH','/Users/taki/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/lib')
sys.path.insert(0,'/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Local_Touchup_20260918/tools')
import xml.etree.ElementTree as E
import numpy as np,pandas as pd,matplotlib as mpl,cairosvg,cairosvg.image,fitz
from cairosvg.surface import cairo
from PIL import Image,ImageDraw
cairosvg.image.IMAGE_RENDERING.update({'pixelated':cairo.FILTER_NEAREST,'crisp-edges':cairo.FILTER_NEAREST})
B=Path('/Users/taki/Desktop/Panvascular_Content_Alignment_ST_Audit_20260919');O=Path('/Users/taki/Desktop/Panvascular_Marker_Spatial_Review_20260919');A=Path('/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Full_Layout_Correction_20260918');N={'s':'http://www.w3.org/2000/svg'};S=N['s'];E.register_namespace('',S);E.register_namespace('xlink','http://www.w3.org/1999/xlink')
num=lambda s:list(map(float,re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?',s)))
checks=[];maps=[];bbrows=[];prior=json.loads((B/'audit/alignment_panels.json').read_text());caps=pd.read_csv(A/'audit/atlas/independent_checks.tsv',sep='\t').dropna(subset=['cohort']).set_index('cohort').cap_Q98.to_dict()

def bbox(root,elements):
 rr=E.Element(root.tag,root.attrib);defs=E.SubElement(rr,'{'+S+'}defs')
 for d in root.findall('.//s:defs',N):
  for c in d:defs.append(copy.deepcopy(c))
 for e in elements:rr.append(copy.deepcopy(e))
 pix=cairosvg.svg2png(bytestring=E.tostring(rr),background_color='white',scale=1.5);ar=np.asarray(Image.open(io.BytesIO(pix)).convert('RGB'));yy,xx=np.where(ar.min(2)<245)
 return [float(xx.min())*.5,float(yy.min())*.5,float(xx.max()+1)*.5,float(yy.max()+1)*.5] if len(xx) else None

def context(root,el):
 parents={c:p for p in root.iter() for c in p};out=copy.deepcopy(el);p=parents[el]
 while p is not root:
  w=E.Element(p.tag,p.attrib);w.append(out);out=w;p=parents[p]
 return out

for no in [2,3,4,5]:
 name=f'ExtendedData_Figure{no}';tree=E.parse(B/'editable_svg'/f'{name}.svg');root=tree.getroot();original=copy.deepcopy(root);cols=[g for g in root.findall('.//s:g',N) if g.get('id','').startswith('PathCollection_')]
 for ci,col in enumerate(cols):
  tilegroup=next(g for g in col if g.get('id','').startswith('FigureYa278_tiles_'));tiles=tilegroup.findall('s:rect',N);paths=col.findall('s:path',N);assert len(tiles)==len(paths)
  if no==2:
   src=A/'panels/ExtendedData_Figure2_b_source.tsv';allv=pd.read_csv(src,sep='\t');groups=list(allv.groupby('L1',sort=False));lineage,v=groups[ci];v=v[np.isfinite(v.mean_log1p)].reset_index(drop=True);area_scale=35;values=v.display_scaled_mean.to_numpy();cmap='RdBu_r';vmin,vmax=-2,2;panel='b';cohort=lineage
  else:
   panel='d' if no==4 and ci==1 else 'b';src=A/f'panels/atlas/{name}_{panel}_full_matrix.tsv';v=pd.read_csv(src,sep='\t');area_scale=17;cohort=v.cohort.iloc[0];values=v.mean_log1p.to_numpy();cmap='Blues';vmin,vmax=0,caps[cohort]
  assert len(v)==len(paths)
  expected=[mpl.colors.to_hex(mpl.colormaps[cmap](mpl.colors.Normalize(vmin,vmax)(val))) for val in values]
  coords=[];widths=[];maxerr=0
  for i,(p,t) in enumerate(zip(paths,tiles)):
   xy=np.array(num(p.get('d'))).reshape(-1,2);x=(xy[:,0].min()+xy[:,0].max())/2;y=(xy[:,1].min()+xy[:,1].max())/2;diam=xy[:,0].max()-xy[:,0].min();err=abs(diam**2-area_scale*v.detection_fraction.iloc[i]);maxerr=max(maxerr,err);assert err<.0001,(name,ci,i,err)
   assert t.get('fill').lower()==expected[i],(name,ci,i,t.get('fill'),expected[i])
   p.set('style',f'fill: {t.get("fill")}; stroke: none')
   coords.append((x,y));widths.append(diam)
   maps.append(dict(figure=name,panel=panel,collection=col.get('id'),source=str(src),source_row=i,mean_log1p=v.mean_log1p.iloc[i],colour_value=values[i],detection_fraction=v.detection_fraction.iloc[i],area_parameter=area_scale*v.detection_fraction.iloc[i],diameter_pt=diam,fill=t.get('fill'),is_zero=bool(v.detection_fraction.iloc[i]==0)))
  col.remove(tilegroup)
  xvals=sorted(set(round(x,4) for x,y in coords));yvals=sorted(set(round(y,4) for x,y in coords));pitchx=min([b-a for a,b in zip(xvals,xvals[1:]) if b-a>.01],default=999);pitchy=min([b-a for a,b in zip(yvals,yvals[1:]) if b-a>.01],default=999)
  if no==2:pitchy=9.5*72/25.4 # frozen repeated facet row pitch in original ED2 code.
  assert max(widths)<min(pitchx,pitchy),(name,ci,max(widths),pitchx,pitchy)
  checks.append(dict(figure=name,panel=panel,collection=col.get('id'),units=len(paths),source=str(src),sha256=hashlib.sha256(src.read_bytes()).hexdigest(),zero_units=int((v.detection_fraction==0).sum()),area_scale=area_scale,max_area_error=maxerr,max_diameter_pt=max(widths),row_pitch_pt=pitchy,col_pitch_pt=pitchx,min_clearance_pt=min(pitchx,pitchy)-max(widths),colour_limits=[vmin,vmax],colour_semantics='within-lineage scaled mean' if no==2 else 'mean log1p(CP10k)',exact_source_colour='PASS'))
 # Only size legends switch from hollow to filled, using a neutral display colour.
 for leg in root.findall('.//s:g',N):
  if not leg.get('id','').startswith('legend_'):continue
  for e in leg.findall('.//s:defs/s:path',N)+leg.findall('.//s:use',N):e.set('style','fill:#7f9ab0;stroke:none')
 # Record exact allocated region plus measured visible bbox; never use bbox to move anything.
 group_by_id={g.get('id'):g for g in root.findall('.//s:g',N)};orig_by_id={g.get('id'):g for g in original.findall('.//s:g',N)}
 for panel in ['b']+(['d'] if no==4 else []):
  r=next(r for r in prior if r['figure']==name and r['panel']==panel);members=[group_by_id[k] for k in r['members']];oldmembers=[orig_by_id[k] for k in r['members']]
  collections=[c for j,c in enumerate(cols) if (no!=4 or (panel=='b' and j==0) or (panel=='d' and j==1))]
  bbrows.append(dict(figure=name,panel=panel,allocated_mm=r['allocated_mm'],baseline_content_bbox_pt=bbox(original,oldmembers),output_content_bbox_pt=bbox(root,members),output_data_bbox_pt=bbox(root,[context(root,c) for c in collections]),panel_translation_pt=[0,0],auto_centering=False))
 out=O/'editable_svg'/f'{name}.svg';tree.write(out,encoding='utf-8',xml_declaration=True);pdf=O/'extended_data'/f'{name}.pdf';cairosvg.svg2pdf(url=str(out),write_to=str(pdf))
 ims=[]
 for path in [B/'extended_data'/f'{name}.pdf',pdf]:
  doc=fitz.open(path);pm=doc[0].get_pixmap(matrix=fitz.Matrix(1.5,1.5));ims.append(Image.frombytes('RGB',[pm.width,pm.height],pm.samples))
 ims[1].save(O/'review'/f'marker_{name}.png');before,after=ims;canvas=Image.new('RGB',(before.width+after.width+20,max(before.height,after.height)+25),'white');canvas.paste(before,(0,25));canvas.paste(after,(before.width+20,25));d=ImageDraw.Draw(canvas);d.text((5,5),'BASELINE',fill='black');d.text((before.width+25,5),'ORDINARY DOTS',fill='black');canvas.save(O/'review'/f'marker_{name}_comparison.png');print(name,flush=True)
(O/'audit/marker_numeric_checks.json').write_text(json.dumps(checks,indent=2));pd.DataFrame(maps).to_csv(O/'audit/marker_source_value_mapping.tsv',sep='\t',index=False);(O/'audit/marker_bboxes.json').write_text(json.dumps(bbrows,indent=2))
