from pathlib import Path
import xml.etree.ElementTree as E,re,json,hashlib
import numpy as np,pandas as pd
O=Path('/Users/taki/Desktop/Panvascular_FigureYa_20260919');B=Path('/Users/taki/Desktop/投稿包_最新版_20260919');S=B/'Source_Data_Figures';A=Path('/Users/taki/Documents/Panvascular_Figure_Archive_20260919');N='{http://www.w3.org/2000/svg}';XL='{http://www.w3.org/1999/xlink}href'
def nums(s):return np.array([float(x) for x in re.findall(r'-?\d+(?:\.\d+)?',s)])
def ax(r,i):return next(e for e in r.iter() if e.get('id')=='axes_'+str(i))
checks=[]
# Read scientific source values; infer only the affine screen transform from retained coordinates.
for name,ai,rel,cols in [('Figure3',5,'Figure3/Fig3d_e_all29.tsv',('meta_logFC','meta_CI_low','meta_CI_high'))]:
 d=pd.read_csv(S/rel,sep='\t');r=E.parse(B/'图件编辑资料/editable_svg'/(name+'.svg')).getroot();a=ax(r,ai);defs={e.get('id'):e for e in r.iter() if e.tag==N+'path' and e.get('id')};points=[];ci=[]
 for u in a.iter():
  if u.tag==N+'use' and 'C' in defs[u.get(XL)[1:]].get('d',''):points.append(float(u.get('x')))
 for g in a:
  if g.get('id','').startswith('LineCollection_'):
   for p in g:
    if p.tag==N+'path':ci.append(nums(p.get('d')).reshape(-1,2)[:,0])
 coef=np.polyfit(d[cols[0]],points,1);err=np.max(np.abs(np.polyval(coef,d[list(cols[1:])])-np.asarray(ci)));assert err<.00001
 checks.append({'figure':name,'n_estimates':len(points),'CI_source_residual_pt':float(err),'original_CI_coordinates_preserved':True})
TOP=pd.read_csv(S/'Figure4/05_top_nodes_reference.tsv',sep='\t')
for name,axes,mods in [('Figure4',[6,7],['CM09','CM03']),('ExtendedData_Figure9',[1],['CM01','CM03','CM09'])]:
 r=E.parse(O/'editable_svg'/(name+'.svg')).getroot();ends=[]
 for ai in axes:
  a=ax(r,ai)
  for g in a:
   if any(e.tag==N+'circle' for e in g):
    p=next(e for e in g if e.tag==N+'path');v=nums(p.get('d')).reshape(-1,2);ends.append(v[1,0]-v[0,0])
 weights=np.concatenate([TOP[(TOP.module==m)&(TOP['rank']<=5)].sort_values('rank').weight for m in mods]);scale=(28 if name=='Figure4' else 87)*72/25.4/.3;err=np.max(np.abs(np.asarray(ends)-weights*scale));assert err<.00001
 checks.append({'figure':name,'n_module_weights':len(ends),'native_weight_endpoint_residual_pt':float(err),'module_order':mods})
# Verify every old CI segment and point use coordinate remains unchanged (only marker definitions differ).
for name,ais in [('Figure3',[5]),('ExtendedData_Figure6',[1,2]),('ExtendedData_Figure10',[5])]:
 before=E.parse(B/'图件编辑资料/editable_svg'/(name+'.svg')).getroot();after=E.parse(O/'editable_svg'/(name+'.svg')).getroot()
 for ai in ais:
  a=ax(before,ai);b=ax(after,ai)
  for g in a:
   if g.get('id','').startswith('LineCollection_'):
    q=next(e for e in b if e.get('id')==g.get('id'));assert E.tostring(g)==E.tostring(q)
  ua=[(u.get(XL),u.get('x'),u.get('y')) for u in a.iter() if u.tag==N+'use'];ub=[(u.get(XL),u.get('x'),u.get('y')) for u in b.iter() if u.tag==N+'use'];assert ua==ub
 checks.append({'figure':name,'all_CI_segments_and_marker_centres_identical':True})
(O/'review/evidence_numeric_identity.json').write_text(json.dumps(checks,indent=2))
print(json.dumps(checks,indent=2))
