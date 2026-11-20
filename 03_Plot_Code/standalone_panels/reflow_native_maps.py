"""Reflow frozen map objects inside their original panel. Images/points are not regenerated."""
from export_native_panels import *

def reflow(name,parts,ncols,widths=None,source_override=None):
 figure=name[:-1];letter=name[-1];source=P/'图件编辑资料/editable_svg'/f'{figure}.svg';original=E.parse(str(source_override or source)).getroot();pieces=[]
 for i,ids in enumerate(parts):
  sub=extract(original,ids);box=content_bbox(sub);pieces.append((sub,box))
 # Each column receives a physical view slot. Isotropic placement only.
 widths=widths or [155]*ncols;gap=15;pad=8
 if name=='Figure1b':
  common=min(widths[i]/(b[2]-b[0]) for i,(_,b) in enumerate(pieces));widths=[(b[2]-b[0])*common for _,b in pieces]
 heights=[]
 for i,(sub,b) in enumerate(pieces):heights.append(widths[i%ncols]*(b[3]-b[1])/(b[2]-b[0]))
 rowh=[max(heights[i:i+ncols]) for i in range(0,len(parts),ncols)];W=sum(widths)+gap*(ncols-1)+pad*2;H=sum(rowh)+gap*(len(rowh)-1)+pad*2
 root=E.Element(Q('svg'),nsmap={None:NS,'xlink':'http://www.w3.org/1999/xlink'},width=f'{W}pt',height=f'{H}pt',viewBox=f'0 0 {W} {H}');defs=E.SubElement(root,Q('defs'));seen=set()
 for i,(sub,b) in enumerate(pieces):
  for de in list(sub):
   if de.tag==Q('defs'):
    for ch in de:
     id=ch.get('id') or E.tostring(ch)
     if id not in seen:defs.append(ch);seen.add(id)
  col=i%ncols;row=i//ncols;scale=widths[col]/(b[2]-b[0]);x=pad+sum(widths[:col])+gap*col;y=pad+sum(rowh[:row])+gap*row
  # Normalise label size after the isotropic content placement, not the data coordinates.
  for tx in sub.iter(Q('text')):
   text=''.join(tx.itertext()).strip();size=8 if text in ['Macrophage','Monocyte','Neutrophil','Dendritic','Fibroblast','Fibromyocyte','Pericyte','Smooth muscle','T','Natural killer','Tissue'] else 7
   tx.set('style',tx.get('style','')+f';font-size:{size/scale}px')
  group=E.SubElement(root,Q('g'),transform=f'translate({x} {y}) scale({scale}) translate({-b[0]} {-b[1]})')
  for ch in list(sub):
   if ch.tag!=Q('defs'):group.append(ch)
 rec=native_save(root,figure,letter,[v for p in parts for v in p]);rec['mode']='reflowed_frozen_map_objects';rec['change']='existing map objects given compact individual slots; isotropic placement; embedded images and number anchors unchanged';return rec

if __name__=='__main__':
 records=[]
 records.append(reflow('Figure1b',[['axes_2','axes_4'],['axes_3','text_148']],2,[255,195]))
 records.append(reflow('Figure1e',[['f1e_axes_7'],['f1e_axes_8','f1e_text_186'],['f1e_axes_9','f1e_text_187'],['f1e_axes_10','f1e_text_188']],4,[75,135,135,135]))
 records.append(reflow('ExtendedData_Figure3a',[[f'axes_{i}',f'text_{145+i}'] for i in range(1,5)],2,[155,155]))
 records.append(reflow('ExtendedData_Figure4c',[['axes_8','text_127'],['axes_10','text_128']],2,[170,170]))
 records.append(reflow('ExtendedData_Figure5a',[['axes_5','text_140'],['axes_1','text_137'],['axes_2','text_138'],['axes_3','text_139']],2,[195,145]))
 # Restore simple object name in EC map; remove unrelated rotated marker-group label.
 original=E.parse(str(P/'图件编辑资料/editable_svg/ExtendedData_Figure4.svg')).getroot();r=extract(original,['axes_2']);x0,y0,x1,y1=content_bbox(r);t=E.SubElement(r,Q('text'),x=str((x0+x1)/2),y=str(y0-7),style='font:7px Arial;fill:#20252A;text-anchor:middle');t.text='Endothelial';rec=native_save(r,'ExtendedData_Figure4','a',['axes_2']);rec['change']='Endothelial object label placed horizontally above frozen9-state map';records.append(rec)
 (O/'review/map_reflow_manifest.json').write_text(json.dumps(records,ensure_ascii=False,indent=2))
