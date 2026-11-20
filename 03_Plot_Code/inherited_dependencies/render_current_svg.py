#!/usr/bin/env python3
"""Export frozen editable SVG masters without re-running scientific analyses."""
from pathlib import Path
import argparse, json, hashlib
import cairosvg, cairosvg.image
from cairosvg.surface import cairo
cairosvg.image.IMAGE_RENDERING.update({"pixelated":cairo.FILTER_NEAREST,"crisp-edges":cairo.FILTER_NEAREST})
import fitz
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);args=p.parse_args()
root=Path(__file__).resolve().parent;args.out.mkdir(parents=True,exist_ok=True)
order=[f'Figure{i}' for i in range(1,6)]+[f'ExtendedData_Figure{i}' for i in range(1,11)]
report=[];merged=fitz.open()
for name in order:
 svg=root/'frozen_svg'/f'{name}.svg';pdf=args.out/f'{name}.pdf'
 cairosvg.svg2pdf(url=str(svg),write_to=str(pdf))
 with fitz.open(pdf) as d:
  assert len(d)==1
  w,h=d[0].rect.width*25.4/72,d[0].rect.height*25.4/72
  merged.insert_pdf(d)
  typ,val=d.xref_get_key(d[0].xref,'Group')
  if typ=='dict':
   assert ' R' not in val;merged.xref_set_key(merged[-1].xref,'Group',val)
  report.append({'name':name,'width_mm':w,'height_mm':h,'svg_sha256':hashlib.sha256(svg.read_bytes()).hexdigest(),'pdf_sha256':hashlib.sha256(pdf.read_bytes()).hexdigest()})
merged.save(args.out/'All_15_figures.pdf',garbage=4,deflate=True)
(args.out/'EXPORT_VALIDATION.json').write_text(json.dumps(report,indent=2))
print('Exported 15 one-page figures at native dimensions. No scientific calculation performed.')
