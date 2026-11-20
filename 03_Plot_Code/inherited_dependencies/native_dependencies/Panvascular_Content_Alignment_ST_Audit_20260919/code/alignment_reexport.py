"""Re-export existing SVG verbatim: no layout edits, nearest-neighbour matrix rendering."""
from pathlib import Path
import os,sys,json,hashlib
os.environ.setdefault('DYLD_FALLBACK_LIBRARY_PATH','/Users/taki/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/lib')
sys.path.insert(0,'/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Local_Touchup_20260918/tools')
import cairosvg,cairosvg.image,fitz
from cairosvg.surface import cairo
from PIL import Image,ImageDraw
cairosvg.image.IMAGE_RENDERING.update({'pixelated':cairo.FILTER_NEAREST,'crisp-edges':cairo.FILTER_NEAREST})
B=Path('/Users/taki/Desktop/Panvascular_FigureYa_20260919');O=Path('/Users/taki/Desktop/Panvascular_Content_Alignment_ST_Audit_20260919');audit=[]
for name in ['Figure1','Figure2','Figure4']+[f'ExtendedData_Figure{i}' for i in range(1,11)]:
 svg=O/'editable_svg'/f'{name}.svg';before=hashlib.sha256(svg.read_bytes()).hexdigest();folder='main' if name.startswith('Figure') else 'extended_data';pdf=O/folder/f'{name}.pdf'
 cairosvg.svg2pdf(url=str(svg),write_to=str(pdf))
 d=fitz.open(pdf);pm=d[0].get_pixmap(matrix=fitz.Matrix(1.3,1.3));im=Image.frombytes('RGB',[pm.width,pm.height],pm.samples);im.save(O/'review'/f'alignment_{name}.png')
 d0=fitz.open(B/folder/f'{name}.pdf');pm=d0[0].get_pixmap(matrix=fitz.Matrix(1.3,1.3));im0=Image.frombytes('RGB',[pm.width,pm.height],pm.samples);out=Image.new('RGB',(im.width+im0.width+20,max(im.height,im0.height)+25),'white');out.paste(im0,(0,25));out.paste(im,(im0.width+20,25));dr=ImageDraw.Draw(out);dr.text((5,5),'BEFORE',fill='black');dr.text((im0.width+25,5),'AFTER',fill='black');out.save(O/'review'/f'alignment_{name}_comparison.png')
 assert before==hashlib.sha256(svg.read_bytes()).hexdigest()
 audit.append(dict(figure=name,svg_sha256=before,layout='VERBATIM_SVG_UNCHANGED',pixelated='FILTER_NEAREST',crisp_edges='FILTER_NEAREST',pdf_images=[{'width':x[2],'height':x[3]} for x in d[0].get_images()]))
 print(name,flush=True)
(O/'audit/alignment_reexport_nearest.json').write_text(json.dumps(audit,indent=2))
