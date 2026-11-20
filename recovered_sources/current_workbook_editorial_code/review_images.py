import fitz
from pathlib import Path
from PIL import Image,ImageOps,ImageDraw
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');O=P/'图件编辑资料/补充表与方法同步_20260921';R=O/'视觉核查';R.mkdir(exist_ok=True)
d=fitz.open(O/'64表逐字段抽样预览.pdf')
for start in range(0,64,4):
 ims=[]
 for p in list(d)[start:start+4]:
  pix=p.get_pixmap(matrix=fitz.Matrix(.82,.82));ims.append(Image.frombytes('RGB',[pix.width,pix.height],pix.samples))
 ww=max(i.width for i in ims);hh=max(i.height for i in ims);out=Image.new('RGB',(ww*2,hh*2),(204,204,204))
 for j,im in enumerate(ims):out.paste(im,((j%2)*ww,(j//2)*hh))
 out.save(R/f'Tables_{start+1:02}_{start+4:02}.png')
for name,terms in [('Manuscript',['Cell-level filters','susieR 0.16.6','Extended Data Fig. 9','Fig. 5 Conditional']),('Supplementary_Information',['Fine-mapping implementation','Single-cell datasets','Spatial and single-cell']),('Figure_Legends',['Fig. 5 Conditional','Extended Data Fig. 9'])]:
 doc=fitz.open(O/'PDF预览'/f'{name}.pdf');idx=[]
 for i,p in enumerate(doc):
  if any(t in p.get_text() for t in terms):idx.append(i)
 for i in sorted(set(idx)):
  doc[i].get_pixmap(matrix=fitz.Matrix(1.35,1.35)).save(R/f'{name}_p{i+1:02}.png')
 print(name,len(doc),'reviewpages',[i+1 for i in idx])
