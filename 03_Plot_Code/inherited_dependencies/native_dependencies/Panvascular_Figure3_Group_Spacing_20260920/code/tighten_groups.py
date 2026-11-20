from common import *
from PIL import Image, ImageChops
BASE=Path('/Users/taki/Desktop/Panvascular_Figure3_Genetic_Transcriptomic_20260920')
root=ET.parse(BASE/'editable_svg/Figure3.svg').getroot()
original=copy.deepcopy(root)
doc=fitz.open(BASE/'figures/Figure3.pdf'); page=doc[0]
spans=[s for b in page.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans']]
records=[]; changed=set()
for bid in ['c_magma','c_pseudobulk']:
 block=next(e for e in root.iter() if e.get('id')==bid)
 tx,ty=map(float,re.findall(r'[-\d.]+',block.get('transform')))
 texts=list(block.iter(q('text')))
 # Exactly nine unchanged 6pt pathway row labels, selected by their original right anchor.
 anchor=174.906709 if bid=='c_magma' else 152.229543
 rows=[t for t in texts if "font: 6px" in t.get('style','') and 'text-anchor: end' in t.get('style','')]
 assert len(rows)==9,[(t.text,t.attrib) for t in rows]
 lines=[p for p in block.iter(q('path')) if 'stroke: #aaaaaa' in p.get('style','')]
 assert len(lines)==3
 for i,label in enumerate(['GO BP','KEGG','Reactome']):
  names=[t.text for t in rows[i*3:i*3+3]]
  bb=[]
  for name in names:
   matches=[s for s in spans if s['text']==name and 145<s['bbox'][1]<222]
   assert len(matches)==1,(name,matches)
   bb.append(matches[0]['bbox'])
  left=min(b[0] for b in bb)
  line=lines[i]; nums=list(map(float,re.findall(r'[-\d.]+',line.get('d')))); assert len(nums)==4
  x=left-tx-1.2*mm
  before=(left-(nums[0]+tx))/mm
  line.set('d',f'M {x:.9f} {nums[1]} L {x:.9f} {nums[3]}')
  lab=next(t for t in texts if t.text==label)
  lab.set('x',str(x-1.0*mm)); lab.attrib.pop('transform',None)
  records.append({'block':bid,'group':label,'pathways':names,'old_gap_mm':before,'new_gap_mm':1.2,'label_to_line_gap_mm':1.0})
  changed.update([id(line),id(lab)])
# Prove only twelve SVG elements changed; all other objects/coordinates/styles preserved.
a=list(original.iter());b=list(root.iter());assert len(a)==len(b)
changed_elements=[]
for old,new in zip(a,b):
 assert old.tag==new.tag and old.text==new.text and old.tail==new.tail
 if old.attrib!=new.attrib:
  assert id(new) in changed
  changed_elements.append(new.tag)
assert len(changed_elements)==12
export(root,'Figure3')
pd.DataFrame(records).to_csv(O/'group_spacing.tsv',sep='\t',index=False)
# Crop previews only; final mother PDF retains its native 183 x 170 mm page.
clip=fitz.Rect(2*mm,49.5*mm,181.8*mm,79.5*mm)
for name,fp in [('before',BASE/'figures/Figure3.pdf'),('after',O/'figures/Figure3.pdf')]:
 d=fitz.open(fp); d[0].get_pixmap(matrix=fitz.Matrix(3,3),clip=clip).save(O/'review'/f'Figure3c_{name}.png')
im1=Image.open(O/'review/Figure3c_before.png');im2=Image.open(O/'review/Figure3c_after.png')
out=Image.new('RGB',(im1.width,im1.height*2+20),'white');out.paste(im1,(0,0));out.paste(im2,(0,im1.height+20));out.save(O/'review/Figure3c_before_after.png')
versions=pd.read_csv(BASE/'FIGURE_VERSIONS.tsv',sep='\t'); row=versions.figure=='Figure3'
versions.loc[row,'source_pdf']=str(O/'figures/Figure3.pdf');versions.loc[row,'source_svg']=str(O/'editable_svg/Figure3.svg');versions.loc[row,'version']=O.name;versions.loc[row,'status']='GROUP_SPACING_REVIEW_PENDING'
versions.loc[row,'sha256']=sha(O/'figures/Figure3.pdf');versions.loc[row,'svg_sha256']=sha(O/'editable_svg/Figure3.svg')
merged=fitz.open()
for _,r in versions.iterrows():
 d=fitz.open(r.source_pdf);assert len(d)==1
 merged.insert_pdf(d)
 typ,val=d.xref_get_key(d[0].xref,'Group')
 if typ=='dict':
  assert not re.search(r'\d+ \d+ R',val)
  merged.xref_set_key(merged[-1].xref,'Group',val)
merged.save(O/'全部图件_15页_分组间距修正版.pdf',garbage=4,deflate=True);merged.close()
merged=fitz.open(O/'全部图件_15页_分组间距修正版.pdf')
for i,(_,r) in enumerate(versions.iterrows()):
 d=fitz.open(r.source_pdf)
 assert merged[i].rect==d[0].rect
 assert merged[i].get_pixmap().samples==d[0].get_pixmap().samples
versions.to_csv(O/'FIGURE_VERSIONS.tsv',sep='\t',index=False)
(O/'修改说明.md').write_text('仅调整 Figure3c 的六处分组竖线与 GO BP / KEGG / Reactome 标签。竖线距各组三条通路名称的最左边界为 1.2 mm，标签距竖线为 1.0 mm。通路名称、数据列、字号及其他图形对象不动。\n\n脚本核验：仅12个SVG元素属性改变；全部科学数值与其他元素保持；15页合并稿与逐页母版尺寸及渲染一致。来源数据和统计解释继续沿用 '+str(BASE)+'。ORA待办不变。\n\nFigure3母版维持183×170 mm；其余14幅沿用登记版本，未重绘。\n')
print(pd.DataFrame(records)[['block','group','old_gap_mm','new_gap_mm']].to_string(index=False))
print('PASS: 12 SVG objects only; 15 merged pages identical to native sources')
