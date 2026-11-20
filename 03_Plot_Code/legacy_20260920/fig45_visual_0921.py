from common_main import *
from PIL import Image
BACK=O/'merge_records/before_20260921_adjustments/editable_svg'
find=lambda r,i:next(e for e in r.iter() if e.get('id')==i)
def move(e,dx=0,dy=0):e.set('transform',f'translate({dx:.8f} {dy:.8f}) '+e.get('transform',''))
def bbox(r,e):
 p={c:a for a in r.iter() for c in a};n=ET.Element(r.tag,r.attrib)
 for d in r.iter(q('defs')):n.append(copy.deepcopy(d))
 chain=[];a=p.get(e)
 while a is not None and a is not r:
  if a.get('transform'):chain.append(a.get('transform'))
  a=p.get(a)
 g=ET.SubElement(n,q('g'),{'transform':' '.join(reversed(chain))});g.append(copy.deepcopy(e))
 b=cairosvg.svg2png(bytestring=ET.tostring(n),output_width=int(float(r.get('width').replace('pt',''))*4),output_height=int(float(r.get('height').replace('pt',''))*4));im=Image.open(io.BytesIO(b));bb=im.getchannel('A').getbbox();return np.array(bb)/4

def center(r,e,x):
 bb=bbox(r,e);dx=x-(bb[0]+bb[2])/2;move(e,dx);return {'id':e.get('id'),'dx_pt':dx,'center_x_pt':x}
checks=[]
r=ET.parse(BACK/'Figure4.svg').getroot()
# 4b follows the actual table data region, not full-page labels.
checks.append(center(r,find(r,'legend_1'),(107.716535+501.732283)/2-12.12992))
for ai,tid,cx in [('axes_6','text_110',130),('axes_7','text_111',389)]:
 a=find(r,ai);rec=center(r,a,cx);checks.append(rec)
 # Translate title to centre over the common 0–0.30 numerical bar region.
 orig=(167.244094+246.614173)/2-7.14862 if ai=='axes_6' else (425.19685+504.566929)/2-5.77854
 t=find(r,tid).find(q('text'));t.set('x',str(orig+rec['dx_pt']));t.set('y','174');t.set('style',t.get('style').replace('text-anchor: start','text-anchor: middle'));t.attrib.pop('transform',None)
# Validate module-to-colour mapping directly from existing 14 row paths, before adding a key.
a=find(r,'axes_10');paths=[p for p in a.iter(q('path')) if 'fill: #' in p.get('style','') and 'stroke: #' in p.get('style','') and '#222222' not in p.get('style','') and '#ffffff' not in p.get('style','')];assert len(paths)==14
cols=[re.search(r'fill: (#[0-9a-f]+)',p.get('style')).group(1) for p in paths]
assert cols[0]=='#657c92' and cols[2]=='#8d7a9f' and cols[8]=='#398f8b' and all(c=='#929fa9' for i,c in enumerate(cols) if i not in [0,2,8])
g=ET.SubElement(r,q('g'),{'id':'Figure4f_module_key_0921'})
for x,label,c in [(0,'CM01',cols[0]),(41,'CM03',cols[2]),(82,'CM09',cols[8]),(123,'Other modules','#929fa9')]:
 ET.SubElement(g,q('circle'),{'cx':str(x),'cy':'354','r':'1.8','fill':c});text(g,x+5,356,label,size=6)
checks.append(center(r,g,(70.866142+490.393701)/2-13.25492))
export(r,'Figure4')
# Final native PDF glyph centres include the title ancestor translations.
def title_centres():
 doc=fitz.open(O/'figures/Figure4.pdf');sp=[s for b in doc[0].get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans'] if 150<s['bbox'][1]<275]
 xs=lambda z:sorted([(s['bbox'][0]+s['bbox'][2])/2 for s in sp if s['text']==z])
 return {lab:{'data_center_pt':(xs('0.00')[j]+xs('0.30')[j])/2,'title_center_pt':xs(lab)[0]} for j,lab in enumerate(['CM09','CM03'])}
before=title_centres()
for tid,lab in [('text_110','CM09'),('text_111','CM03')]:
 t=find(r,tid).find(q('text'));delta=before[lab]['data_center_pt']-before[lab]['title_center_pt'];t.set('x',str(float(t.get('x'))+delta))
export(r,'Figure4')
after=title_centres()
for lab in after:
 after[lab]['error_pt']=after[lab]['title_center_pt']-after[lab]['data_center_pt'];assert abs(after[lab]['error_pt'])<.01
(O/'notes/Figure4_title_alignment_0921.json').write_text(json.dumps({'before':before,'after':after,'only_title_x_adjusted':True},indent=2))
# Figure5: stored artists only; no arrays, colours, values or masks recalculated.
r=ET.parse(BACK/'Figure5.svg').getroot();orig=copy.deepcopy(r)
remove_ids(r,{'r2_text_23','r2_text_27','r2_legend_4','r2_legend_5'})
for i in ['r2_axes_6','r2_axes_7','r2_axes_8','r2_axes_9','r2_axes_10','r2_text_28','r2_text_29','r2_text_30','r2_text_31','r2_text_32']:
 move(find(r,i),0,-4*mm)
move(find(r,'r2_legend_3'),0,-4*mm);checks.append(center(r,find(r,'r2_legend_3'),17*mm))
# 5b: identity keys centred in corresponding sides of the zero line.
checks.append(center(r,find(r,'r2_legend_1'),(348.661417+419.98476)/2));checks.append(center(r,find(r,'r2_legend_2'),(419.98476+459.212598)/2))
# Pooled-symbol key is centred over the numerical forest extent. Package its symbol and text together.
p={c:a for a in r.iter() for c in a};g=ET.SubElement(r,q('g'),{'id':'Figure5b_pooled_key_0921'})
for i in ['r2_axes_12','r2_text_39']:
 e=find(r,i);p[e].remove(e);g.append(e)
checks.append(center(r,g,(348.661417+459.212598)/2))
# 5c centre over the actual forest numerical data interval (preserve all 29 row centres).
forest=find(r,'FigureYa6_native_forest')
# Native svglite table column 85.87–221.93 is the unchanged numerical forest.
cx=(85.87+221.93)/2+2.81693+2.834645669291339
checks.append(center(r,find(r,'legend_2'),cx))
checks.append(center(r,find(r,'legend_3'),(425.19685+507.401575)/2-25.46014))
checks.append(center(r,find(r,'text_156'),(342.992126+507.401575)/2-13.08514))
checks.append(center(r,find(r,'legend_4'),(396.850394+484.724409)/2-20.77264))
# Assert all image payloads and original data paths are unchanged apart from ancestor translations.
images=lambda rr:[e.attrib for e in rr.iter(q('image'))]
assert images(r)==images(orig)
for i in ['FigureYa6_native_forest','axes_8','axes_9','axes_10']:
 assert ET.tostring(find(r,i))==ET.tostring(find(orig,i)),i
# Fixed outer panel letters must remain above image artists.
p={c:a for a in r.iter() for c in a}
for i in ['text_131','text_132','text_136','text_137','text_148','text_149']:
 e=find(r,i);p[e].remove(e);r.append(e)
export(r,'Figure5')
(O/'notes/Figure45_0921_layout.json').write_text(json.dumps({'changes':checks,'Figure5_images_identical':True,'Figure5_c_f_data_groups_identical':True,'module_colours_verified_from_row_paths':dict(zip(['CM01','CM03','CM09','Other modules'],[cols[0],cols[2],cols[8],'#929fa9'])),'models_recomputed':False,'ORA':'PENDING','donor_identity':'PENDING'},indent=2))
(O/'notes/Figure45_0921_legend_suggestions.md').write_text('Figure4f新增颜色图例仅解释原有模块颜色：CM01、CM03、CM09及其余模块，不增加证据等级。\n\nFigure5a图面删除切片身份长标题；图注和Source Data仍注明：上行Anterior tibial artery / CA1，1779个既有入图spots；下行Coronary artery / T1077，1363个既有入图spots。三视图均为H&E、各Fibroblast状态边际q05丰度之和、FN1 log1p(CP10k)。只平移第二组三视图，没有裁切组织、改变坐标/图像、mask、色标和原值。\n\nFigure5b仍为60个gene–unit estimates、34个recorded unit IDs，不称独立供者。FN1和SCARB1图例各自居中于零线左、右对应图区；其余符号图例居中于各自数据区域。\n\nORA及供者身份核查继续PENDING；本轮不作生物学结论变化。\n')
print('exported Figure4/Figure5')
