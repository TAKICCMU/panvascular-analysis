from common import *
from matplotlib.patches import Rectangle
import shutil
from PIL import Image
BASE=Path('/Users/taki/Desktop/Panvascular_Marker_Communication_Enrichment_20260920')
m=pd.read_csv(O/'source_data/MAGMA_display9.tsv',sep='\t');p=pd.read_csv(O/'source_data/pseudobulk_display9_all7_lineages.tsv',sep='\t');rank=pd.read_csv(O/'source_data/pseudobulk_display9_selection.tsv',sep='\t')
labels={
'GOBP_STEROL_HOMEOSTASIS':'Sterol homeostasis','GOBP_LIPID_HOMEOSTASIS':'Lipid homeostasis','GOBP_STEROL_TRANSPORT':'Sterol transport',
'KEGG_PYRIMIDINE_METABOLISM':'Pyrimidine metabolism','KEGG_SELENOAMINO_ACID_METABOLISM':'Selenoamino acid metabolism','KEGG_SMALL_CELL_LUNG_CANCER':'Small cell lung cancer gene set',
'REACTOME_PLASMA_LIPOPROTEIN_ASSEMBLY_REMODELING_AND_CLEARANCE':'Lipoprotein assembly / remodelling / clearance','REACTOME_PLASMA_LIPOPROTEIN_CLEARANCE':'Plasma lipoprotein clearance','REACTOME_LDL_CLEARANCE':'LDL clearance',
'GOBP_NEGATIVE_REGULATION_OF_SYNCYTIUM_FORMATION_BY_PLASMA_MEMBRANE_FUSION':'Syncytium formation (negative regulation)',
'GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS':'Ribonucleoprotein complex biogenesis','GOBP_RIBOSOME_BIOGENESIS':'Ribosome biogenesis','KEGG_RIBOSOME':'Ribosome','KEGG_SPLICEOSOME':'Spliceosome','KEGG_OXIDATIVE_PHOSPHORYLATION':'Oxidative phosphorylation',
'REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE':'SRP-dependent targeting to membrane',
'REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION':'Eukaryotic translation elongation','REACTOME_EUKARYOTIC_TRANSLATION_INITIATION':'Eukaryotic translation initiation'}
assert set(m.pathway)|set(rank.pathway)<=set(labels)
labelmap=pd.DataFrame([{'pathway_id':k,'display_label':v} for k,v in labels.items()]);labelmap.to_csv(O/'source_data/pathway_display_labels.tsv',sep='\t',index=False)
# Construct in physical points; labels and legends are included in the cropped export.
# No whole-panel stretching: crop removes only external white/transparent canvas.
records=[]
def canvas(w,h):
 fig=plt.figure(figsize=(w/25.4,h/25.4));ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,w);ax.set_ylim(h,0);ax.axis('off');return fig,ax

def crop(fig,name):
 fig.canvas.draw();bounds=fig.get_tightbbox(fig.canvas.get_renderer());svg=O/'companion'/f'{name}.svg';fig.savefig(svg,bbox_inches='tight',pad_inches=.008,transparent=True);fig.savefig(O/'companion'/f'{name}.pdf',bbox_inches='tight',pad_inches=.008,transparent=True);plt.close(fig)
 r=ET.parse(svg).getroot();vb=list(map(float,r.get('viewBox').split()));oldvb=vb.copy()
 rgba=Image.open(io.BytesIO(cairosvg.svg2png(url=str(svg),scale=4))).convert('RGBA');ink=rgba.getchannel('A').getbbox();px=oldvb[2]/rgba.width;py=oldvb[3]/rgba.height;pad=.65
 cropbox=[max(0,ink[0]*px-pad),max(0,ink[1]*py-pad),min(oldvb[2],ink[2]*px+pad),min(oldvb[3],ink[3]*py+pad)]
 vb=[cropbox[0],cropbox[1],cropbox[2]-cropbox[0],cropbox[3]-cropbox[1]]
 r.set('viewBox',' '.join(map(str,vb)));r.set('width',str(vb[2])+'pt');r.set('height',str(vb[3])+'pt');ET.ElementTree(r).write(svg,encoding='utf-8',xml_declaration=True);cairosvg.svg2pdf(url=str(svg),write_to=str(O/'companion'/f'{name}.pdf'))
 records.append({'name':name,'cropped_width_mm':vb[2]/mm,'cropped_height_mm':vb[3]/mm,'artist_content_bbox_inches':list(bounds.bounds),'scale':1.0,'actual_ink_bbox_px':ink,'original_viewbox':oldvb,'cropped_viewbox':vb});return r,vb
# MAGMA: 3 pathways per database, FDR selection rule fixed in extract.py.
f,a=canvas(88,29);a.text(1,1,'MAGMA · genetic pathways',fontsize=6.5,fontweight='bold',va='top');a.text(68,1,'NES',ha='center',va='top',fontsize=6);a.text(84,1,'FDR',ha='center',va='top',fontsize=6)
ys=np.linspace(6,25.3,9);cx0,cx1=63,78.5
for i,(_,row) in enumerate(m.iterrows()):
 y=ys[i];a.text(61.5,y,labels[row.pathway],ha='right',va='center',fontsize=6.0)
 a.add_patch(Rectangle((cx0,y-.55),(cx1-cx0)*row.NES/2,.0+1.1,facecolor='#356b8c',edgecolor='none'))
 a.text(84,y,f'{row.padj:.3f}',ha='center',va='center',fontsize=5.8)
for x,val in zip([cx0,(cx0+cx1)/2,cx1],[0,1,2]):a.plot([x,x],[26.2,26.65],c='#555555',lw=.4);a.text(x,27.1,str(val),ha='center',va='top',fontsize=5.7)
a.plot([cx0,cx1],[26.2,26.2],c='#555555',lw=.4)
for start,label in [(0,'GO BP'),(3,'KEGG'),(6,'Reactome')]:
 a.plot([10,10],[ys[start]-.8,ys[start+2]+.8],c='#aaaaaa',lw=.45);a.text(8.8,np.mean(ys[start:start+3]),label,ha='right',va='center',fontsize=5.4)
left,lv=crop(f,'Figure3c_MAGMA_cropped')
# Pseudobulk matrix: original signed NES for every one of seven primary L1s.
states=[s for s in P.L1 if s in set(p.cell_state)];short={'Endothelial cell':'EC','Fibroblast':'Fib','Smooth muscle cell':'SMC','Macrophage':'Mac','Monocyte':'Mo','T cell':'T','Natural killer cell':'NK'}
mat=p.pivot(index='pathway',columns='cell_state',values='NES').reindex(index=rank.pathway,columns=states);qmat=p.pivot(index='pathway',columns='cell_state',values='padj').reindex(index=rank.pathway,columns=states);vmax=float(np.nanmax(np.abs(mat.to_numpy())))
f,a=canvas(88,29);a.text(1,1,'Pseudobulk GSEA',fontsize=6.5,fontweight='bold',va='top')
xx=np.linspace(57,85,7);yy=np.linspace(7,26.3,9);norm=colors.Normalize(-vmax,vmax);cmap=plt.get_cmap('RdBu_r')
for j,st in enumerate(states):a.text(xx[j],5.1,short[st],ha='center',va='bottom',fontsize=5.8)
for i,path in enumerate(rank.pathway):
 a.text(53.5,yy[i],labels[path],ha='right',va='center',fontsize=6)
 for j,st in enumerate(states):
  val=mat.iloc[i,j];qval=qmat.iloc[i,j];a.add_patch(Rectangle((xx[j]-1.9,yy[i]-.86),3.8,1.72,facecolor=cmap(norm(val)) if np.isfinite(val) else '#eeeeee',edgecolor='none'))
  if not np.isfinite(val):a.text(xx[j],yy[i],'×',ha='center',va='center',fontsize=5,color='#555555')
  elif np.isfinite(qval) and qval<.05:a.plot(xx[j],yy[i],marker='o',ms=1.5,color='white' if abs(val)>vmax*.6 else '#333333',markeredgewidth=0)
for start,label in [(0,'GO BP'),(3,'KEGG'),(6,'Reactome')]:
 a.plot([10,10],[yy[start]-.8,yy[start+2]+.8],c='#aaaaaa',lw=.45);a.text(8.8,np.mean(yy[start:start+3]),label,ha='right',va='center',fontsize=5.4)
# Compact key in the header, separated from pathway names and data rows.
for j,val in enumerate(np.linspace(-vmax,vmax,60)):a.add_patch(Rectangle((24+j*.13,1.1),.14,.8,facecolor=cmap(norm(val)),ec='none'))
a.text(33,1.5,f'NES: −{vmax:.2f} to +{vmax:.2f}',fontsize=5.2,va='center');a.plot(65,1.5,'o',ms=1.5,c='#333333');a.text(66.3,1.5,'FDR < 0.05',fontsize=5.2,va='center')
right,rv=crop(f,'Figure3c_pseudobulk_cropped')
root=ET.parse(BASE/'editable_svg/Figure3.svg').getroot()
# Remove the previous c plot only, including its obsolete database key.
remove_ids(root,{'enrichc_figure_1','text_104','text_105','patch_35','patch_36'})
parents={c:g for g in root.iter() for c in g}
for t in list(root.iter(q('text'))):
 if (t.text in ['■','▪','●']) and 200<float(t.get('y','-1'))<230:
  # Only the original c database symbols have these literal positions and group transforms.
  g=parents[t]
  if g.get('id') in ['text_102','text_103']:parents[g].remove(g)
# Protected b, d/e and f retain their full transformations and geometry.

def add_native(block,name,x,y):
 ids={e.get('id'):name+'_'+e.get('id') for e in block.iter() if e.get('id')}
 for e in block.iter():
  for k,v in list(e.attrib.items()):
   if k=='id':e.set(k,ids[v]);continue
   for old,new in ids.items():
    v=v.replace('url(#'+old+')','url(#'+new+')')
    if v=='#'+old:v='#'+new
   e.set(k,v)
 vb=list(map(float,block.get('viewBox').split()));g=ET.SubElement(root,q('g'),{'id':name,'transform':f'translate({x*mm-vb[0]} {y*mm-vb[1]})'})
 for child in block:
  if child.tag in [q('g'),q('defs')]:g.append(child)
add_native(left,'c_magma',3,50.4);add_native(right,'c_pseudobulk',94,50.4)
# c label remains at its original outer anchor; no repositioning of d/e shared rows.
export(root,'Figure3');(O/'notes/cropped_block_geometry.json').write_text(json.dumps(records,indent=2))
# Keep original GSEA/ORA and primary L2 displays as companions, without repainting them.
for name in ['Genetic_enrichment_seven_pathways.pdf','Genetic_enrichment_seven_pathways.svg','Pseudobulk_primary_enrichment.pdf']:
 shutil.copy2(BASE/'companion'/name,O/'companion'/name)
print(records)
