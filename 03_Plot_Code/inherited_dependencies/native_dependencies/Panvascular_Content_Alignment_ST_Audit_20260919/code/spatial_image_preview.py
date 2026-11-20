from pathlib import Path
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
import pandas as pd,json
O=Path(__file__).resolve().parents[1];A=Path('/Users/taki/Desktop/Codex任务1/spatial_visium_geo_local/extracted');S=Path('/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Nature_Reconstruction_20260918/SourceData/spatial')
for gene,gse,sec,bounds in [('LIPA','GSE248608','GSM7917895_CA2',(300,1100,200,920)),('STAT3','GSE314851','GSM9414344_LCA',(450,1300,850,1550))]:
 img=Image.open(A/gse/(sec+'_tissue_hires_image.png'));sf=json.load(open(A/gse/(sec+'_scalefactors_json.json')))['tissue_hires_scalef'];df=pd.read_csv(S/f'Figure5a_{gene}_spots.tsv',sep='\t');x=df.x*sf;y=df.y*sf
 f,axs=plt.subplots(1,2,figsize=(14,7))
 for ax in axs:
  ax.imshow(img,interpolation='nearest');ax.set_xlim(bounds[:2]);ax.set_ylim(bounds[3],bounds[2]);ax.set_xlabel('Hires pixel column (x)');ax.set_ylabel('Hires pixel row (y)');ax.grid(alpha=.2)
 axs[1].scatter(x,y,s=30,facecolors='none',edgecolors='cyan',lw=.5)
 f.tight_layout();f.savefig(O/'review'/f'{gene}_raw_coordinates.png',dpi=160)
