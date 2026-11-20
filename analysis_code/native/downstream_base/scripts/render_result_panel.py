"""Render a reviewed result table, preserving source precision and missingness."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('--spec',required=True);p.add_argument('--out',required=True);a=p.parse_args()
s=json.loads(Path(a.spec).read_text());out=Path(a.out)
if out.exists():raise SystemExit('New output directory required')
source=Path(s['input_table']);digest=hashlib.sha256(source.read_bytes()).hexdigest()
if digest!=s['input_sha256']:raise SystemExit('Source table hash differs from reviewed panel spec')
d=pd.read_csv(source,sep='\t');out.mkdir(parents=True);fig,ax=plt.subplots(figsize=s.get('figsize',[8,6]))
kind=s['kind']
if kind=='forest':
 x=pd.to_numeric(d[s['effect']],errors='coerce');se=pd.to_numeric(d[s['se']],errors='coerce');keep=x.notna()&se.notna();y=np.arange(keep.sum());ax.errorbar(x[keep],y,xerr=1.96*se[keep],fmt='o',capsize=2);ax.set_yticks(y,d.loc[keep,s['label']]);ax.axvline(0,c='grey',ls='--')
elif kind=='heatmap':
 v=d.pivot(index=s['row'],columns=s['column'],values=s['value']);im=ax.imshow(np.ma.masked_invalid(v.to_numpy()),aspect='auto',cmap=s.get('cmap','viridis'));ax.set_xticks(range(v.shape[1]),v.columns,rotation=90);ax.set_yticks(range(v.shape[0]),v.index);fig.colorbar(im,ax=ax)
elif kind=='scatter':ax.scatter(d[s['x']],d[s['y']],s=8,alpha=.7)
elif kind=='bar':ax.barh(d[s['label']],d[s['value']])
else:raise SystemExit('Unsupported reviewed panel type')
ax.set_title(s['title']);fig.tight_layout()
for ext in ['png','pdf','svg']:fig.savefig(out/(s['panel_id']+'.'+ext),dpi=300,bbox_inches='tight')
d.to_csv(out/'source_data.tsv',sep='\t',index=False)
s.update(rendering_only=True,computed_new_analysis=False,source_missing_count=d.isna().sum().to_dict());(out/'panel_manifest.json').write_text(json.dumps(s,indent=2,ensure_ascii=False))
