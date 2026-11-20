from pathlib import Path
import json,hashlib,zipfile
import pandas as pd
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908/production_launches')
RUN=R/'phase1_/20260913T183720_920345Z/outputs/runs/20260913T183912_636388_077099775370';OUT=R/'execution_control_20260914/review_packages';SUMMARY=OUT/'marker_review';SUMMARY.mkdir(exist_ok=True)
status=json.loads((RUN/'RUN_STATUS.json').read_text());assert all(x['status']=='COMPLETE' for x in status.values())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
records=[]
for p in sorted((RUN/'markers').glob('*.tsv.gz')):
 parts=[];rows=0
 for c in pd.read_csv(p,sep='\t',chunksize=100000):
  rows+=len(c);assert c.detection_fraction.between(0,1).all();assert c.mean_log1p.ge(0).all()
  parts.append(c.loc[c.descriptive_rank.le(30)])
 dest=SUMMARY/(p.name.replace('.tsv.gz','_display_top30.tsv.gz'));pd.concat(parts).to_csv(dest,sep='\t',index=False)
 records.append({'server_full_table':str(p),'sha256':sha(p),'bytes':p.stat().st_size,'rows':rows,'display_subset':dest.name,'display_rule':'descriptive_rank <=30 including ties; display only, no model selection or annotation changes'})
(SUMMARY/'FULL_MARKER_TABLE_MANIFEST.json').write_text(json.dumps(records,indent=2))
groups={'phase1_marker_review':list(SUMMARY.iterdir())+[RUN/'markers/MARKER_METHOD.json'], 'phase1_dotplot_review':list((RUN/'figures').glob('*DotPlot*'))+[RUN/'dotplot.receipt.json',RUN/'dotplot.log']}
entries=[]
for name,files in groups.items():
 z=OUT/(name+'.zip');assert not z.exists()
 with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as a:
  d={}
  for p in files:a.write(p,p.name);d[p.name]=sha(p)
  a.writestr('PACKAGE_SHA256.json',json.dumps(d,indent=2))
 with zipfile.ZipFile(z) as a:assert a.testzip() is None
 assert z.stat().st_size<35*1024**2
 entries.append({'path':str(z),'bytes':z.stat().st_size,'sha256':sha(z)})
existing=json.loads((OUT/'AVAILABLE_PACKAGES.json').read_text());(OUT/'AVAILABLE_PACKAGES.json').write_text(json.dumps(existing+entries,indent=2));print(json.dumps(entries),flush=True)
