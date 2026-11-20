import subprocess,json
from pathlib import Path
import pandas as pd,numpy as np
O=Path('/Users/taki/Desktop/Panvascular_Spatial_Evidence_Assessment_20260920');script='''import h5py,json,numpy as np
for gse,stem in [('GSE248608','GSM7917895_CA2'),('GSE314851','GSM9414344_LCA')]:
 p='/home/anzhen_vhd/scRNAseq/codex/spatial_validation_v5_locked71_gpu_final/02_cell2location/'+gse+'__final_atlas_cell2location.h5ad'
 with h5py.File(p,'r') as f:
  def read(obj):
   if isinstance(obj,h5py.Group):
    c=obj['codes'][:];v=obj['categories'].asstr()[:];return np.array([v[i] if i>=0 else '' for i in c])
   return obj.asstr()[:] if obj.dtype.kind in ['O','S'] else obj[:]
  section=read(f['obs/section_id']);inds=np.where(section==gse+'::'+stem)[0];keys=read(f['obs/_index'])[inds];coords=f['obsm/spatial'][:][inds];print(json.dumps({'gse':gse,'section':gse+'::'+stem,'file':p,'keys':keys.tolist(),'coords':coords.tolist(),'obs_fields':list(f['obs'])}))
'''
out=subprocess.check_output(['ssh','-o','BatchMode=yes','codex-shpc-gpu','/home/anzhen_vhd/miniconda3/envs/cell2location/bin/python','-'],input=script.encode()).decode();checks=[]
for line in out.splitlines():
 r=json.loads(line);g='LIPA' if r['gse']=='GSE248608' else 'STAT3';p=pd.read_csv('/Users/taki/Desktop/Panvascular_Content_Alignment_ST_Audit_20260919/audit/spatial_counts_'+g+'_plotted_spots.tsv',sep='\t');h=pd.DataFrame({'historical_obs_id':r['keys'],'historical_x':np.array(r['coords'])[:,0],'historical_y':np.array(r['coords'])[:,1]});h['barcode']=h.historical_obs_id.str.rsplit('::',n=1).str[-1];h.to_csv(O/f'tables/historical_{g}_spot_ids_and_coords.tsv',sep='\t',index=False);j=p.merge(h,on='barcode',validate='one_to_one');checks.append({'gene':g,'section_id':r['section'],'historical_file':r['file'],'historical_n':len(h),'current_plot_n':len(p),'barcode_set_identical':set(h.barcode)==set(p.barcode),'matched_n':len(j),'coordinate_max_abs_diff':float(max((j.historical_x-j.plot_source_x).abs().max(),(j.historical_y-j.plot_source_y).abs().max())),'inference_boundary':'Same barcode/coordinates does not equate count-gene universe, model context, residual transform, or neighborhood test to current raw expression display'})
(O/'tables/historical_target_section_spot_comparison.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))
