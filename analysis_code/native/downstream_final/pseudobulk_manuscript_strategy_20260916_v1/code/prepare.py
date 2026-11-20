import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,hashlib,gzip,shutil,csv
import numpy as np,pandas as pd
from scipy import sparse,io
from threadpoolctl import threadpool_limits
B=Path(__file__).resolve().parents[1];W=B.parent
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908')
OLD=Path('/home/anzhen_vhd/scRNAseq/codex/pseudobulk_target_discovery_v1')
PB=R/'pseudobulk/descriptive_final_atlas_v1/outputs'
S=W/'disease_background_3groups_review_v2/sample_background_sidecar_v2.tsv'
MAP=R/'production_launches/magma_recovery_roundtrip_fix_20260914T0930/outputs/runs/20260914T013548_608658_d62bd9f120de/magma/primary_r0.1_L2_cell_weighted_model_mapping.tsv'
MASTER=Path('/home/anzhen_vhd/更新版MAGMA_scDRS/final_atlas_MAGMA_scDRS_v1/GWAS_gene_priority_master.tsv')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def dump(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n')
assert not (B/'PREPARATION.json').exists()
for d in ['inputs','tables','logs','results','tests']: (B/d).mkdir(exist_ok=True)
bounds={str(S):sha(S),str(MAP):sha(MAP),str(MASTER):sha(MASTER)}
old=pd.read_csv(OLD/'sample_disease_map.tsv',sep='\t',keep_default_na=False)
shutil.copy2(OLD/'sample_disease_map.tsv',B/'evidence/original_sample_disease_map.tsv');shutil.copy2(MASTER,B/'evidence/original_GWAS_master.tsv')
oldrows=[]
for _,r in old.iterrows():
 for sid in r.technical_sample_uids.split(';'):
  oldrows.append(dict(sample_uid=sid,original_context=r.disease_context,original_include=str(r.include_primary).lower()=='true',original_condition=r.condition_for_pseudobulk,original_donor=r.donor_uid_v2,original_tissue=r.tissue_v2))
o=pd.DataFrame(oldrows);assert o.sample_uid.is_unique
m=pd.read_csv(S,sep='\t',keep_default_na=False);m=m.merge(o,on='sample_uid',how='left',validate='one_to_one')
assert len(m)==247 and m.n_cells.sum()==1175115
m['selected']=m.original_include.fillna(False).astype(bool)
sel=m[m.selected].copy();assert len(sel)==53 and sel.study.nunique()==6
assert (sel.original_condition==sel.group_binary).all()
assert (sel.original_tissue==sel.tissue).all()
assert sel.donor_uid.ne('').all() and sel.donor_resolved.astype(str).str.lower().eq('true').all()
m.to_csv(B/'tables/sample_scope_audit.tsv',sep='\t',index=False)
m[m.study.str.contains('HCA')].to_csv(B/'tables/HCA_original_scope.tsv',sep='\t',index=False)
sel[sel.original_donor.ne(sel.donor_uid)].to_csv(B/'tables/corrected_donor_differences.tsv',sep='\t',index=False)
sel.groupby(['study','original_context','group_binary']).agg(samples=('sample_uid','size'),donors=('donor_uid','nunique'),cells=('n_cells','sum')).reset_index().to_csv(B/'tables/included_study_summary.tsv',sep='\t',index=False)
g=pd.read_csv(PB/'genes.tsv',sep='\t',keep_default_na=False);valid=g.stable_ensembl_id.str.fullmatch(r'ENSG\d{11}')
genes=pd.DataFrame({'gene_id':g.loc[valid,'stable_ensembl_id'].values,'gene_symbol':g.loc[valid,'gene_symbol_std'].values});assert len(genes)==36601 and genes.gene_id.is_unique;genes.to_csv(B/'genes.tsv',sep='\t',index=False)
elig=set(pd.read_csv(MAP,sep='\t')['group']);assert len(elig)==103
manifest=json.loads((PB/'FILES.json').read_text());roster=[];lineage=[];conservation=[];gate=[]
for level in ['L1','L2']:
 for fn in ['genes.tsv',f'sample_{level}/units.tsv.gz',f'sample_{level}/counts_units_by_genes.npz']:
  h=sha(PB/fn);assert h==manifest[fn]['sha256'];bounds[str(PB/fn)]=h
 u=pd.read_csv(PB/f'sample_{level}/units.tsv.gz',sep='\t',keep_default_na=False);x=sparse.load_npz(PB/f'sample_{level}/counts_units_by_genes.npz').tocsr();assert x.shape==(len(u),len(g)) and u.n_cells.sum()==1175115
 assert (x.data>=0).all() and np.equal(x.data,np.floor(x.data)).all()
 u['source_row']=np.arange(len(u));z=u.merge(sel[['sample_uid','study','donor_uid','tissue','group_binary','age','sex','original_context']],on='sample_uid',how='inner',validate='many_to_one',suffixes=('','_curated'))
 # Existing PB units carry only frozen annotations; explicitly take sidecar metadata.
 for k in ['study','donor_uid','tissue','group_binary','age','sex']:
  if k+'_curated' in z:z[k]=z[k+'_curated']
 z['cell_state']=z.level1_annotation if level=='L1' else z.level1_annotation+' | '+z.level2_annotation
 z=z[~z.level1_annotation.isin(['Schwann cell','Epithelial cell','Unknown'])]
 if level=='L2':z=z[z.cell_state.isin(elig)]
 cols=['original_context','cell_state','study','donor_uid','tissue','group_binary'];units=[];rr=[];cc=[]
 for i,(keys,v) in enumerate(z.groupby(cols,sort=True)):
  d=dict(zip(cols,keys));d.update(pseudobulk_id=level+'_'+hashlib.sha256(json.dumps(keys).encode()).hexdigest()[:16],n_cells=int(v.n_cells.sum()),sample_uids=';'.join(sorted(v.sample_uid.unique())))
  for k in ['age','sex']:
   a=[str(q) for q in v[k].unique() if str(q) not in ['','NA','NR','nan']];assert len(a)<=1;d[k]=a[0] if a else ''
  units.append(d);rr.extend([i]*len(v));cc.extend(v.source_row)
  for _,q in v.iterrows():lineage.append(dict(level=level,aggregate=d['pseudobulk_id'],sample_uid=q.sample_uid,source_pseudobulk_id=q.pseudobulk_id,n_cells=int(q.n_cells)))
 agg=sparse.csr_matrix((np.ones(len(rr),dtype='int64'),(rr,cc)),shape=(len(units),len(u)));y=(agg@x[:,valid.values]).tocsr();un=pd.DataFrame(units)
 assert int(y.sum())==int(x[cc,:][:,valid.values].sum());assert int(un.n_cells.sum())==int(z.n_cells.sum())
 conservation.append(dict(level=level,cells=int(un.n_cells.sum()),counts=int(y.sum()),units=len(un)))
 un.to_csv(B/f'tables/aggregated_{level}.tsv',sep='\t',index=False)
 for (ctx,state),v in un.groupby(['original_context','cell_state'],sort=True):
  v=v[v.n_cells>=20].copy();tab=pd.crosstab(v.study,v.group_binary);studies=[s for s in tab.index if tab.loc[s].get('healthy',0)>=2 and tab.loc[s].get('disease',0)>=2]
  for s,row in tab.iterrows():gate.append(dict(level=level,context=ctx,state=state,study=s,case=int(row.get('disease',0)),reference=int(row.get('healthy',0)),included=s in studies))
  v=v[v.study.isin(studies)];nc=sum(v.group_binary=='disease');nr=sum(v.group_binary=='healthy');name=level+'_'+ctx+'_'+hashlib.sha256(state.encode()).hexdigest()[:12]
  rec=dict(model_id=name,level=level,context=ctx,cell_state=state,n_case=int(nc),n_reference=int(nr),n_studies=v.study.nunique(),status='READY' if min(nc,nr)>=5 else 'INSUFFICIENT_UNITS');roster.append(rec)
  if rec['status']!='READY':continue
  p=B/'inputs'/name;p.mkdir();sub=y[v.index,:].T.tocsc()
  with gzip.open(p/'counts.mtx.gz','wb') as f:
   with threadpool_limits(limits=1):io.mmwrite(f,sub,field='integer')
  v=v.rename(columns={'study':'study_canonical','donor_uid':'donor_uid_v2','tissue':'tissue_v2','age':'age_v2','sex':'sex_v2','original_context':'disease_context','group_binary':'condition'})
  v.to_csv(p/'units.tsv',sep='\t',index=False)
 print(level,'READY',sum(r['level']==level and r['status']=='READY' for r in roster),flush=True)
 del x,y,agg
pd.DataFrame(roster).to_csv(B/'tables/model_roster.tsv',sep='\t',index=False);pd.DataFrame(gate).to_csv(B/'tables/within_study_gate.tsv',sep='\t',index=False);pd.DataFrame(lineage).to_csv(B/'tables/aggregation_lineage.tsv',sep='\t',index=False)
for fn in ['all_disease_context_edgeR_results.tsv.gz','GWAS_prioritised_cross_disease_meta_analysis.tsv.gz']:
 p=OLD/'edgeR_meta_results_v4'/fn;bounds[str(p)]=sha(p);pieces=[]
 for a in pd.read_csv(p,sep='\t',chunksize=100000,keep_default_na=False):pieces.append(a[a.gene_symbol.eq('CKB')])
 pd.concat(pieces).to_csv(B/'evidence'/('old_CKB_'+fn.replace('.gz','')),sep='\t',index=False)
for p in list((B/'code').glob('*'))+list((B/'evidence').glob('*'))+list((B/'inputs').rglob('*'))+[B/'genes.tsv',B/'tables/model_roster.tsv']:
 if p.is_file():bounds[str(p)]=sha(p)
base=json.loads((W/'RECEIPT.json').read_text());atlas=R/'panvascular_CPU_theta2_L1L2_final_v1.h5ad';st=atlas.stat();assert [st.st_dev,st.st_ino,st.st_size,st.st_mtime_ns]==base['atlas_stat']
dump(B/'RECEIPT.json',dict(authorization='User requested original manuscript thresholds/grouping/strategy replication and CKB audit',atlas_sha256=base['atlas_sha256'],atlas_stat=base['atlas_stat'],bound_files=bounds,scope='Corrected current counts and frozen L1/L2; original six-study map, 20 cells, within-study2/2 and context5/5 units, exact original R functions, full-genome DE, original2199 priority-set REML/BH and leave-one-study; corrected donor metadata retained; original L2 eligibility mechanism uses current103 labels; no atlas changes',conservation=conservation))
dump(B/'PREPARATION.json',dict(status='PASS',models=sum(r['status']=='READY' for r in roster),technical_samples=len(sel),studies=6,conservation=conservation))
