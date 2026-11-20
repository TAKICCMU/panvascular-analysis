import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,gzip,sys
from scipy import sparse,io
from threadpoolctl import threadpool_limits
B=Path(__file__).resolve().parents[1];W=B.parent
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908')
PB=R/'pseudobulk/descriptive_final_atlas_v1/outputs';S=W/'disease_background_3groups_review_v2/sample_background_sidecar_v2.tsv'
MAP=R/'production_launches/magma_recovery_roundtrip_fix_20260914T0930/outputs/runs/20260914T013548_608658_d62bd9f120de/magma/primary_r0.1_L2_cell_weighted_model_mapping.tsv'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def dump(p,d):Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2,default=str))
if (B/'PREPARATION.json').exists():raise RuntimeError('Existing preparation: do not overwrite')
(B/'inputs').mkdir(exist_ok=True);(B/'tables').mkdir(exist_ok=True);(B/'logs').mkdir(exist_ok=True)
m=pd.read_csv(S,sep='\t',keep_default_na=False).set_index('sample_uid',drop=False)
assert len(m)==247 and m.index.is_unique and m.n_cells.sum()==1175115
assert sum(m.main_analysis_scope_candidate.astype(str).str.lower().eq('true'))==241
assert sha(S)=='852c08353c53a1b4d32b0b5928c9d1ac6a0278e365748912adfaf2ff30cf2ce4'
manifest=json.loads((PB/'FILES.json').read_text());bounds={str(S):sha(S),str(MAP):sha(MAP),str(PB/'FILES.json'):sha(PB/'FILES.json')}
files=['genes.tsv','sample_L1/units.tsv.gz','sample_L1/counts_units_by_genes.npz','sample_L2/units.tsv.gz','sample_L2/counts_units_by_genes.npz']
for f in files:
 h=sha(PB/f);assert h==manifest[f]['sha256'],f;bounds[str(PB/f)]=h
print('INPUT_HASHES_PASS',flush=True)
genes=pd.read_csv(PB/'genes.tsv',sep='\t',keep_default_na=False);valid=genes.stable_ensembl_id.str.fullmatch(r'ENSG\d{11}');g=genes[valid].copy();g['original_gene_index']=np.where(valid)[0];g.to_csv(B/'tables/gene_input_map.tsv',sep='\t',index=False)
assert not g.stable_ensembl_id.duplicated().any(),'Need explicit duplicate Ensembl aggregation'
gg=pd.DataFrame({'gene_id':g.stable_ensembl_id.values,'gene_symbol':g.gene_symbol_std.values});gg.to_csv(B/'genes.tsv',sep='\t',index=False)
assert len(gg)==36601 and gg.loc[gg.gene_symbol.eq('CKB'),'gene_id'].tolist()==['ENSG00000166165']
elig=set(pd.read_csv(MAP,sep='\t')['group']);assert len(elig)==103
roster=[];lineage=[];excluded=[];conservation=[]
for level in ['L1','L2']:
 u=pd.read_csv(PB/('sample_'+level)/'units.tsv.gz',sep='\t',keep_default_na=False)
 x=sparse.load_npz(PB/('sample_'+level)/'counts_units_by_genes.npz').tocsr()
 assert x.shape==(len(u),len(genes)) and u.pseudobulk_id.is_unique and u.n_cells.sum()==1175115
 assert (x.data>=0).all() and np.equal(x.data,np.floor(x.data)).all()
 for sample,z in u.groupby('sample_uid'):assert z.n_cells.sum()==m.loc[sample,'n_cells']
 assert set(u.sample_uid)==set(m.index)
 z=u[['pseudobulk_id','sample_uid','n_cells','level1_annotation']+(['level2_annotation'] if level=='L2' else [])].copy();z['source_row']=z.index
 for c in ['study','donor_uid','tissue','group_binary','age','sex','assay_type','proposed_background_id','main_analysis_scope_candidate','donor_resolved']:z[c]=z.sample_uid.map(m[c])
 z['cell_state']=z.level1_annotation if level=='L1' else z.level1_annotation+' | '+z.level2_annotation
 z['eligible_state']=~z.level1_annotation.isin(['Schwann cell','Epithelial cell','Unknown'])
 if level=='L2':z['eligible_state'] &=z.cell_state.isin(elig)
 z['in_scope']=z.main_analysis_scope_candidate.astype(str).str.lower().eq('true')
 selected=z[z.in_scope&z.eligible_state].copy()
 for a in z[~(z.in_scope&z.eligible_state)].itertuples():excluded.append({'level':level,'source_pseudobulk_id':a.pseudobulk_id,'cell_state':a.cell_state,'sample_uid':a.sample_uid,'reason':'unmatched_tissue' if not a.in_scope else 'original_L1_exclusion_or_not_in_frozen_L2_eligibility'})
 cols=['proposed_background_id','cell_state','study','donor_uid','tissue','group_binary']
 units=[];rr=[];cc=[]
 for i,(keys,v) in enumerate(selected.groupby(cols,sort=True,dropna=False)):
  d=dict(zip(cols,keys));d['pseudobulk_id']=level+'_'+hashlib.sha256(json.dumps(keys).encode()).hexdigest()[:16];d['n_cells']=int(v.n_cells.sum());d['sample_uids']=';'.join(sorted(v.sample_uid.unique()))
  for a in ['age','sex','assay_type']:
   vals=[str(q) for q in v[a].unique() if str(q) not in ['','nan','NA','NR']]
   if len(vals)>1:raise ValueError('Conflicting metadata '+a+' '+str(keys)+' '+str(vals))
   d[a]=vals[0] if vals else ''
  d['study_block']=d['study']
  if d['study']=='PRJNA562645':
   branches={s.split('__')[0] for s in v.sample_uid};assert len(branches)==1;d['study_block']=branches.pop()
  d['donor_resolved']=all(v.donor_resolved.astype(str).str.lower().eq('true'))
  units.append(d);rr.extend([i]*len(v));cc.extend(v.source_row.tolist())
  for q in v.itertuples():lineage.append({'level':level,'aggregate_pseudobulk_id':d['pseudobulk_id'],'source_pseudobulk_id':q.pseudobulk_id,'sample_uid':q.sample_uid,'source_row':q.source_row,'n_cells':q.n_cells})
 agg=sparse.csr_matrix((np.ones(len(rr),dtype=np.int64),(rr,cc)),shape=(len(units),len(u)))
 y=(agg@x[:,valid.values]).tocsr();un=pd.DataFrame(units)
 assert y.sum()==x[np.array(cc),:][:,valid.values].sum()
 conservation.append({'level':level,'input_selected_count_sum':int(x[np.array(cc),:][:,valid.values].sum()),'aggregate_count_sum':int(y.sum()),'selected_cells':int(selected.n_cells.sum()),'aggregate_cells':int(un.n_cells.sum()),'aggregate_units':len(un)})
 un.to_csv(B/'tables'/('all_aggregated_'+level+'_units.tsv'),sep='\t',index=False)
 ckbidx=int(np.where(gg.gene_symbol.eq('CKB'))[0][0]);ckb=un.copy();ckb['CKB_raw_count']=y[:,ckbidx].toarray().ravel();ckb['library_size_mapped_genes']=np.asarray(y.sum(1)).ravel();ckb['CKB_unadjusted_CPM']=1e6*ckb.CKB_raw_count/ckb.library_size_mapped_genes.replace(0,np.nan);ckb.to_csv(B/'tables'/('CKB_all_aggregated_'+level+'.tsv'),sep='\t',index=False)
 for (ctx,state),v in un.groupby(['proposed_background_id','cell_state'],sort=True):
  v=v[v.n_cells.ge(20)&v.donor_resolved].copy();v['comparison_group']=v.group_binary.map({'healthy':'reference','disease':'case'})
  arms=v.groupby('comparison_group').donor_uid.nunique();nc=int(arms.get('case',0));nr=int(arms.get('reference',0));name=level+'_'+ctx+'_'+hashlib.sha256(state.encode()).hexdigest()[:12]
  rec={'model_id':name,'level':level,'context':ctx,'cell_state':state,'n_units':len(v),'case_donors':nc,'reference_donors':nr,'status':'READY_FOR_DESIGN_CHECK' if min(nc,nr)>=5 else 'INSUFFICIENT_DONORS','input':str(B/'inputs'/name)};roster.append(rec)
  if min(nc,nr)<5:continue
  p=B/'inputs'/name;p.mkdir();sub=y[v.index,:].T.tocsc()
  with gzip.open(p/'counts.mtx.gz','wb') as f:
   with threadpool_limits(limits=1):io.mmwrite(f,sub,field='integer')
  v.to_csv(p/'units.tsv',sep='\t',index=False)
 print(level,'prepared',sum(q['level']==level and q['status']=='READY_FOR_DESIGN_CHECK' for q in roster),'models',flush=True)
 del x,y,agg
pd.DataFrame(roster).to_csv(B/'tables/model_roster.tsv',sep='\t',index=False);pd.DataFrame(lineage).to_csv(B/'tables/aggregation_lineage.tsv',sep='\t',index=False);pd.DataFrame(excluded).to_csv(B/'tables/excluded_partitions.tsv',sep='\t',index=False)
dump(B/'tables/COUNT_CONSERVATION.json',conservation)
for p in (B/'inputs').rglob('*'):
 if p.is_file():bounds[str(p)]=sha(p)
for p in [B/'genes.tsv',B/'tables/model_roster.tsv',B/'tests/native_v1/PASS.json']+list((B/'code').glob('*')):
 if p.is_file():bounds[str(p)]=sha(p)
atlas=R/'panvascular_CPU_theta2_L1L2_final_v1.h5ad';st=atlas.stat();base=json.loads((W/'RECEIPT.json').read_text());assert [st.st_dev,st.st_ino,st.st_size,st.st_mtime_ns]==base['atlas_stat']
dump(B/'RECEIPT.json',{'authorization':'User: 开始分析吧，帮我看一下CKB的变化; approved 3 background v2 and healthy controls','created_utc':__import__('datetime').datetime.utcnow().isoformat(),'scope':'pseudobulk L1 and eligible L2 three approved contexts; CKB focused review; no atlas or annotation mutations','atlas_sha256':base['atlas_sha256'],'atlas_stat':base['atlas_stat'],'bound_files':bounds,'gene_count':len(gg),'unmapped_features_excluded':int((~valid).sum()),'L2_eligibility':'103 frozen MAGMA eligibility labels; frequency qualification is not independent biological validation','parameters':{'min_cells_per_donor_tissue_condition_state':20,'min_donors_per_arm':5,'filterByExpr_min_count':10,'min_tested_genes':100,'edgeR':'4.0.16','normalization':'TMM','robust_dispersion':True,'robust_QL':True,'seed':0,'no_within_study_both_arm_gate':True,'max_CPU_heavy_project_jobs':4,'max_project_threads':40}})
dump(B/'PREPARATION.json',{'status':'PASS_READY_FOR_MODEL_EXECUTION','n_models':len(roster),'ready_models':sum(q['status']=='READY_FOR_DESIGN_CHECK' for q in roster),'input_hashes_verified':True,'new_aggregation_conservation':conservation,'atlas_unchanged':True})
print('PREPARATION_COMPLETE',flush=True)
