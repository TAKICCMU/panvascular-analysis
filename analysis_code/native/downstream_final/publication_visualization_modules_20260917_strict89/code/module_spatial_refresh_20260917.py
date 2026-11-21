import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='8'
os.environ['CUDA_VISIBLE_DEVICES']=''
from pathlib import Path
import numpy as np,pandas as pd,ast,json,hashlib,collections,itertools,time,math
from scipy.stats import spearmanr,binomtest
from sklearn.neighbors import NearestNeighbors
from threadpoolctl import threadpool_limits
threadpool_limits(8)
W=Path('/data/anzhen_vhd/scRNAseq/codex/panvascular_downstream_storage_v1/expanded_parallel_20260916_v1');J=W/'overnight_manuscript_integration_20260917/spatial_target_excluded_joint_v1';M=W/'overnight_manuscript_integration_20260917/spatial_target_excluded_v1/outputs';O=W/'publication_visualization_modules_20260917_strict89';O.mkdir(exist_ok=True)
legacy=W/'overnight_manuscript_integration_20260917/spatial_joint_evidence_v1/code/legacy_spatial.py';K_NEIGHBORS=6
funcs=['zscore','matrix_correlations','upper_mean','spatial_cross_matrix','module_scores','lineage_from_state','bh_fdr','strip_cell2location_prefix'];tree=ast.parse(legacy.read_text());exec(compile(ast.Module(body=[x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name in funcs],type_ignores=[]),str(legacy),'exec'),globals())
def neighbour_lag(v,co,k=6):
 ids=NearestNeighbors(n_neighbors=min(k+1,len(co)),n_jobs=8).fit(co).kneighbors(co,return_distance=False);ids=np.array([[j for j in row if j!=i][:k] for i,row in enumerate(ids)]);return v[ids].mean(axis=1)
weights_path=W/'covarnet/stable_NMF/module_node_weights_K14.tsv';node_path=W/'covarnet/CoVarNet_node_parent_map.tsv';weights=pd.read_csv(weights_path,sep='\t').set_index('cell_state');nm=pd.read_csv(node_path,sep='\t');nm['label']=nm.parent_lineage+' | '+nm.level2_annotation;weights.index=weights.index.map(nm.set_index('cell_state').label);assert weights.shape==(89,14) and weights.index.notna().all()
factors=pd.read_csv(M/'reference_factors.tsv',sep='\t');fmap=dict(zip(factors.reference_state,factors.reference_group));spots=pd.read_csv(J/'spot_joint_source.tsv.gz',sep='\t').set_index('spot_id');audit=pd.read_csv(J/'section_input_audit.tsv',sep='\t');assert len(audit)==49 and len(spots)==45348
programs=[c for c in spots if c.startswith('program__')];assert len(programs)==6
B=1000;SEED=20260917;rng=np.random.default_rng(SEED);allstates=list(factors.reference_group);assert len(allstates)==128;parents=np.array([lineage_from_state(x) for x in allstates]); tops={m:list(weights[m].sort_values(ascending=False).head(10).index) for m in weights};nulls={};drawrows=[]
for mod,top in tops.items():
 counts=collections.Counter(lineage_from_state(x) for x in top)
 for method in ['legacy_fallback','strict_parent']:
  draw=[]
  for b in range(B):
   selected=[];fallback=0
   if method=='strict_parent':
    for par,n in counts.items():selected.extend(rng.choice(np.array(weights.index)[np.array([lineage_from_state(x) for x in weights.index])==par],n,replace=False))
   else:
    for par in [lineage_from_state(x) for x in top]:
     pool=[x for x in allstates if lineage_from_state(x)==par and x not in top and x not in selected]
     if not pool:pool=[x for x in allstates if x not in top and x not in selected];fallback+=1
     selected.append(str(rng.choice(pool)))
   ix=[allstates.index(x) for x in selected];draw.append(ix);drawrows.append({'module':mod,'method':method,'draw':b,'fallback_nodes':fallback,'states':';'.join(selected),'observed_top10_weights':';'.join(map(str,weights.loc[top,mod])),'strict_parent_combination_count':str(math.prod(math.comb(sum(lineage_from_state(x)==par for x in weights.index),n) for par,n in counts.items())) if method=='strict_parent' else 'NA'})
  nulls[(mod,method)]=np.array(draw)
pd.DataFrame(drawrows).to_csv(O/'null_draw_identity.tsv.gz',sep='\t',index=False)
rows=[];nullrecords=[];hashes={};tri=np.triu_indices(10,1)
for gse,au in audit.groupby('GSE',sort=True):
 for post,prefix in [('q05','q05'),('mean','means')]:
  path=M/f'{gse}_{prefix}_named.tsv.gz';hashes[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest();ab=pd.read_csv(path,sep='\t',index_col=0);nuisance=M/f'{gse}_{prefix}_nuisance.tsv.gz';hashes[str(nuisance)]=hashlib.sha256(nuisance.read_bytes()).hexdigest();ab=pd.concat([ab,pd.read_csv(nuisance,sep='\t',index_col=0)],axis=1);ab.columns=[fmap[strip_cell2location_prefix(x)] for x in ab.columns];ab=ab.reindex(columns=allstates);assert ab.shape[1]==128 and not ab.isna().any().any()
  for mr in au.itertuples():
   ss=spots[spots.section_id.eq(mr.section_id)];z=ab.loc[ss.index];assert len(z)==mr.n_spots;v=np.log1p(z.to_numpy());corr=matrix_correlations(v);cross=spatial_cross_matrix(v,ss[['x','y']].to_numpy());scores=module_scores(z,weights)
   for mod in weights:
    top=tops[mod];ix=[allstates.index(x) for x in top];observed=[upper_mean(corr,ix),upper_mean(cross,ix)];base={'GSE':gse,'section_id':mr.section_id,'biological_unit_id':mr.biological_unit_id,'posterior':post,'module':mod,'QC_fallback':mr.QC_fallback,'n_spots':len(z),'within_spot_coherence':observed[0],'spatial_lag_coherence':observed[1]}
    for prog in programs:base['rho__'+prog.removeprefix('program__')]=float(spearmanr(scores[mod],ss[prog]).statistic)
    for method in ['legacy_fallback','strict_parent']:
     ni=nulls[(mod,method)];nv=np.column_stack([mat[ni[:,tri[0]],ni[:,tri[1]]].mean(axis=1) for mat in [corr,cross]])
     for j,key in enumerate(['within','lag']):base[f'{method}__{key}__p']=(1+np.sum(nv[:,j]>=observed[j]))/(B+1);base[f'{method}__{key}__nullmean']=nv[:,j].mean()
     nullrecords.append((mr.section_id,post,mod,method,nv))
    rows.append(base)
   print('DONE',post,mr.section_id,flush=True)
 pd.DataFrame(rows).to_csv(O/'module_section_stats.tsv',sep='\t',index=False)
s=pd.DataFrame(rows);units=s.groupby(['posterior','module','GSE','biological_unit_id']).median(numeric_only=True).reset_index();units.to_csv(O/'module_biological_unit_stats.tsv',sep='\t',index=False);datasets=units.groupby(['posterior','module','GSE']).median(numeric_only=True).reset_index();datasets.to_csv(O/'module_dataset_stats.tsv',sep='\t',index=False)
# Preserve common draw identities through section->biological unit->dataset aggregation.
lookup={(a,b,c,d):v for a,b,c,d,v in nullrecords};coherence=[]
for post in ['q05','mean']:
 for mod in weights:
  obs=datasets.query('posterior==@post and module==@mod'); med=obs[['within_spot_coherence','spatial_lag_coherence']].median().to_numpy()
  for method in ['legacy_fallback','strict_parent']:
   bydataset=[]
   for gse,au in audit.groupby('GSE'):
    byunit=[]
    for uid,x in au.groupby('biological_unit_id'):byunit.append(np.median(np.stack([lookup[(sec,post,mod,method)] for sec in x.section_id]),axis=0))
    bydataset.append(np.median(np.stack(byunit),axis=0))
   nv=np.median(np.stack(bydataset),axis=0)
   for j,metric in enumerate(['within_spot','spatial_lag']):coherence.append({'posterior':post,'module':mod,'null_method':method,'metric':metric,'observed_median':med[j],'null_median':float(np.median(nv[:,j])),'observed_minus_null':med[j]-float(np.median(nv[:,j])),'P':float((1+sum(nv[:,j]>=med[j]))/(B+1)),'datasets':7,'biological_units':35,'draws':B})
c=pd.DataFrame(coherence);c['FDR_BH_28']=np.nan
for key,x in c.groupby(['posterior','null_method']):c.loc[x.index,'FDR_BH_28']=bh_fdr(x.P)
c.to_csv(O/'module_coherence_meta.tsv',sep='\t',index=False)
pr=[]
for (post,mod),x in datasets.groupby(['posterior','module']):
 for prog in programs:
  col='rho__'+prog.removeprefix('program__');v=x[col].dropna();n=len(v);pos=sum(v>0);pr.append({'posterior':post,'module':mod,'program':prog.removeprefix('program__'),'n_datasets':n,'positive_datasets':pos,'median_rho':float(v.median()),'min_rho':float(v.min()),'max_rho':float(v.max()),'P_direction_positive':binomtest(pos,n,.5,alternative='greater').pvalue if n else np.nan})
pr=pd.DataFrame(pr);pr['FDR_BH_84']=np.nan
for post,x in pr.groupby('posterior'):pr.loc[x.index,'FDR_BH_84']=bh_fdr(x.P_direction_positive)
pr.to_csv(O/'module_program_direction_meta.tsv',sep='\t',index=False)
for p in [weights_path,node_path,legacy,J/'spot_joint_source.tsv.gz',J/'section_input_audit.tsv',M/'reference_factors.tsv',Path(__file__)]:hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
(O/'SOURCE_SHA256.json').write_text(json.dumps(hashes,indent=2));(O/'COMPLETE.json').write_text(json.dumps({'status':'COMPLETE','modules':14,'nodes':89,'sections':49,'spots':45348,'biological_units':35,'datasets':7,'posterior':['q05','mean'],'seed':SEED,'null_draws':B,'threads':8,'coherence_BH_family':28,'program_BH_family':84,'method_note':'legacy_fallback excludes observedtop10 and maycrossparent; strict_parent samples sameparentcounts withoutreplacement frommapped89eligiblepool allowingobservednodes; top10weights attached in nodeorder; originalcoherence statistic unweighted. Common nulldrawstates acrosssections. Coherence aggregatedmediansectionswithinunit thenunitswithindataset thendatasets. Programs mediansectionrho withinunit thenunitswithindataset; one-sidedpositive binomial acrossdatasets. No confidenceinterval inferred from7datasets.'},indent=2));print('COMPLETE',flush=True)
