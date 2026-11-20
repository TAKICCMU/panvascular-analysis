#!/usr/bin/env python3
"""Audited CD-NMF sensitivity against a frozen nsNMF reference, not same-method CV."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]=os.environ.get('HELDOUT_THREADS','2')
import argparse,hashlib,importlib.util,json,sys,time,warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment,nnls
from scipy.stats import random_table,spearmanr
from sklearn.decomposition import NMF
from sklearn.exceptions import ConvergenceWarning


def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(p,d):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,allow_nan=False,default=str)+'\n')
def finite_matrix(x):
 z=np.asarray(x,dtype=float)
 if z.ndim!=2 or not np.isfinite(z).all():raise ValueError('Expected finite 2D matrix')
 return z

def validate_counts(x):
 z=finite_matrix(x)
 if (z<0).any() or (z!=np.floor(z)).any() or (z>2**53-1).any():raise ValueError('Counts must be exact nonnegative integers before conversion')
 if sum(map(int,z.ravel()))>np.iinfo(np.int64).max:raise ValueError('Integer margin overflow')
 return z.astype(np.int64)
def minmax_fit_transform(values):
 values=finite_matrix(values);low=np.min(values,axis=1);high=np.max(values,axis=1);scale=np.maximum(high-low,1e-12)
 return np.clip((values-low[:,None])/scale[:,None],0,1),low,high

def minmax_apply(values,low,high):return np.clip((finite_matrix(values)-low[:,None])/np.maximum(high-low,1e-12)[:,None],0,1)
def normalize_columns(x):return x/np.maximum(np.linalg.norm(x,axis=0,keepdims=True),1e-12)
def cosine_matrix(a,b):return normalize_columns(a).T@normalize_columns(b)
def top_jaccard(a,b,top_n=10):
 aa=set(np.argsort(a)[-top_n:]);bb=set(np.argsort(b)[-top_n:]);return len(aa&bb)/max(len(aa|bb),1)
def fit_nmf(matrix,k,seed,max_iter=1000):
 matrix=finite_matrix(matrix)
 if (matrix<0).any() or min(matrix.shape)<k:raise ValueError('Invalid NMF nonnegative dimension')
 model=NMF(n_components=k,init='nndsvda',solver='cd',beta_loss='frobenius',max_iter=max_iter,tol=1e-5,random_state=seed)
 with warnings.catch_warnings(record=True) as ws:
  warnings.simplefilter('always');w=model.fit_transform(matrix)
 h=model.components_;err=float(np.linalg.norm(matrix-w@h)/max(np.linalg.norm(matrix),1e-12))
 return w,h,err,dict(n_iter=int(model.n_iter_),solver_converged=not any(issubclass(q.category,ConvergenceWarning) for q in ws),warnings=[str(q.message) for q in ws])
def module_coherence(proportions,indices):
 if proportions.shape[1]<3 or len(indices)<2:return dict(median_pair_rho=np.nan,positive_pair_fraction=np.nan,n_pairs=0)
 rr=[]
 for i in range(len(indices)):
  for j in range(i+1,len(indices)):
   with warnings.catch_warnings():
    warnings.simplefilter('ignore');rho=spearmanr(proportions[indices[i]],proportions[indices[j]],nan_policy='omit').statistic
   if np.isfinite(rho):rr.append(float(rho))
 return dict(median_pair_rho=float(np.median(rr)) if rr else np.nan,positive_pair_fraction=float(np.mean(np.asarray(rr)>0)) if rr else np.nan,n_pairs=len(rr))
def empirical_p(observed,null):
 null=np.asarray(null,float);null=null[np.isfinite(null)]
 return (float((1+(null>=observed).sum())/(1+len(null))) if np.isfinite(observed) and len(null) else np.nan,len(null))
def bh(values):
 x=np.asarray(values,float);out=np.full(len(x),np.nan);ii=np.flatnonzero(np.isfinite(x));order=ii[np.argsort(x[ii])];m=len(order)
 if m:out[order]=np.minimum.accumulate((x[order]*m/np.arange(1,m+1))[::-1])[::-1].clip(0,1)
 return out

def matrix_file(p):
 x=pd.read_csv(p,sep='\t',index_col=0)
 if not x.index.is_unique or not x.columns.is_unique:raise ValueError('Duplicate matrix IDs '+str(p))
 finite_matrix(x);return x

def compare_matrix(actual,expected,*,atol=0):
 if not actual.index.equals(expected.index) or not actual.columns.equals(expected.columns):raise ValueError('Exact ordered row/profile IDs differ')
 a=finite_matrix(actual);b=finite_matrix(expected)
 if not np.allclose(a,b,rtol=0,atol=atol):raise ValueError('Numeric mismatch beyond absolute tolerance '+str(atol))
 return float(np.max(np.abs(a-b))) if a.size else 0.
def verify_manifest(path):
 d=json.loads(Path(path).read_text())
 if d.get('status')!='completed':raise ValueError('Incomplete dependency '+str(path))
 for key in ['dependency_sha256','output_sha256']:
  if not d.get(key):raise ValueError('Missing manifest SHA records')
  for p,h in d[key].items():
   if sha(p)!=h:raise ValueError('Dependency/output SHA mismatch '+str(p))
 return d

def load_bridge(bridge):
 bridge=Path(bridge).resolve();done=verify_manifest(bridge/'BRIDGE_COMPLETE.json');k=int(done['selected_K']);txt=(bridge/'selected_K.txt').read_text().splitlines();kk=[int(q.split('=',1)[1]) for q in txt if q.startswith('K_FINAL=')]
 if kk!=[k]:raise ValueError('Selected-K mismatch')
 c=matrix_file(bridge/'CoVarNet_node_counts.tsv');counts=validate_counts(c);p=matrix_file(bridge/'CoVarNet_node_proportion_of_all_atlas_cells.tsv');x=matrix_file(bridge/'NMF_input_minmax.tsv');w=matrix_file(bridge/f'module_node_weights_K{k:02d}.tsv');h=matrix_file(bridge/f'sample_module_abundance_K{k:02d}.tsv');md=pd.read_csv(bridge/'CoVarNet_sample_metadata.tsv',sep='\t',keep_default_na=False).set_index('covarnet_sample_uid')
 if not md.index.is_unique or not md.index.equals(c.columns):raise ValueError('Metadata profile order/identity mismatch')
 if md.study_canonical.eq('').any():raise ValueError('Missing canonical study')
 for q in [p,x]:
  if not q.index.equals(c.index) or not q.columns.equals(c.columns):raise ValueError('Bridge matrix exact IDs mismatch')
 if not w.index.equals(c.index) or w.shape[1]!=k or not h.index.equals(c.columns) or not h.columns.equals(w.columns):raise ValueError('W/H/module identity mismatch')
 if (w.to_numpy()<0).any() or not np.allclose(w.sum(axis=0),1,rtol=0,atol=1e-12):raise ValueError('Reference basis invalid or not column normalized')
 den=pd.to_numeric(md.n_total_atlas_cells).to_numpy(float)
 if not np.isfinite(den).all() or (den<=0).any() or (den!=np.floor(den)).any() or (counts.sum(0)>den).any():raise ValueError('Invalid all-atlas denominator')
 prop=counts/den;pe=compare_matrix(p,pd.DataFrame(prop,index=c.index,columns=c.columns),atol=1e-12);recomputed,_,_=minmax_fit_transform(prop);xe=compare_matrix(x,pd.DataFrame(recomputed,index=c.index,columns=c.columns),atol=1e-12)
 return dict(bridge=bridge,k=k,modules=list(w.columns),counts_df=c,counts=counts,props=prop,reference_W=w.to_numpy(),metadata=md,all_atlas_denominator=den,input_numeric_error=xe,proportion_numeric_error=pe,bridge_manifest_sha256=sha(bridge/'BRIDGE_COMPLETE.json'))

def production_gate(root,annotation_dir,bridge,d):
 root=Path(root).resolve();run=bridge.parent
 for name in ['INPUTS_COMPLETE.json','STABLE_COMPLETE.json']:verify_manifest(run/name)
 spec=importlib.util.spec_from_file_location('shared_annotation_validator',root/'pseudobulk/code/validate_frozen_annotations.py');v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
 ov=pd.read_csv(root/'metadata/cell_metadata_overlay.tsv.gz',sep='\t',usecols=['cell_id','sample_uid'],dtype=str);a,g=v.validate_frozen_annotations(annotation_dir,['Level1','Level2'],expected_cell_ids=ov.cell_id,review_only=False);fr=json.loads(Path(g['freeze_file']).read_text())
 if fr.get('is_partial') is True or not g['publication_ready'] or g['review_only']:raise ValueError('Explicit production complete L2 required')
 a['Level2_CoVarNet_eligible']=v.strict_bool(a.Level2_CoVarNet_eligible)
 deps=json.loads((run/'INPUTS_COMPLETE.json').read_text())['dependency_sha256'];required=[Path(g['annotation_file']),Path(g['freeze_file']),Path(g['L1_freeze_file']),root/'metadata/cell_metadata_overlay.tsv.gz',root/'metadata/sample_metadata_contract.tsv',Path(fr['roster_file'])]
 for p in required:
  if deps.get(str(p))!=sha(p):raise ValueError('Bridge input chain does not bind current freeze/metadata '+str(p))
 if sha(fr['roster_file'])!=fr.get('roster_sha256'):raise ValueError('Actual roster SHA mismatch')
 sys.path.insert(0,str(root/'covarnet/code'));from build_covarnet_inputs import build
 samples=pd.read_csv(root/'metadata/sample_metadata_contract.tsv',sep='\t');r=build(a,ov,samples);counts=r['counts'].loc[d['counts_df'].index,d['counts_df'].columns];err=compare_matrix(d['counts_df'],counts,atol=0)
 md=r['profiles'].set_index('covarnet_sample_uid').loc[d['metadata'].index]
 for col in ['study_canonical','covarnet_donor_key','condition_v2','tissue_v2','donor_uid_analysis','donor_resolution','n_total_atlas_cells']:
  if not md[col].astype(str).equals(d['metadata'][col].astype(str)):raise ValueError('Authoritative donor/profile metadata changed '+col)
 return dict(annotation_gate=g,exact_current_cell_to_count_rebuild=True,count_max_abs_error=err,n_all_atlas_cells=len(a),n_profiles=len(md),dependency_sha256={str(p):sha(p) for p in required+[root/'pseudobulk/code/validate_frozen_annotations.py',root/'covarnet/code/build_covarnet_inputs.py',root/'composition/code/composition_inputs.py']})

def load_evidence(path,d):
 if path is None:return None,{'status':'not_supplied_similarity_P_not_estimable'}
 c=json.loads(Path(path).read_text());p=Path(c['evidence_file'])
 if c.get('status')!='completed' or c.get('bridge_complete_sha256')!=d['bridge_manifest_sha256'] or sha(p)!=c.get('evidence_sha256'):raise ValueError('Evidence not bound to this exact bridge')
 e=pd.read_csv(p,sep='\t').set_index('module')
 if not e.index.is_unique or set(e.index)!=set(d['modules']) or 'leave_one_study' not in e:raise ValueError('Dynamic module evidence identity/schema mismatch')
 return e.loc[d['modules']],dict(contract_sha256=sha(path),evidence_sha256=sha(p))

def analyze(d,out,*,n_perm=200,min_heldout=5,seed=20260824,denominator='all_atlas',evidence=None):
 out=Path(out);out.mkdir(parents=True,exist_ok=True)
 if any(out.iterdir()):raise ValueError('Analysis output must be new/empty; do not overwrite')
 c=d['counts'];k=d['k'];mods=d['modules'];wref=d['reference_W'];md=d['metadata'];den=d['all_atlas_denominator'] if denominator=='all_atlas' else np.maximum(c.sum(0),1);props=c/den;rng=np.random.default_rng(seed);null=[];fitlog=[]
 for perm in range(1,n_perm+1):
  nc=random_table.rvs(c.sum(1),c.sum(0),random_state=rng);assert np.array_equal(nc.sum(1),c.sum(1)) and np.array_equal(nc.sum(0),c.sum(0));np_=nc/den;xx,_,_=minmax_fit_transform(np_);w,_,err,fit=fit_nmf(xx,k,seed+perm,500);fitlog.append(dict(stage='fixed_margin_null',iteration=perm,**fit));cos=cosine_matrix(wref,w)
  for i,mod in enumerate(mods):
   j=int(np.argmax(cos[i]));co=module_coherence(np_,np.argsort(wref[:,i])[-5:]);null.append(dict(permutation=perm,module=mod,best_null_component=j+1,maximum_cosine_to_reference=cos[i,j],top10_jaccard_to_reference=top_jaccard(wref[:,i],w[:,j]),reference_top5_null_median_pair_rho=co['median_pair_rho'],reference_top5_null_positive_pair_fraction=co['positive_pair_fraction'],null_relative_reconstruction_error=err,solver_converged=fit['solver_converged']))
 nt=pd.DataFrame(null);nt.to_csv(out/'01_fixed_margin_compositional_null.tsv.gz',sep='\t',index=False);summary=[]
 for i,mod in enumerate(mods):
  t=nt[nt.module==mod];co=module_coherence(props,np.argsort(wref[:,i])[-5:]);obs=float(evidence.loc[mod,'leave_one_study']) if evidence is not None else np.nan;ps,ns=empirical_p(obs,t.maximum_cosine_to_reference);pc,nc=empirical_p(co['median_pair_rho'],t.reference_top5_null_median_pair_rho);summary.append(dict(module=mod,n_permutations=n_perm,observed_leave_one_study_similarity=obs,empirical_P_similarity=ps,n_valid_null_similarity=ns,observed_top5_median_pair_rho=co['median_pair_rho'],empirical_P_coherence=pc,n_valid_null_coherence=nc,n_null_solver_converged=int(t.solver_converged.sum()),null_max_cosine_median=t.maximum_cosine_to_reference.median(),null_max_cosine_95pct=t.maximum_cosine_to_reference.quantile(.95),null_top5_rho_median=t.reference_top5_null_median_pair_rho.median(),null_top5_rho_95pct=t.reference_top5_null_median_pair_rho.quantile(.95)))
 ns=pd.DataFrame(summary)
 for metric in ['similarity','coherence']:ns['BH_FDR_'+metric]=bh(ns['empirical_P_'+metric])
 ns.to_csv(out/'02_fixed_margin_null_summary.tsv',sep='\t',index=False)
 held=[];ineligible=[];seedlog=[];eligible_index=0
 for study,holdmd in md.groupby('study_canonical',sort=True):
  mask=d['counts_df'].columns.isin(holdmd.index);nh=int(mask.sum());ntr=int((~mask).sum())
  if nh<min_heldout or ntr<k+2:ineligible.append(dict(study=study,n_hold=nh,n_train=ntr,reason='heldout_profiles_below_minimum' if nh<min_heldout else 'training_profiles_below_Kplus2'));continue
  xx,lo,hi=minmax_fit_transform(props[:,~mask]);hx=minmax_apply(props[:,mask],lo,hi);fitseed=seed+10000+2*eligible_index;eligible_index+=1;w,_,err,fit=fit_nmf(xx,k,fitseed,1500);fitlog.append(dict(stage='heldout',study=study,**fit));co=cosine_matrix(wref,w);ri,ci=linear_sum_assignment(-co);alignment=dict(zip(ri,ci));hh=np.column_stack([nnls(w,hx[:,j])[0] for j in range(nh)]);herr=float(np.linalg.norm(hx-w@hh)/max(np.linalg.norm(hx),1e-12));sid=f'study{eligible_index:03d}';np.savez_compressed(out/(sid+'_train_transform.npz'),low=lo,high=hi,train_profile_ids=np.asarray(d['counts_df'].columns[~mask],str),heldout_profile_ids=np.asarray(d['counts_df'].columns[mask],str),W_train=w,H_heldout=hh);seedlog.append(dict(study=study,fit_seed=fitseed,transform_file=sid+'_train_transform.npz'))
  for i,mod in enumerate(mods):
   j=alignment[i];ch=module_coherence(props[:,mask],np.argsort(w[:,j])[-5:]);valid=np.isfinite(ch['median_pair_rho']);struct=bool(co[i,j]>=.70 and top_jaccard(wref[:,i],w[:,j])>=.30);rep=bool(ch['positive_pair_fraction']>=.60 and ch['median_pair_rho']>0) if valid else pd.NA;held.append(dict(heldout_study=study,module=mod,n_heldout_samples=nh,n_training_samples=ntr,aligned_training_component=j+1,training_to_reference_cosine=co[i,j],training_to_reference_top10_jaccard=top_jaccard(wref[:,i],w[:,j]),heldout_top5_median_pair_rho=ch['median_pair_rho'],heldout_top5_positive_pair_fraction=ch['positive_pair_fraction'],heldout_top5_n_pairs=ch['n_pairs'],training_relative_reconstruction_error=err,heldout_relative_projection_error=herr,heldout_module_score_SD=float(np.std(hh[j])),solver_converged=fit['solver_converged'],structural_replication=struct,within_study_coherence_replication=rep,joint_replication=(struct and rep) if valid else pd.NA,coherence_estimable=valid))
 pd.DataFrame(ineligible,columns=['study','n_hold','n_train','reason']).to_csv(out/'ineligible_studies.tsv',sep='\t',index=False);ht=pd.DataFrame(held);ht.to_csv(out/'03_strict_heldout_study_projection.tsv',sep='\t',index=False);hr=[]
 for mod in mods:
  t=ht[ht.module==mod] if len(ht) else pd.DataFrame();valid=t[t.coherence_estimable] if len(t) else t
  hr.append(dict(module=mod,n_eligible_heldout_studies=len(t),n_coherence_estimable_studies=len(valid),n_nonconverged_training_fits=int((~t.solver_converged).sum()) if len(t) else 0,structural_replication_fraction=float(t.structural_replication.mean()) if len(t) else np.nan,coherence_replication_fraction=float(valid.within_study_coherence_replication.mean()) if len(valid) else np.nan,joint_replication_fraction=float(valid.joint_replication.mean()) if len(valid) else np.nan,median_training_to_reference_cosine=float(t.training_to_reference_cosine.median()) if len(t) else np.nan,median_heldout_pair_rho=float(valid.heldout_top5_median_pair_rho.median()) if len(valid) else np.nan,minimum_training_to_reference_cosine=float(t.training_to_reference_cosine.min()) if len(t) else np.nan,median_heldout_projection_error=float(t.heldout_relative_projection_error.median()) if len(t) else np.nan))
 pd.DataFrame(hr).to_csv(out/'04_heldout_study_summary.tsv',sep='\t',index=False);save(out/'SOLVER_AUDIT.json',fitlog);save(out/'TRAIN_ONLY_TRANSFORM_AUDIT.json',seedlog)
 return dict(status='completed_CD_algorithm_sensitivity',reference_algorithm='frozen_R_nsNMF',refit_algorithm='sklearn_CD_NMF_nndsvda_frobenius',same_algorithm_validation=False,selected_K=k,modules=mods,n_perm=n_perm,seed=seed,min_heldout_profiles=min_heldout,denominator=denominator,n_eligible_heldout_studies=eligible_index,n_solver_nonconvergence=sum(not z['solver_converged'] for z in fitlog),all_fit_convergence_verified=all(z['solver_converged'] for z in fitlog),multiple_testing='BH across all finite module tests separately for similarity and coherence; families never pooled',similarity_P_interpretation='historical cross-algorithm calibration, not same-method hypothesis test',heldout_fit_seed='seed+10000+2*zero_based_eligible_study_index; preserves original two-target seed sequence independently of current module count',leakage_scope='Training scaling/refit/top5 selection use training only; frozen full-data node universe/K/reference used for conditional naming and structural comparison, not nested independent discovery',missing_observed_similarity='NA if current bound followup evidence absent',figures_generated=False)


def main():
 p=argparse.ArgumentParser();p.add_argument('--takeover-root',type=Path,required=True);p.add_argument('--annotation-dir',type=Path,required=True);p.add_argument('--bridge',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--evidence-contract',type=Path);p.add_argument('--denominator',choices=['all_atlas','legacy_stable_node_closure'],default='all_atlas');p.add_argument('--n-perm',type=int,default=200);p.add_argument('--min-heldout-profiles',type=int,default=5);p.add_argument('--seed',type=int,default=20260824);p.add_argument('--run',action='store_true');a=p.parse_args()
 if min(a.n_perm,a.min_heldout_profiles)<1:raise ValueError('Positive parameters required')
 
 if json.loads((a.bridge/'BRIDGE_COMPLETE.json').read_text()).get('fixture_only'):raise ValueError('Synthetic fixture bridge cannot enter production CLI')
 if not a.out.resolve().is_relative_to(Path(__file__).resolve().parents[1]):raise ValueError('Output must stay in exclusive heldout_validation_v2 directory')
 d=load_bridge(a.bridge);gate=production_gate(a.takeover_root,a.annotation_dir,a.bridge.resolve(),d);ev,eg=load_evidence(a.evidence_contract,d)
 info=dict(status='validated_production_inputs_models_not_run',gate=gate,evidence=eg,bridge_sha256=d['bridge_manifest_sha256'],input_numeric_max_abs_error=d['input_numeric_error'],proportion_numeric_max_abs_error=d['proportion_numeric_error'],code_sha256=sha(__file__))
 if not a.run:save(a.out/'INPUT_VALIDATION_ONLY.json',info);print(json.dumps(info,indent=2));return
 result=analyze(d,a.out,n_perm=a.n_perm,min_heldout=a.min_heldout_profiles,seed=a.seed,denominator=a.denominator,evidence=ev);result['inputs']=info;import scipy,sklearn;result['versions']=dict(python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,sklearn=sklearn.__version__);result['parameters']=dict(init='nndsvda',solver='cd',beta_loss='frobenius',tol=1e-5,null_max_iter=500,heldout_max_iter=1500,top_coherence=5,top_jaccard=10,structural_cosine=.70,structural_jaccard=.30,coherence_positive_pair_fraction=.60,coherence_rho_strictly_positive=True);result['output_sha256']={str(q.name):sha(q) for q in sorted(a.out.iterdir()) if q.is_file()};save(a.out/'RUN_COMPLETE.json',result)
if __name__=='__main__':main()
