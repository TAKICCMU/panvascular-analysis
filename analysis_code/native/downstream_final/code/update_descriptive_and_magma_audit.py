from common import *
import pandas as pd,numpy as np,anndata as ad,ast
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=W/'descriptive_updated';out.mkdir(exist_ok=False);(out/'figures').mkdir()
m=pd.read_csv(AUDIT/'sample_sidecar_updated.tsv',sep='\t',keep_default_na=False).set_index('sample_uid');m['healthy_definition']='study reference group; does not imply absence of all diseases'
source=R/'production_launches/phase1_/20260913T183720_920345Z/outputs/runs/20260913T183912_636388_077099775370/counts/sample_state_counts.tsv'
x=pd.read_csv(source,sep='\t');assert set(x.sample_uid)==set(m.index)
a=ad.read_h5ad(ATLAS,backed='r');obs=a.obs[['sample_uid','study','level1_annotation','level2_annotation']].copy();a.file.close();assert len(obs)==1175115 and obs.index.is_unique
assert np.array_equal(obs.study.astype(str).to_numpy(),obs.sample_uid.astype(str).map(m.study).to_numpy()),'Study assignment changed: reassess study-balanced MAGMA'
actual=obs.groupby(['sample_uid','level1_annotation','level2_annotation'],observed=True).size();decl=x.set_index(['sample_uid','level1_annotation','level2_annotation']).n_cells
assert not decl.index.duplicated().any() and decl.sum()==len(obs)
assert np.array_equal(actual.reindex(decl.index,fill_value=0).to_numpy(),decl.to_numpy()),'Frozen counts/state membership mismatch'
assert m.n_cells.sum()==len(obs)
m.to_csv(out/'sample_study_donor_group_all247.tsv',sep='\t');l1=x.groupby(['sample_uid','level1_annotation'],observed=True).n_cells.sum().reset_index();l1['state']=l1.level1_annotation
l2=x.copy();l2['state']=l2.level1_annotation+' | '+l2.level2_annotation
reports={}
for level,df in [('L1',l1),('L2',l2)]:
 z=df.merge(m.reset_index(),on='sample_uid',validate='many_to_one',suffixes=('','_metadata'));z['fraction_all_cells']=z.n_cells/z.n_cells_metadata
 z.to_csv(out/(level+'_sample_composition.tsv.gz'),sep='\t',index=False)
 globalc=z.groupby('state').n_cells.sum().sort_values(ascending=False);pd.DataFrame({'cell_number':globalc,'percentage':globalc/len(obs)*100}).to_csv(out/(level+'_global_composition.tsv'),sep='\t')
 for field in ['study','tissue','group_binary','cohort_disease_context']:
  c=z.groupby([field,'state'],observed=True).n_cells.sum().unstack(fill_value=0);den=m.groupby(field).n_cells.sum();f=c.div(den,axis=0)
  c.to_csv(out/(level+'_'+field+'_cell_numbers.tsv'),sep='\t');f.to_csv(out/(level+'_'+field+'_cell_weighted_fraction.tsv'),sep='\t')
  mean=z.groupby([field,'state'],observed=True).fraction_all_cells.mean().unstack();mean.to_csv(out/(level+'_'+field+'_sample_mean_fraction.tsv'),sep='\t')
  donor=z.groupby([field,'donor_uid','state'],observed=True).fraction_all_cells.mean().groupby([field,'state']).mean().unstack();donor.to_csv(out/(level+'_'+field+'_donor_mean_fraction.tsv'),sep='\t')
  if level=='L1':
   fig,ax=plt.subplots(figsize=(max(9,len(f)*.38),6));f.plot.bar(stacked=True,ax=ax,width=.85,colormap='tab20');ax.set_ylabel('Fraction of all QC-passed cells');ax.set_xlabel(field);ax.legend(bbox_to_anchor=(1.01,1),fontsize=7);fig.tight_layout();fig.savefig(out/'figures'/('L1_'+field+'.png'),dpi=200);fig.savefig(out/'figures'/('L1_'+field+'.pdf'));plt.close(fig)
  else:
   fig,ax=plt.subplots(figsize=(max(10,len(f)*.35),max(9,len(f.columns)*.12)));im=ax.imshow(f.T,aspect='auto',interpolation='nearest');ax.set_yticks(range(len(f.columns)),f.columns,fontsize=4);ax.set_xticks(range(len(f)),f.index,rotation=90,fontsize=6);fig.colorbar(im,ax=ax,label='Fraction of all QC-passed cells');fig.tight_layout();fig.savefig(out/'figures'/('L2_'+field+'_heatmap.png'),dpi=200);plt.close(fig)
 cov=z[z.covarnet_include.astype(str).str.lower().eq('true')];cov.groupby('state').n_cells.sum().to_csv(out/(level+'_CoVarNet_scope_counts_ONLY.tsv'),sep='\t')
 reports[level]=len(globalc)
rec=l2[l2.n_cells>0].merge(m.reset_index()[['sample_uid','study','donor_uid','tissue']],on='sample_uid',validate='many_to_one').groupby('state').agg(n_cells=('n_cells','sum'),samples=('sample_uid','nunique'),studies=('study','nunique'),donors=('donor_uid','nunique'),tissues=('tissue','nunique'))
rec.to_csv(out/'L2_cross_study_recurrence.tsv',sep='\t')
l1[l1.level1_annotation.eq('Unknown')].merge(m.reset_index(),on='sample_uid').to_csv(out/'Unknown_origin.tsv',sep='\t',index=False)
ma=R/'production_launches/magma_recovery_roundtrip_fix_20260914T0930/magma_adapter_v3/src/external_modules.py';text=ma.read_text();tree=ast.parse(text);fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='specificity_from_counts');src=ast.get_source_segment(text,fn)
(out/'actual_MAGMA_specificity_function.py.txt').write_text(src)
assert 'study_key' in src and 'nc[i]>=min_study_cells' in src
assert not any(q in src for q in ['donor_uid','group_binary','biological_sample_uid','covarnet_include'])
mp=R/'production_launches/magma_recovery_roundtrip_fix_20260914T0930/outputs/runs/20260914T013548_608658_d62bd9f120de/magma'
complete=json.loads((mp/'COMPLETE.json').read_text());dump(out/'MAGMA_completion_source.json',complete)
dump(out/'MAGMA_METADATA_IMPACT.json',dict(status='NO_CHANGE_TO_SPECIFICITY_INPUTS',actual_production_source=str(ma),source_sha256=sha(ma),cell_counts_match_phase1=True,study_per_cell_identical=True,L1_L2_unchanged=True,cell_weighted='mean CP1e6/log1p across all cells per state',study_balanced='equal mean of state-specific study means; >=20 cells in each contributing study',unused_fields=['donor','sample aggregation','healthy/disease group'],CoVarNet_exclusions_not_applied=True,qualification_thresholds_unchanged=True,rerun_MAGMA_required=False))
(out/'README_CN.md').write_text('本轮只更新样本/供者/分组连接及描述图表。全量细胞和正式L1/L2未改变，原全局细胞总数可复用；旧分组、样本/供者汇总图由本目录替代。比例明确区分细胞加权、样本均值和供者均值；L2保留原within-parent比例。健康组表示本研究参照组，不代表供者没有其他疾病。CoVarNet排除规则仅用于专门标记的范围表，不扩展到图谱描述、MAGMA或通讯。\n原稿对应Figure1d、Figure1f及Supplementary Figure1；不继承旧稿错误细胞数、固定71状态或14模块。无疾病显著性检验。\n')
dump(out/'COMPLETE.json',dict(status='COMPLETE_DESCRIPTIVE_UPDATE',n_cells=len(obs),samples=len(m),donors=m.donor_uid.nunique(),studies=m.study.nunique(),categories=reports,source_counts_sha256=sha(source),metadata_sha256=sha(AUDIT/'sample_sidecar_updated.tsv')))
print('Descriptive update and live MAGMA weighting audit complete',flush=True)
