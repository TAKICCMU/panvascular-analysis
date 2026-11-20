from common import *
import pandas as pd,numpy as np,base64,html,zipfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=W/'descriptive_updated_v2';x=pd.read_csv(out/'L1_sample_composition.tsv.gz',sep='\t');m=pd.read_csv(out/'sample_study_donor_group_all247.tsv',sep='\t',keep_default_na=False)
f=x.pivot(index='sample_uid',columns='state',values='fraction_all_cells');order=m.sort_values(['study','sample_uid']).sample_uid;f=f.loc[order]
for i in range(0,len(f),65):
 fig,ax=plt.subplots(figsize=(22,8));f.iloc[i:i+65].plot.bar(stacked=True,ax=ax,width=.9,colormap='tab20');ax.tick_params(axis='x',labelsize=6);ax.legend(bbox_to_anchor=(1.01,1),fontsize=8);ax.set_ylabel('Fraction of all QC-passed cells');fig.tight_layout();fig.savefig(out/'figures'/f'L1_individual_samples_{i+1:03d}.png',dpi=180);plt.close(fig)
rare=[]
for state,z in x.groupby('state'):
 total=int(z.n_cells.sum());r=dict(state=state,n_cells=total,percentage_all=100*total/1175115)
 for key in ['sample_uid','study','donor_uid']:
  n=z.groupby(key).n_cells.sum().sort_values(ascending=False);r[key+'_number']=int((n>0).sum());r[key+'_largest']=n.index[0];r[key+'_largest_fraction']=float(n.iloc[0]/total) if total else None
 r['majority_one_sample']=r['sample_uid_largest_fraction']>.5;r['majority_one_study']=r['study_largest_fraction']>.5;rare.append(r)
pd.DataFrame(rare).to_csv(out/'L1_population_concentration_audit.tsv',sep='\t',index=False)
rows=[];completion=json.loads((W/'pseudobulk/COMPLETE.json').read_text())
for r in pd.read_csv(W/'pseudobulk/model_roster.tsv',sep='\t').to_dict('records'):
 p=Path(r['input'])/'results';assert sha(p/'DE.tsv')==completion['model_results'][r['model_id']];z=pd.read_csv(p/'DE.tsv',sep='\t');model=json.loads((p/'MODEL.json').read_text());assert model['status']=='COMPLETE'
 assert np.isfinite(z[['PValue','FDR','logFC']]).all().all() and z.PValue.between(0,1).all() and z.FDR.between(0,1).all() and not z.gene_id.duplicated().any()
 rows.append({**r,'tested_genes':len(z),'FDR_below_0_05':int((z.FDR<.05).sum()),'formula':model['formula'],'independent_donors':model['n_donors']})
pd.DataFrame(rows).to_csv(W/'pseudobulk/RESULTS_SUMMARY.tsv',sep='\t',index=False)
dump(W/'pseudobulk/RESULTS_VALIDATION.json',dict(status='PASS',models=len(rows),all_result_hashes_match=True,finite_valid_P_FDR=True,unique_gene_rows=True,cross_disease_meta=False))
body=['<!doctype html><html lang="zh"><meta charset="utf-8"><title>新版样本信息与描述图表审核</title><style>body{font:16px sans-serif;max-width:1500px;margin:30px auto}img{max-width:100%}table{border-collapse:collapse;font-size:12px}td,th{padding:5px;border:1px solid #ddd}h2{margin-top:45px}</style><h1>更新后的样本信息与描述性分析</h1><p>全量1,175,115细胞，247独立样本，172名供者，26个study。全部正式标签保持不变。健康组是本研究参照组，不等于供者没有其他疾病。图表只描述组成，不作跨疾病因果或显著性推断。</p><p>全局细胞总数及每个细胞标签不变；样本／供者汇总及健康／疾病分组采用用户已确认信息。CoVarNet排除规则只应用于专门的网络范围表。</p><p>对应原稿Figure 1d、Figure 1f及Supplementary Figure 1。下图为按细胞数计算的组成，逐样本均值和供者均值另存为源数据TSV，不能混作同一分母。</p>']
for p in sorted((out/'figures').glob('*.png')):
 body+=['<h2>'+html.escape(p.stem)+'</h2><img src="data:image/png;base64,'+base64.b64encode(p.read_bytes()).decode()+'">']
body+=['<h2>L1来源集中度审核</h2>',pd.DataFrame(rare).to_html(index=False),'<h2>样本—供者—分组</h2>',m[['study','sample_uid','donor_uid','group_binary','n_cells','source_diagnosis']].to_html(index=False),'</html>']
(out/'REVIEW_CN.html').write_text(''.join(body))
dump(out/'FILES_SHA256.json',{str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file() and p.name!='FILES_SHA256.json'})
with zipfile.ZipFile(W/'updated_metadata_descriptive_and_pseudobulk_review.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in out.rglob('*'):
  if p.is_file():z.write(p,'descriptive_updated/'+str(p.relative_to(out)))
 for n in ['RESULTS_SUMMARY.tsv','RESULTS_VALIDATION.json','all_contrast_readiness.tsv','model_roster.tsv','PREPARATION.json','COMPLETE.json']:z.write(W/'pseudobulk'/n,'pseudobulk/'+n)
print('Review package and 25 native DE result validations complete')
