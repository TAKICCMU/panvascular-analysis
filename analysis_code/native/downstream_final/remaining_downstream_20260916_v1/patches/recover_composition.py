import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
from scipy.stats import false_discovery_control
N=Path(__file__).resolve().parents[1]
src=N/'composition/production';out=N/'composition/aggregation_patch_v2';out.mkdir(exist_ok=False)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
audit=pd.read_csv(src/'model_audit.tsv',sep='\t');assert len(audit)==816 and audit.model_id.is_unique
bound={str(src/'model_audit.tsv'):sha(src/'model_audit.tsv')};rows=[]
for i,row in audit.iterrows():
    if row.status!='COMPLETE':continue
    p=src/row.model_id/'effects.tsv';m=p.with_name('MODEL.json')
    bound[str(p)]=sha(p);bound[str(m)]=sha(m)
    d=pd.read_csv(p,sep='\t');assert len(d)==1 and d.model_id.eq(row.model_id).all()
    assert json.loads(m.read_text())['status']=='COMPLETE'
    valid=np.isfinite(d[['estimate','SE','z','P']]).all(axis=1)&d.SE.gt(0)&d.P.between(0,1)
    d['inference_status']=np.where(valid,'VALID_ESTIMATE','NOT_ESTIMABLE_NONFINITE_COEFFICIENT_UNCERTAINTY')
    if not valid.all():audit.loc[i,'status']='NOT_ESTIMABLE_NONFINITE_COEFFICIENT_UNCERTAINTY';audit.loc[i,'reason']='Native optimizer reported success but tested coefficient SE/z/P is nonfinite; no P imputation and no model refit.'
    rows.append(d)
e=pd.concat(rows,ignore_index=True);e['FDR_BH']=np.nan
for family,z in e[e.inference_status.eq('VALID_ESTIMATE')].groupby('FDR_family'):
    fdr=false_discovery_control(z.P.to_numpy());e.loc[z.index,'FDR_BH']=fdr
    # Independent explicit BH formula, retaining original family definition.
    p=z.P.to_numpy();order=np.argsort(p);s=p[order];expected=np.minimum.accumulate((s*len(s)/np.arange(1,len(s)+1))[::-1])[::-1].clip(0,1)
    assert np.allclose(fdr[order],expected)
e.to_csv(out/'effects_with_estimability_and_BH.tsv',sep='\t',index=False)
audit.to_csv(out/'model_audit_corrected.tsv',sep='\t',index=False)
e[e.inference_status.ne('VALID_ESTIMATE')].to_csv(out/'nonfinite_native_effects.tsv',sep='\t',index=False)
e.groupby('FDR_family').agg(native_rows=('P','size'),testable_rows=('P','count')).to_csv(out/'BH_family_sizes.tsv',sep='\t')
assert e.inference_status.eq('VALID_ESTIMATE').sum()==577 and e.inference_status.ne('VALID_ESTIMATE').sum()==5
fixture=np.array([.01,.03,np.nan,1.]);ok=np.isfinite(fixture);ff=np.full(4,np.nan);ff[ok]=false_discovery_control(fixture[ok]);assert np.isnan(ff[2]) and np.allclose(ff[ok],[.03,.045,1.])
(out/'METHODS_CN.md').write_text('只汇总已有原生模型，不重新拟合。5个模型虽返回optimizer成功，但SE/z/P缺失，改列不可估计并保留NA；另外52原生拟合失败和182资格不足继续列明。577有效结果按既定FDR_family进行BH，不把缺失P填0或1；无可计算P的模型不进入受检家族。保留全部原始结果和日志。\n')
(out/'RECEIPT.json').write_text(json.dumps({'sources':bound,'code_sha256':sha(Path(__file__))},indent=2))
(out/'COMPLETE.json').write_text(json.dumps({'status':'COMPLETE_AGGREGATION_WITH_EXPLICIT_UNESTIMABLE_MODELS','valid_models':577,'nonfinite_uncertainty':5,'native_failed':52,'eligibility_excluded':182,'BH_validated':True,'outputs':{str(p):sha(p) for p in out.iterdir() if p.is_file()}},indent=2))
print('COMPOSITION_RECOVERED 577 valid; 57 native non-estimable; 182 eligibility excluded',flush=True)
