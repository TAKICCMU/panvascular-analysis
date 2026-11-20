import run,pandas as pd,numpy as np
rng=np.random.default_rng(20260916);out=run.N/'native_test';out.mkdir(exist_ok=False)
d=pd.DataFrame({'group_binary':['healthy','disease']*40,'donor_uid':[f'd{i}' for i in range(80)],'biological_sample_uid':[f'b{i}' for i in range(80)],'donor_resolved':True,'study':['S1']*40+['S2']*40,'tissue':'artery','assay_type':'scRNA','n_parent':1000})
d['n_state']=rng.binomial(1000,rng.beta(5+3*(d.group_binary=='disease'),25));eligible,reason=run.gate(d);assert reason is None
s=run.spec(eligible,'synthetic','L1','TEST','TEST');assert 'stratum' in s['formula'];assert run.fit(eligible,s,out/'fit')==0
bad=d.copy();bad['study']=bad.group_binary;_,r=run.gate(bad);assert r=='NO_COMPARABLE_STUDY_TISSUE_ASSAY_STRATUM'
one=d[d.study=='S1'].copy();o,r=run.gate(one);assert r is None and 'stratum' not in run.spec(o,'one','L1','TEST','TEST')['formula']
rep=d.copy();rep['donor_uid']=[f'd{i//2}' for i in range(80)];r,_=run.gate(rep);assert '(1|donor_uid)' in run.spec(r,'repeat','L1','TEST','TEST')['formula']
missing=d.copy();missing['donor_uid']=None;_,r=run.gate(missing);assert r
run.dump(out/'receipt.json',{'status':'PASS','native_beta_binomial':True,'confounded_strata_blocked':True,'single_homogeneous_stratum_allowed':True,'repeated_donor_random_effect_required':True,'missing_donor_blocked':True,'seed':20260916})
