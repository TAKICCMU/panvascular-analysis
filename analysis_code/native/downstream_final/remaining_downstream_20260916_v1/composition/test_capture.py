import run,pandas as pd,numpy as np
rng=np.random.default_rng(19)
d=pd.DataFrame({'group_binary':['healthy','disease']*40,'sample_uid':[f's{i}' for i in range(80)],'analysis_profile_uid':[f's{i}' for i in range(80)],'donor_uid':[f'd{i}' for i in range(80)],'biological_sample_uid':[f'b{i}' for i in range(80)],'replicate_relationship':'independent_sample_independent_donor','donor_resolved':True,'study':['S1']*40+['S2']*40,'tissue':'artery','assay_type':'scRNA','n_parent':1000})
d['n_state']=rng.binomial(1000,rng.beta(5+3*(d.group_binary=='disease'),25))
for i,name in [(0,'GSE210248__GSM6424655_rep1'),(2,'GSE210248__GSM6424655_rep2')]:
 d.loc[i,['sample_uid','analysis_profile_uid','donor_uid','biological_sample_uid','replicate_relationship']]=[name,name,'capture_donor','capture_bio','independent_capture_same_donor']
e,r=run.gate(d);assert r is None and e.donor_uid.nunique()==79 and len(e)==80
s=run.spec(e,'capture','L1','TEST','TEST');assert '(1|donor_uid)' in s['formula']
out=run.N/'native_capture_test';out.mkdir(exist_ok=False);assert run.fit(e,s,out/'fit')==0
bad=d.copy();bad.loc[0,'sample_uid']='unknown';_,r=run.gate(bad);assert r=='UNKNOWN_DUPLICATE_BIOLOGICAL_UNIT'
run.dump(out/'receipt.json',dict(status='PASS',observations=80,independent_donors=79,native_repeated_capture_fit=True,unknown_duplicate_blocked=True))
