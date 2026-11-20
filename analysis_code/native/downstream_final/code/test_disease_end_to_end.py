from common import *
import sys,importlib.util,pandas as pd,numpy as np
sys.path.insert(0,str(C/'src'));s=importlib.util.spec_from_file_location('native_meta',C/'scripts/external_communication_meta.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
out=W/'tests/communication_disease_native_end_to_end';out.parent.mkdir(exist_ok=True)
rng=np.random.default_rng(741);dd=[];xx=[]
for study in ['studyA','studyB']:
 for donor in range(8):
  base=rng.uniform(.1,.2)
  for group in ['control','disease']:
   sample=f'{study}_{donor}_{group}';dd.append(dict(sample_unit=sample,donor_uid=f'{study}_d{donor}',donor_verified=True,study_canonical=study,tissue_v2='artery',disease_context='synthetic_test',condition_class=group))
   xx.append(dict(sample_unit=sample,interaction_key='A||B||L||R',source='A',target='B',ligand='L',receptor='R',score=base+.15*(group=='disease')+rng.normal(0,.02),threshold_scheme='synthetic',method='LIANA'))
os.environ['G08_R_LIBRARY']=str(W/'Rlib');d=pd.DataFrame(dd);x=pd.DataFrame(xx);m.calculate(x,d,out,'/usr/bin/Rscript')
z=pd.read_csv(out/'cohort_effects_with_CR2.tsv',sep='\t');assert len(z)==2 and z.donor_model.eq('CR2_Satterthwaite').all() and z.se.gt(0).all()
meta=pd.read_csv(out/'context_DL_meta.tsv',sep='\t');assert len(meta)==1
from external_statistics import communication_effects_allow_repeated
try:communication_effects_allow_repeated(pd.concat([x,x.iloc[[0]]]),d)
except ValueError as e:assert 'Repeated edge rows' in str(e)
else:raise AssertionError('Duplicate edge not rejected')
dump(out/'PASS.json',dict(status='PASS',native_CR2_cohorts=2,context_meta_rows=1,duplicate_edge_rejected=True,source_sha256=sha(C/'scripts/external_communication_meta.py')))
print('Native disease communication complete pipeline synthetic test PASS')
