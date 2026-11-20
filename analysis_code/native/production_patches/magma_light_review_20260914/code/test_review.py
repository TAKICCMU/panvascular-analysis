import copy,json
from pathlib import Path
from review_core import bh,validate_family
W=Path(__file__).resolve().parents[1];family='synthetic_L1_cell_weighted';passed=[]
rows=[];mapping=[];variables=[];native=[]
for i,p in enumerate([.01,.2,.3]):
 k=f'S{i:04}';model=str(i+1);group=f'Test{i}';var={'variable':k,'group':group};variables.append(var);mapping.append({'VARIABLE':k,'MODEL':model,'group':group,'condition_variable':'Average'})
 x={'VARIABLE':k,'TYPE':'COVAR','MODEL':model,'NGENES':'100','BETA':'0.1','BETA_STD':'0.01','SE':'0.2','P':str(p)};native.append(dict(x,VARIABLE='Average'));native.append(x.copy());rows.append(dict(x,**var,FDR_BH=str(bh([.01,.2,.3])[i]),BH_family=family,n_tests_in_family='3'))
assert validate_family(rows,mapping,variables,native,family)['n_tested']==3;passed.append('valid_family_and_BH')
def fail(name,fn):
 try:fn()
 except (AssertionError,ValueError,KeyError):passed.append(name);return
 raise AssertionError(name)
r=copy.deepcopy(rows);r[1]['VARIABLE']=r[0]['VARIABLE'];fail('duplicate_variable',lambda:validate_family(r,mapping,variables,native,family))
r=copy.deepcopy(rows);r[0]['FDR_BH']='0.001';fail('changed_FDR',lambda:validate_family(r,mapping,variables,native,family))
r=copy.deepcopy(rows);r[0]['P']='nan';fail('nonfinite_P',lambda:validate_family(r,mapping,variables,native,family))
r=copy.deepcopy(rows);r[0]['group']='wrong';fail('mapping_mismatch',lambda:validate_family(r,mapping,variables,native,family))
r=copy.deepcopy(rows);r[0]['BH_family']='other';fail('family_changed',lambda:validate_family(r,mapping,variables,native,family))
r=copy.deepcopy(rows);r[0]['BETA']='0.2';fail('native_value_changed',lambda:validate_family(r,mapping,variables,native,family))
fail('missing_target',lambda:validate_family(rows[:-1],mapping,variables,native,family))
report={'status':'PASS','synthetic_tests':passed,'native_scientific_calls':0};(W/'TEST_RESULTS.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))
