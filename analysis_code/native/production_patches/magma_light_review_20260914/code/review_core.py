import csv,math
from pathlib import Path

def read_tsv(path):
 with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))
def native_rows(path):
 lines=[s for s in Path(path).read_text().splitlines() if s.strip() and not s.startswith('#')]
 header=lines[0].split();return [dict(zip(header,s.split())) for s in lines[1:]]
def bh(p):
 order=sorted(range(len(p)),key=lambda i:p[i]);q=[0.]*len(p);running=1.
 for rank in range(len(p),0,-1):
  i=order[rank-1];running=min(running,p[i]*len(p)/rank);q[i]=running
 return q

def validate_family(rows,mapping,variables,native,family):
 n=len(rows);assert n>0
 vm={x['variable']:x['group'] for x in variables};assert len(vm)==len(variables) and len(set(vm.values()))==len(vm)
 keys=[x['VARIABLE'] for x in rows];assert len(set(keys))==n and set(keys)==set(vm)
 mm={x['VARIABLE']:x for x in mapping};assert len(mm)==len(mapping)==n and set(mm)==set(vm)
 assert all(x['condition_variable']=='Average' for x in mapping)
 native_targets=[x for x in native if x['VARIABLE']!='Average'];nn={x['VARIABLE']:x for x in native_targets}
 assert len(nn)==len(native_targets)==n and set(nn)==set(vm)
 nuisance=[x for x in native if x['VARIABLE']=='Average'];assert len(nuisance)==n
 for row in rows:
  k=row['VARIABLE'];assert k!='Average' and row['variable']==k and row['group']==vm[k]==mm[k]['group']
  assert row['MODEL']==mm[k]['MODEL']==nn[k]['MODEL'];assert row['BH_family']==family and int(row['n_tests_in_family'])==n
  assert row['TYPE']==nn[k]['TYPE']=='COVAR'
  assert sum(x['MODEL']==row['MODEL'] for x in nuisance)==1
  for field in ['NGENES','BETA','BETA_STD','SE','P']:
   assert math.isfinite(float(row[field])) and float(row[field])==float(nn[k][field]),(family,k,field)
  assert int(row['NGENES'])>0 and float(row['SE'])>0
  assert 0<=float(row['P'])<=1 and 0<=float(row['FDR_BH'])<=1
 p=[float(x['P']) for x in rows];q=bh(p);err=max(abs(z-float(x['FDR_BH'])) for z,x in zip(q,rows));assert err<=1e-12,(family,err)
 return {'n_tested':n,'n_FDR_lt_0_05':sum(float(x['FDR_BH'])<.05 for x in rows),'n_P_lt_0_05':sum(float(x['P'])<.05 for x in rows),'min_P':min(p),'min_original_FDR':min(float(x['FDR_BH']) for x in rows),'BH_check_max_abs_error':err}
