import tempfile,json,copy,ast
from pathlib import Path
import resume_scdrs as r
passed=[]
def reject(name,fn):
 try:fn()
 except (ValueError,FileExistsError):passed.append(name);return
 raise AssertionError('Did not reject '+name)
with tempfile.TemporaryDirectory() as t:
 p=Path(t);f=p/'input';f.write_text('fixture');s=f.stat()
 cert=dict(path=str(f),sha256=r.sha(f),full_sha256_verified=True,stat=dict(size=s.st_size,mtime_ns=s.st_mtime_ns,inode=s.st_ino,device=s.st_dev,nlink=s.st_nlink))
 assert r.check_stat(cert,f,r.sha(f));passed.append('new_fullhash_stat_certificate_pass')
 for key in ['size','mtime_ns','inode','device','nlink']:
  bad=copy.deepcopy(cert);bad['stat'][key]+=1;reject('reject_changed_'+key,lambda:r.check_stat(bad,f,r.sha(f)))
 reject('reject_wrong_SHA',lambda:r.check_stat(cert,f,'wrong'))
 bad=copy.deepcopy(cert);bad['observed_sha256']='wrong';reject('reject_observed_SHA_mismatch',lambda:r.check_stat(bad,f,r.sha(f)))
 if Path('/proc/sys/kernel/random/boot_id').exists():
  bad=copy.deepcopy(cert);bad['boot_id']='old-boot';reject('reject_previous_boot',lambda:r.check_stat(bad,f,r.sha(f)))
 reject('reject_wrong_path',lambda:r.check_stat(cert,p/'other',r.sha(f)))
 bad=copy.deepcopy(cert);bad['full_sha256_verified']=False;reject('reject_unverified_hash',lambda:r.check_stat(bad,f,r.sha(f)))
 old=p/'old.json';sc=[dict(name='s',geneset='g',sha256='h')];old.write_text(json.dumps(dict(scenarios=sc,scdrs_executable='/approved/scdrs')))
 r.W=p;c=dict(original_config=str(old),scenarios=sc,native_threads=12,n_ctrl=1000,graph_policy='frozen_global_n30',scdrs_executable='/approved/scdrs',rebuild_derived=False,group_keys=['level1_annotation_scdrs_test','level2_annotation_scdrs_test'],output=str(p/'outputs'))
 r.validate_parameters(c);passed.append('same_scientific_parameters_pass')
 for key,value in [('n_ctrl',999),('native_threads',13),('graph_policy','rebuild'),('rebuild_derived',True),('group_keys',['oldlabels']),('scenarios',[])]:
  b=copy.deepcopy(c);b[key]=value;reject('reject_changed_'+key,lambda:r.validate_parameters(b))
 (p/'STARTED.json').write_text('{}');reject('reject_duplicate_started_attempt',lambda:r.validate_parameters(c));(p/'STARTED.json').unlink()
 (p/'outputs').mkdir();reject('reject_existing_output',lambda:r.validate_parameters(c))
 q=dict(communication_root='/comm',gpu_root='/gpu')
 a=r.process_budget(q,[],set());assert a['projected']==28 and a['pass_gate'];passed.append('future_comm10_GPU4_reserved')
 a=r.process_budget(q,[dict(pid=1,state='S',threads=10,cmd='/comm/worker'),dict(pid=2,state='R',threads=4,cmd='/gpu/worker')],set());assert a['projected']==28;passed.append('live_comm_GPU_not_double_counted')
 a=r.process_budget(q,[dict(pid=1,state='T',threads=36,cmd='/old'),dict(pid=2,state='R',threads=20,cmd='/other')],set());assert a['projected']==48 and not a['pass_gate'];passed.append('stopped_excluded_live_budget_rejected')
 source=Path(r.__file__).read_text();tree=ast.parse(source)
 assert not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='build_derived' for n in ast.walk(tree));passed.append('no_build_derived_call')
print(json.dumps(dict(status='PASS',n_tests=len(passed),tests=passed),indent=2))
