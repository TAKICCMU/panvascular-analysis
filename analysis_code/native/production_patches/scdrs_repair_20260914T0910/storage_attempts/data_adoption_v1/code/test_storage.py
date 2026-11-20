from storage_utils import *
import tempfile,sys
passed=[]
with tempfile.TemporaryDirectory(prefix='scdrs_storage_fixture_') as t:
 p=Path(t);old=p/'scores';old.mkdir();target=p/'target';target.mkdir();link=p/'alias';os.symlink(target,link)
 ino=old.stat().st_ino;exchange(old,link)
 assert old.is_symlink() and old.resolve()==target.resolve() and link.is_dir() and link.stat().st_ino==ino
 (old/'future_score').write_text('unchanged computation');assert (target/'future_score').read_text()=='unchanged computation';assert not list(link.iterdir());passed+=['atomic_directory_symlink_exchange','old_empty_inode_retained','future_output_resolves_data_alias']
 try:assert_no_score_writer([dict(target=str(old/'open.gz'),flags=os.O_WRONLY)],old)
 except ValueError:passed.append('open_score_writer_rejected')
 else:raise AssertionError('writer accepted')
 me=proc(os.getpid());assert same(me);bad=dict(me,start_ticks='wrong');assert not same(bad);passed+=['live_process_identity','PID_reuse_rejected']
 assert decode_exit(0)==0 and decode_exit(256)==1;passed.append('native_exit_wait_status_decoding')
print(json.dumps(dict(status='PASS',n_tests=len(passed),tests=passed),indent=2))
