"""One-time RAM-preserving score-directory exchange; requires approved sealed receipt."""
from pathlib import Path
import os,json,time,signal,fcntl
from storage_utils import *
W=Path(__file__).resolve().parents[1]
def main():
 c=json.loads((W/'config.json').read_text());receipt=json.loads((W/'RECEIPT.json').read_text())
 for x in receipt['bound_files']:
  if sha(x['path'])!=x['sha256']:raise ValueError('Changed receipt binding')
 if receipt['scope']!='RAM_preserving_scDRS_storage_adoption':raise ValueError('Wrong receipt')
 if (W/'MIGRATION_COMPLETE.json').exists():raise FileExistsError('Already migrated')
 budget=mount_budget(c)
 if budget['reasons']:raise ValueError('Mount gate failed '+repr(budget))
 controller=c['old_controller'];native=c['native'];guard=c['old_guard']
 if not all(same(x) for x in [controller,native,guard]):raise ValueError('Process identity changed')
 if proc(controller['pid'])['state'] in ['T','t','Z'] or proc(native['pid'])['state'] in ['T','t','Z']:raise ValueError('Preexisting pause/death requires review')
 gs=json.loads(Path(c['old_guard_status']).read_text())
 if gs['status']!='WATCHING_CURRENT_SCORE':raise ValueError('Old guard owns a pause; refuse migration')
 scores=Path(c['old_scores']);target=Path(c['data_scores']);backup=scores.with_name('scores.home_empty_before_data');alias=scores.with_name('scores.data_exchange_pending')
 if scores.is_symlink() or not scores.is_dir() or list(scores.iterdir()) or backup.exists() or alias.exists() or target.exists():raise ValueError('Expected untouched empty old scores and absent targets')
 before=fds(native);assert_no_score_writer(before,scores)
 record=dict(time=time.time(),controller=controller,native=native,old_guard=guard,old_guard_status=gs,native_fds_before=before,budget=budget,old_scores=str(scores),data_scores=str(target),backup=str(backup),old_scores_stat=dict(inode=scores.stat().st_ino,device=scores.stat().st_dev))
 write(W/'PRE_MIGRATION.json',record)
 # Retire only WATCHING guard; it owns no resource pause. Preserve its code/status/log.
 send(guard,signal.SIGTERM)
 send(controller,signal.SIGSTOP);wait_state(controller,{'T','t'})
 target.mkdir(parents=True,exist_ok=False);os.symlink(str(target),alias)
 paused=False;switched=False;t0=time.monotonic()
 try:
  send(native,signal.SIGSTOP);paused=True;wait_state(native,{'T','t'})
  frozen=fds(native);assert_no_score_writer(frozen,scores)
  if list(scores.iterdir()):raise ValueError('Scores appeared before stop')
  exchange(scores,alias);switched=True;alias.rename(backup)
  if not scores.is_symlink() or scores.resolve()!=target.resolve() or list(backup.iterdir()):raise ValueError('Bad alias exchange')
  after=fds(native)
  if [(x['fd'],x['target']) for x in before if x['fd'] in [1,2]]!=[(x['fd'],x['target']) for x in after if x['fd'] in [1,2]]:raise ValueError('Log descriptors changed')
 finally:
  if paused and same(native):send(native,signal.SIGCONT)
  write(W/'MIGRATION_WINDOW.json',dict(paused_seconds=time.monotonic()-t0,alias_switched=switched,native=proc(native['pid']),controller=proc(controller['pid'])))
 record.update(status='COMPLETE_ALIAS_SWITCH_NATIVE_CONTINUES',alias_target=os.readlink(scores),native_fds_at_stop=frozen,native_fds_after=after,paused_seconds=time.monotonic()-t0,receipt_sha256=sha(W/'RECEIPT.json'),sealed_old_code_changed=False,raw_or_derived_moved=False)
 write(W/'MIGRATION_COMPLETE.json',record);print(json.dumps(record))
if __name__=='__main__':main()
