import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from production_runner import remaining_disk
G=1024**3
base={'SC':dict(allocation_gib=32,baseline_bytes=0,current_bytes=0),'GPU':dict(allocation_gib=16,baseline_bytes=2*G,current_bytes=2*G),'MAGMA':dict(allocation_gib=4,baseline_bytes=G,current_bytes=G)}
a=remaining_disk(74,base);assert a['disk_pass'] and a['required_free_gib']==72
base['GPU']['current_bytes']+=8*G;a=remaining_disk(66,base);assert a['disk_pass'] and a['required_free_gib']==64
assert not remaining_disk(63,base)['disk_pass'] # unknown unrelated consumption does not get credit
base['SC']['current_bytes']=33*G;a=remaining_disk(41,base);assert not a['disk_pass'] and a['over_allocation']==['SC']
base['SC']['current_bytes']=0;base['GPU']['current_bytes']=G;assert remaining_disk(74,base)['required_free_gib']==72 # no credit for pre-baseline outputs
print('PASS 5 remaining-disk accounting cases')
