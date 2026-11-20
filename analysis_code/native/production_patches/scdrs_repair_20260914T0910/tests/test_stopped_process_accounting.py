import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from production_runner_r2 import active_thread_counts
p=[dict(pid=i,state=s,threads=t) for i,s,t in [(1,'T',36),(2,'t',1),(3,'R',12),(4,'S',13),(5,'D',2)]]
n,stopped=active_thread_counts(p);assert n==27 and sum(x['threads'] for x in stopped)==37
# A stopped ancestor does not exempt a separately runnable descendant.
n,stopped=active_thread_counts([dict(pid=1,ppid=0,state='T',threads=36),dict(pid=2,ppid=1,state='S',threads=13)]);assert n==13 and len(stopped)==1
print('PASS T/t excluded; R/S/D counted; each descendant state counted independently')
