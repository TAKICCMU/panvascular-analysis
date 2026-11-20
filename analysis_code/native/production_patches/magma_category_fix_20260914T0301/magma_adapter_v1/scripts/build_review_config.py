"""Copy the current reviewed template to a new file; never rebuild from stale defaults."""
from pathlib import Path
import argparse,json
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args()
c=json.loads((root/'configs/analysis_config.yaml').read_text())
c['execution']['run_full_analysis']=False;c['execution']['mode']='design_and_test_only';c['execution']['approval_receipt']=None
with Path(a.output).open('x') as f:json.dump(c,f,ensure_ascii=False,indent=2)
print(a.output)
