#!/usr/bin/env python3
import argparse,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from external_gate import load_config,approved_preflight,check_receipt,authorize
p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out');p.add_argument('--check-receipt');a=p.parse_args()
if a.check_receipt:check_receipt(a.config,a.check_receipt)
elif a.out:approved_preflight(a.config,a.out)
else:authorize(a.config)
