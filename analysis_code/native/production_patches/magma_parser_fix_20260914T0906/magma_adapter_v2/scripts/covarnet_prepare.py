#!/usr/bin/env python3
"""Prepare-only CLI. Full nsNMF is available only via the root approval gate."""
import argparse,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import pandas as pd
from src.covarnet import prepare,write_prepared,sha
p=argparse.ArgumentParser();p.add_argument('--counts',required=True);p.add_argument('--cohort',required=True);p.add_argument('--ontology',required=True);p.add_argument('--out',required=True);p.add_argument('--approve-named-state-frequency-rule',action='store_true');a=p.parse_args()
result=prepare(pd.read_csv(a.counts,sep='\t'),pd.read_csv(a.cohort,sep='\t',keep_default_na=False),pd.read_csv(a.ontology,sep='\t',keep_default_na=False),{'approve_named_state_frequency_rule':a.approve_named_state_frequency_rule})
print(json.dumps(write_prepared(result,a.out,{p:sha(p) for p in [a.counts,a.cohort,a.ontology]}),indent=2))
