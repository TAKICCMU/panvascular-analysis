#!/usr/bin/env python3
"""Prepare dependency review evidence. This never grants execution approval."""
import argparse,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.review_dependencies import reviewed_dependencies
p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--output',required=True);a=p.parse_args()
import yaml
cfg=yaml.safe_load(Path(a.config).read_text());digest,manifest=reviewed_dependencies(cfg,a.config)
out=Path(a.output)
with out.open('x') as f:json.dump({'reviewed_dependency_digest':digest,'reviewed_dependency_manifest':manifest,'approved':False,'notice':'Review evidence only; subsequent user authorization required'},f,indent=2)
print(out)
