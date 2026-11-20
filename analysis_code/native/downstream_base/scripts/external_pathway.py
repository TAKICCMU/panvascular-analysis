#!/usr/bin/env python3
"""Common-authorized launcher for inherited fgsea/ORA native R adapter."""
import argparse,sys,subprocess
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from external_gate import bootstrap_cli,approved_preflight
bootstrap_cli()
p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--ranks',required=True);p.add_argument('--genesets',required=True);p.add_argument('--out',required=True);p.add_argument('--rscript',required=True);a=p.parse_args()
approved_preflight(a.config,Path(a.out).with_name(Path(a.out).name+'.authorization'))
subprocess.run([a.rscript,str(Path(__file__).with_suffix('.R')),a.ranks,a.genesets,a.out],check=True)
