#!/usr/bin/env python3
"""Adapt inherited complete-cis-region coloc.abf to new target/QTL manifest.
Default render-only. No network fetch, no automatic target list, no binary-GWAS coercion.
"""
import argparse,json,re,sys,subprocess,difflib
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from external_gate import bootstrap_cli,approved_preflight
if __name__=='__main__':bootstrap_cli()
import pandas as pd,yaml
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from external_modules import sha256
p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);p.add_argument('--execute',action='store_true');a=p.parse_args()
c=yaml.safe_load(Path(a.config).read_text());cfg=c['external']['coloc']
if a.execute:approved_preflight(a.config,Path(a.out).with_name(Path(a.out).name+'.authorization'))
if not cfg.get('qtl_manifest') or not cfg.get('gwas_scenarios'):raise RuntimeError('BLOCKED_MISSING_QTL_GWAS_RESOURCE: provide audited complete-region inputs')
out=Path(a.out);out.mkdir(parents=True,exist_ok=False)
source=Path(__file__).resolve().parents[1]/'external_legacy/colocalisation/run_full_coloc.R';s=source.read_text();before=s
if cfg.get('gwas_effect_type')!='standardized_quantitative_latent_G':raise ValueError('Inherited sdY=1 assumes standardized latent G; disease case/control GWAS needs separately reviewed model')
def rq(x):return json.dumps(str(x),ensure_ascii=False)
s=s.replace('ROOT <- "/path/to/panvas_work/Tier1_full_colocalisation_v1"','ROOT <- '+rq(out/'analysis'))
s=s.replace('MANIFEST <- file.path(ROOT, "00_metadata", "qtl_download_manifest.tsv")','MANIFEST <- '+rq(cfg['qtl_manifest']))
block='GWAS <- c(\n'+',\n'.join('  '+k+' = '+rq(v['path']) for k,v in cfg['gwas_scenarios'].items())+'\n)'
s,n=re.subn(r'GWAS <- c\(.*?\n\)',lambda m:block,s,count=1,flags=re.S)
if n!=1 or '/path/to/' in s:raise ValueError('Inherited source mismatch/unresolved resources')
# Keep full original priors, prior sensitivity grid, SNP orientation, minimum100, and reports.
s='source('+rq(Path(__file__).with_name('external_native_gate.R').resolve())+')\n'+s
(out/'coloc_complete_cis.R').write_text(s)
(out/'source.diff').write_text(''.join(difflib.unified_diff(before.splitlines(True),s.splitlines(True))))
(out/'PROVENANCE.json').write_text(json.dumps({'source_sha256':sha256(source),'rendered_sha256':sha256(out/'coloc_complete_cis.R'),'executed':False},indent=2))
if a.execute:

 q=pd.read_csv(cfg['qtl_manifest'],sep='\t')
 for col in ['gene','context','source','sample_size','output','n_qtl_rows','nominated_lineage']:
  if col not in q:raise ValueError('QTL manifest lacks '+col)
 if not cfg.get('genome_build_and_alleles_reviewed'):raise ValueError('Build and allele audit required')
 for g in cfg['gwas_scenarios'].values():
  if sha256(g['path'])!=g['sha256']:raise ValueError('GWAS source hash mismatch')
 if sha256(cfg['qtl_manifest'])!=cfg['qtl_manifest_sha256']:raise ValueError('QTL manifest hash mismatch')
 subprocess.run([cfg['rscript'],str(out/'coloc_complete_cis.R')],check=True)
