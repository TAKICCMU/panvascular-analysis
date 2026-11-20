from pathlib import Path
import ast,csv,hashlib,json,re,sys,collections,subprocess,shutil
import pandas as pd,numpy as np
B=Path(__file__).resolve().parents[2];A=B/'audit';S=B/'evidence/source_data';T=B/'evidence/supplementary_tables'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def write(name,rows):pd.DataFrame(rows).to_csv(A/name,sep='\t',index=False)
def read(p,**kw):return pd.read_csv(p,sep='\t',**kw)
checks=[]
def check(name,passed,detail):checks.append(dict(check=name,status='PASS' if passed else 'REVIEW_REQUIRED',detail=str(detail)))
def norm(s):return re.sub('[^a-z0-9]','',str(s).lower())
# original-machine column spellings preserved in exported TSV, display headers mapped.
colmap=read(T/'column_name_map.tsv',keep_default_na=False);struct=json.loads((T/'table_structure.json').read_text())
def table(n):
 sh=f'Table {n:02d}'; st=next(v for v in struct if v['sheet']==sh)
 d=read(T/(sh+'.tsv.gz'),header=2,nrows=st['rows']);m=colmap[colmap.Worksheet==sh]
 check('table_header_'+sh,list(d.columns)==list(m['Original column']),'Preserved machine headers match mapping; workbook uses mapped display headers')
 d.columns=m['Original column'].tolist();return d.rename(columns={c:norm(c) for c in d.columns})
def cmp(a,b,keys,fields,name):
 a=a.rename(columns={c:norm(c) for c in a});b=b.rename(columns={c:norm(c) for c in b});ks=list(map(norm,keys));fs=list(map(norm,fields))
 try:
  m=a[ks+fs].merge(b[ks+fs],on=ks,suffixes=('_a','_b'),validate='one_to_one',how='outer',indicator=True)
  allok=m['_merge'].eq('both').all();errors={}
  for f in fs:
   x=pd.to_numeric(m[f+'_a'],errors='coerce');y=pd.to_numeric(m[f+'_b'],errors='coerce')
   if (x.notna()|y.notna()).any():good=np.isclose(x,y,rtol=1e-11,atol=1e-12,equal_nan=True)
   else:good=m[f+'_a'].fillna('<NA>').astype(str).eq(m[f+'_b'].fillna('<NA>').astype(str))
   errors[f]=int((~good).sum());allok=allok and bool(np.all(good))
  check(name,allok,{'rows':len(m),'unmatched':int(m['_merge'].ne('both').sum()),'mismatches':errors})
 except Exception as e:check(name,False,type(e).__name__+': '+str(e))
# All 64 exported table headers and row counts
for n in range(1,65):
 d=table(n);st=next(x for x in struct if x['sheet']==f'Table {n:02d}')
 check(f'table_rows_{n:02}',len(d)==st['rows'],len(d))
# Scientific figure/table equivalence
cmp(read(S/'Figure3/Fig3d_e_all29.tsv'),table(13),['gene_symbol','cell_state'],['meta_logFC','meta_SE','meta_P','meta_CI_low','meta_CI_high','target_set_meta_FDR','direction_consistency'],'Figure3de_Table13')
cmp(read(S/'Figure5_final_companion/P86_Table32_all87_original_entries.tsv'),table(32),['gene','cell_context','term'],['beta','se','p','CI_low','CI_high','FDR_BH_all_29_pairs_x_3_terms'],'Figure5_all87_Table32')
f5=read(S/'Figure5_final_companion/Fig5c_all29.tsv');meta=read(S/'Figure5_final_companion/all29_original_meta.tsv');sub=meta[(meta.posterior=='q05')&(meta.sensitivity=='all_sections')&(meta.term=='cell_abundance')]
cmp(f5,sub,['gene','cell_context','term'],['beta','CI_low','CI_high','FDR_BH_all_29_pairs_x_3_terms'],'Figure5c_original_meta29')
u=read(S/'Figure5_final_companion/Figure5b_60_gene_unit_estimates.tsv'); fullu=read(S/'Figure5_final_companion/all29_original_unit_models.tsv');keys=['gene','biological_unit_id','posterior','term','sensitivity']
subu=fullu.merge(u[keys],on=keys,validate='many_to_one')
cmp(u,subu,keys,['beta','se','n_sections'],'Figure5b_original_unit_estimates')
check('Figure5b_count',len(u)==60 and u.biological_unit_id.nunique()==34 and u.groupby('gene').size().to_dict()=={'FN1':34,'SCARB1':26},{'estimates':len(u),'recorded_IDs':u.biological_unit_id.nunique(),'by_gene':u.groupby('gene').size().to_dict()})
# Stored marker display transform only, no raw expression recomputation
mroot=S/'Current_Display_20260920/marker_audit_20260921'
m=read(mroot/'original_mean_log1p.tsv',index_col=0);d=read(mroot/'original_detection_fraction.tsv',index_col=0);z=read(mroot/'per_gene_scaled_mean.tsv',index_col=0)
check('Figure1g_geometry',m.shape==d.shape==z.shape==(14,28) and m.index.equals(z.index) and m.columns.equals(z.columns),'14 x 28 stored matrix, row/column identity')
w=(m-m.min())/(m.max()-m.min());check('Figure1g_scaled_storage',np.allclose(w,z,atol=1e-12,equal_nan=True),'Stored B vs min-max of stored A; not new expression computation')
check('Figure1g_detection_range',bool(((d>=0)&(d<=1)).all().all()),'Fractions 0–1')
ma=read(mroot/'all392_full_audit.tsv');check('Figure1g_prior_raw_recheck',ma.mean_abs_error.max()<1e-6 and ma.detection_abs_error.max()<1e-12,{'max_mean_error':ma.mean_abs_error.max(),'max_detection_error':ma.detection_abs_error.max(),'scope':'Earlier raw recheck record read only; not repeated now','coverage_limit':ma.measurement_coverage_status.unique().tolist()})
# GSEA displays vs original outputs
mg=read(S/'Current_Display_20260920/Figure3c_current/MAGMA_display9.tsv');allmg=read(S/'Current_Display_20260920/Figure3c_current/MAGMA_complete_GSEA.tsv')
sel=allmg.merge(mg[['pathway','database']],on=['pathway','database'],validate='one_to_one')
cmp(mg,sel,['pathway','database'],['pval','padj','NES','size','n_ranked_genes'],'Figure3c_MAGMA_original_GSEA')
check('MAGMA_display_scope',len(mg)==9 and mg.database.nunique()==3 and mg.n_ranked_genes.eq(17885).all(),mg.groupby('database').size().to_dict())
pb=read(S/'Current_Display_20260920/Figure3c_current/pseudobulk_display9_all7_lineages.tsv');pbparts=[]
for aid,g in pb.groupby('analysis_id'):
 f=B/'evidence/table20_full_results'/aid/'fgsea_all.tsv.gz';raw=read(f);sel=raw.merge(g[['pathway','database']],on=['pathway','database'],validate='one_to_one');sel['analysis_id']=aid;pbparts.append(sel)
cmp(pb,pd.concat(pbparts),['analysis_id','pathway','database'],['NES','pval','padj','size','n_ranked_genes'],'Figure3c_pseudobulk_full_DL_GSEA')
check('pseudobulk_display_scope',len(pb)==63 and pb.cell_state.nunique()==7 and pb.pathway.nunique()==9 and pb.analysis_id.str.startswith('PRIMARY_META_DL_L1').all(),'9 pathways x 7 L1 groups; signed DL rank')
# Atlas coordinate and named-state identity
coords=read(S/'Figure1/full_cells_frozen_umap.tsv.gz'); manifest=json.loads((B/'analysis_code/native/atlas_finalization/PREPARED_RELEASE_MANIFEST.json').read_text())
(A/'coordinate_columns.json').write_text(json.dumps({'columns':coords.columns.tolist(),'rows':len(coords)},indent=2))
check('atlas_n_cells',len(coords)==1175115,len(coords))
check('named_state_table29',len(table(29))==116,len(table(29)))
check('full_recurrence_groups',len(read(S/'Figure1/L2_cross_study_recurrence.tsv'))==128,len(read(S/'Figure1/L2_cross_study_recurrence.tsv')))
# Each adopted master and export code copy is bound to real final figure hashes
versions=read(B/'evidence/图件采用版本.tsv');figrows=[]
import fitz
for _,r in versions.iterrows():
 pdf=B/'evidence/current_figures'/(r.figure+'.pdf');svg=B/'figure_code/frozen_svg'/(r.figure+'.svg')
 with fitz.open(pdf) as doc:w,h=doc[0].rect.width*25.4/72,doc[0].rect.height*25.4/72
 ok=sha(pdf)==r.package_pdf_sha256 and sha(svg)==r.package_svg_sha256
 figrows.append(dict(figure=r.figure,pdf='evidence/current_figures/'+pdf.name,svg=str(svg.relative_to(B)),pdf_sha256=sha(pdf),svg_sha256=sha(svg),width_mm=w,height_mm=h,adopted_version=r.version,hash_match=ok))
 check('adopted_'+r.figure,ok,'PDF/SVG vs adopted manifest; this is identity check, not new visual approval')
write('FIGURE_MASTER_BINDINGS.tsv',figrows)
# Code inventory, non-executing static syntax, duplicate groups and literal input/output paths
files=[];syntax=[];paths=[];imports=[];secrets=[]
for root in [B/'analysis_code',B/'figure_code',B/'recovered_sources']:
 for f in sorted(root.rglob('*')):
  if not f.is_file():continue
  rel=str(f.relative_to(B));digest=sha(f);files.append(dict(path=rel,bytes=f.stat().st_size,sha256=digest,suffix=f.suffix))
  if f.suffix not in ['.py','.R','.r','.sh']:continue
  text=f.read_text(errors='replace')
  if f.suffix=='.py':
   try:
    tree=ast.parse(text);status='PASS_PARSE';detail='AST only; not imported'
    for node in ast.walk(tree):
     if isinstance(node,(ast.Import,ast.ImportFrom)):
      mods=[a.name for a in node.names] if isinstance(node,ast.Import) else [node.module or '']
      for mod in mods:imports.append(dict(path=rel,module=mod))
     if isinstance(node,ast.Constant) and isinstance(node.value,str) and node.value.startswith(('/Users/','/home/','/data/','/usr/','/tmp/')) and '\n' not in node.value and len(node.value)<1000:
      paths.append(dict(code=rel,line=getattr(node,'lineno',None),literal_path=node.value,local_exists=Path(node.value).exists(),role='literal may be input/output/history; not inferred',host_scope='server' if node.value.startswith(('/home/','/data/')) else 'local'))
   except SyntaxError as e:status='SYNTAX_REVIEW';detail=f'{e.msg}; line {e.lineno}'
   syntax.append(dict(path=rel,language='Python',status=status,detail=detail,sha256=digest))
  for i,line in enumerate(text.splitlines(),1):
   if re.search(r'(?:password|api_key|access_token|secret_key)\s*=\s*[\'\"][^\'\"]{16,}[\'\"]',line,re.I):secrets.append({'file':rel,'line':i,'status':'REVIEW_LITERAL_NO_VALUE_EXPOSED'})
write('CODE_FILE_INVENTORY.tsv',files);write('PYTHON_SYNTAX.tsv',syntax);write('ABSOLUTE_PATH_REFERENCES.tsv',paths);write('IMPORTED_MODULES.tsv',imports)
write('CREDENTIAL_PATTERN_SCAN.tsv',secrets or [{'file':'NONE','line':'','status':'NO_LITERAL_PATTERN_HIT; heuristic not security certification'}])
# Confirm known packaging portability defect without executing historical launchers
resolver=B/'figure_code/native_plot_code/current/resolve_paths.py';expected=resolver.parents[1]/'merge_records/PATH_MAP.json'
check('packaged_native_resolver_dependency',expected.exists(),str(expected.relative_to(B))+' missing in inherited package; original manifest recovered separately, not silently patched')
# Adopted parameters and historical notebook execution indices
params=json.loads((B/'analysis_code/native/upstream/CPU_full_atlas_theta2_20260911_160002_388642/pipeline_parameters.json').read_text());(A/'ADOPTED_CPU_THETA2_PARAMETERS.json').write_text(json.dumps(params,ensure_ascii=False,indent=2))
notebook=B/'analysis_code/native/upstream/CPU_theta2_author_executed_notebook.py';lines=notebook.read_text().splitlines();cells=[]
for i,l in enumerate(lines,1):
 if re.match(r'#.*(?:[Cc]ell|execution_count)',l):cells.append({'line':i,'marker':l})
write('UPSTREAM_NOTEBOOK_CELL_MARKERS.tsv',cells)
write('RESULT_CONSISTENCY_CHECKS.tsv',checks)
(A/'AUDIT_COUNTS.json').write_text(json.dumps({'files':len(files),'code_python':len(syntax),'python_parse_ok':sum(r['status']=='PASS_PARSE' for r in syntax),'checks':len(checks),'check_pass':sum(c['status']=='PASS' for c in checks),'review_required':[c for c in checks if c['status']!='PASS']},ensure_ascii=False,indent=2))
print((A/'AUDIT_COUNTS.json').read_text())
