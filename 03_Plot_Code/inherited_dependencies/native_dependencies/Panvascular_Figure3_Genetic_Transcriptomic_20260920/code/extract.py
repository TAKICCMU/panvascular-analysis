from common import *
import zipfile,shutil
base=Path('/Users/taki/Desktop/投稿包_最新版_20260919');z=zipfile.ZipFile(base/'Supplementary_Tables.xlsx');ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'};ss=[]
if 'xl/sharedStrings.xml' in z.namelist():ss=[''.join(t.text or '' for t in v.iter('{'+ns['s']+'}t')) for v in ET.fromstring(z.read('xl/sharedStrings.xml'))]
records=[]
for row in ET.fromstring(z.read('xl/worksheets/sheet21.xml')).findall('.//s:row',ns)[3:]:
 d={}
 for c in row.findall('s:c',ns):
  v=c.find('s:v',ns);s=v.text if v is not None else ''.join(t.text or '' for t in c.findall('.//s:t',ns));s=ss[int(s)] if c.get('t')=='s' else s;d[re.sub('[0-9]','',c.get('r'))]=s
 records.append(d)
m=pd.DataFrame(records);m=m[m.M.eq('preranked_GSEA')].rename(columns={'A':'pathway','B':'pval','C':'padj','F':'NES','G':'size','H':'database','L':'n_ranked_genes'});m=m[['pathway','pval','padj','NES','size','database','n_ranked_genes']]
for c in ['pval','padj','NES','size','n_ranked_genes']:m[c]=pd.to_numeric(m[c],errors='coerce')
assert not m.duplicated(['database','pathway']).any();m.to_csv(O/'source_data/MAGMA_complete_GSEA.tsv',sep='\t',index=False)
dbs=['GO_BP','KEGG_LEGACY','REACTOME'];sel=[]
for db in dbs:sel.append(m[m.database.eq(db)&m.padj.lt(.05)].sort_values(['padj','pathway']).head(3))
m=pd.concat(sel);m.to_csv(O/'source_data/MAGMA_display9.tsv',sep='\t',index=False)
sumfile=Path('/Users/taki/Desktop/Codex任务1/pseudobulk_manuscript_strategy_enrichment_v1/completed_review_v1/analysis_summary.tsv');summary=pd.read_csv(sumfile,sep='\t');primary=summary[summary.role.eq('primary_original_pathway_DL')&summary.level.eq('L1')];assert len(primary)==7
ps=[];provenance=[]
for _,row in primary.iterrows():
 path=base/'Supplementary_Table20_full_results'/row.analysis_id/'fgsea_all.tsv.gz';shutil.copy2(path,O/'source_data'/f'{row.analysis_id}_fgsea.tsv.gz');d=pd.read_csv(path,sep='\t');assert d.analysis_id.eq(row.analysis_id).all();d['cell_state']=row.cell_state;ps.append(d);provenance.append({'analysis_id':row.analysis_id,'state':row.cell_state,'source':str(path),'sha256':sha(path)})
p=pd.concat(ps);assert not p.duplicated(['cell_state','database','pathway']).any()
# Display selection only: recurrent significant pathways first, then median absolute NES and pathway ID.
# No new hypothesis tests, FDR adjustment, gene filtering or state selection.
p['significant']=p.padj.lt(.05);p['abs_NES']=p.NES.abs();rank=p.groupby(['database','pathway']).agg(significant_states=('significant','sum'),evaluable_states=('NES','count'),median_abs_NES=('abs_NES','median')).reset_index();chosen=[]
for db in dbs:chosen.append(rank[rank.database.eq(db)&rank.significant_states.gt(0)].sort_values(['significant_states','median_abs_NES','pathway'],ascending=[False,False,True]).head(3))
rank=pd.concat(chosen);rank.to_csv(O/'source_data/pseudobulk_display9_selection.tsv',sep='\t',index=False)
states=[s for s in P.L1 if s in set(primary.cell_state)];keys=rank[['database','pathway']].merge(pd.DataFrame({'cell_state':states}),how='cross');display=keys.merge(p,on=['database','pathway','cell_state'],how='left',validate='one_to_one');display.to_csv(O/'source_data/pseudobulk_display9_all7_lineages.tsv',sep='\t',index=False)
pd.DataFrame(provenance).to_csv(O/'source_data/pseudobulk_input_manifest.tsv',sep='\t',index=False)
print('MAGMA',m[['pathway','NES','padj']].to_string(index=False));print('PB',rank.to_string(index=False));print(states)
