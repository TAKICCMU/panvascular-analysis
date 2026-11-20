from common import *
import csv
rows=list(csv.DictReader(open(W/'pseudobulk/model_roster.tsv'),delimiter='\t'))
for r in rows:
 p=Path(r['input']);run_job('pseudobulk_'+r['model_id'],['/usr/bin/Rscript','--vanilla',W/'code/pseudobulk_edgeR.R',p/'counts.mtx.gz',p/'units.tsv',W/'pseudobulk/genes.tsv',p/'selected_design.json',p/'results'],W/'logs'/('pseudobulk_'+r['model_id']+'.log'),reserve=3)
 assert json.loads((p/'results/MODEL.json').read_text())['status']=='COMPLETE'
dump(W/'pseudobulk/COMPLETE.json',dict(status='COMPLETE_ELIGIBLE_WITHIN_STUDY_DE',models=len(rows),cross_disease_meta_run=False,model_results={r['model_id']:sha(Path(r['input'])/'results/DE.tsv') for r in rows}))
