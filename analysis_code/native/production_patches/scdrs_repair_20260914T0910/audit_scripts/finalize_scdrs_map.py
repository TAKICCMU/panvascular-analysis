from pathlib import Path
import json,hashlib,gzip,re
import pandas as pd
import h5py,anndata as ad
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908')
W=R/'production_launches/scdrs_repair_20260914T0910'
I=W/'inputs';I.mkdir(parents=True,exist_ok=True)
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()
source=R.parent/'01_corrected_counts/audit/full_gene_identifier_mapping.tsv.gz'
old=R/'production_launches/traits_20260914T023900Z/inputs/gene_map_MAGMA_READY.tsv'
gm=pd.read_csv(old,sep='\t',keep_default_na=False)
audit=pd.read_csv(source,sep='\t',keep_default_na=False)
assert gm.feature_id.tolist()==audit.feature_name.tolist()
ap=R/'panvascular_CPU_theta2_L1L2_final_v1.h5ad'
with h5py.File(ap,'r') as h:
    var=ad.io.read_elem(h['var'])
assert gm.feature_id.tolist()==var.index.tolist()
assert len(gm)==37336 and gm.feature_id.is_unique
missing=gm.gene_symbol.eq('')
assert missing.sum()==735
assert (audit.loc[missing,'gene_mapping_status']=='symbol_fallback').all()
assert (audit.loc[missing,'stable_ensembl_id']=='').all()
assert (audit.loc[missing,'display_gene_symbol'].to_numpy()==gm.loc[missing,'feature_id'].to_numpy()).all()
assert (audit.loc[missing,'stable_feature_id'].to_numpy()==('SYMBOL_FALLBACK::'+gm.loc[missing,'feature_id']).to_numpy()).all()
gm['scdrs_gene_symbol']=gm.gene_symbol.where(~missing,gm.feature_id)
assert gm.scdrs_gene_symbol.is_unique and gm.scdrs_gene_symbol.str.len().gt(0).all()
gm['scdrs_mapping_provenance']='frozen_var_symbol_with_stable_ensembl'
gm.loc[missing,'scdrs_mapping_provenance']='exact_archived_source_symbol_only_no_ensembl_inference'
gm['ensembl_mapping_changed']=False
gm['count_aggregation_required']=False
gm.to_csv(I/'gene_map_scDRS_READY.tsv',sep='\t',index=False)
hits=json.loads((I/'local_gtf_exact_symbol_hits.json').read_text())
v43=next(v for k,v in hits.items() if 'v43' in k)
hgnc=json.loads((I/'HGNC_24509.json').read_text())['response']['docs']
assert len(hgnc)==1 and hgnc[0]['hgnc_id']=='HGNC:24509' and 'TBC1D29' in hgnc[0]['prev_symbol']
resources=json.loads((R/'production_launches/traits_20260914T023900Z/evidence/TRAIT_RESOURCES_VERIFIED.json').read_text())
traitrows=[];stats=[]
for rr in resources:
    if rr.get('module')!='scdrs' or 'path' not in rr:continue
    assert sha(rr['path'])==rr['sha256']
    gs=pd.read_csv(rr['path'],sep='\t');assert len(gs)==1
    genes=[g.split(':')[0] for g in gs.GENESET.iloc[0].split(',')]
    assert len(genes)==len(set(genes))==1000
    for g in sorted(set(genes)&set(gm.loc[missing,'feature_id'])):
        if g=='TBC1D29':proof='HGNC:24509 prev_symbol; retain legacy name, no alias rewrite';ids=hgnc[0]['ensembl_gene_id']
        else:
            assert len(v43.get(g,[]))==1,('Unresolved trait source symbol',g)
            proof='GENCODE v43 exact gene_name singleton';ids=v43[g][0]
        traitrows.append({'scenario':rr['trait'],'source_symbol':g,'evidence':proof,'external_reference_id_NOT_assigned_to_atlas':ids})
    matched=set(genes)&set(gm.scdrs_gene_symbol)
    stats.append({'scenario':rr['trait'],'geneset_size':1000,'pre_native_filter_exact_overlap':len(matched),'native_filter_overlap_pending':True,'source_symbol_only_trait_hits':len(set(genes)&set(gm.loc[missing,'feature_id']))})
    pd.DataFrame({'gene':genes,'present_in_input':[g in matched for g in genes]}).to_csv(I/(rr['trait']+'_geneset_overlap.tsv'),sep='\t',index=False)
pd.DataFrame(traitrows).to_csv(I/'source_symbol_trait_evidence.tsv',sep='\t',index=False)
gm.loc[missing].to_csv(I/'source_symbol_only_features.tsv',sep='\t',index=False)
proof={'status':'PASS_SOURCE_SYMBOL_MAPPING_FOR_NATIVE_scDRS','n_features':len(gm),'confirmed_stable_id_rows':int((~missing).sum()),'source_symbol_only_rows':int(missing.sum()),'final_var_exact_order':True,'archived_source_exact_order':True,'symbol_unique':True,'counts_columns_dropped':0,'counts_columns_merged':0,'gene_set_files_changed':False,'ensembl_ids_inferred':0,'canonical_symbol_claim_for_fallback':False,'mapping_policy':'Use already archived exact source symbols as native feature names. Ensembl remains unresolved for 735; no cross-gene alias collapse, no uppercasing, no modern nomenclature substitution. All three trait-set intersections with fallback names have independent annotation evidence; other fallback features retain original background role. This is feature identity verification, not recovery of canonical Ensembl mapping.','scenarios':stats,'sources':[{ 'path':str(p),'sha256':sha(p)} for p in [source,old,I/'HGNC_24509.json',I/'local_gtf_exact_symbol_hits.json',Path('/home/anzhen_vhd/software/scE2G/resources/genome_annotations/gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz')]],'ready_map_sha256':sha(I/'gene_map_scDRS_READY.tsv'),'hgnc_url':'https://rest.genenames.org/fetch/hgnc_id/HGNC:24509'}
(I/'SCDRS_GENE_MAPPING_VALIDATION.json').write_text(json.dumps(proof,indent=2))
print(json.dumps(proof,indent=2))
