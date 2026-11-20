from pathlib import Path
import json,hashlib,shutil
from zipfile import ZipFile
from xml.etree import ElementTree as E
import numpy as np,pandas as pd
from scipy.io import mmread
O=Path(__file__).resolve().parents[1]; P=Path('/Users/taki/Desktop/投稿包_最新版_20260919'); S=P/'Source_Data_Figures'; B=Path('/Users/taki/Desktop/Panvascular_Marker_Spatial_Review_20260919'); R=Path('/Users/taki/Desktop/Codex任务1/spatial_visium_geo_local/extracted/GSE248608'); stem='GSM7917894_CA1';section='GSE248608::'+stem
hashes={}
def track(p):
 p=Path(p);hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();return p
for p in list((B/'main').glob('*.pdf'))+list((B/'extended_data').glob('*.pdf'))+[P/'Manuscript.docx',P/'Supplementary_Information.docx'] :track(p)
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
r=E.fromstring(ZipFile(P/'Manuscript.docx').read('word/document.xml'));pars={f'P{i+1}':''.join(t.text or '' for t in e.findall('.//w:t',ns)) for i,e in enumerate(r.findall('./w:body/w:p',ns))}
(O/'source_data/manuscript_spatial_verbatim.json').write_text(json.dumps({k:v for k,v in pars.items() if k in ['P32','P59','P84','P85','P86','P87','P95','P96','P137','P138','P140','P250','P262']},ensure_ascii=False,indent=2))
m=pd.read_csv(track(S/'Figure5/target_meta_descriptive.tsv'),sep='\t');u=pd.read_csv(track(S/'ExtendedData10/target_biological_unit_models.tsv'),sep='\t');t=pd.read_csv(track(S/'ExtendedData10/target_section_models.tsv'),sep='\t');roster=pd.read_csv(track(S/'ExtendedData10/target_model_roster.tsv'),sep='\t')
for name,df in [('all29_original_meta',m),('all29_original_unit_models',u),('all29_original_section_models',t),('all29_original_model_roster',roster)]:df.to_csv(O/'source_data'/f'{name}.tsv',sep='\t',index=False)
main=m[(m.posterior=='q05')&(m.sensitivity=='all_sections')&(m.term=='cell_abundance')];assert len(main)==29
selected=['FN1','SCARB1'];pool=['CCNL1','C17ORF58','STAT3','SCARB1','BAG5','SBF2','FN1']
main[main.gene.isin(pool)].to_csv(O/'source_data/candidate_pool_primary.tsv',sep='\t',index=False)
un=u[u.gene.isin(selected)&(u.posterior=='q05')&(u.sensitivity=='all_sections')&(u.term=='cell_abundance')].copy();assert len(un)==60
un['unit_CI_status']='NOT_SAVED_IN_ORIGINAL_TABLE';un['display_interval']='beta +/- original SE; NOT 95% CI';un.to_csv(O/'source_data/focus_all_contributing_units.tsv',sep='\t',index=False)
m[m.gene.isin(selected)].to_csv(O/'source_data/focus_all_terms_and_sensitivity.tsv',sep='\t',index=False)
rr=roster[roster.gene.isin(selected)&(roster.posterior=='q05')].copy();rr.to_csv(O/'source_data/focus_all_section_eligibility.tsv',sep='\t',index=False)
t[t.gene.isin(selected)&(t.posterior=='q05')].to_csv(O/'source_data/focus_all_section_coefficients.tsv',sep='\t',index=False)
j=pd.read_csv(track(S/'Figure5/spot_joint_source.tsv.gz'),sep='\t');z=j[j.section_id.eq(section)].copy();assert len(z)==1779 and z.spot_id.is_unique
pos=pd.read_csv(track(R/f'{stem}_tissue_positions_list.csv'),header=None,names=['barcode','in_tissue','array_row','array_col','fullres_y','fullres_x']);z['barcode']=z.spot_id.str.split('::').str[-1]
assert set(z.barcode)<=set(pos.barcode)
z=z.merge(pos,on='barcode',how='left',validate='one_to_one');assert z.in_tissue.eq(1).all() and z.library_size_original.ge(500).all();assert np.array_equal(z.x,z.fullres_x) and np.array_equal(z.y,z.fullres_y)
sf=json.loads(track(R/f'{stem}_scalefactors_json.json').read_text());z['hires_x']=z.x*sf['tissue_hires_scalef'];z['hires_y']=z.y*sf['tissue_hires_scalef'];track(R/f'{stem}_tissue_hires_image.png')
feat=pd.read_csv(track(R/f'{stem}_features.tsv'),sep='\t',header=None);bc=pd.read_csv(track(R/f'{stem}_barcodes.tsv'),sep='\t',header=None)[0];X=mmread(track(R/f'{stem}_matrix.mtx')).tocsr();idx=pd.Index(bc).get_indexer(z.barcode);assert (idx>=0).all();lib=np.asarray(X[:,idx].sum(0)).ravel();assert np.array_equal(lib,z.library_size_original)
for g in selected:
 ix=np.flatnonzero(feat[1].eq(g));assert len(ix)==1
 z['raw__'+g]=X[ix[0],idx].toarray().ravel();assert np.allclose(np.log1p(z['raw__'+g]*1e4/lib),z['expr__'+g],atol=1e-12)
z.to_csv(O/'source_data/CA1_existing_spot_values.tsv',sep='\t',index=False)
# New examples have no inherited cutoffs. Freeze the existing figure's positive 2nd–98th percentile display convention once; no inference.
limits={g:list(map(float,np.quantile(z.loc[z['expr__'+g]>0,'expr__'+g],[.02,.98]))) for g in selected}
limits['Fib_abundance']=[float(z['abundance__Fibroblast'].min()),float(z['abundance__Fibroblast'].max())]
summary=[]
for g in selected:
 q=main[main.gene.eq(g)].iloc[0];uu=un[un.gene.eq(g)];ss=t[t.gene.eq(g)&t.section_id.eq(section)&t.posterior.eq('q05')&t.term.eq('cell_abundance')].iloc[0]
 summary.append({'gene':g,'overall_beta':q.beta,'overall_CI_low':q.CI_low,'overall_CI_high':q.CI_high,'overall_FDR':q.FDR_BH_all_29_pairs_x_3_terms,'I2':q.I2,'units':len(uu),'studies':uu.GSE.nunique(),'unit_positive':int(uu.beta.gt(0).sum()),'unit_negative':int(uu.beta.lt(0).sum()),'section':section,'section_beta':ss.beta,'section_SE':ss.se,'section_nominal_P':ss.p,'section_FDR':'NOT_PROVIDED','n_spots':len(z),'n_positive':int(z['raw__'+g].gt(0).sum()),'total_target_UMI':int(z['raw__'+g].sum()),'section_biological_unit':ss.biological_unit_id,'metadata_status':ss.metadata_status})
pd.DataFrame(summary).to_csv(O/'source_data/focus_summary.tsv',sep='\t',index=False)
(O/'notes/display_parameters.json').write_text(json.dumps({'section':section,'section_choice':'Results P86 names FN1 and SCARB1; one non-QC-fallback calcified arterial section with known paired identity and both targets estimable. Both section slopes align with their overall sign; SCARB1 is not section-significant. Chosen after examining existing results; not blinded and not a maximum signal or minimum P selection.','shared_field':'entire original hires image; no new crop, ROI or mask','n_spots':1779,'limits':limits,'new_example_scale_rule':'Positive-expression 2nd–98th percentiles, the pre-existing plotting convention; zero separately grey; no numeric expression changes. Abundance full min–max, not expression scale.','image':str(R/f'{stem}_tissue_hires_image.png'),'scalefactors':sf,'histology':'NOT_CONFIRMED; no named layers; no new per-spot pathology screening','model':'primary q05/all_sections/cell_abundance; conditional at section program mean','unit_CI':'not stored; show original SE explicitly in detailed quantitative prototype'},indent=2))
(O/'notes/readonly_inputs_sha256.json').write_text(json.dumps(hashes,indent=2))
print(pd.DataFrame(summary).to_string(index=False))
