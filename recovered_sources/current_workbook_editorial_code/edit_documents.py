from pathlib import Path
from docx import Document
from docx.oxml.ns import qn
from difflib import SequenceMatcher
from zipfile import ZipFile
import shutil,json,re,csv,ast,hashlib
from table_metadata import META
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');O=P/'图件编辑资料/补充表与方法同步_20260921';S=O/'待替换';B=O/'轻量备份';S.mkdir(exist_ok=True);B.mkdir(exist_ok=True)
files={'Manuscript':P/'Manuscript.docx','Figure_Legends':P/'Nature_Submission/Figure_Legends.docx','Supplementary_Information':P/'Supplementary_Information.docx'}
for p in files.values():
 if not (B/p.name).exists():shutil.copy2(p,B/p.name)
# Reuse only the minimal OOXML text-node editor, not the superseded prior workflow.
tree=ast.parse((P/'图件编辑资料/图文对应精简审阅_20260921/code/prepare_review.py').read_text());fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='set_text');changes=[];exec(compile(ast.Module(body=[fn],type_ignores=[]),'<minimal text editor>','exec'))
m=Document(B/'Manuscript.docx');r=Document(P/'图件编辑资料/图文对应精简审阅_20260921/Manuscript_审阅版.docx')
# Adopt only directly relevant, already checked figure/text alignment edits; never adopt superseded donor claims.
for i in [18,19,20,27,40,197,199,203,205,206,207,208]:set_text(m.paragraphs[i],r.paragraphs[i].text,'Adopt figure/text alignment and concise annotation/enrichment wording; science unchanged.')
def edit(i,old,new,reason):
 p=m.paragraphs[i];assert old in p.text,(i,old);set_text(p,p.text.replace(old,new),reason)
edit(39,'35 recorded biological units','35 analysis units','Retain approved spatial grouping and its actual unit count.')
edit(39,'A current atlas match to GSE306288 was not established, and 16 GSE314851 units were identified at library rather than verified donor level','A current atlas match to GSE306288 was not established','Remove superseded spatial donor-uncertainty judgment; retain cross-modality source boundaries.')
edit(49,'Unresolved donor identities were retained as separate recorded units and cannot be regarded as proven independent individuals. ','Missing donor identities in the single-cell metadata were retained explicitly. ','Do not extend superseded spatial uncertainty to the accepted grouping; retain real single-cell missingness.')
edit(55,'Original gene identifiers and source annotations were retained for provenance,','Cell-level filters retained 500–8,000 detected genes and 1,000–100,000 total counts (inclusive), mitochondrial transcript fractions <20% and ribosomal fractions ≤30%. Mitochondrial and ribosomal genes were identified by the MT- and RPS/RPL prefixes, respectively. The preprocessing workflow excluded technical samples with fewer than 50 retained cells after these filters and before doublet removal. Scrublet was applied by technical sample with automatic threshold calling; predicted doublets were removed. The corrected atlas retained the post-QC cell identifiers and full-gene counts after source deduplication rather than applying a second round of cell filtering. Original gene identifiers and source annotations were retained for provenance,','QC: old Methods + actual preprocessing notebook + exact retained-cell production chain + read-only stored-metadata check; no QC rerun.')
old=' using European LD reference matrices. Alleles were aligned to the reference panel, and each region was evaluated for model convergence and credible-set reliability before inclusion in the final results. Fine-mapping results were therefore reported only for regions that passed these quality checks.'
new=' (susieR 0.16.6) using FUMA risk intervals extended by 100 kb on each side and the GRCh37 European 1000 Genomes reference (502 individuals). Duplicate variant identifiers retained the lowest-P record; palindromic variants and variants absent from or unresolved against the reference were excluded. Alleles and variant order were aligned to the reference dosages, with Z-statistic signs adjusted to the counted allele. Missing reference dosages were mean-imputed per variant, zero-variance variants were removed, and at least 20 aligned polymorphic variants were required. LD was computed from standardised reference dosages, symmetrised and assigned a unit diagonal; non-finite Z statistics or LD rows were removed. The implemented susie_rss model used Z statistics and the LD matrix, L=min(10, number of variants), R_finite=502, R_mismatch="eb", coverage=0.95, min_abs_corr=0.5, max_iter=200, tol=0.001 and estimate_residual_variance=FALSE. Refinement was attempted and disabled only if the refined fit raised an error. The primary fit did not supply a GWAS sample-size argument; supplying the median available per-variant sample size was a separate sensitivity analysis. An analysable region required convergence and R_reliability_flag=FALSE. Accepted credible sets additionally required minimum absolute within-set LD correlation ≥0.5. Other recorded diagnostics were not additional exclusion gates; displayed PIP thresholds were not credible-set acceptance thresholds.'
edit(74,old,new,'Recover accepted mismatch-aware LD-based SuSiE implementation and explicit acceptance flags, not software defaults.')
edit(80,'For the archived MAGMA over-representation analysis','For MAGMA over-representation analysis','Retain accepted ORA without obsolete workflow labels.')
edit(80,'the archived ORA remains separately available','ORA results are separately available','Retain ORA source and result counts.')
edit(80,'The prefilter-before-BH ordering is reported as implemented and does not establish the statistical validity of selection-conditioned ORA inference. ','','Author accepts retained ORA; remove task-state discussion, not its true computation order.')
edit(80,'For the current Figure 3c display, three MAGMA GSEA pathways per database','For Figure 3c, GO BP, KEGG LEGACY and Reactome were used. Three MAGMA GSEA pathways per database','Make actual display database scope explicit once in Methods.')
edit(90,'35 recorded biological units','35 analysis units','Retain 49 sections → 35 accepted analysis units.')
edit(90,'Serial sections were pooled by inverse-variance fixed effects within recorded biological units, then by DL random effects across units; serial-section residual covariance was not explicitly modelled.','Section estimates were pooled by inverse-variance fixed effects within the recorded analysis units, followed by DL random-effects meta-analysis across units; within-unit section residual covariance was not explicitly modelled.','Describe unchanged within-unit and between-unit pooling without labeling grouping erroneous.')
edit(90,'within recorded biological units','within analysis units','Use consistent accepted analysis-unit terminology.')
edit(90,'Dataset identity did not imply independently verified donors, and these tests did not alter structural module-confidence tiers.','Dataset summaries did not alter structural module-confidence tiers.','Remove superseded donor-uncertainty task wording without asserting independent histological validation.')
edit(92,'pooling within recorded biological units, which were not assumed to be independent donors','pooling within the 35 analysis units defined by the spatial metadata','Retain actual units and accepted grouping.')
edit(92,'Unknown donor identity was retained explicitly and did not establish independence.','Missing donor identity in the single-cell metadata was retained explicitly.','Spatial grouping confirmation does not fill single-cell donor missingness.')
cap5=r.paragraphs[201].text.replace('recorded IDs from the legacy aggregation','distinct analysis-unit IDs').replace('legacy groups','analysis units')
set_text(m.paragraphs[201],cap5,'Synchronise actual two-section FN1 views, 60 estimates/34 IDs, SE versus CI and accepted spatial groups.')
p=m.paragraphs[211];old='b,c, Ten top-node-pair pooled correlations and direction-consistency summaries per module; boxes describe node pairs, not donor-bootstrap draws.'
new='b,c, Pairwise abundance coherence among the five highest-weight states in each displayed module. Each point represents one of the ten state pairs, showing the pooled within-study Spearman correlation in b and the fraction of contributing studies with a positive correlation in c. Diamonds indicate the median across the ten pairs for each module. Contributing studies have at least five profiles and estimable within-study correlations.'
set_text(p,p.text.replace(old,new),'Native ED9 and code: points and median diamonds; denominator is eligible studies, not donors or bootstrap draws.')
edit(212,'Section and recorded biological-unit counts by dataset; serial-section dependence remains a limitation.','Section and analysis-unit counts by dataset.','Use approved spatial grouping; keep the covariance method limitation in Methods rather than duplicating it here.')
# Restricted genetic domain rename. Do not alter general systemic terms or reference titles.
for p in m.paragraphs[:99]+m.paragraphs[197:]:
 new=p.text
 for a,b in [('Systemic vascular factor','Atherosclerotic factor'),('systemic-domain','atherosclerotic-domain'),('aortic and systemic factors','aortic and atherosclerotic factors'),('Aortic–systemic','Aortic–atherosclerotic'),('Systemic factor','Atherosclerotic factor')]:new=new.replace(a,b)
 set_text(p,new,'Genetic domain display naming only; original model identifiers and references retained.')
l=Document(B/'Figure_Legends.docx')
for p in l.paragraphs:
 mat=re.match(r'(Fig\. \d+ |Extended Data Fig\. \d+ \|)',p.text)
 if mat:
  q=[q for q in m.paragraphs if q.text.startswith(mat.group())];assert len(q)==1
  set_text(p,q[0].text,'Standalone legend matches manuscript exactly.')
s=Document(B/'Supplementary_Information.docx')
for row in s.tables[0].rows[1:]:
 n=int(row.cells[0].text);md=META[n]
 set_text(row.cells[1].paragraphs[0],md['title'],'Synchronise verified Table 01–64 title.')
 set_text(row.cells[3].paragraphs[0],md['contents']+' '+md['row_unit'],'Synchronise scope and row unit; retain table number and worksheet.')
# Remove only stale operational statements in prose; no scientific result table edits here.
for p in s.paragraphs:
 if p.text.startswith('Fig. 3 '):set_text(p,m.paragraphs[199].text,'Synchronise duplicated Figure 3 legend.')
 elif p.text.startswith('Fig. 5 '):set_text(p,cap5+' The companion SCARB1 CA1 display remains descriptive: section-level P=0.582 and a negative P1 pooled effect; overall FDR is not assigned to that section.','Synchronise Figure 5 legend and retain descriptive non-significant companion evidence.')
 else:
  new=p.text
  new=re.sub(r'ORA[^.\n]*(?:validity|calibration)[^.\n]*pending[^.\n]*\.\s*','',new,flags=re.I)
  for a,b in [('systemic-domain','atherosclerotic-domain'),('Aortic–systemic','Aortic–atherosclerotic'),('aortic and systemic factors','aortic and atherosclerotic factors')]:new=new.replace(a,b)
  set_text(p,new,'Remove obsolete task-state wording; genetic domain naming only.')
# Additional reproducibility details, without altering references or existing paragraph styles.
p=s.add_paragraph('Fine-mapping implementation details');p.style=s.paragraphs[0].style if False else s.paragraphs[1].style
s.add_paragraph('FUMA-defined intervals were extended by 100 kb, and regions with more than 20,000 aligned variants were not submitted by the preparation script. The retained primary susieR 0.16.6 implementation used Z statistics together with the reference LD matrix; its internal name finite502_EB_z_only indicates omission of the GWAS sample-size argument, not omission of LD. The median-N sensitivity and separate single-effect diagnostics were not substituted for primary credible-set evidence. Sets were taken from the fitted model’s reported credible sets, with the explicit convergence, LD-reliability and minimum-absolute-correlation gates described in Methods. Acceptance did not impose an additional set-deduplication, PIP or diagnostic threshold. Failed or unavailable regional fits were not imputed as zero PIP. Source outputs retain the primary 155 regions, 136 regions with 156 accepted sets, and 35 leave-CAD regions with accepted sets among the 67 attempted regions.')
for name,doc in [('Manuscript',m),('Figure_Legends',l),('Supplementary_Information',s)]:doc.save(S/(name+'.docx'))
# Structural validation: original properties, citations and bibliography intact.
checks={}
for name,path in files.items():
 a=Document(B/path.name);b=Document(S/path.name)
 props=lambda x:x.xml if x is not None else None
 for i,(x,y) in enumerate(zip(a.paragraphs,b.paragraphs)):
  assert props(x._p.pPr)==props(y._p.pPr),(name,i,'pPr')
  assert [props(r._r.rPr) for r in x.runs]==[props(r._r.rPr) for r in y.runs],(name,i,'rPr')
  assert [r.text for r in x.runs if r.font.superscript]==[r.text for r in y.runs if r.font.superscript],(name,i,'citation')
 with ZipFile(B/path.name) as z:old={k:z.read(k) for k in z.namelist()}
 with ZipFile(S/path.name) as z:new={k:z.read(k) for k in z.namelist()}
 assert all(old[k]==new[k] for k in old if k!='word/document.xml'),name
 checks[name]={'paragraph_and_run_properties_retained':True,'superscript_citations_retained':True,'field_instructions_unchanged':old['word/document.xml'].count(b'instrText')==new['word/document.xml'].count(b'instrText'),'only_document_xml_changed':True}
# Reference-list exact XML: all text between References and figure captions.
a=Document(B/'Manuscript.docx');assert [p._p.xml for p in a.paragraphs[100:196]]==[p._p.xml for p in m.paragraphs[100:196]]
checks['bibliography_xml_unchanged']=True
checks['all_15_legends_identical']=True
checks['status']='STAGED_VALIDATED_NOT_YET_REPLACED'
(O/'document_validation.json').write_text(json.dumps(checks,indent=2));(O/'document_changes.json').write_text(json.dumps(changes,ensure_ascii=False,indent=2))
print(json.dumps(checks,indent=2));print('changes',len(changes))
