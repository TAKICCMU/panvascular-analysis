# One entry per existing numbered table; no row/column rearrangement.
ROWS='''
Single-cell datasets and publications|Public single-cell sources, assay types and atlas-specific coverage.|One canonical study.
Sample and donor metadata|Technical-library, specimen and donor crosswalk with clinical fields and retained grouping decisions.|One technical sample identifier.
Cell-state annotations and coverage|Named L2 states, parent identities, marker lists and sample/study coverage.|One named L2 state.
GWAS datasets|Phenotypes, case/control counts, ancestry and original GWAS sources.|One GWAS phenotype.
MAGMA and scDRS cell-group associations|MAGMA gene-property results and matched scDRS results across scenarios, levels and weighting schemes.|One cell group, GWAS scenario and weighting scheme; scDRS values can recur across MAGMA weightings.
Pairwise genetic correlations|LDSC genetic correlations between the eight input phenotypes.|One ordered phenotype pair.
Exploratory factor loadings|Promax-rotated exploratory factor loadings for the input phenotypes.|One phenotype, with two loading columns.
Bifactor model loadings|Shared, Aortic and Atherosclerotic domain estimates and standard errors.|One model parameter.
Leave-one-trait factor correlations|Correlated-factor sensitivity estimates after omission of individual traits.|One recorded parameter and omission analysis.
Component-disease locus novelty assessment|Lead-variant annotations and component-disease novelty evidence; novelty is not a causal designation.|One assessed lead variant.
Locus catalogue and literature matches|FUMA loci and catalogue/literature matches under the stated LD and significance rules.|One genomic risk locus.
Regulatory evidence for prioritised genes|FUMA eQTL mapping and current colocalisation summaries across primary and leave-CAD scenarios.|One prioritised gene–cell-population pair.
Prioritised gene–cell-population pairs|The 29 selected pairs with REML expression effects, genetic evidence and regulatory context.|One prioritised gene–cell-population pair.
CoVarNet module evidence|Structural stability, context breadth and within-study coherence for the 14 modules.|One CoVarNet module.
CoVarNet cell-state weights|NMF weights for 89 included states across 14 modules.|One eligible cell state; columns are module weights.
CoVarNet evaluation criteria|Recorded structural gates and contextual evidence used to classify module confidence.|One CoVarNet module.
Section-level module spatial statistics|Spatial coherence, programme correlations and null comparisons by section and posterior summary.|One section, module and posterior setting.
Dataset-level module spatial statistics|Dataset medians after section-to-unit pooling; P columns are descriptive medians, not combined dataset P values.|One dataset, module and posterior setting.
MAGMA pathway enrichment|Complete retained MAGMA GSEA and ORA results. ORA retains overlap ≥3 before testing and BH within database and selected-gene definition.|One pathway and enrichment analysis; GSEA- and ORA-specific fields are not interchangeable.
Pseudobulk pathway enrichment|Analysis index and significant GSEA/ORA records. Complete results are in Supplementary_Table20_full_results/ and the Supplementary_Data upload archive.|Index rows: one analysis; result rows: one pathway and analysis, distinguished by Record type.
Spatial sample metadata|The 49 spatial sections, their 35 analysis units, source files and inclusion metadata.|One spatial section.
Availability of spatial region annotations|Source-reported regional information and availability of reusable spot-level annotations.|One spatial dataset.
Plaque and calcification gene programmes|Source differential-expression evidence used to define the plaque/calcification programmes; these are not new atlas differential tests.|One source gene–cluster result.
Software versions|Recorded software/resource versions and the workflows to which they apply.|One software or resource record.
LDSC intercepts|Cross-trait LDSC intercepts for the eight input GWASs.|One ordered phenotype pair.
Marker expression summaries|L1/L2 marker means, detection fractions and separate canonical-programme summaries, distinguished by Record type.|One marker–group result or programme summary; missing fields depend on Record type.
Integration benchmark|The eight representations and original biology/batch metrics and composite scores.|One integration representation.
Sample-level post-QC cell counts|Retained cell counts and sample/specimen/donor crosswalk. This table does not contain per-cell QC distributions.|One technical sample identifier.
Cell-state study coverage|Coverage of the 116 named L2 states; the complete 128 parent-qualified groups remain in Figure 1 Source Data.|One named L2 state.
Cell-type composition by study|Mean sample-level cell-type fractions within each study.|One study–L1 combination.
Unknown-cell source coverage|Unknown-cell counts and source metadata with the corresponding total sample counts.|One sample's Unknown-cell group.
Candidate spatial associations|All 87 primary q05/all-section candidate model terms, including unsupported or non-estimable terms.|One of 29 gene–cell-population pairs and one of three model terms.
Module spatial coherence|Observed and null module coherence with posterior-specific multiple-testing families.|One module, coherence metric, null method and posterior setting.
Module–programme direction across datasets|Dataset-level direction consistency for six expression programmes and 14 modules.|One module–programme pair and posterior setting.
Disease-associated cell proportions|Estimated disease-associated proportion effects and uncertainty, retaining invalid-uncertainty records.|One fitted state–comparison model.
Cell-type fractions by tissue|Sample-averaged L1 composition within tissue categories.|One tissue category; columns are L1 fractions.
Recurrent pseudobulk pathways|Descriptive recurrence of pathways across significant pseudobulk enrichment analyses.|One recurrent pathway and database.
Cell–cell communication effects|Matched LIANA and CellChat disease effects and cohort support for concordantly supported interactions.|One disease-context interaction with paired method-specific results.
Communication consensus by disease context|Agreement between communication methods within each tissue and disease background.|One tissue–disease-background combination.
GWAS–eQTL colocalisation|Complete ABF colocalisation tests, including non-supporting results, for prioritised genes.|One gene, GWAS scenario and QTL context.
QTL coverage by gene and GWAS scenario|Requested, available, missing and tested QTL contexts with evidence summaries.|One gene and GWAS scenario.
Candidate MAGMA evidence across scenarios|Gene-level MAGMA evidence for the candidate pairs across three GWAS scenarios.|One candidate pair and GWAS scenario.
Pseudobulk meta-analysis of prioritised genes|All 15,639 evaluable gene–state pairs with REML effects and target-set FDR, not only selected pairs.|One evaluable priority-gene–cell-population pair.
CoVarNet rank selection|NMF quality metrics and the recorded automatic rank-selection criteria for ranks 2–25.|One tested NMF rank.
Leave-one-study CoVarNet refits|Module matching under study omission, with similarity, overlap and matching diagnostics.|One reference-module match in a study-omission refit.
Multi-signal colocalisation sensitivity|Historical-scope multi-signal results and fit statuses; this is not a complete matched-LD screen of the current 29 candidates.|One historical gene–QTL-context sensitivity test.
Analysis sample sizes and units|The distinct cell, sample, specimen, profile and spatial analysis units used by each analysis.|One analysis-specific unit/count definition.
Spatial posterior and QC sensitivity|Candidate model terms under four posterior/QC settings, each retaining its 87-test BH family.|One candidate pair, model term and posterior/QC setting.
Cohort contributions to communication effects|Pooled interactions with their recorded cohort-effect contributions; the contribution field contains multiple cohorts per row.|One method–context interaction, not one independent cohort.
Relationships among evidence sources|Shared inputs and interpretation boundaries between evidence streams.|One relationship between analysis/evidence sources.
Multiple-testing families|Prespecified test families and correction boundaries across analyses.|One multiple-testing family definition.
Genetic cell-group enrichment summaries|Counts of tests and significant groups by genetic scenario, method, annotation level and weighting.|One scenario–method–level–weighting summary.
Target-exclusion spatial sensitivity|Paired original-reference and target-excluded spatial results, retaining the original model terms and settings.|One candidate pair, model term and posterior/QC setting.
cell2location model diagnostics|Recorded numerical/output diagnostics by spatial dataset; these do not alone establish convergence.|One spatial dataset.
Parent-lineage-matched CoVarNet null results|Matched-null comparisons for the original module-structure statistics.|One module and evaluated statistic.
Spatial and single-cell source overlap|Dataset-level source-overlap records; spatial grouping acceptance does not resolve cross-modality donor overlap.|One spatial source-overlap record.
Unit-level spatial coefficients|Within-unit pooled coefficients and SE used in the candidate meta-analysis; SE is not a confidence interval.|One candidate pair, model term, analysis unit and posterior/QC setting.
Section-level spatial coefficients|Section model coefficients, uncertainty, eligibility and expression coverage.|One candidate pair, model term, section and posterior setting.
Spatial model eligibility|Spot/positive-spot coverage and fit eligibility, retaining missing or failed results.|One candidate pair, section and posterior setting.
Cell-proportion model eligibility|All planned cell-proportion models, including ineligible and failed combinations.|One planned state–disease-background–denominator model.
CoVarNet sample and donor metadata|Profiles and grouping metadata used in CoVarNet disease association, distinct from the structural profile set.|One CoVarNet disease-association profile.
Donor-aware communication sensitivity|Full pooled communication effects and donor-aware cohort contributions across tested interactions.|One method–context interaction.
CoVarNet resampling and refits|Reference-to-refit module matches across structural perturbations and hierarchical resampling.|One reference-module match in a perturbation/refit.
Genetic model fit statistics|Fit statistics, roles and stable identifiers for the four compared genetic models; source underflow P=0 is retained.|One confirmatory genetic model.
'''.strip().splitlines()
META={i+1:dict(zip(['title','contents','row_unit'],r.split('|'))) for i,r in enumerate(ROWS)}
assert len(META)==64
