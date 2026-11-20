# Supplementary Tables：当前Excel实际表头与定义
来源：Supplementary_Tables.xlsx。只读导出；未改写工作簿。表头来自各表第3行，定义和单位原样取自现有Data dictionary；缺少定义时明确注明。完整列序与记录类型保留。

## Supplementary Table 1 — Table 01
最终标题：Supplementary Table 1 | Single-cell datasets and publications
实际数据：26行；Excel第4–29行；13列。
表内说明：Public single-cell sources, assay types and atlas-specific coverage. One canonical study. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One canonical study.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Canonical study accession**
   - 定义：Recorded Canonical study accession; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Technical samples (n)**
   - 定义：Technical samples represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Biological samples (n)**
   - 定义：Biological samples represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Resolved donor identifiers (n)**
   - 定义：Resolved donor identifiers represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Cells (n)**
   - 定义：Cells represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Tissue**
   - 定义：Recorded tissue for one canonical study.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Condition**
   - 定义：Recorded condition for one canonical study.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Type**
   - 定义：Recorded type for one canonical study.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Publication**
   - 定义：Recorded publication for one canonical study.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — DOI**
   - 定义：Recorded dOI for one canonical study.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — PMID**
   - 定义：Recorded PMID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Journal or repository**
   - 定义：Recorded journal or repository for one canonical study.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Year**
   - 定义：Recorded year for one canonical study.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 2 — Table 02
最终标题：Supplementary Table 2 | Sample and donor metadata
实际数据：247行；Excel第4–250行；37列。
表内说明：Technical-library, specimen and donor crosswalk with clinical fields and retained grouping decisions. One technical sample identifier. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One technical sample identifier.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Sample ID**
   - 定义：Original technical sample identifier.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Study ID**
   - 定义：Source-study accession or its canonical identifier.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Analysis sample ID**
   - 定义：Recorded Analysis sample ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Specimen ID**
   - 定义：Biological specimen identifier, distinct from a technical library and from donor identity.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Donor ID**
   - 定义：Recorded donor identifier; NR/blank is not a new independent donor.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Tissue**
   - 定义：Recorded tissue for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Tissue detail**
   - 定义：Recorded tissue detail for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Cells (n)**
   - 定义：Cells represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Replicate relationship**
   - 定义：Recorded replicate relationship for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Evidence source**
   - 定义：Recorded evidence source for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Decision**
   - 定义：Recorded decision for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Treat same-donor samples as independent**
   - 定义：Recorded treat same-donor samples as independent for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Automatic specimen pooling**
   - 定义：Recorded automatic specimen pooling for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Design approval record**
   - 定义：Recorded design approval record for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Previous specimen ID**
   - 定义：Recorded Previous specimen ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — Previous donor ID**
   - 定义：Recorded Previous donor ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — Original clinical group**
   - 定义：Recorded original clinical group for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
18. **R — Covarnet group**
   - 定义：Recorded covarnet group for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
19. **S — Group evidence**
   - 定义：Recorded group evidence for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
20. **T — Covarnet include**
   - 定义：Recorded covarnet include for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
21. **U — Grouping-review flag**
   - 定义：Recorded grouping-review flag for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
22. **V — Source diagnosis**
   - 定义：Recorded source diagnosis for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
23. **W — Pooling**
   - 定义：Recorded pooling for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
24. **X — Analysis approval record**
   - 定义：Recorded analysis approval record for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
25. **Y — Group binary**
   - 定义：Recorded group binary for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
26. **Z — Sample name**
   - 定义：Recorded sample name for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
27. **AA — Age**
   - 定义：Recorded age for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
28. **AB — Sex**
   - 定义：Recorded sex for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
29. **AC — Assay type**
   - 定义：Recorded assay type for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
30. **AD — Vascular bed**
   - 定义：Recorded vascular bed for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
31. **AE — Cohort disease context**
   - 定义：Recorded cohort disease context for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
32. **AF — Donor disease context**
   - 定义：Recorded donor disease context for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
33. **AG — Tissue disease status**
   - 定义：Recorded tissue disease status for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
34. **AH — Donor resolved**
   - 定义：Recorded donor resolved for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
35. **AI — Donor block**
   - 定义：Recorded donor block for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
36. **AJ — Profile ID**
   - 定义：Retained analysis profile identifier; not interchangeable with donor count.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
37. **AK — Healthy definition**
   - 定义：Recorded healthy definition for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 3 — Table 03
最终标题：Supplementary Table 3 | Cell-state annotations and coverage
实际数据：116行；Excel第4–119行；15列。
表内说明：Named L2 states, parent identities, marker lists and sample/study coverage. One named L2 state. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One named L2 state.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — L1 category**
   - 定义：Final Level 1 annotation. Unknown is a valid annotation, not missing data.
   - 单位/尺度：category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — L2 state**
   - 定义：Final Level 2 annotation within the stated parent lineage.
   - 单位/尺度：category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — State**
   - 定义：Recorded state for one named L2 state.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — State ID**
   - 定义：Recorded State ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Parent-qualified state**
   - 定义：Recorded parent-qualified state for one named L2 state.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Safe state label**
   - 定义：Recorded safe state label for one named L2 state.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Spatial factor label**
   - 定义：Recorded spatial factor label for one named L2 state.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Reference state name**
   - 定义：Recorded reference state name for one named L2 state.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Reference safe label**
   - 定义：Recorded reference safe label for one named L2 state.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Cells (n)**
   - 定义：Cells represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Samples (n)**
   - 定义：Samples represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Studies (n)**
   - 定义：Studies represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Donors (n)**
   - 定义：Donors represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Tissues (n)**
   - 定义：Tissues represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Marker genes**
   - 定义：Recorded marker genes for one named L2 state.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 4 — Table 04
最终标题：Supplementary Table 4 | GWAS datasets
实际数据：8行；Excel第4–11行；8列。
表内说明：Phenotypes, case/control counts, ancestry and original GWAS sources. One GWAS phenotype. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One GWAS phenotype.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Phenotype**
   - 定义：Recorded phenotype for one GWAS phenotype.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Cases (n)**
   - 定义：Cases represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Controls (n)**
   - 定义：Controls represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Total sample size**
   - 定义：Recorded total sample size for one GWAS phenotype.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Source**
   - 定义：Recorded source for one GWAS phenotype.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — PMID**
   - 定义：Recorded PMID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — DOI**
   - 定义：Recorded dOI for one GWAS phenotype.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Ancestry**
   - 定义：Recorded ancestry for one GWAS phenotype.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 5 — Table 05
最终标题：Supplementary Table 5 | MAGMA and scDRS cell-group associations
实际数据：708行；Excel第4–711行；20列。
表内说明：MAGMA gene-property results and matched scDRS results across scenarios, levels and weighting schemes. One cell group, GWAS scenario and weighting scheme; scDRS values can recur across MAGMA weightings. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One cell group, GWAS scenario and weighting scheme; scDRS values can recur across MAGMA weightings.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Annotation level**
   - 定义：Recorded annotation level for one cell group, GWAS scenario and weighting scheme; scDRS values can recur across MAGMA weightings.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — GWAS scenario**
   - 定义：Recorded gWAS scenario for one cell group, GWAS scenario and weighting scheme; scDRS values can recur across MAGMA weightings.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Weighting**
   - 定义：Recorded weighting for one cell group, GWAS scenario and weighting scheme; scDRS values can recur across MAGMA weightings.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — MAGMA variable ID**
   - 定义：Recorded MAGMA variable ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — TYPE**
   - 定义：Recorded tYPE for one cell group, GWAS scenario and weighting scheme; scDRS values can recur across MAGMA weightings.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — MODEL**
   - 定义：Recorded mODEL for one cell group, GWAS scenario and weighting scheme; scDRS values can recur across MAGMA weightings.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — MAGMA genes (n)**
   - 定义：MAGMA genes represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — MAGMA β**
   - 定义：Recorded mAGMA β for one cell group, GWAS scenario and weighting scheme; scDRS values can recur across MAGMA weightings.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — MAGMA standardised β**
   - 定义：Recorded mAGMA standardised β for one cell group, GWAS scenario and weighting scheme; scDRS values can recur across MAGMA weightings.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — MAGMA SE**
   - 定义：Standard error on the corresponding effect/estimate scale; not a confidence interval.
   - 单位/尺度：effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — MAGMA P**
   - 定义：Original P for MAGMA in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Mapped variable ID**
   - 定义：Recorded Mapped variable ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Group**
   - 定义：Recorded group for one cell group, GWAS scenario and weighting scheme; scDRS values can recur across MAGMA weightings.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — MAGMA BH FDR**
   - 定义：Original adjusted P for the stated testing family; no correction has been recomputed.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — BH family**
   - 定义：Recorded bH family for one cell group, GWAS scenario and weighting scheme; scDRS values can recur across MAGMA weightings.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — Family tests (n)**
   - 定义：Family tests represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — Cells (n)**
   - 定义：Cells represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
18. **R — ScDRS Z**
   - 定义：Recorded scDRS Z for one cell group, GWAS scenario and weighting scheme; scDRS values can recur across MAGMA weightings.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
19. **S — ScDRS P**
   - 定义：Original P for ScDRS in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
20. **T — ScDRS FDR**
   - 定义：Original adjusted P for the stated testing family; no correction has been recomputed.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 6 — Table 06
最终标题：Supplementary Table 6 | Pairwise genetic correlations
实际数据：64行；Excel第4–67行；3列。
表内说明：LDSC genetic correlations between the eight input phenotypes. One ordered phenotype pair. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One ordered phenotype pair.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Row phenotype**
   - 定义：Recorded row phenotype for one ordered phenotype pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Column phenotype**
   - 定义：Recorded column phenotype for one ordered phenotype pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Genetic correlation**
   - 定义：Recorded genetic correlation for one ordered phenotype pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 7 — Table 07
最终标题：Supplementary Table 7 | Exploratory factor loadings
实际数据：8行；Excel第4–11行；3列。
表内说明：Promax-rotated exploratory factor loadings for the input phenotypes. One phenotype, with two loading columns. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One phenotype, with two loading columns.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Phenotype**
   - 定义：Recorded phenotype for one phenotype, with two loading columns.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Exploratory factor 1 loading**
   - 定义：Recorded exploratory factor 1 loading for one phenotype, with two loading columns.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Exploratory factor 2 loading**
   - 定义：Recorded exploratory factor 2 loading for one phenotype, with two loading columns.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 8 — Table 08
最终标题：Supplementary Table 8 | Bifactor model loadings
实际数据：16行；Excel第4–19行；10列。
表内说明：Shared, Aortic and Atherosclerotic domain estimates and standard errors. One model parameter. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One model parameter.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Latent factor**
   - 定义：Recorded latent factor for one model parameter.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Latent factor description**
   - 定义：Recorded latent factor description for one model parameter.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Parameter**
   - 定义：Recorded parameter for one model parameter.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Observed phenotype**
   - 定义：Recorded observed phenotype for one model parameter.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Unstandardised estimate**
   - 定义：Recorded unstandardised estimate for one model parameter.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Unstandardised standard error**
   - 定义：Standard error on the corresponding effect/estimate scale; not a confidence interval.
   - 单位/尺度：effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Standardised estimate (genotype scale)**
   - 定义：Recorded standardised estimate (genotype scale) for one model parameter.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Standard error (genotype scale)**
   - 定义：Standard error on the corresponding effect/estimate scale; not a confidence interval.
   - 单位/尺度：effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Fully standardised estimate**
   - 定义：Recorded fully standardised estimate for one model parameter.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — P**
   - 定义：Original P for P in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 9 — Table 09
最终标题：Supplementary Table 9 | Leave-one-trait factor correlations
实际数据：9行；Excel第4–12行；13列。
表内说明：Correlated-factor sensitivity estimates after omission of individual traits. One recorded parameter and omission analysis. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One recorded parameter and omission analysis.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Latent factor**
   - 定义：Recorded latent factor for one recorded parameter and omission analysis. Fixed model symbols are retained; AORTIC denotes Aortic and ATHERO denotes Atherosclerotic.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Parameter**
   - 定义：Recorded parameter for one recorded parameter and omission analysis.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Observed phenotype**
   - 定义：Recorded observed phenotype for one recorded parameter and omission analysis.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Unstandardised estimate**
   - 定义：Recorded unstandardised estimate for one recorded parameter and omission analysis.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Unstandardised standard error**
   - 定义：Standard error on the corresponding effect/estimate scale; not a confidence interval.
   - 单位/尺度：effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Standardised estimate (genotype scale)**
   - 定义：Recorded standardised estimate (genotype scale) for one recorded parameter and omission analysis.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Standard error (genotype scale)**
   - 定义：Standard error on the corresponding effect/estimate scale; not a confidence interval.
   - 单位/尺度：effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Fully standardised estimate**
   - 定义：Recorded fully standardised estimate for one recorded parameter and omission analysis.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — P**
   - 定义：Original P for P in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Analysis**
   - 定义：Recorded analysis for one recorded parameter and omission analysis.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Omitted phenotype abbreviation**
   - 定义：Recorded omitted phenotype abbreviation for one recorded parameter and omission analysis.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Omitted phenotype**
   - 定义：Recorded omitted phenotype for one recorded parameter and omission analysis.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Model**
   - 定义：Recorded model for one recorded parameter and omission analysis.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 10 — Table 10
最终标题：Supplementary Table 10 | Component-disease locus novelty assessment
实际数据：5行；Excel第4–8行；13列。
表内说明：Lead-variant annotations and component-disease novelty evidence; novelty is not a causal designation. One assessed lead variant. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One assessed lead variant.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Lead variant**
   - 定义：Recorded lead variant for one assessed lead variant.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Chromosome**
   - 定义：Recorded chromosome for one assessed lead variant.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Position (GRCh37)**
   - 定义：Recorded position (GRCh37) for one assessed lead variant.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — P**
   - 定义：Original P for P in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Regional or nearest genes**
   - 定义：Recorded regional or nearest genes for one assessed lead variant.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — VEP consequence**
   - 定义：Recorded vEP consequence for one assessed lead variant.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — VEP impact**
   - 定义：Recorded vEP impact for one assessed lead variant.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — CADD**
   - 定义：Recorded cADD for one assessed lead variant.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Primary posterior inclusion probability (PIP)**
   - 定义：Recorded primary posterior inclusion probability (PIP) for one assessed lead variant.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Included in accepted credible set**
   - 定义：Recorded included in accepted credible set for one assessed lead variant.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Prior evidence context**
   - 定义：Recorded prior evidence context for one assessed lead variant.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Conservative interpretation**
   - 定义：Recorded conservative interpretation for one assessed lead variant.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Novelty definition**
   - 定义：Recorded novelty definition for one assessed lead variant.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 11 — Table 11
最终标题：Supplementary Table 11 | Locus catalogue and literature matches
实际数据：10行；Excel第4–13行；35列。
表内说明：FUMA loci and catalogue/literature matches under the stated LD and significance rules. One genomic risk locus. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One genomic risk locus.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Genomic locus**
   - 定义：Recorded genomic locus for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Unique variant identifier**
   - 定义：Recorded unique variant identifier for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — DbSNP identifier**
   - 定义：Recorded dbSNP identifier for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Chromosome**
   - 定义：Recorded chromosome for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Position (GRCh37)**
   - 定义：Recorded position (GRCh37) for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — P**
   - 定义：Original P for P in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Start**
   - 定义：Recorded start for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — End**
   - 定义：Recorded end for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — SNPs in locus (n)**
   - 定义：SNPs in locus represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — GWAS SNPs (n)**
   - 定义：GWAS SNPs represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Independent significant SNPs (n)**
   - 定义：Independent significant SNPs represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Independent significant SNPs**
   - 定义：Recorded independent significant SNPs for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Lead SNPs (n)**
   - 定义：Lead SNPs represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Lead SNPs**
   - 定义：Recorded lead SNPs for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — No. current catalog GWS matches**
   - 定义：Recorded no. current catalog GWS matches for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — Component disease known**
   - 定义：Recorded component disease known for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — Broad vascular known**
   - 定义：Recorded broad vascular known for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
18. **R — Cardiovascular metabolic known**
   - 定义：Recorded cardiovascular metabolic known for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
19. **S — Strongest prior p**
   - 定义：Recorded strongest prior p for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
20. **T — Max proxy r²**
   - 定义：Recorded max proxy r² for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
21. **U — Prior catalog SNPs**
   - 定义：Recorded prior catalog SNPs for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
22. **V — Prior traits**
   - 定义：Recorded prior traits for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
23. **W — Prior mapped traits**
   - 定义：Recorded prior mapped traits for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
24. **X — Prior pmids**
   - 定义：Recorded prior pmids for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
25. **Y — Prior study accessions**
   - 定义：Recorded Prior study accessions; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
26. **Z — Mapped genes**
   - 定义：Recorded mapped genes for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
27. **AA — Candidate novel vs component diseases**
   - 定义：Recorded candidate novel vs component diseases for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
28. **AB — Candidate novel vs broad vascular traits**
   - 定义：Recorded candidate novel vs broad vascular traits for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
29. **AC — Candidate novel vs cardiovascular metabolic traits**
   - 定义：Recorded candidate novel vs cardiovascular metabolic traits for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
30. **AD — Candidate globally novel in catalog**
   - 定义：Recorded candidate globally novel in catalog for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
31. **AE — Novelty status**
   - 定义：Recorded novelty status for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
32. **AF — Publication novelty class**
   - 定义：Recorded publication novelty class for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
33. **AG — Linkage-disequilibrium reference**
   - 定义：Recorded linkage-disequilibrium reference for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
34. **AH — Linkage-disequilibrium matching rule**
   - 定义：Recorded linkage-disequilibrium matching rule for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
35. **AI — GWAS Catalog significance rule**
   - 定义：Recorded gWAS Catalog significance rule for one genomic risk locus.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 12 — Table 12
最终标题：Supplementary Table 12 | Regulatory evidence for prioritised genes
实际数据：29行；Excel第4–32行；15列。
表内说明：FUMA eQTL mapping and current colocalisation summaries across primary and leave-CAD scenarios. One prioritised gene–cell-population pair. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One prioritised gene–cell-population pair.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Gene**
   - 定义：Reported gene symbol; original gene identities are unchanged.
   - 单位/尺度：gene symbol
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Cell state**
   - 定义：Recorded cell state for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — QTL status**
   - 定义：Recorded qTL status for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Primary FUMA eQTL rows (n)**
   - 定义：Primary FUMA eQTL rows represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Leave-CAD FUMA eQTL rows (n)**
   - 定义：Leave-CAD FUMA eQTL rows represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Primary tested contexts (n)**
   - 定义：Primary tested contexts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Primary strong contexts (n)**
   - 定义：Primary strong contexts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Primary maximum PP.H4**
   - 定义：Recorded primary maximum PP.H4 for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Primary best context**
   - 定义：Recorded primary best context for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Primary best evidence class**
   - 定义：Recorded primary best evidence class for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Leave-CAD tested contexts (n)**
   - 定义：Leave-CAD tested contexts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Leave-CAD strong contexts (n)**
   - 定义：Leave-CAD strong contexts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Leave-CAD maximum PP.H4**
   - 定义：Recorded leave-CAD maximum PP.H4 for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Leave-CAD best context**
   - 定义：Recorded leave-CAD best context for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Leave-CAD best evidence class**
   - 定义：Recorded leave-CAD best evidence class for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 13 — Table 13
最终标题：Supplementary Table 13 | Prioritised gene–cell-population pairs
实际数据：29行；Excel第4–32行；55列。
表内说明：The 29 selected pairs with REML expression effects, genetic evidence and regulatory context. One prioritised gene–cell-population pair. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One prioritised gene–cell-population pair.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Gene**
   - 定义：Reported gene symbol; original gene identities are unchanged.
   - 单位/尺度：gene symbol
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Annotation level**
   - 定义：Recorded annotation level for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Cell state**
   - 定义：Recorded cell state for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Disease contexts (n)**
   - 定义：Disease contexts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — GWAS-aligned contexts (n)**
   - 定义：GWAS-aligned contexts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Disease contexts**
   - 定义：Recorded disease contexts for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Meta log2FC**
   - 定义：REML random-effects meta-analysis log2 disease-versus-reference fold change.
   - 单位/尺度：log2 fold change
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Meta SE**
   - 定义：Standard error of the REML meta-analysis log2 fold change.
   - 单位/尺度：log2 fold change
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Meta Z**
   - 定义：Recorded meta Z for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Meta P**
   - 定义：Original P for Meta in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Meta 95% CI lower**
   - 定义：Recorded confidence bound on the corresponding effect scale; no bounds were reconstructed.
   - 单位/尺度：effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Meta 95% CI upper**
   - 定义：Recorded confidence bound on the corresponding effect scale; no bounds were reconstructed.
   - 单位/尺度：effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — I² (%)**
   - 定义：Between-unit heterogeneity statistic on its original percentage scale.
   - 单位/尺度：percent (0–100)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Direction agreement**
   - 定义：Fraction of contributing effects with the recorded consensus direction; the qualifying units are defined by this table.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Consensus direction**
   - 定义：Recorded consensus direction for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — Significant contexts (n)**
   - 定义：Significant contexts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — Median context log2FC**
   - 定义：Recorded median context log2FC for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
18. **R — Ensembl gene ID**
   - 定义：Stable source gene identifier; not a row number or expression value.
   - 单位/尺度：Ensembl identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
19. **S — Target-set BH FDR**
   - 定义：BH-adjusted meta-analysis P within cell state over the prespecified priority-gene set; distinct from all-gene FDR.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
20. **T — SYMBOL**
   - 定义：Recorded sYMBOL for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
21. **U — Bonferroni 0.05 leave one CAD r0.1**
   - 定义：Recorded bonferroni 0.05 leave one CAD r0.1 for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
22. **V — Bonferroni 0.05 primary r0.1**
   - 定义：Recorded bonferroni 0.05 primary r0.1 for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
23. **W — Bonferroni 0.05 primary r0.6**
   - 定义：Recorded bonferroni 0.05 primary r0.6 for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
24. **X — FDR BH leave one CAD r0.1**
   - 定义：Original adjusted P for the stated testing family; no correction has been recomputed.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
25. **Y — FDR BH primary r0.1**
   - 定义：Original adjusted P for the stated testing family; no correction has been recomputed.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
26. **Z — FDR BH primary r0.6**
   - 定义：Original adjusted P for the stated testing family; no correction has been recomputed.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
27. **AA — P leave one CAD r0.1**
   - 定义：Recorded p leave one CAD r0.1 for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
28. **AB — P primary r0.1**
   - 定义：Recorded p primary r0.1 for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
29. **AC — P primary r0.6**
   - 定义：Recorded p primary r0.6 for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
30. **AD — Primary Bonferroni**
   - 定义：Recorded primary Bonferroni for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
31. **AE — Primary FDR**
   - 定义：Original adjusted P for the stated testing family; no correction has been recomputed.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
32. **AF — Retained leave one CAD nominal**
   - 定义：Recorded retained leave one CAD nominal for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
33. **AG — Retained leave one CAD FDR**
   - 定义：Original adjusted P for the stated testing family; no correction has been recomputed.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
34. **AH — Tier**
   - 定义：Recorded tier for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
35. **AI — Strict three context shared**
   - 定义：Recorded strict three context shared, with its analysis-specific denominator; no conversion to percent.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
36. **AJ — N manifest contexts**
   - 定义：Recorded n manifest contexts for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
37. **AK — N existing ge100 qtl contexts**
   - 定义：Recorded n existing ge100 qtl contexts for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
38. **AL — N existing coloc rows**
   - 定义：Recorded n existing coloc rows for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
39. **AM — QTL status**
   - 定义：Recorded qTL status for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
40. **AN — Max historical ABF H4**
   - 定义：Recorded max historical ABF H4 for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
41. **AO — Regulatory interpretation**
   - 定义：Recorded regulatory interpretation for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
42. **AP — Primary FUMA eQTL rows (n)**
   - 定义：Primary FUMA eQTL rows represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
43. **AQ — Primary FUMA eQTL tissues**
   - 定义：Recorded primary FUMA eQTL tissues for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
44. **AR — Leave-CAD FUMA eQTL rows (n)**
   - 定义：Leave-CAD FUMA eQTL rows represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
45. **AS — Leave-CAD FUMA eQTL tissues**
   - 定义：Recorded leave-CAD FUMA eQTL tissues for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
46. **AT — Primary tested contexts (n)**
   - 定义：Primary tested contexts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
47. **AU — Primary strong contexts (n)**
   - 定义：Primary strong contexts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
48. **AV — Primary maximum PP.H4**
   - 定义：Recorded primary maximum PP.H4 for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
49. **AW — Primary best context**
   - 定义：Recorded primary best context for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
50. **AX — Primary best evidence class**
   - 定义：Recorded primary best evidence class for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
51. **AY — Leave-CAD tested contexts (n)**
   - 定义：Leave-CAD tested contexts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
52. **AZ — Leave-CAD strong contexts (n)**
   - 定义：Leave-CAD strong contexts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
53. **BA — Leave-CAD maximum PP.H4**
   - 定义：Recorded leave-CAD maximum PP.H4 for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
54. **BB — Leave-CAD best context**
   - 定义：Recorded leave-CAD best context for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
55. **BC — Leave-CAD best evidence class**
   - 定义：Recorded leave-CAD best evidence class for one prioritised gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 14 — Table 14
最终标题：Supplementary Table 14 | CoVarNet module evidence
实际数据：14行；Excel第4–17行；35列。
表内说明：Structural stability, context breadth and within-study coherence for the 14 modules. One CoVarNet module. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One CoVarNet module.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Module**
   - 定义：Recorded module for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Cell depth downsample**
   - 定义：Recorded cell depth downsample for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Condition only**
   - 定义：Recorded condition only for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Hierarchical study bootstrap**
   - 定义：Recorded hierarchical study bootstrap for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Leave one tissue**
   - 定义：Recorded leave one tissue for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Transformation**
   - 定义：Recorded transformation for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Adjacent rank median**
   - 定义：Recorded adjacent rank median for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Leave one study**
   - 定义：Recorded leave one study for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Adjacent rank median pass**
   - 定义：Recorded adjacent rank median pass for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Leave one study pass**
   - 定义：Recorded leave one study pass for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Leave one tissue pass**
   - 定义：Recorded leave one tissue pass for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Condition only pass**
   - 定义：Recorded condition only pass for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Hierarchical study bootstrap pass**
   - 定义：Recorded hierarchical study bootstrap pass for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Transformation pass**
   - 定义：Recorded transformation pass for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Stable structure**
   - 定义：Recorded stable structure for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — Failed gates**
   - 定义：Recorded failed gates for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — Top study**
   - 定义：Recorded top study for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
18. **R — Top study share**
   - 定义：Recorded top study share, with its analysis-specific denominator; no conversion to percent.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
19. **S — Effective n studies**
   - 定义：Recorded effective n studies for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
20. **T — N studies gt 1pct**
   - 定义：Recorded n studies gt 1pct for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
21. **U — Top tissue**
   - 定义：Recorded top tissue for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
22. **V — Top tissue share**
   - 定义：Recorded top tissue share, with its analysis-specific denominator; no conversion to percent.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
23. **W — Effective n tissues**
   - 定义：Recorded effective n tissues for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
24. **X — N tissues gt 1pct**
   - 定义：Recorded n tissues gt 1pct for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
25. **Y — N top5 pairs meta analyzed**
   - 定义：Recorded n top5 pairs meta analyzed for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
26. **Z — Median pooled rho**
   - 定义：Recorded median pooled rho for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
27. **AA — Minimum pooled rho**
   - 定义：Recorded minimum pooled rho for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
28. **AB — Positive pair fraction**
   - 定义：Fraction of evaluated top-five-node pairs with positive pooled correlation.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
29. **AC — Positive FDR pair fraction**
   - 定义：Recorded positive fdr pair fraction, with its analysis-specific denominator; no conversion to percent.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
30. **AD — Median direction consistency**
   - 定义：Median fraction of contributing studies with positive within-study correlation across the top-five-node pairs.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
31. **AE — Median I² (%)**
   - 定义：Between-unit heterogeneity statistic on its original percentage scale.
   - 单位/尺度：percent (0–100)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
32. **AF — Broad context**
   - 定义：Recorded broad context for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
33. **AG — Coherent within study**
   - 定义：Recorded coherent within study for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
34. **AH — Module tier**
   - 定义：Recorded module tier for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
35. **AI — Module tier (Chinese)**
   - 定义：Recorded module tier (Chinese) for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 15 — Table 15
最终标题：Supplementary Table 15 | CoVarNet cell-state weights
实际数据：89行；Excel第4–92行；15列。
表内说明：NMF weights for 89 included states across 14 modules. One eligible cell state; columns are module weights. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One eligible cell state; columns are module weights.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Cell state**
   - 定义：Recorded cell state for one eligible cell state; columns are module weights.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — CM01**
   - 定义：Recorded cM01 for one eligible cell state; columns are module weights.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — CM02**
   - 定义：Recorded cM02 for one eligible cell state; columns are module weights.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — CM03**
   - 定义：Recorded cM03 for one eligible cell state; columns are module weights.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — CM04**
   - 定义：Recorded cM04 for one eligible cell state; columns are module weights.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — CM05**
   - 定义：Recorded cM05 for one eligible cell state; columns are module weights.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — CM06**
   - 定义：Recorded cM06 for one eligible cell state; columns are module weights.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — CM07**
   - 定义：Recorded cM07 for one eligible cell state; columns are module weights.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — CM08**
   - 定义：Recorded cM08 for one eligible cell state; columns are module weights.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — CM09**
   - 定义：Recorded cM09 for one eligible cell state; columns are module weights.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — CM10**
   - 定义：Recorded cM10 for one eligible cell state; columns are module weights.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — CM11**
   - 定义：Recorded cM11 for one eligible cell state; columns are module weights.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — CM12**
   - 定义：Recorded cM12 for one eligible cell state; columns are module weights.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — CM13**
   - 定义：Recorded cM13 for one eligible cell state; columns are module weights.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — CM14**
   - 定义：Recorded cM14 for one eligible cell state; columns are module weights.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 16 — Table 16
最终标题：Supplementary Table 16 | CoVarNet evaluation criteria
实际数据：14行；Excel第4–17行；12列。
表内说明：Recorded structural gates and contextual evidence used to classify module confidence. One CoVarNet module. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One CoVarNet module.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Module**
   - 定义：Recorded module for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Adjacent rank median pass**
   - 定义：Recorded adjacent rank median pass for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Leave one study pass**
   - 定义：Recorded leave one study pass for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Leave one tissue pass**
   - 定义：Recorded leave one tissue pass for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Condition only pass**
   - 定义：Recorded condition only pass for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Hierarchical study bootstrap pass**
   - 定义：Recorded hierarchical study bootstrap pass for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Transformation pass**
   - 定义：Recorded transformation pass for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Stable structure**
   - 定义：Recorded stable structure for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Failed gates**
   - 定义：Recorded failed gates for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Broad context**
   - 定义：Recorded broad context for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Coherent within study**
   - 定义：Recorded coherent within study for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Module tier**
   - 定义：Recorded module tier for one CoVarNet module.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 17 — Table 17
最终标题：Supplementary Table 17 | Section-level module spatial statistics
实际数据：1,372行；Excel第4–1375行；23列。
表内说明：Spatial coherence, programme correlations and null comparisons by section and posterior summary. One section, module and posterior setting. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One section, module and posterior setting.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Study ID**
   - 定义：Source-study accession or its canonical identifier.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Section ID**
   - 定义：Recorded Section ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Analysis unit ID**
   - 定义：Recorded grouping used to combine section estimates before between-unit meta-analysis; the grouping is retained as approved.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Posterior summary**
   - 定义：cell2location posterior summary used for abundance: q05 or posterior mean.
   - 单位/尺度：summary label
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Module**
   - 定义：Recorded module for one section, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — QC fallback**
   - 定义：Recorded qC fallback for one section, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Spots (n)**
   - 定义：Spots represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Within spot coherence**
   - 定义：Recorded within spot coherence for one section, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Spatial lag coherence**
   - 定义：Recorded spatial lag coherence for one section, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — ρ: fibrous-cap core**
   - 定义：Recorded ρ: fibrous-cap core for one section, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — ρ: fibrous-cap nonstructural sensitivity**
   - 定义：Recorded ρ: fibrous-cap nonstructural sensitivity for one section, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — ρ: calcification**
   - 定义：Recorded ρ: calcification for one section, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — ρ: lipid foam**
   - 定义：Recorded ρ: lipid foam for one section, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — ρ: inflammation**
   - 定义：Recorded ρ: inflammation for one section, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — ρ: endothelial activation**
   - 定义：Recorded ρ: endothelial activation for one section, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — Fallback within-spot P**
   - 定义：Original P for Fallback within-spot in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — Fallback within-spot null mean**
   - 定义：Recorded fallback within-spot null mean for one section, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
18. **R — Fallback lag P**
   - 定义：Original P for Fallback lag in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
19. **S — Fallback lag null mean**
   - 定义：Recorded fallback lag null mean for one section, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
20. **T — Parent-matched within-spot P**
   - 定义：Original P for Parent-matched within-spot in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
21. **U — Parent-matched within-spot null mean**
   - 定义：Recorded parent-matched within-spot null mean for one section, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
22. **V — Parent-matched lag P**
   - 定义：Original P for Parent-matched lag in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
23. **W — Parent-matched lag null mean**
   - 定义：Recorded parent-matched lag null mean for one section, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 18 — Table 18
最终标题：Supplementary Table 18 | Dataset-level module spatial statistics
实际数据：196行；Excel第4–199行；21列。
表内说明：Dataset medians after section-to-unit pooling; P columns are descriptive medians, not combined dataset P values. One dataset, module and posterior setting. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One dataset, module and posterior setting.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Posterior summary**
   - 定义：cell2location posterior summary used for abundance: q05 or posterior mean.
   - 单位/尺度：summary label
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Module**
   - 定义：Recorded module for one dataset, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Study ID**
   - 定义：Source-study accession or its canonical identifier.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Median QC-fallback indicator**
   - 定义：Median of the recorded section-to-unit-to-dataset summary; P medians are descriptive, not a combined significance test.
   - 单位/尺度：original scale, median summary
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Median spots per unit**
   - 定义：Median of the recorded section-to-unit-to-dataset summary; P medians are descriptive, not a combined significance test.
   - 单位/尺度：original scale, median summary
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Within spot coherence**
   - 定义：Recorded within spot coherence for one dataset, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Spatial lag coherence**
   - 定义：Recorded spatial lag coherence for one dataset, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — ρ: fibrous-cap core**
   - 定义：Recorded ρ: fibrous-cap core for one dataset, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — ρ: fibrous-cap nonstructural sensitivity**
   - 定义：Recorded ρ: fibrous-cap nonstructural sensitivity for one dataset, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — ρ: calcification**
   - 定义：Recorded ρ: calcification for one dataset, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — ρ: lipid foam**
   - 定义：Recorded ρ: lipid foam for one dataset, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — ρ: inflammation**
   - 定义：Recorded ρ: inflammation for one dataset, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — ρ: endothelial activation**
   - 定义：Recorded ρ: endothelial activation for one dataset, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Median Fallback within-spot P (descriptive)**
   - 定义：Median of the recorded section-to-unit-to-dataset summary; P medians are descriptive, not a combined significance test.
   - 单位/尺度：original scale, median summary
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Fallback within-spot null mean**
   - 定义：Recorded fallback within-spot null mean for one dataset, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — Median Fallback lag P (descriptive)**
   - 定义：Median of the recorded section-to-unit-to-dataset summary; P medians are descriptive, not a combined significance test.
   - 单位/尺度：original scale, median summary
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — Fallback lag null mean**
   - 定义：Recorded fallback lag null mean for one dataset, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
18. **R — Median Parent-matched within-spot P (descriptive)**
   - 定义：Median of the recorded section-to-unit-to-dataset summary; P medians are descriptive, not a combined significance test.
   - 单位/尺度：original scale, median summary
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
19. **S — Parent-matched within-spot null mean**
   - 定义：Recorded parent-matched within-spot null mean for one dataset, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
20. **T — Median Parent-matched lag P (descriptive)**
   - 定义：Median of the recorded section-to-unit-to-dataset summary; P medians are descriptive, not a combined significance test.
   - 单位/尺度：original scale, median summary
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
21. **U — Parent-matched lag null mean**
   - 定义：Recorded parent-matched lag null mean for one dataset, module and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 19 — Table 19
最终标题：Supplementary Table 19 | MAGMA pathway enrichment
实际数据：16,194行；Excel第4–16197行；19列。
表内说明：Complete retained MAGMA GSEA and ORA results. ORA retains overlap ≥3 before testing and BH within database and selected-gene definition. One pathway and enrichment analysis; GSEA- and ORA-specific fields are not interchangeable. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One pathway and enrichment analysis; GSEA- and ORA-specific fields are not interchangeable.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Pathway**
   - 定义：Recorded pathway for one pathway and enrichment analysis; GSEA- and ORA-specific fields are not interchangeable.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — P**
   - 定义：Original P for P in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Adjusted P**
   - 定义：Original adjusted P for the stated enrichment analysis; adjustment families follow the Methods.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — log2 error estimate**
   - 定义：Recorded log2 error estimate for one pathway and enrichment analysis; GSEA- and ORA-specific fields are not interchangeable.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — ES**
   - 定义：Recorded eS for one pathway and enrichment analysis; GSEA- and ORA-specific fields are not interchangeable.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — NES**
   - 定义：Normalised enrichment score from the indicated preranked GSEA; positive and negative values preserve direction.
   - 单位/尺度：signed enrichment score
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Gene-set size (n)**
   - 定义：Number of pathway genes represented in the ranked GSEA input.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Database**
   - 定义：Recorded database for one pathway and enrichment analysis; GSEA- and ORA-specific fields are not interchangeable.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Analysis ID**
   - 定义：Recorded Analysis ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Annotation level**
   - 定义：Recorded annotation level for one pathway and enrichment analysis; GSEA- and ORA-specific fields are not interchangeable.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Cell state**
   - 定义：Recorded cell state for one pathway and enrichment analysis; GSEA- and ORA-specific fields are not interchangeable.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Ranked genes (n)**
   - 定义：Ranked genes represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Analysis**
   - 定义：Recorded analysis for one pathway and enrichment analysis; GSEA- and ORA-specific fields are not interchangeable.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Overlap (n)**
   - 定义：Number of selected genes overlapping the tested pathway.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Pathway genes (n)**
   - 定义：Pathway genes represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — OR**
   - 定义：MAGMA ORA odds ratio with 0.5 continuity correction.
   - 单位/尺度：odds ratio
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — Direction**
   - 定义：Recorded direction for one pathway and enrichment analysis; GSEA- and ORA-specific fields are not interchangeable.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
18. **R — Selected genes (n)**
   - 定义：Selected genes represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
19. **S — Background genes (n)**
   - 定义：Background genes represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 20 — Table 20
最终标题：Supplementary Table 20 | Pseudobulk pathway enrichment
实际数据：36,752行；Excel第4–36755行；27列。
表内说明：Analysis index and significant GSEA/ORA records. Complete results are in Supplementary_Table20_full_results/ and the Supplementary_Data upload archive. Index rows: one analysis; result rows: one pathway and analysis, distinguished by Record type. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：Index rows: one analysis; result rows: one pathway and analysis, distinguished by Record type.
混合记录类型（同一数据区，不另建区块）：Record type；"ANALYSIS_INDEX": 93行; "BH_FDR_lt_0.05": 36659行。不适用于某类型的列保留原空白，不能当0。

1. **A — Annotation level**
   - 定义：Recorded annotation level for index rows: one analysis; result rows: one pathway and analysis, distinguished by Record type.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Cell state**
   - 定义：Recorded cell state for index rows: one analysis; result rows: one pathway and analysis, distinguished by Record type.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Analysis ID**
   - 定义：Recorded Analysis ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Record type**
   - 定义：Recorded record type for index rows: one analysis; result rows: one pathway and analysis, distinguished by Record type.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Role**
   - 定义：Recorded role for index rows: one analysis; result rows: one pathway and analysis, distinguished by Record type.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — N GSEA**
   - 定义：N GSEA represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — GSEA FDR05**
   - 定义：Number of GSEA pathways meeting native adjusted P <0.05 for the indexed analysis.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Native NA**
   - 定义：Count of unavailable native enrichment outputs for this analysis index; not a zero P value.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — N ORA**
   - 定义：N ORA represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — ORA FDR05 overlap3**
   - 定义：Number of ORA pathways meeting native adjusted P <0.05 and recorded overlap ≥3 for the indexed analysis.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Pathway**
   - 定义：Recorded pathway for index rows: one analysis; result rows: one pathway and analysis, distinguished by Record type.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — P**
   - 定义：Original P for P in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Adjusted P**
   - 定义：Original adjusted P for the stated enrichment analysis; adjustment families follow the Methods.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — log2 error estimate**
   - 定义：Recorded log2 error estimate for index rows: one analysis; result rows: one pathway and analysis, distinguished by Record type.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — ES**
   - 定义：Recorded eS for index rows: one analysis; result rows: one pathway and analysis, distinguished by Record type.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — NES**
   - 定义：Normalised enrichment score from the indicated preranked GSEA; positive and negative values preserve direction.
   - 单位/尺度：signed enrichment score
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — Gene-set size (n)**
   - 定义：Number of pathway genes represented in the ranked GSEA input.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
18. **R — Leading-edge genes**
   - 定义：GSEA leading-edge subset supplied by fgsea, not the full pathway membership.
   - 单位/尺度：gene list
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
19. **S — Database**
   - 定义：Recorded database for index rows: one analysis; result rows: one pathway and analysis, distinguished by Record type.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
20. **T — Ranked genes (n)**
   - 定义：Ranked genes represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
21. **U — Analysis kind**
   - 定义：Recorded analysis kind for index rows: one analysis; result rows: one pathway and analysis, distinguished by Record type.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
22. **V — Pathway genes (n)**
   - 定义：Pathway genes represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
23. **W — Overlap (n)**
   - 定义：Number of selected genes overlapping the tested pathway.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
24. **X — Direction**
   - 定义：Recorded direction for index rows: one analysis; result rows: one pathway and analysis, distinguished by Record type.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
25. **Y — Background genes (n)**
   - 定义：Background genes represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
26. **Z — Selected genes (n)**
   - 定义：Selected genes represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
27. **AA — Display overlap ge3**
   - 定义：Recorded display overlap ge3 for index rows: one analysis; result rows: one pathway and analysis, distinguished by Record type.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 21 — Table 21
最终标题：Supplementary Table 21 | Spatial sample metadata
实际数据：49行；Excel第4–52行；21列。
表内说明：The 49 spatial sections, their 35 analysis units, source files and inclusion metadata. One spatial section. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One spatial section.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Study ID**
   - 定义：Source-study accession or its canonical identifier.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Sample ID**
   - 定义：Original technical sample identifier.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Section ID**
   - 定义：Recorded Section ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Analysis unit ID**
   - 定义：Recorded grouping used to combine section estimates before between-unit meta-analysis; the grouping is retained as approved.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Section label**
   - 定义：Recorded section label for one spatial section.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Sex**
   - 定义：Recorded sex for one spatial section.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Age**
   - 定义：Recorded age for one spatial section.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Metadata status**
   - 定义：Source metadata status retained from model execution. It does not supersede author confirmation of the current analysis-unit grouping.
   - 单位/尺度：source status
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Condition**
   - 定义：Recorded condition for one spatial section.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Tissue**
   - 定义：Recorded tissue for one spatial section.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Dataset role**
   - 定义：Recorded dataset role for one spatial section.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Include mapping**
   - 定义：Recorded include mapping for one spatial section.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Include within section**
   - 定义：Recorded include within section for one spatial section.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Include group contrast**
   - 定义：Recorded include group contrast for one spatial section.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Matrix**
   - 定义：Recorded matrix for one spatial section.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — Barcodes**
   - 定义：Recorded barcodes for one spatial section.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — Features**
   - 定义：Recorded features for one spatial section.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
18. **R — Positions**
   - 定义：Recorded positions for one spatial section.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
19. **S — Scalefactors**
   - 定义：Recorded scalefactors for one spatial section.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
20. **T — Hires image**
   - 定义：Recorded hires image for one spatial section.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
21. **U — Metadata note**
   - 定义：Recorded metadata note for one spatial section.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 22 — Table 22
最终标题：Supplementary Table 22 | Availability of spatial region annotations
实际数据：7行；Excel第4–10行；8列。
表内说明：Source-reported regional information and availability of reusable spot-level annotations. One spatial dataset. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One spatial dataset.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — GEO accession**
   - 定义：Recorded GEO accession; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Tissue or disease context**
   - 定义：Recorded tissue or disease context for one spatial dataset.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Sections (n)**
   - 定义：Sections represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Analysis units (n)**
   - 定义：Analysis units represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Reusable section-level label**
   - 定义：Recorded reusable section-level label for one spatial dataset.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Author-reported regional information**
   - 定义：Recorded author-reported regional information for one spatial dataset.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Machine-readable spot-level region of interest**
   - 定义：Recorded machine-readable spot-level region of interest for one spatial dataset.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Role in this study**
   - 定义：Recorded role in this study for one spatial dataset.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 23 — Table 23
最终标题：Supplementary Table 23 | Plaque and calcification gene programmes
实际数据：50行；Excel第4–53行；8列。
表内说明：Source differential-expression evidence used to define the plaque/calcification programmes; these are not new atlas differential tests. One source gene–cluster result. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One source gene–cluster result.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Gene symbol (source column)**
   - 定义：Recorded gene symbol (source column) for one source gene–cluster result.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — P**
   - 定义：Original P for P in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Average log2 fold change**
   - 定义：Recorded average log2 fold change for one source gene–cluster result.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Within-cluster detected fraction**
   - 定义：Source fraction of cells with detectable expression within the source cluster; includes zero-expression cells in the denominator.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Outside-cluster detected fraction**
   - 定义：Source fraction of cells with detectable expression outside the source cluster.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Adjusted P**
   - 定义：Original adjusted P for the stated enrichment analysis; adjustment families follow the Methods.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Cluster**
   - 定义：Recorded cluster for one source gene–cluster result.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Gene symbol (harmonised)**
   - 定义：Recorded gene symbol (harmonised) for one source gene–cluster result.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 24 — Table 24
最终标题：Supplementary Table 24 | Software versions
实际数据：15行；Excel第4–18行；3列。
表内说明：Recorded software/resource versions and the workflows to which they apply. One software or resource record. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One software or resource record.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Software or resource**
   - 定义：Recorded software or resource for one software or resource record.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Version**
   - 定义：Recorded version for one software or resource record.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Recorded scope**
   - 定义：Recorded recorded scope for one software or resource record.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 25 — Table 25
最终标题：Supplementary Table 25 | LDSC intercepts
实际数据：64行；Excel第4–67行；3列。
表内说明：Cross-trait LDSC intercepts for the eight input GWASs. One ordered phenotype pair. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One ordered phenotype pair.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Phenotype**
   - 定义：Recorded phenotype for one ordered phenotype pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Column phenotype**
   - 定义：Recorded column phenotype for one ordered phenotype pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — LDSC intercept**
   - 定义：Recorded lDSC intercept for one ordered phenotype pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 26 — Table 26
最终标题：Supplementary Table 26 | Marker expression summaries
实际数据：6,114行；Excel第4–6117行；16列。
表内说明：L1/L2 marker means, detection fractions and separate canonical-programme summaries, distinguished by Record type. One marker–group result or programme summary; missing fields depend on Record type. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One marker–group result or programme summary; missing fields depend on Record type.
混合记录类型（同一数据区，不另建区块）：Record type；"L1 canonical-union descriptive expression": 1584行; "L2 displayed-panel descriptive expression": 4290行; "Canonical-panel mean of gene-wise across-L1 z scores": 240行。不适用于某类型的列保留原空白，不能当0。

1. **A — L1 category**
   - 定义：Final Level 1 annotation. Unknown is a valid annotation, not missing data.
   - 单位/尺度：category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — L2 state**
   - 定义：Final Level 2 annotation within the stated parent lineage.
   - 单位/尺度：category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Record type**
   - 定义：Recorded record type for one marker–group result or programme summary; missing fields depend on Record type.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Gene identifier**
   - 定义：Gene symbol or source gene identifier present in the marker export; not uniformly an Ensembl ID.
   - 单位/尺度：source identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Mean log1p expression**
   - 定义：Mean full-library log1p(CP10k) expression in the indicated group, including zero-expression cells.
   - 单位/尺度：log1p(CP10k)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Detected fraction**
   - 定义：Within-group fraction of cells with raw counts >0. Missing measurement is not converted to zero.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Other-cell mean log1p**
   - 定义：Recorded other-cell mean log1p for one marker–group result or programme summary; missing fields depend on Record type.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Other-cell detected fraction**
   - 定义：Recorded other-cell detected fraction, with its analysis-specific denominator; no conversion to percent.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Mean expression difference**
   - 定义：Recorded mean expression difference for one marker–group result or programme summary; missing fields depend on Record type.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Cells (n)**
   - 定义：Cells represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Other cells (n)**
   - 定义：Other cells represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Reference scope**
   - 定义：Recorded reference scope for one marker–group result or programme summary; missing fields depend on Record type.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Descriptive rank**
   - 定义：Recorded descriptive rank for one marker–group result or programme summary; missing fields depend on Record type.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Cohort**
   - 定义：Recorded cohort for one marker–group result or programme summary; missing fields depend on Record type.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Canonical programme**
   - 定义：Recorded canonical programme for one marker–group result or programme summary; missing fields depend on Record type.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — Mean gene Z score**
   - 定义：Canonical programme summary on the recorded gene-standardised scale; not raw mean expression or Figure 1g min–max colour.
   - 单位/尺度：gene Z score
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 27 — Table 27
最终标题：Supplementary Table 27 | Integration benchmark
实际数据：8行；Excel第4–11行；15列。
表内说明：The eight representations and original biology/batch metrics and composite scores. One integration representation. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One integration representation.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Method**
   - 定义：Recorded method for one integration representation.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — BRAS**
   - 定义：Batch-removal adjusted silhouette metric on its original scIB scale.
   - 单位/尺度：score
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — cLISI**
   - 定义：Recorded cLISI for one integration representation.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Graph connectivity**
   - 定义：Recorded graph connectivity for one integration representation.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — iLISI**
   - 定义：Recorded iLISI for one integration representation.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Isolated-label score**
   - 定义：Recorded isolated-label score for one integration representation.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — kBET**
   - 定义：Recorded kBET for one integration representation.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — KMeans ARI**
   - 定义：Recorded kMeans ARI for one integration representation.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — KMeans NMI**
   - 定义：Recorded kMeans NMI for one integration representation.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — PCR comparison**
   - 定义：Recorded pCR comparison for one integration representation.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Cell-type ASW**
   - 定义：Recorded cell-type ASW for one integration representation.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Batch score**
   - 定义：Recorded batch score for one integration representation.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Biology score**
   - 定义：Recorded biology score for one integration representation.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Overall score**
   - 定义：0.6 × biology score + 0.4 × batch score on the implemented metric scales.
   - 单位/尺度：score
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Rank**
   - 定义：Recorded rank for one integration representation.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 28 — Table 28
最终标题：Supplementary Table 28 | Sample-level post-QC cell counts
实际数据：247行；Excel第4–250行；6列。
表内说明：Retained cell counts and sample/specimen/donor crosswalk. This table does not contain per-cell QC distributions. One technical sample identifier. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One technical sample identifier.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Sample ID**
   - 定义：Original technical sample identifier.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Study ID**
   - 定义：Source-study accession or its canonical identifier.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Cells (n)**
   - 定义：Cells represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Specimen ID**
   - 定义：Biological specimen identifier, distinct from a technical library and from donor identity.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Donor ID**
   - 定义：Recorded donor identifier; NR/blank is not a new independent donor.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Tissue**
   - 定义：Recorded tissue for one technical sample identifier.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 29 — Table 29
最终标题：Supplementary Table 29 | Cell-state study coverage
实际数据：116行；Excel第4–119行；14列。
表内说明：Coverage of the 116 named L2 states; the complete 128 parent-qualified groups remain in Figure 1 Source Data. One named L2 state. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One named L2 state.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — L1 category**
   - 定义：Final Level 1 annotation. Unknown is a valid annotation, not missing data.
   - 单位/尺度：category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — L2 state**
   - 定义：Final Level 2 annotation within the stated parent lineage.
   - 单位/尺度：category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — State**
   - 定义：Recorded state for one named L2 state.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — State ID**
   - 定义：Recorded State ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Parent-qualified state**
   - 定义：Recorded parent-qualified state for one named L2 state.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Safe state label**
   - 定义：Recorded safe state label for one named L2 state.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Spatial factor label**
   - 定义：Recorded spatial factor label for one named L2 state.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Reference state name**
   - 定义：Recorded reference state name for one named L2 state.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Reference safe label**
   - 定义：Recorded reference safe label for one named L2 state.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Cells (n)**
   - 定义：Cells represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Samples (n)**
   - 定义：Samples represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Studies (n)**
   - 定义：Studies represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Donors (n)**
   - 定义：Donors represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Tissues (n)**
   - 定义：Tissues represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 30 — Table 30
最终标题：Supplementary Table 30 | Cell-type composition by study
实际数据：416行；Excel第4–419行；3列。
表内说明：Mean sample-level cell-type fractions within each study. One study–L1 combination. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One study–L1 combination.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Study ID**
   - 定义：Source-study accession or its canonical identifier.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — L1 lineage**
   - 定义：Recorded l1 lineage for one study–L1 combination.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Mean fraction across samples**
   - 定义：Recorded mean fraction across samples, with its analysis-specific denominator; no conversion to percent.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 31 — Table 31
最终标题：Supplementary Table 31 | Unknown-cell source coverage
实际数据：247行；Excel第4–250行；40列。
表内说明：Unknown-cell counts and source metadata with the corresponding total sample counts. One sample's Unknown-cell group. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One sample's Unknown-cell group.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — L1 category**
   - 定义：Final Level 1 annotation. Unknown is a valid annotation, not missing data.
   - 单位/尺度：category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — State**
   - 定义：Recorded state for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Sample ID**
   - 定义：Original technical sample identifier.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Unknown cells (n)**
   - 定义：Unknown cells represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Study ID**
   - 定义：Source-study accession or its canonical identifier.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Analysis sample ID**
   - 定义：Recorded Analysis sample ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Specimen ID**
   - 定义：Biological specimen identifier, distinct from a technical library and from donor identity.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Donor ID**
   - 定义：Recorded donor identifier; NR/blank is not a new independent donor.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Tissue**
   - 定义：Recorded tissue for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Tissue detail**
   - 定义：Recorded tissue detail for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Total sample cells (n)**
   - 定义：Total sample cells represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Replicate relationship**
   - 定义：Recorded replicate relationship for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Evidence source**
   - 定义：Recorded evidence source for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Decision**
   - 定义：Recorded decision for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Treat same-donor samples as independent**
   - 定义：Recorded treat same-donor samples as independent for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — Automatic specimen pooling**
   - 定义：Recorded automatic specimen pooling for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — Design approval record**
   - 定义：Recorded design approval record for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
18. **R — Previous specimen ID**
   - 定义：Recorded Previous specimen ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
19. **S — Previous donor ID**
   - 定义：Recorded Previous donor ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
20. **T — Original clinical group**
   - 定义：Recorded original clinical group for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
21. **U — Covarnet group**
   - 定义：Recorded covarnet group for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
22. **V — Group evidence**
   - 定义：Recorded group evidence for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
23. **W — Covarnet include**
   - 定义：Recorded covarnet include for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
24. **X — Grouping-review flag**
   - 定义：Recorded grouping-review flag for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
25. **Y — Source diagnosis**
   - 定义：Recorded source diagnosis for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
26. **Z — Pooling**
   - 定义：Recorded pooling for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
27. **AA — Analysis approval record**
   - 定义：Recorded analysis approval record for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
28. **AB — Group binary**
   - 定义：Recorded group binary for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
29. **AC — Sample name**
   - 定义：Recorded sample name for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
30. **AD — Age**
   - 定义：Recorded age for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
31. **AE — Sex**
   - 定义：Recorded sex for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
32. **AF — Assay type**
   - 定义：Recorded assay type for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
33. **AG — Vascular bed**
   - 定义：Recorded vascular bed for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
34. **AH — Cohort disease context**
   - 定义：Recorded cohort disease context for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
35. **AI — Donor disease context**
   - 定义：Recorded donor disease context for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
36. **AJ — Tissue disease status**
   - 定义：Recorded tissue disease status for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
37. **AK — Donor resolved**
   - 定义：Recorded donor resolved for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
38. **AL — Donor block**
   - 定义：Recorded donor block for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
39. **AM — Profile ID**
   - 定义：Retained analysis profile identifier; not interchangeable with donor count.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
40. **AN — Healthy definition**
   - 定义：Recorded healthy definition for one sample's Unknown-cell group.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 32 — Table 32
最终标题：Supplementary Table 32 | Candidate spatial associations
实际数据：87行；Excel第4–90行；17列。
表内说明：All 87 primary q05/all-section candidate model terms, including unsupported or non-estimable terms. One of 29 gene–cell-population pairs and one of three model terms. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One of 29 gene–cell-population pairs and one of three model terms.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Gene**
   - 定义：Reported gene symbol; original gene identities are unchanged.
   - 单位/尺度：gene symbol
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Cell context**
   - 定义：Recorded cell context for one of 29 gene–cell-population pairs and one of three model terms.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Units (n)**
   - 定义：Number of contributing analysis units in the spatial meta-analysis, not number of spots or gene–unit estimates.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — β**
   - 定义：Conditional spatial log-link coefficient for the stated model term.
   - 单位/尺度：natural-log link coefficient
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — SE**
   - 定义：Standard error of the conditional spatial coefficient; not a confidence interval.
   - 单位/尺度：natural-log link coefficient
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Z**
   - 定义：Recorded z for one of 29 gene–cell-population pairs and one of three model terms.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — P**
   - 定义：Original P for P in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — τ²**
   - 定义：Estimated between-unit effect variance.
   - 单位/尺度：squared effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — I² (%)**
   - 定义：Between-unit heterogeneity statistic on its original percentage scale.
   - 单位/尺度：percent (0–100)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Posterior summary**
   - 定义：cell2location posterior summary used for abundance: q05 or posterior mean.
   - 单位/尺度：summary label
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Expression programme**
   - 定义：Molecular expression signature used in the spatial model; not an independently assigned histopathological region.
   - 单位/尺度：programme name
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Model term**
   - 定义：cell_abundance: conditional abundance coefficient; pathology: expression-programme coefficient; interaction: abundance × programme coefficient.
   - 单位/尺度：model term
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Section set**
   - 定义：all_sections or the recorded exclusion of QC-fallback sections; not an expression-based exclusion.
   - 单位/尺度：setting
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Datasets (n)**
   - 定义：Datasets represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — 95% CI lower**
   - 定义：Lower endpoint of the recorded 95% confidence interval.
   - 单位/尺度：natural-log link coefficient
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — 95% CI upper**
   - 定义：Upper endpoint of the recorded 95% confidence interval.
   - 单位/尺度：natural-log link coefficient
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — BH FDR (87 tests)**
   - 定义：BH-adjusted P across 29 candidate pairs × 3 model terms within the stated posterior/QC setting. Non-estimable terms retain their recorded P=1 accounting.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 33 — Table 33
最终标题：Supplementary Table 33 | Module spatial coherence
实际数据：112行；Excel第4–115行；12列。
表内说明：Observed and null module coherence with posterior-specific multiple-testing families. One module, coherence metric, null method and posterior setting. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One module, coherence metric, null method and posterior setting.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Posterior summary**
   - 定义：cell2location posterior summary used for abundance: q05 or posterior mean.
   - 单位/尺度：summary label
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Module**
   - 定义：Recorded module for one module, coherence metric, null method and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Null method**
   - 定义：Recorded null method for one module, coherence metric, null method and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Metric**
   - 定义：Recorded metric for one module, coherence metric, null method and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Observed median**
   - 定义：Recorded observed median for one module, coherence metric, null method and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Null median**
   - 定义：Recorded null median for one module, coherence metric, null method and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Observed minus null**
   - 定义：Recorded observed minus null for one module, coherence metric, null method and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — P**
   - 定义：Original P for P in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Datasets (n)**
   - 定义：Datasets represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Analysis units (n)**
   - 定义：Analysis units represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Draws**
   - 定义：Recorded draws for one module, coherence metric, null method and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — BH FDR (28 tests)**
   - 定义：Original adjusted P for the stated testing family; no correction has been recomputed.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 34 — Table 34
最终标题：Supplementary Table 34 | Module–programme direction across datasets
实际数据：168行；Excel第4–171行；10列。
表内说明：Dataset-level direction consistency for six expression programmes and 14 modules. One module–programme pair and posterior setting. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One module–programme pair and posterior setting.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Posterior summary**
   - 定义：cell2location posterior summary used for abundance: q05 or posterior mean.
   - 单位/尺度：summary label
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Module**
   - 定义：Recorded module for one module–programme pair and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Program**
   - 定义：Recorded program for one module–programme pair and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Datasets (n)**
   - 定义：Datasets represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Positive datasets**
   - 定义：Recorded positive datasets for one module–programme pair and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Median rho**
   - 定义：Recorded median rho for one module–programme pair and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Min rho**
   - 定义：Recorded min rho for one module–programme pair and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Max rho**
   - 定义：Recorded max rho for one module–programme pair and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — P direction positive**
   - 定义：Recorded p direction positive for one module–programme pair and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — BH FDR (84 tests)**
   - 定义：Original adjusted P for the stated testing family; no correction has been recomputed.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 35 — Table 35
最终标题：Supplementary Table 35 | Disease-associated cell proportions
实际数据：582行；Excel第4–585行；11列。
表内说明：Estimated disease-associated proportion effects and uncertainty, retaining invalid-uncertainty records. One fitted state–comparison model. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One fitted state–comparison model.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Annotation level**
   - 定义：Recorded annotation level for one fitted state–comparison model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — State**
   - 定义：Recorded state for one fitted state–comparison model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Model id**
   - 定义：Recorded model id for one fitted state–comparison model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Estimate**
   - 定义：Recorded estimate for one fitted state–comparison model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — SE**
   - 定义：Standard error on the corresponding effect/estimate scale; not a confidence interval.
   - 单位/尺度：effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Z**
   - 定义：Recorded z for one fitted state–comparison model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — P**
   - 定义：Original P for P in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Model term**
   - 定义：Disease-versus-reference coefficient in the stated beta-binomial proportion model.
   - 单位/尺度：model term
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — FDR family**
   - 定义：Recorded fDR family for one fitted state–comparison model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Inference status**
   - 定义：Recorded inference status for one fitted state–comparison model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — BH FDR**
   - 定义：Original adjusted P for the stated testing family; no correction has been recomputed.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 36 — Table 36
最终标题：Supplementary Table 36 | Cell-type fractions by tissue
实际数据：12行；Excel第4–15行；17列。
表内说明：Sample-averaged L1 composition within tissue categories. One tissue category; columns are L1 fractions. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One tissue category; columns are L1 fractions.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Tissue**
   - 定义：Recorded tissue for one tissue category; columns are L1 fractions.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — B cell**
   - 定义：Recorded b cell for one tissue category; columns are L1 fractions.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Dendritic cell**
   - 定义：Recorded dendritic cell for one tissue category; columns are L1 fractions.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Endothelial cell**
   - 定义：Recorded endothelial cell for one tissue category; columns are L1 fractions.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Fibroblast**
   - 定义：Recorded fibroblast for one tissue category; columns are L1 fractions.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Fibromyocyte**
   - 定义：Recorded fibromyocyte for one tissue category; columns are L1 fractions.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Macrophage**
   - 定义：Recorded macrophage for one tissue category; columns are L1 fractions.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Mast cell**
   - 定义：Recorded mast cell for one tissue category; columns are L1 fractions.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Monocyte**
   - 定义：Recorded monocyte for one tissue category; columns are L1 fractions.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Natural killer cell**
   - 定义：Recorded natural killer cell for one tissue category; columns are L1 fractions.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Neutrophil**
   - 定义：Recorded neutrophil for one tissue category; columns are L1 fractions.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Pericyte**
   - 定义：Recorded pericyte for one tissue category; columns are L1 fractions.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Plasma cell**
   - 定义：Recorded plasma cell for one tissue category; columns are L1 fractions.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Schwann cell**
   - 定义：Recorded schwann cell for one tissue category; columns are L1 fractions.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Smooth muscle cell**
   - 定义：Recorded smooth muscle cell for one tissue category; columns are L1 fractions.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — T cell**
   - 定义：Recorded t cell for one tissue category; columns are L1 fractions.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — Unknown**
   - 定义：Recorded unknown for one tissue category; columns are L1 fractions.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 37 — Table 37
最终标题：Supplementary Table 37 | Recurrent pseudobulk pathways
实际数据：3,728行；Excel第4–3731行；5列。
表内说明：Descriptive recurrence of pathways across significant pseudobulk enrichment analyses. One recurrent pathway and database. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One recurrent pathway and database.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Pathway**
   - 定义：Recorded pathway for one recurrent pathway and database.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Database**
   - 定义：Recorded database for one recurrent pathway and database.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Significant analysis units**
   - 定义：Recorded significant analysis units for one recurrent pathway and database.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Median NES**
   - 定义：Recorded median NES for one recurrent pathway and database.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Minimum native FDR**
   - 定义：Original adjusted P for the stated testing family; no correction has been recomputed.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 38 — Table 38
最终标题：Supplementary Table 38 | Cell–cell communication effects
实际数据：176行；Excel第4–179行；43列。
表内说明：Matched LIANA and CellChat disease effects and cohort support for concordantly supported interactions. One disease-context interaction with paired method-specific results. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One disease-context interaction with paired method-specific results.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — LIANA method**
   - 定义：LIANA: Recorded method for one disease-context interaction with paired method-specific results.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Disease context**
   - 定义：Recorded disease context for one disease-context interaction with paired method-specific results.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Eligibility scheme**
   - 定义：Recorded eligibility scheme for one disease-context interaction with paired method-specific results.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Interaction ID**
   - 定义：Recorded Interaction ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Sender state**
   - 定义：Recorded sender state for one disease-context interaction with paired method-specific results.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Receiver state**
   - 定义：Recorded receiver state for one disease-context interaction with paired method-specific results.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Ligand**
   - 定义：Recorded ligand for one disease-context interaction with paired method-specific results.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Receptor**
   - 定义：Recorded receptor for one disease-context interaction with paired method-specific results.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — LIANA Signalling category**
   - 定义：LIANA: Recorded signalling category for one disease-context interaction with paired method-specific results.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — LIANA Effect (disease − control)**
   - 定义：LIANA: Recorded effect (disease − control) for one disease-context interaction with paired method-specific results.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — LIANA Pooled SE**
   - 定义：LIANA: Standard error on the corresponding effect/estimate scale; not a confidence interval.
   - 单位/尺度：effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — LIANA 95% CI lower**
   - 定义：LIANA: Recorded confidence bound on the corresponding effect scale; no bounds were reconstructed.
   - 单位/尺度：effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — LIANA 95% CI upper**
   - 定义：LIANA: Recorded confidence bound on the corresponding effect scale; no bounds were reconstructed.
   - 单位/尺度：effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — LIANA Pooled P**
   - 定义：Original P for LIANAooled in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — LIANA Cohorts (n)**
   - 定义：LIANA Cohorts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — LIANA Studies (n)**
   - 定义：LIANA Studies represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — LIANA Tissues (n)**
   - 定义：LIANA Tissues represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
18. **R — LIANA Direction agreement**
   - 定义：LIANA: Fraction of contributing effects with the recorded consensus direction; the qualifying units are defined by this table.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
19. **S — LIANA Largest cohort effect share**
   - 定义：Recorded liana largest cohort effect share, with its analysis-specific denominator; no conversion to percent.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
20. **T — LIANA I² (%)**
   - 定义：Between-unit heterogeneity statistic on its original percentage scale.
   - 单位/尺度：percent (0–100)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
21. **U — LIANA τ²**
   - 定义：Estimated between-unit effect variance.
   - 单位/尺度：squared effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
22. **V — LIANA Cohort contributions**
   - 定义：LIANA: Recorded cohort contributions for one disease-context interaction with paired method-specific results.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
23. **W — LIANA BH FDR**
   - 定义：Original adjusted P for the stated testing family; no correction has been recomputed.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
24. **X — LIANA Context support**
   - 定义：LIANA: Recorded context support for one disease-context interaction with paired method-specific results.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
25. **Y — LIANA evidence class**
   - 定义：LIANA: Recorded evidence class for one disease-context interaction with paired method-specific results.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
26. **Z — CellChat method**
   - 定义：CellChat: Recorded method for one disease-context interaction with paired method-specific results.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
27. **AA — CellChat Signalling category**
   - 定义：CellChat: Recorded signalling category for one disease-context interaction with paired method-specific results.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
28. **AB — CellChat Effect (disease − control)**
   - 定义：CellChat: Recorded effect (disease − control) for one disease-context interaction with paired method-specific results.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
29. **AC — CellChat Pooled SE**
   - 定义：CellChat: Standard error on the corresponding effect/estimate scale; not a confidence interval.
   - 单位/尺度：effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
30. **AD — CellChat 95% CI lower**
   - 定义：CellChat: Recorded confidence bound on the corresponding effect scale; no bounds were reconstructed.
   - 单位/尺度：effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
31. **AE — CellChat 95% CI upper**
   - 定义：CellChat: Recorded confidence bound on the corresponding effect scale; no bounds were reconstructed.
   - 单位/尺度：effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
32. **AF — CellChat Pooled P**
   - 定义：Original P for CellChatooled in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
33. **AG — CellChat Cohorts (n)**
   - 定义：CellChat Cohorts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
34. **AH — CellChat Studies (n)**
   - 定义：CellChat Studies represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
35. **AI — CellChat Tissues (n)**
   - 定义：CellChat Tissues represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
36. **AJ — CellChat Direction agreement**
   - 定义：CellChat: Fraction of contributing effects with the recorded consensus direction; the qualifying units are defined by this table.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
37. **AK — CellChat Largest cohort effect share**
   - 定义：Recorded cellchat largest cohort effect share, with its analysis-specific denominator; no conversion to percent.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
38. **AL — CellChat I² (%)**
   - 定义：Between-unit heterogeneity statistic on its original percentage scale.
   - 单位/尺度：percent (0–100)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
39. **AM — CellChat τ²**
   - 定义：Estimated between-unit effect variance.
   - 单位/尺度：squared effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
40. **AN — CellChat Cohort contributions**
   - 定义：CellChat: Recorded cohort contributions for one disease-context interaction with paired method-specific results.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
41. **AO — CellChat BH FDR**
   - 定义：Original adjusted P for the stated testing family; no correction has been recomputed.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
42. **AP — CellChat Context support**
   - 定义：CellChat: Recorded context support for one disease-context interaction with paired method-specific results.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
43. **AQ — CellChat evidence class**
   - 定义：CellChat: Recorded evidence class for one disease-context interaction with paired method-specific results.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 39 — Table 39
最终标题：Supplementary Table 39 | Communication consensus by disease context
实际数据：6行；Excel第4–9行；5列。
表内说明：Agreement between communication methods within each tissue and disease background. One tissue–disease-background combination. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One tissue–disease-background combination.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Background**
   - 定义：Recorded background for one tissue–disease-background combination.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Tissue**
   - 定义：Recorded tissue for one tissue–disease-background combination.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Both measured**
   - 定义：Recorded both measured for one tissue–disease-background combination.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Same direction**
   - 定义：Recorded same direction for one tissue–disease-background combination.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Both supported same direction**
   - 定义：Recorded both supported same direction for one tissue–disease-background combination.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 40 — Table 40
最终标题：Supplementary Table 40 | GWAS–eQTL colocalisation
实际数据：418行；Excel第4–421行；26列。
表内说明：Complete ABF colocalisation tests, including non-supporting results, for prioritised genes. One gene, GWAS scenario and QTL context. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One gene, GWAS scenario and QTL context.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Gene**
   - 定义：Reported gene symbol; original gene identities are unchanged.
   - 单位/尺度：gene symbol
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — GWAS scenario**
   - 定义：Recorded gWAS scenario for one gene, GWAS scenario and QTL context.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Ensembl gene ID**
   - 定义：Stable source gene identifier; not a row number or expression value.
   - 单位/尺度：Ensembl identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Nominated lineage**
   - 定义：Recorded nominated lineage for one gene, GWAS scenario and QTL context.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Context**
   - 定义：Recorded context for one gene, GWAS scenario and QTL context.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — QTL source**
   - 定义：Recorded qTL source for one gene, GWAS scenario and QTL context.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — QTL sample size**
   - 定义：Recorded qTL sample size for one gene, GWAS scenario and QTL context.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Shared variants (n)**
   - 定义：Shared variants represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — GWAS effective N**
   - 定义：Recorded gWAS effective N for one gene, GWAS scenario and QTL context.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Minimum GWAS P**
   - 定义：Original P for Minimum GWAS in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Minimum eQTL P**
   - 定义：Original P for Minimum eQTL in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — PP.H0**
   - 定义：Colocalisation hypothesis posterior; H0 neither, H1 GWAS only, H2 QTL only, H3 distinct signals, H4 shared signal.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — PP.H1**
   - 定义：Colocalisation hypothesis posterior; H0 neither, H1 GWAS only, H2 QTL only, H3 distinct signals, H4 shared signal.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — PP.H2**
   - 定义：Colocalisation hypothesis posterior; H0 neither, H1 GWAS only, H2 QTL only, H3 distinct signals, H4 shared signal.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — PP.H3**
   - 定义：Colocalisation hypothesis posterior; H0 neither, H1 GWAS only, H2 QTL only, H3 distinct signals, H4 shared signal.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — PP.H4**
   - 定义：Colocalisation hypothesis posterior; H0 neither, H1 GWAS only, H2 QTL only, H3 distinct signals, H4 shared signal.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — PP.H4 / (PP.H3+PP.H4)**
   - 定义：Posterior probability of a shared signal conditional on either a distinct or shared two-trait signal.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
18. **R — Top shared variant**
   - 定义：Recorded top shared variant, with its analysis-specific denominator; no conversion to percent.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
19. **S — Top-variant SNP.PP.H4**
   - 定义：Variant-level posterior conditional on H4; not the region-level PP.H4.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
20. **T — Top gwas beta**
   - 定义：Recorded top gwas beta for one gene, GWAS scenario and QTL context.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
21. **U — Top eqtl beta aligned**
   - 定义：Recorded top eqtl beta aligned for one gene, GWAS scenario and QTL context.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
22. **V — Top direction interpretable**
   - 定义：Recorded top direction interpretable for one gene, GWAS scenario and QTL context.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
23. **W — Top effect direction**
   - 定义：Recorded top effect direction for one gene, GWAS scenario and QTL context.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
24. **X — Evidence class**
   - 定义：Recorded evidence class for one gene, GWAS scenario and QTL context.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
25. **Y — Primary prior p12**
   - 定义：Recorded primary prior p12 for one gene, GWAS scenario and QTL context.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
26. **Z — Method**
   - 定义：Recorded method for one gene, GWAS scenario and QTL context.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 41 — Table 41
最终标题：Supplementary Table 41 | QTL coverage by gene and GWAS scenario
实际数据：58行；Excel第4–61行；15列。
表内说明：Requested, available, missing and tested QTL contexts with evidence summaries. One gene and GWAS scenario. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One gene and GWAS scenario.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Gene**
   - 定义：Reported gene symbol; original gene identities are unchanged.
   - 单位/尺度：gene symbol
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — GWAS scenario**
   - 定义：Recorded gWAS scenario for one gene and GWAS scenario.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Requested contexts (n)**
   - 定义：Requested contexts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Available QTL contexts (n)**
   - 定义：Available QTL contexts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Missing QTL contexts (n)**
   - 定义：Missing QTL contexts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Retrieval failures (n)**
   - 定义：Retrieval failures represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Tested contexts (n)**
   - 定义：Tested contexts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Strong contexts (n)**
   - 定义：Strong contexts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Suggestive contexts (n)**
   - 定义：Suggestive contexts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Distinct-signal contexts (n)**
   - 定义：Distinct-signal contexts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Maximum PP.H4**
   - 定义：Recorded maximum PP.H4 for one gene and GWAS scenario.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Maximum H4 given H3H4**
   - 定义：Recorded maximum H4 given H3H4 for one gene and GWAS scenario.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Best context**
   - 定义：Recorded best context for one gene and GWAS scenario.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Best qtl source**
   - 定义：Recorded best qtl source for one gene and GWAS scenario.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Best classification**
   - 定义：Recorded best classification for one gene and GWAS scenario.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 42 — Table 42
最终标题：Supplementary Table 42 | Candidate MAGMA evidence across scenarios
实际数据：87行；Excel第4–90行；7列。
表内说明：Gene-level MAGMA evidence for the candidate pairs across three GWAS scenarios. One candidate pair and GWAS scenario. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One candidate pair and GWAS scenario.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Gene**
   - 定义：Reported gene symbol; original gene identities are unchanged.
   - 单位/尺度：gene symbol
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Level**
   - 定义：Recorded level for one candidate pair and GWAS scenario.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — State**
   - 定义：Recorded state for one candidate pair and GWAS scenario.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — GWAS scenario**
   - 定义：Recorded gWAS scenario for one candidate pair and GWAS scenario.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — P**
   - 定义：Original P for P in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — FDR**
   - 定义：Original adjusted P for the stated testing family; no correction has been recomputed.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Bonferroni significant**
   - 定义：Recorded bonferroni significant for one candidate pair and GWAS scenario.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 43 — Table 43
最终标题：Supplementary Table 43 | Pseudobulk meta-analysis of prioritised genes
实际数据：15,639行；Excel第4–15642行；19列。
表内说明：All 15,639 evaluable gene–state pairs with REML effects and target-set FDR, not only selected pairs. One evaluable priority-gene–cell-population pair. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One evaluable priority-gene–cell-population pair.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Gene**
   - 定义：Reported gene symbol; original gene identities are unchanged.
   - 单位/尺度：gene symbol
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Annotation level**
   - 定义：Recorded annotation level for one evaluable priority-gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Cell state**
   - 定义：Recorded cell state for one evaluable priority-gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Disease contexts (n)**
   - 定义：Disease contexts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — GWAS-aligned contexts (n)**
   - 定义：GWAS-aligned contexts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Disease contexts**
   - 定义：Recorded disease contexts for one evaluable priority-gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Meta log2FC**
   - 定义：REML random-effects meta-analysis log2 disease-versus-reference fold change.
   - 单位/尺度：log2 fold change
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Meta SE**
   - 定义：Standard error of the REML meta-analysis log2 fold change.
   - 单位/尺度：log2 fold change
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Meta Z**
   - 定义：Recorded meta Z for one evaluable priority-gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Meta P**
   - 定义：Original P for Meta in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Meta 95% CI lower**
   - 定义：Recorded confidence bound on the corresponding effect scale; no bounds were reconstructed.
   - 单位/尺度：effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Meta 95% CI upper**
   - 定义：Recorded confidence bound on the corresponding effect scale; no bounds were reconstructed.
   - 单位/尺度：effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — I² (%)**
   - 定义：Between-unit heterogeneity statistic on its original percentage scale.
   - 单位/尺度：percent (0–100)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Direction agreement**
   - 定义：Fraction of contributing effects with the recorded consensus direction; the qualifying units are defined by this table.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Consensus direction**
   - 定义：Recorded consensus direction for one evaluable priority-gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — Significant contexts (n)**
   - 定义：Significant contexts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — Median context log2FC**
   - 定义：Recorded median context log2FC for one evaluable priority-gene–cell-population pair.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
18. **R — Ensembl gene ID**
   - 定义：Stable source gene identifier; not a row number or expression value.
   - 单位/尺度：Ensembl identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
19. **S — Target-set BH FDR**
   - 定义：BH-adjusted meta-analysis P within cell state over the prespecified priority-gene set; distinct from all-gene FDR.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 44 — Table 44
最终标题：Supplementary Table 44 | CoVarNet rank selection
实际数据：24行；Excel第4–27行；31列。
表内说明：NMF quality metrics and the recorded automatic rank-selection criteria for ranks 2–25. One tested NMF rank. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One tested NMF rank.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Rank label**
   - 定义：Recorded rank label for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Method**
   - 定义：Recorded method for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Seed**
   - 定义：Recorded seed for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Rng**
   - 定义：Recorded rng for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Metric**
   - 定义：Recorded metric for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Rank**
   - 定义：Recorded rank for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Sparseness.basis**
   - 定义：Recorded sparseness.basis for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Sparseness.coef**
   - 定义：Recorded sparseness.coef for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Rss**
   - 定义：Recorded rss for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Evar**
   - 定义：Recorded evar for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Silhouette.coef**
   - 定义：Recorded silhouette.coef for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Silhouette.basis**
   - 定义：Recorded silhouette.basis for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Residuals**
   - 定义：Recorded residuals for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Niter**
   - 定义：Recorded niter for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Cpu**
   - 定义：Recorded cpu for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — Cpu.all**
   - 定义：Recorded cpu.all for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — Nrun**
   - 定义：Recorded nrun for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
18. **R — Cophenetic**
   - 定义：Recorded cophenetic for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
19. **S — Dispersion**
   - 定义：Recorded dispersion for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
20. **T — Silhouette.consensus**
   - 定义：Recorded silhouette.consensus for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
21. **U — Scaled cophenetic**
   - 定义：Recorded scaled cophenetic for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
22. **V — Scaled dispersion**
   - 定义：Recorded scaled dispersion for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
23. **W — Scaled silhouette consensus**
   - 定义：Recorded scaled silhouette consensus for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
24. **X — Composite stability quality**
   - 定义：Recorded composite stability quality for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
25. **Y — Rss elbow rank**
   - 定义：Recorded rss elbow rank for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
26. **Z — Within rss elbow window**
   - 定义：Recorded within rss elbow window for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
27. **AA — Stable by cophenetic**
   - 定义：Recorded stable by cophenetic for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
28. **AB — Automatic candidate**
   - 定义：Recorded automatic candidate for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
29. **AC — Selected automatically**
   - 定义：Recorded selected automatically for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
30. **AD — Selected for primary analysis**
   - 定义：Recorded selected for primary analysis for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
31. **AE — Selection mode**
   - 定义：Recorded selection mode for one tested NMF rank.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 45 — Table 45
最终标题：Supplementary Table 45 | Leave-one-study CoVarNet refits
实际数据：700行；Excel第4–703行；9列。
表内说明：Module matching under study omission, with similarity, overlap and matching diagnostics. One reference-module match in a study-omission refit. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One reference-module match in a study-omission refit.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Reference module**
   - 定义：Recorded reference module for one reference-module match in a study-omission refit.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Run type**
   - 定义：Recorded run type for one reference-module match in a study-omission refit.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Omitted or bootstrap**
   - 定义：Recorded omitted or bootstrap for one reference-module match in a study-omission refit.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Matching mode**
   - 定义：Recorded matching mode for one reference-module match in a study-omission refit.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Matched module**
   - 定义：Recorded matched module for one reference-module match in a study-omission refit.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Cosine similarity**
   - 定义：Recorded cosine similarity for one reference-module match in a study-omission refit.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Top node jaccard**
   - 定义：Recorded top node jaccard for one reference-module match in a study-omission refit.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — N reference modules matching candidate**
   - 定义：Recorded n reference modules matching candidate for one reference-module match in a study-omission refit.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Duplicate match**
   - 定义：Recorded duplicate match for one reference-module match in a study-omission refit.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 46 — Table 46
最终标题：Supplementary Table 46 | Multi-signal colocalisation sensitivity
实际数据：8行；Excel第4–11行；15列。
表内说明：Historical-scope multi-signal results and fit statuses; this is not a complete matched-LD screen of the current 29 candidates. One historical gene–QTL-context sensitivity test. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One historical gene–QTL-context sensitivity test.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Analysis scenario**
   - 定义：Recorded analysis scenario for one historical gene–QTL-context sensitivity test.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Gene**
   - 定义：Reported gene symbol; original gene identities are unchanged.
   - 单位/尺度：gene symbol
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Context**
   - 定义：Recorded context for one historical gene–QTL-context sensitivity test.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Primary approximate-Bayes-factor PP.H4**
   - 定义：Recorded primary approximate-Bayes-factor PP.H4 for one historical gene–QTL-context sensitivity test.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Primary approximate-Bayes-factor classification**
   - 定义：Recorded primary approximate-Bayes-factor classification for one historical gene–QTL-context sensitivity test.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Variants (n)**
   - 定义：Variants represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — GWAS credible sets (n)**
   - 定义：GWAS credible sets represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — QTL credible sets (n)**
   - 定义：QTL credible sets represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Maximum SuSiE PP.H4**
   - 定义：Recorded maximum SuSiE PP.H4 for one historical gene–QTL-context sensitivity test.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Best GWAS signal**
   - 定义：Recorded best GWAS signal for one historical gene–QTL-context sensitivity test.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Best QTL signal**
   - 定义：Recorded best QTL signal for one historical gene–QTL-context sensitivity test.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Analysis status**
   - 定义：Recorded analysis status for one historical gene–QTL-context sensitivity test.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Error or exclusion reason**
   - 定义：Recorded error or exclusion reason for one historical gene–QTL-context sensitivity test.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Conservative interpretation**
   - 定义：Recorded conservative interpretation for one historical gene–QTL-context sensitivity test.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Current scope note**
   - 定义：Recorded current scope note for one historical gene–QTL-context sensitivity test.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 47 — Table 47
最终标题：Supplementary Table 47 | Analysis sample sizes and units
实际数据：9行；Excel第4–12行；4列。
表内说明：The distinct cell, sample, specimen, profile and spatial analysis units used by each analysis. One analysis-specific unit/count definition. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One analysis-specific unit/count definition.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Analysis**
   - 定义：Recorded analysis for one analysis-specific unit/count definition.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Unit**
   - 定义：Recorded unit for one analysis-specific unit/count definition.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Count**
   - 定义：Recorded count for one analysis-specific unit/count definition.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Scope**
   - 定义：Analysis-specific scope of the associated count or result.
   - 单位/尺度：text
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 48 — Table 48
最终标题：Supplementary Table 48 | Spatial posterior and QC sensitivity
实际数据：348行；Excel第4–351行；17列。
表内说明：Candidate model terms under four posterior/QC settings, each retaining its 87-test BH family. One candidate pair, model term and posterior/QC setting. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One candidate pair, model term and posterior/QC setting.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Gene**
   - 定义：Reported gene symbol; original gene identities are unchanged.
   - 单位/尺度：gene symbol
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Cell context**
   - 定义：Recorded cell context for one candidate pair, model term and posterior/QC setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Units (n)**
   - 定义：Number of contributing analysis units in the spatial meta-analysis, not number of spots or gene–unit estimates.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — β**
   - 定义：Conditional spatial log-link coefficient for the stated model term.
   - 单位/尺度：natural-log link coefficient
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — SE**
   - 定义：Standard error of the conditional spatial coefficient; not a confidence interval.
   - 单位/尺度：natural-log link coefficient
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Z**
   - 定义：Recorded z for one candidate pair, model term and posterior/QC setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — P**
   - 定义：Original P for P in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — τ²**
   - 定义：Estimated between-unit effect variance.
   - 单位/尺度：squared effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — I² (%)**
   - 定义：Between-unit heterogeneity statistic on its original percentage scale.
   - 单位/尺度：percent (0–100)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Posterior summary**
   - 定义：cell2location posterior summary used for abundance: q05 or posterior mean.
   - 单位/尺度：summary label
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Expression programme**
   - 定义：Molecular expression signature used in the spatial model; not an independently assigned histopathological region.
   - 单位/尺度：programme name
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Model term**
   - 定义：cell_abundance: conditional abundance coefficient; pathology: expression-programme coefficient; interaction: abundance × programme coefficient.
   - 单位/尺度：model term
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Section set**
   - 定义：all_sections or the recorded exclusion of QC-fallback sections; not an expression-based exclusion.
   - 单位/尺度：setting
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Datasets (n)**
   - 定义：Datasets represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — 95% CI lower**
   - 定义：Lower endpoint of the recorded 95% confidence interval.
   - 单位/尺度：natural-log link coefficient
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — 95% CI upper**
   - 定义：Upper endpoint of the recorded 95% confidence interval.
   - 单位/尺度：natural-log link coefficient
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — BH FDR (87 tests)**
   - 定义：BH-adjusted P across 29 candidate pairs × 3 model terms within the stated posterior/QC setting. Non-estimable terms retain their recorded P=1 accounting.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 49 — Table 49
最终标题：Supplementary Table 49 | Cohort contributions to communication effects
实际数据：5,506行；Excel第4–5509行；15列。
表内说明：Pooled interactions with their recorded cohort-effect contributions; the contribution field contains multiple cohorts per row. One method–context interaction, not one independent cohort. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One method–context interaction, not one independent cohort.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Method**
   - 定义：Recorded method for one method–context interaction, not one independent cohort.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Disease context**
   - 定义：Recorded disease context for one method–context interaction, not one independent cohort.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Eligibility scheme**
   - 定义：Recorded eligibility scheme for one method–context interaction, not one independent cohort.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Interaction ID**
   - 定义：Recorded Interaction ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Sender state**
   - 定义：Recorded sender state for one method–context interaction, not one independent cohort.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Receiver state**
   - 定义：Recorded receiver state for one method–context interaction, not one independent cohort.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Ligand**
   - 定义：Recorded ligand for one method–context interaction, not one independent cohort.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Receptor**
   - 定义：Recorded receptor for one method–context interaction, not one independent cohort.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Effect (disease − control)**
   - 定义：Recorded effect (disease − control) for one method–context interaction, not one independent cohort.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Cohorts (n)**
   - 定义：Cohorts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Direction agreement**
   - 定义：Fraction of contributing effects with the recorded consensus direction; the qualifying units are defined by this table.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Largest cohort effect share**
   - 定义：Recorded largest cohort effect share, with its analysis-specific denominator; no conversion to percent.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Cohort contributions**
   - 定义：Recorded cohort contributions for one method–context interaction, not one independent cohort.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — BH FDR**
   - 定义：Original adjusted P for the stated testing family; no correction has been recomputed.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Context support**
   - 定义：Recorded context support for one method–context interaction, not one independent cohort.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 50 — Table 50
最终标题：Supplementary Table 50 | Relationships among evidence sources
实际数据：7行；Excel第4–10行；3列。
表内说明：Shared inputs and interpretation boundaries between evidence streams. One relationship between analysis/evidence sources. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One relationship between analysis/evidence sources.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Analysis**
   - 定义：Recorded analysis for one relationship between analysis/evidence sources.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Shared evidence**
   - 定义：Recorded shared evidence, with its analysis-specific denominator; no conversion to percent.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Interpretation boundary**
   - 定义：Recorded interpretation boundary for one relationship between analysis/evidence sources.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 51 — Table 51
最终标题：Supplementary Table 51 | Multiple-testing families
实际数据：6行；Excel第4–9行；4列。
表内说明：Prespecified test families and correction boundaries across analyses. One multiple-testing family definition. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One multiple-testing family definition.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Analysis**
   - 定义：Recorded analysis for one multiple-testing family definition.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Multiplicity family**
   - 定义：Recorded multiplicity family for one multiple-testing family definition.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Correction**
   - 定义：Recorded correction for one multiple-testing family definition.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Unit / boundary**
   - 定义：Recorded unit / boundary for one multiple-testing family definition.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 52 — Table 52
最终标题：Supplementary Table 52 | Genetic cell-group enrichment summaries
实际数据：18行；Excel第4–21行；6列。
表内说明：Counts of tests and significant groups by genetic scenario, method, annotation level and weighting. One scenario–method–level–weighting summary. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One scenario–method–level–weighting summary.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Annotation level**
   - 定义：Recorded annotation level for one scenario–method–level–weighting summary.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — GWAS scenario**
   - 定义：Recorded gWAS scenario for one scenario–method–level–weighting summary.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Weighting**
   - 定义：Recorded weighting for one scenario–method–level–weighting summary.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Tests (n)**
   - 定义：Tests represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — FDR <0.05 (n)**
   - 定义：FDR <0.05 represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Method**
   - 定义：Recorded method for one scenario–method–level–weighting summary.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 53 — Table 53
最终标题：Supplementary Table 53 | Target-exclusion spatial sensitivity
实际数据：348行；Excel第4–351行；33列。
表内说明：Paired original-reference and target-excluded spatial results, retaining the original model terms and settings. One candidate pair, model term and posterior/QC setting. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One candidate pair, model term and posterior/QC setting.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Gene**
   - 定义：Reported gene symbol; original gene identities are unchanged.
   - 单位/尺度：gene symbol
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Cell context**
   - 定义：Recorded cell context for one candidate pair, model term and posterior/QC setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Original-reference Units (n)**
   - 定义：Original-reference Units represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Original-reference β**
   - 定义：Original-reference model. Conditional spatial log-link coefficient for the stated model term.
   - 单位/尺度：natural-log link coefficient
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Original-reference SE**
   - 定义：Original-reference model. Standard error of the conditional spatial coefficient; not a confidence interval.
   - 单位/尺度：natural-log link coefficient
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Original-reference Z**
   - 定义：Original-reference model. Recorded z for one candidate pair, model term and posterior/QC setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Original-reference P**
   - 定义：Original P for Original-reference in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Original-reference τ²**
   - 定义：Estimated between-unit effect variance.
   - 单位/尺度：squared effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Original-reference I² (%)**
   - 定义：Between-unit heterogeneity statistic on its original percentage scale.
   - 单位/尺度：percent (0–100)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Posterior summary**
   - 定义：cell2location posterior summary used for abundance: q05 or posterior mean.
   - 单位/尺度：summary label
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Expression programme**
   - 定义：Molecular expression signature used in the spatial model; not an independently assigned histopathological region.
   - 单位/尺度：programme name
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Model term**
   - 定义：cell_abundance: conditional abundance coefficient; pathology: expression-programme coefficient; interaction: abundance × programme coefficient.
   - 单位/尺度：model term
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Section set**
   - 定义：all_sections or the recorded exclusion of QC-fallback sections; not an expression-based exclusion.
   - 单位/尺度：setting
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Original-reference Datasets (n)**
   - 定义：Original-reference Datasets represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Original-reference 95% CI lower**
   - 定义：Original-reference model. Lower endpoint of the recorded 95% confidence interval.
   - 单位/尺度：natural-log link coefficient
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — Original-reference 95% CI upper**
   - 定义：Original-reference model. Upper endpoint of the recorded 95% confidence interval.
   - 单位/尺度：natural-log link coefficient
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — Original-reference BH FDR (87 tests)**
   - 定义：Original adjusted P for the stated testing family; no correction has been recomputed.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
18. **R — Target-excluded Units (n)**
   - 定义：Target-excluded Units represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
19. **S — Target-excluded β**
   - 定义：Target-excluded model. Conditional spatial log-link coefficient for the stated model term.
   - 单位/尺度：natural-log link coefficient
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
20. **T — Target-excluded SE**
   - 定义：Target-excluded model. Standard error of the conditional spatial coefficient; not a confidence interval.
   - 单位/尺度：natural-log link coefficient
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
21. **U — Target-excluded Z**
   - 定义：Target-excluded model. Recorded z for one candidate pair, model term and posterior/QC setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
22. **V — Target-excluded P**
   - 定义：Original P for Target-excluded in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
23. **W — Target-excluded τ²**
   - 定义：Estimated between-unit effect variance.
   - 单位/尺度：squared effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
24. **X — Target-excluded I² (%)**
   - 定义：Between-unit heterogeneity statistic on its original percentage scale.
   - 单位/尺度：percent (0–100)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
25. **Y — Target-excluded Datasets (n)**
   - 定义：Target-excluded Datasets represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
26. **Z — Target-excluded 95% CI lower**
   - 定义：Target-excluded model. Lower endpoint of the recorded 95% confidence interval.
   - 单位/尺度：natural-log link coefficient
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
27. **AA — Target-excluded 95% CI upper**
   - 定义：Target-excluded model. Upper endpoint of the recorded 95% confidence interval.
   - 单位/尺度：natural-log link coefficient
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
28. **AB — Target-excluded BH FDR (87 tests)**
   - 定义：Original adjusted P for the stated testing family; no correction has been recomputed.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
29. **AC — Merge**
   - 定义：Recorded merge for one candidate pair, model term and posterior/QC setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
30. **AD — Beta difference**
   - 定义：Recorded beta difference for one candidate pair, model term and posterior/QC setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
31. **AE — Sign changed**
   - 定义：Recorded sign changed for one candidate pair, model term and posterior/QC setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
32. **AF — Original-reference FDR <0.05**
   - 定义：Boolean indicator of BH FDR <0.05 for the stated reference implementation, not a count.
   - 单位/尺度：boolean
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
33. **AG — Target-excluded FDR <0.05**
   - 定义：Boolean indicator of BH FDR <0.05 for the stated reference implementation, not a count.
   - 单位/尺度：boolean
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 54 — Table 54
最终标题：Supplementary Table 54 | cell2location model diagnostics
实际数据：7行；Excel第4–10行；10列。
表内说明：Recorded numerical/output diagnostics by spatial dataset; these do not alone establish convergence. One spatial dataset. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One spatial dataset.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Study ID**
   - 定义：Source-study accession or its canonical identifier.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Sections (n)**
   - 定义：Sections represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Spots (n)**
   - 定义：Spots represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — N model genes**
   - 定义：Recorded n model genes for one spatial dataset.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — N factors**
   - 定义：Recorded n factors for one spatial dataset.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Quantile-order check**
   - 定义：Recorded quantile-order check for one spatial dataset.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Finite/nonnegative check**
   - 定义：Recorded finite/nonnegative check for one spatial dataset.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Mean outside q05–q95 (n)**
   - 定义：Mean outside q05–q95 represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Spot-identity check**
   - 定义：Recorded spot-identity check for one spatial dataset.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Source path**
   - 定义：Recorded source path for one spatial dataset.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 55 — Table 55
最终标题：Supplementary Table 55 | Parent-lineage-matched CoVarNet null results
实际数据：70行；Excel第4–73行；6列。
表内说明：Matched-null comparisons for the original module-structure statistics. One module and evaluated statistic. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One module and evaluated statistic.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Module**
   - 定义：Recorded module for one module and evaluated statistic.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Metric**
   - 定义：Recorded metric for one module and evaluated statistic.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Observed**
   - 定义：Recorded observed for one module and evaluated statistic.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Empirical P**
   - 定义：Original P for Empirical in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — N finite null**
   - 定义：Recorded n finite null for one module and evaluated statistic.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — FDR**
   - 定义：Original adjusted P for the stated testing family; no correction has been recomputed.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 56 — Table 56
最终标题：Supplementary Table 56 | Spatial and single-cell source overlap
实际数据：7行；Excel第4–10行；8列。
表内说明：Dataset-level source-overlap records; spatial grouping acceptance does not resolve cross-modality donor overlap. One spatial source-overlap record. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One spatial source-overlap record.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — GEO accession**
   - 定义：Recorded GEO accession; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Source-overlap category**
   - 定义：Recorded source-overlap category for one spatial source-overlap record.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Related atlas accession**
   - 定义：Recorded Related atlas accession; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Donor-overlap status**
   - 定义：Recorded donor-overlap status for one spatial source-overlap record.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Analysis group**
   - 定义：Recorded analysis group for one spatial source-overlap record.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — No. datasets**
   - 定义：Recorded no. datasets for one spatial source-overlap record.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — No. sections**
   - 定义：Recorded no. sections for one spatial source-overlap record.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Analysis units (n)**
   - 定义：Analysis units represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 57 — Table 57
最终标题：Supplementary Table 57 | Unit-level spatial coefficients
实际数据：8,574行；Excel第4–8577行；12列。
表内说明：Within-unit pooled coefficients and SE used in the candidate meta-analysis; SE is not a confidence interval. One candidate pair, model term, analysis unit and posterior/QC setting. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One candidate pair, model term, analysis unit and posterior/QC setting.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Gene**
   - 定义：Reported gene symbol; original gene identities are unchanged.
   - 单位/尺度：gene symbol
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Cell context**
   - 定义：Recorded cell context for one candidate pair, model term, analysis unit and posterior/QC setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Posterior summary**
   - 定义：cell2location posterior summary used for abundance: q05 or posterior mean.
   - 单位/尺度：summary label
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Expression programme**
   - 定义：Molecular expression signature used in the spatial model; not an independently assigned histopathological region.
   - 单位/尺度：programme name
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Model term**
   - 定义：cell_abundance: conditional abundance coefficient; pathology: expression-programme coefficient; interaction: abundance × programme coefficient.
   - 单位/尺度：model term
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Study ID**
   - 定义：Source-study accession or its canonical identifier.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Analysis unit ID**
   - 定义：Recorded grouping used to combine section estimates before between-unit meta-analysis; the grouping is retained as approved.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Section set**
   - 定义：all_sections or the recorded exclusion of QC-fallback sections; not an expression-based exclusion.
   - 单位/尺度：setting
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — β**
   - 定义：Conditional spatial log-link coefficient for the stated model term.
   - 单位/尺度：natural-log link coefficient
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — SE**
   - 定义：Standard error of the conditional spatial coefficient; not a confidence interval.
   - 单位/尺度：natural-log link coefficient
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Sections (n)**
   - 定义：Sections represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Section covariance unmodelled**
   - 定义：Recorded flag that within-unit section residual covariance was not explicitly estimated. It is not a grouping-error flag.
   - 单位/尺度：boolean
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 58 — Table 58
最终标题：Supplementary Table 58 | Section-level spatial coefficients
实际数据：6,132行；Excel第4–6135行；22列。
表内说明：Section model coefficients, uncertainty, eligibility and expression coverage. One candidate pair, model term, section and posterior setting. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One candidate pair, model term, section and posterior setting.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Gene**
   - 定义：Reported gene symbol; original gene identities are unchanged.
   - 单位/尺度：gene symbol
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Cell context**
   - 定义：Recorded cell context for one candidate pair, model term, section and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Study ID**
   - 定义：Source-study accession or its canonical identifier.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Section ID**
   - 定义：Recorded Section ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Analysis unit ID**
   - 定义：Recorded grouping used to combine section estimates before between-unit meta-analysis; the grouping is retained as approved.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Metadata status**
   - 定义：Source metadata status retained from model execution. It does not supersede author confirmation of the current analysis-unit grouping.
   - 单位/尺度：source status
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — QC fallback**
   - 定义：Recorded qC fallback for one candidate pair, model term, section and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Posterior summary**
   - 定义：cell2location posterior summary used for abundance: q05 or posterior mean.
   - 单位/尺度：summary label
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Expression programme**
   - 定义：Molecular expression signature used in the spatial model; not an independently assigned histopathological region.
   - 单位/尺度：programme name
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Expected direction**
   - 定义：Recorded expected direction for one candidate pair, model term, section and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Target removed from program**
   - 定义：Recorded target removed from program for one candidate pair, model term, section and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — N program genes**
   - 定义：Recorded n program genes for one candidate pair, model term, section and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Spots (n)**
   - 定义：Spots represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Positive spots (n)**
   - 定义：Positive spots represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Gene present in source**
   - 定义：Recorded gene present in source for one candidate pair, model term, section and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — Model status**
   - 定义：Recorded model status for one candidate pair, model term, section and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — Model term**
   - 定义：cell_abundance: conditional abundance coefficient; pathology: expression-programme coefficient; interaction: abundance × programme coefficient.
   - 单位/尺度：model term
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
18. **R — β**
   - 定义：Conditional spatial log-link coefficient for the stated model term.
   - 单位/尺度：natural-log link coefficient
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
19. **S — SE**
   - 定义：Standard error of the conditional spatial coefficient; not a confidence interval.
   - 单位/尺度：natural-log link coefficient
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
20. **T — Z**
   - 定义：Recorded z for one candidate pair, model term, section and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
21. **U — P**
   - 定义：Original P for P in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
22. **V — Model**
   - 定义：Recorded model for one candidate pair, model term, section and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 59 — Table 59
最终标题：Supplementary Table 59 | Spatial model eligibility
实际数据：2,842行；Excel第4–2845行；16列。
表内说明：Spot/positive-spot coverage and fit eligibility, retaining missing or failed results. One candidate pair, section and posterior setting. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One candidate pair, section and posterior setting.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Gene**
   - 定义：Reported gene symbol; original gene identities are unchanged.
   - 单位/尺度：gene symbol
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Cell context**
   - 定义：Recorded cell context for one candidate pair, section and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Study ID**
   - 定义：Source-study accession or its canonical identifier.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Section ID**
   - 定义：Recorded Section ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Analysis unit ID**
   - 定义：Recorded grouping used to combine section estimates before between-unit meta-analysis; the grouping is retained as approved.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Metadata status**
   - 定义：Source metadata status retained from model execution. It does not supersede author confirmation of the current analysis-unit grouping.
   - 单位/尺度：source status
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — QC fallback**
   - 定义：Recorded qC fallback for one candidate pair, section and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Posterior summary**
   - 定义：cell2location posterior summary used for abundance: q05 or posterior mean.
   - 单位/尺度：summary label
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Expression programme**
   - 定义：Molecular expression signature used in the spatial model; not an independently assigned histopathological region.
   - 单位/尺度：programme name
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Expected direction**
   - 定义：Recorded expected direction for one candidate pair, section and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Target removed from program**
   - 定义：Recorded target removed from program for one candidate pair, section and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — N program genes**
   - 定义：Recorded n program genes for one candidate pair, section and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Spots (n)**
   - 定义：Spots represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Positive spots (n)**
   - 定义：Positive spots represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Gene present in source**
   - 定义：Recorded gene present in source for one candidate pair, section and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — Model status**
   - 定义：Recorded model status for one candidate pair, section and posterior setting.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 60 — Table 60
最终标题：Supplementary Table 60 | Cell-proportion model eligibility
实际数据：816行；Excel第4–819行；12列。
表内说明：All planned cell-proportion models, including ineligible and failed combinations. One planned state–disease-background–denominator model. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One planned state–disease-background–denominator model.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Annotation level**
   - 定义：Recorded annotation level for one planned state–disease-background–denominator model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — State**
   - 定义：Recorded state for one planned state–disease-background–denominator model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Model id**
   - 定义：Recorded model id for one planned state–disease-background–denominator model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Background**
   - 定义：Recorded background for one planned state–disease-background–denominator model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Denominator**
   - 定义：Recorded denominator for one planned state–disease-background–denominator model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Samples (n)**
   - 定义：Samples represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Donors (n)**
   - 定义：Donors represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Status**
   - 定义：Recorded status for one planned state–disease-background–denominator model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Reason**
   - 定义：Recorded reason for one planned state–disease-background–denominator model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Formula**
   - 定义：Recorded formula for one planned state–disease-background–denominator model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Scope kind**
   - 定义：Recorded scope kind for one planned state–disease-background–denominator model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Exit code**
   - 定义：Recorded exit code for one planned state–disease-background–denominator model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 61 — Table 61
最终标题：Supplementary Table 61 | CoVarNet sample and donor metadata
实际数据：217行；Excel第4–220行；42列。
表内说明：Profiles and grouping metadata used in CoVarNet disease association, distinct from the structural profile set. One CoVarNet disease-association profile. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One CoVarNet disease-association profile.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — CoVarNet profile ID**
   - 定义：Recorded CoVarNet profile ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Study ID**
   - 定义：Source-study accession or its canonical identifier.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Analysis sample ID**
   - 定义：Recorded Analysis sample ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Specimen ID**
   - 定义：Biological specimen identifier, distinct from a technical library and from donor identity.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Donor ID**
   - 定义：Recorded donor identifier; NR/blank is not a new independent donor.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Tissue**
   - 定义：Recorded tissue for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Tissue detail**
   - 定义：Recorded tissue detail for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Cells (n)**
   - 定义：Cells represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Replicate relationship**
   - 定义：Recorded replicate relationship for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Evidence source**
   - 定义：Recorded evidence source for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Decision**
   - 定义：Recorded decision for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — Treat same-donor samples as independent**
   - 定义：Recorded treat same-donor samples as independent for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — Automatic specimen pooling**
   - 定义：Recorded automatic specimen pooling for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Design approval record**
   - 定义：Recorded design approval record for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Previous specimen ID**
   - 定义：Recorded Previous specimen ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — Previous donor ID**
   - 定义：Recorded Previous donor ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — Original clinical group**
   - 定义：Recorded original clinical group for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
18. **R — Covarnet group**
   - 定义：Recorded covarnet group for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
19. **S — Group evidence**
   - 定义：Recorded group evidence for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
20. **T — Covarnet include**
   - 定义：Recorded covarnet include for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
21. **U — Grouping-review flag**
   - 定义：Recorded grouping-review flag for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
22. **V — Source diagnosis**
   - 定义：Recorded source diagnosis for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
23. **W — Pooling**
   - 定义：Recorded pooling for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
24. **X — Analysis approval record**
   - 定义：Recorded analysis approval record for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
25. **Y — Condition**
   - 定义：Recorded condition for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
26. **Z — Sample name**
   - 定义：Recorded sample name for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
27. **AA — Age**
   - 定义：Recorded age for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
28. **AB — Sex**
   - 定义：Recorded sex for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
29. **AC — Assay type**
   - 定义：Recorded assay type for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
30. **AD — Vascular bed**
   - 定义：Recorded vascular bed for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
31. **AE — Cohort disease context**
   - 定义：Recorded cohort disease context for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
32. **AF — Donor disease context**
   - 定义：Recorded donor disease context for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
33. **AG — Tissue disease status**
   - 定义：Recorded tissue disease status for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
34. **AH — Donor resolved**
   - 定义：Recorded donor resolved for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
35. **AI — CoVarNet donor key**
   - 定义：Recorded coVarNet donor key for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
36. **AJ — Profile ID**
   - 定义：Retained analysis profile identifier; not interchangeable with donor count.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
37. **AK — Total atlas cells (n)**
   - 定义：Total atlas cells represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
38. **AL — Technical samples (n)**
   - 定义：Technical samples represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
39. **AM — Background**
   - 定义：Recorded background for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
40. **AN — Group**
   - 定义：Recorded group for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
41. **AO — Association donor ID**
   - 定义：Recorded Association donor ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
42. **AP — Association scope eligible**
   - 定义：Recorded association scope eligible for one CoVarNet disease-association profile.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 62 — Table 62
最终标题：Supplementary Table 62 | Donor-aware communication sensitivity
实际数据：5,506行；Excel第4–5509行；26列。
表内说明：Full pooled communication effects and donor-aware cohort contributions across tested interactions. One method–context interaction. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One method–context interaction.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Method**
   - 定义：Recorded method for one method–context interaction.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Disease context**
   - 定义：Recorded disease context for one method–context interaction.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Eligibility scheme**
   - 定义：Recorded eligibility scheme for one method–context interaction.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Interaction ID**
   - 定义：Recorded Interaction ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Sender state**
   - 定义：Recorded sender state for one method–context interaction.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Receiver state**
   - 定义：Recorded receiver state for one method–context interaction.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Ligand**
   - 定义：Recorded ligand for one method–context interaction.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Receptor**
   - 定义：Recorded receptor for one method–context interaction.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Signalling category**
   - 定义：Recorded signalling category for one method–context interaction.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Effect (disease − control)**
   - 定义：Recorded effect (disease − control) for one method–context interaction.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
11. **K — Pooled SE**
   - 定义：Standard error on the corresponding effect/estimate scale; not a confidence interval.
   - 单位/尺度：effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
12. **L — 95% CI lower**
   - 定义：Recorded confidence bound on the corresponding effect scale; no bounds were reconstructed.
   - 单位/尺度：effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
13. **M — 95% CI upper**
   - 定义：Recorded confidence bound on the corresponding effect scale; no bounds were reconstructed.
   - 单位/尺度：effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
14. **N — Pooled P**
   - 定义：Original P for Pooled in the stated model; source numeric underflow zeros and eligibility flags are retained.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
15. **O — Cohorts (n)**
   - 定义：Cohorts represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
16. **P — Studies (n)**
   - 定义：Studies represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
17. **Q — Tissues (n)**
   - 定义：Tissues represented by this row.
   - 单位/尺度：count
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
18. **R — Direction agreement**
   - 定义：Fraction of contributing effects with the recorded consensus direction; the qualifying units are defined by this table.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
19. **S — Largest cohort effect share**
   - 定义：Recorded largest cohort effect share, with its analysis-specific denominator; no conversion to percent.
   - 单位/尺度：fraction (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
20. **T — I² (%)**
   - 定义：Between-unit heterogeneity statistic on its original percentage scale.
   - 单位/尺度：percent (0–100)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
21. **U — τ²**
   - 定义：Estimated between-unit effect variance.
   - 单位/尺度：squared effect units
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
22. **V — Cohort contributions**
   - 定义：Recorded cohort contributions for one method–context interaction.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
23. **W — BH FDR**
   - 定义：Original adjusted P for the stated testing family; no correction has been recomputed.
   - 单位/尺度：probability (0–1)
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
24. **X — Context support**
   - 定义：Recorded context support for one method–context interaction.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
25. **Y — Evidence class**
   - 定义：Recorded evidence class for one method–context interaction.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
26. **Z — Tissue**
   - 定义：Recorded tissue for one method–context interaction.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 63 — Table 63
最终标题：Supplementary Table 63 | CoVarNet resampling and refits
实际数据：1,568行；Excel第4–1571行；9列。
表内说明：Reference-to-refit module matches across structural perturbations and hierarchical resampling. One reference-module match in a perturbation/refit. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One reference-module match in a perturbation/refit.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Reference module**
   - 定义：Recorded reference module for one reference-module match in a perturbation/refit.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Run type**
   - 定义：Recorded run type for one reference-module match in a perturbation/refit.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Omitted or bootstrap**
   - 定义：Recorded omitted or bootstrap for one reference-module match in a perturbation/refit.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Matching mode**
   - 定义：Recorded matching mode for one reference-module match in a perturbation/refit.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Matched module**
   - 定义：Recorded matched module for one reference-module match in a perturbation/refit.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — Cosine similarity**
   - 定义：Recorded cosine similarity for one reference-module match in a perturbation/refit.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Top node jaccard**
   - 定义：Recorded top node jaccard for one reference-module match in a perturbation/refit.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — N reference modules matching candidate**
   - 定义：Recorded n reference modules matching candidate for one reference-module match in a perturbation/refit.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Duplicate match**
   - 定义：Recorded duplicate match for one reference-module match in a perturbation/refit.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.

## Supplementary Table 64 — Table 64
最终标题：Supplementary Table 64 | Genetic model fit statistics
实际数据：4行；Excel第4–7行；10列。
表内说明：Fit statistics, roles and stable identifiers for the four compared genetic models; source underflow P=0 is retained. One confirmatory genetic model. Blank/NA values retain missing or inapplicable status, never numeric zero.
行单位：One confirmatory genetic model.
数据区：第3行单一表头；数据行共用该列序。未发现独立Record type列；不同模型/场景由相应身份字段区分。

1. **A — Model**
   - 定义：Recorded model for one confirmatory genetic model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
2. **B — Prespecified role**
   - 定义：Recorded prespecified role for one confirmatory genetic model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
3. **C — Chi-square**
   - 定义：Recorded chi-square for one confirmatory genetic model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
4. **D — Degrees of freedom**
   - 定义：Recorded degrees of freedom for one confirmatory genetic model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
5. **E — Raw software P value**
   - 定义：Recorded raw software P value for one confirmatory genetic model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
6. **F — P-value interpretation**
   - 定义：Recorded p-value interpretation for one confirmatory genetic model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
7. **G — Akaike information criterion**
   - 定义：Recorded akaike information criterion for one confirmatory genetic model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
8. **H — Comparative fit index**
   - 定义：Recorded comparative fit index for one confirmatory genetic model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
9. **I — Standardised root mean square residual**
   - 定义：Recorded standardised root mean square residual for one confirmatory genetic model.
   - 单位/尺度：source scale / category
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
10. **J — Stable model ID**
   - 定义：Recorded Stable model ID; preserved as an identifier, not recoded into another unit.
   - 单位/尺度：identifier
   - 缺失含义：Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.
