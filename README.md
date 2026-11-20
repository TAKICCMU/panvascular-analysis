# Panvascular manuscript release — 2026.11.20

This release corresponds to 4 main figures, 10 Extended Data figures and 60 Supplementary Tables in the current manuscript. Use `06_Maps/CURRENT_PANEL_DATA_CODE_MAP.tsv` and `06_Maps/ALL60_TABLE_CODE_MAP.tsv` as the current navigation. Historical filenames retain their scientific identities and do not define current panel numbers.

The final atlas contains 1,175,115 cells and 37,336 genes from 26 studies and 247 sample identifiers. It uses the CPU Harmony theta=2 branch. There are 16 L1 annotation values and 116 named L2 states; Unknown and Not_assessed are additional annotation values. Neither is a newly identified biological lineage. The complete expression matrix is not included.

Current figure identities: Figure1f is clinical percentage summaries; Figure1g is L2 coverage; Figure1h is the 14-category × 28-marker display with per-gene 0–1 colour scaling and the original detection fractions. Figure3c shows separate MAGMA and pseudobulk GSEA analyses (Tables21/22). Figure4 combines the adopted CoVarNet and spatial results; former Figure5 source filenames remain provenance, not a fifth current main figure. Cell communication results and the four associated former tables are not part of the current manuscript.

Scientific definitions remain unchanged. Candidate meta uses REML, whereas pathway ranks use DL. The five joint candidate conditions retain 29 gene–state pairs from 15,639 evaluable pairs; the 1,658-gene P/E/M union is not an upstream count for the 2,199-gene prioritisation universe. Spatial analyses retain the accepted 49-section → 35-unit grouping. The displayed FN1/SCARB1 results contain 60 gene–unit estimates spanning 34 recorded unit IDs. Table8 is a donor-identifier mean; Figure1e uses a sample mean. These are different summaries.

The packages collect existing source files and adopted outputs. They do not replace any input data, scientific statistics or model fits. Existing execution logs describe their original runs. A packaging hash check is not a new biological validation or an end-to-end rerun.

## Licences
Author analysis/plotting code: MIT (see LICENSE). Author-derived data: CC BY 4.0. Original public inputs, external software, illustrations and BioRender assets retain their original rights and terms; their inclusion does not relicense them. Full component GWAS files, individual genotypes, full count matrices and spatial model weights are not redistributed.

The code archive contains `CODE_GUIDE.md`; the data archive contains `DATA_GUIDE.md`. Unpack both archives into the same parent directory. Their shared navigation files are identical. Checksums are separately named CODE_MANIFEST_SHA256.tsv and DATA_MANIFEST_SHA256.tsv.
