# CoVarNet remaining work prepared 2026-09-16
Production was NOT launched by this preparation agent.

Commands (each reserves 1 CPU; may run independently, at most 2 workers total):
`python3 code/run.py associations`
`python3 code/run.py refits`
Use absolute path to run.py from any directory. Runtime R: refits communication_r/bin/Rscript; associations /usr/bin/Rscript with W/Rlib clubSandwich. BLAS and OpenMP capped to one; NMF p1.

Inputs: completed W/covarnet primary K14 (223 profiles,89 nodes), latest three-background sidecar SHA256 852c08353c53a1b4d32b0b5928c9d1ac6a0278e365748912adfaf2ff30cf2ce4. Scope exclusions retained. All outputs only here. Primary model not rerun.

Refits: inherited nsNMF K12,13,15,16 x100 starts; remaining sensitivities x30 starts: robust minmax, Hellinger, leave-tissue, healthy-only/disease-only, 30 hierarchical study then donor-block bootstraps, 10 multinomial depth500 replicates (stable nodes + other denominator). Exact Hungarian one-to-one with zero similarity for unmatched reference modules at smaller rank. Reuse already-completed denominator, donor bootstrap, LOS and parent-null from primary directory; no duplicate computation. Original multinomial sampling is a cell-depth perturbation with replacement, not literal without-replacement subsampling. Study/tissue concentration is descriptive profile-weighted relative H coefficient mass, not cell abundance.

Disease: latest metadata join by analysis_profile_uid. CLR H with pseudocount1e-6 and1e-4 sensitivity. Each background x tissue x study separately; donor CR2 variance and Satterthwaite df, requires >=3 donors/group and >=6 total. DL random-effects meta only across studies with identical background+tissue; sparse/no-control designs explicitly NOT_ESTIMABLE. Spearman top-five node pairs computed within study+tissue on donor means, pooled only within tissue using approximate Fisher-z variance. These approximations and exploratory status are written to METHODS. No all-bed pooled P values.

run.py binds all inputs and code by SHA256 before resumption; per-fit R caches additionally bind exact matrix, rank, starts, seed and code. STATUS.json records completion only after zero return code and hashes all outputs. Existing historical outputs are not modified.

Validation: test_cr2_final/TEST_PASS.txt: actual nsNMF small fit, exact equal/lower-rank matching, repeated-donor synthetic disease and node-correlation end-to-end fixture. This is executable validation, not production results. Full production fit robustness assessed only after completion.

Final disease adapter adds explicit eligible G1/G2/G3 scope, original-H descriptive summaries, cross-study donor-overlap meta blocking, and exploratory CLR interpretation. It is not a verbatim reproduction of the original manuscript.
