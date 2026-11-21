# CPO diagnostics, version 2.4.2 (2026-11-21)

This code supplements, rather than replaces, the corrected R2 principal analysis. No original observations are removed from the principal analysis: the frozen totals remain 4,801 cases and 220 outbreaks. The archive contains code and documentation only; input data and cached fitted objects must be supplied separately under the conditions described in the R2 README.

Set R2_CODE_ROOT to the absolute revision_R2 directory. First run the corrected-moisture pipeline to generate output/fits/*_corrected__M4__cases__1.rds. Preserve input/analysis_input.rds, the adjacency graph and the original row ordering. Use the recorded INLA environment (sessionInfo.txt); numerical approximations and warning locations may vary by software version/platform.

Run in this order:

1. Rscript cpo_diagnostics/gcpo.R
2. Rscript cpo_diagnostics/loo_all.R
3. Rscript cpo_diagnostics/loo_extra.R
4. Rscript cpo_diagnostics/precision_Namwon.R
5. Rscript cpo_diagnostics/independent_numerics/check_grid.R
6. Rscript cpo_diagnostics/export_gcpo.R
7. python3 cpo_diagnostics/summarize_diagnostics.py

The first direct-refit script addresses all flags in the first principal fit and adds a simultaneous-exclusion stress test in each stratum. The extra script includes national indices 455 and 578, additional locations flagged among the 24 repeated offset-comparison fits. These indices refer to the frozen, sorted district-year input, not arbitrary new datasets. The expected union is 17 national, 11 urban and one rural observation (17 distinct district-years), and there are 32 direct fits including three stress tests. Cached direct fits are reused; remove the relevant output cache only when intentionally recomputing.

The response alone is set to NA, preserving covariates, offset, graph and latent indices. Predictive masses average over joint posterior draws of the held-out predictor and negative-binomial size. All single-observation calculations use 20,000 draws except national Namwon (100,000). The separate Daegu Seo-gu calculation uses full Laplace/grid integration. MCSE measures simulation variability conditional on the INLA approximation; it does not quantify approximation bias. Persistent local-search warnings must be reported, not suppressed. GCPO is an additional comparison, not a universally equivalent replacement for direct refitting.

The companion ../pool_offset_runs.R takes four paths: principal fit directory, first offset-repeat fit directory, second offset-repeat fit directory, and output directory. Run the offset sensitivity pipeline independently in two output roots to obtain the two repeated sets; do not overwrite the first set.

Expected publication summaries: all national/urban moisture intervals include one in the single-exclusion checks; livestock intervals exclude one in 14/17 national and 9/11 urban fits. The simultaneous-exclusion test is selective and is not an alternative primary analysis. The primary model was not chosen using CPO.
