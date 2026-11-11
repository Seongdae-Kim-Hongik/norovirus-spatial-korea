# Norovirus spatial Korea — R2 correction (v2.4.0)

This is the author-authorised R2 moisture-correction release, version 2.4.0. The existing v2.3.0 DOI covers R1 only; the repository concept DOI resolves to the versioned Zenodo archive.

## Scope
The frozen outcomes are 4,801 cases and 220 outbreaks in 1,115 district-years (223 districts, 2020–2024). This update corrects plant/subtotal aggregation and zero codes in sludge moisture, retains the other R1 covariates, and repeats principal, paediatric, outbreak, graph, prior, influence and availability analyses. No individual patient information is included.

## Required inputs
Place `analysis_input.rds` in `input/`: the archived R1 analysis object containing `cor_merged`, `shp_main`, `g_main`, `nb_obj`, `FIX_FORMS` and `FAMILY`. It must be the R1 object underlying version 2.3.0; new outcome data are not interchangeable. This object and the official annual sewerage workbooks are **not bundled in this code-only archive**. Obtain them from the corresponding author and the Ministry of Environment source. Workbook names are `공식_2020_하수도통계.xlsx` through `공식_2024_하수도통계.xlsx`, using the sludge-generation/treatment worksheet. Expected workbook hashes are in validated_results/source_correction_summary.json. The two small R1 panel/region CSVs, fixed covariate list and verified outbreak panel are included.

## Run
Set the working directory to this release directory (or set R2_CODE_ROOT to its absolute path for R). Python needs pandas, numpy and openpyxl; R needs INLA, dplyr, spdep, sf, glmmTMB and car. See sessionInfo_final.txt for versions.

```sh
python3 correct_moisture_r2.py
Rscript reanalyse_corrected_moisture.R
Rscript reanalyse_additional_r2.R
Rscript reanalyse_final_checks_r2.R
```

Outputs are written only to `output/`. Do not run the R scripts concurrently: later stages extend the result CSV. Cached fits are reused. Raw physical moisture remains missing when not observed; a centred design entry and availability indicator are used to retain the outcome panel. Four source-reported 0.1% values are kept and tested separately without an unsupported unit conversion.

## Numerical qualification
The manuscript reports first valid corrected M4 fits, their independent repeats, and the finite-start M5 refit (repeat_id=2). The initial M5 run with WAIC 30,958.6 is retained as a documented unstable result, not a reported model-comparison estimate. Failed/non-finite/non-positive-DIC fits were rejected; spatial fits were compared with corresponding non-spatial fits. Repeatability is not proof of global optimality. Rural estimates are especially uncertain. Full-precision validated results are provided for comparison; rounding is half up only at document generation.

## Provenance and limitations
The outcome source hash and all district-level R1 comparisons are recorded in outcome_freeze_verification.json. The reconstruction corrects moisture only; it does not independently recreate every administrative covariate from all raw sources. Reproducibility requires the specified R1 object, not merely the small CSVs included here. R1 publication: https://doi.org/10.5281/zenodo.22821175 ; concept DOI: https://doi.org/10.5281/zenodo.20725142 .
