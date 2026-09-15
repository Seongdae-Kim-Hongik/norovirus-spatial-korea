# Foodborne norovirus — district-level Bayesian spatial analysis, Republic of Korea, 2020–2024

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20725142.svg)](https://doi.org/10.5281/zenodo.20725142)


Reproducible R analysis code for the manuscript:
**"Divergent urban and rural environmental drivers of foodborne norovirus infections and their dissociation from paediatric disease burden in South Korea"** (under review at *Scientific Reports*).

## Version 2 (revision)
Version 2 corrects the data preparation and model specification of version 1:
- **Covariate-year harmonisation.** Several administrative covariate series are not published every year and
  had been stored as zero, rather than missing, for all districts in those years. Such structurally missing
  years are now set to missing and filled from the nearest available year within the same district.
- **Missing values in binary transforms** are kept as missing (version 1 converted them to zero).
- **District industrial wastewater discharge, 2020–2023,** is read from `ext_wastewater_kosis.csv`
  (Ministry of Environment, *Industrial Wastewater Generation and Treatment*, KOSIS table DT_106N_01_0100069;
  public official statistics; cities with several gu are averaged over gu, as in the original series).
- **On-site sludge self-treatment** is the sum of the primary on-site treatment categories reported in every year.
- **Covariate transformations are fixed a priori** from each covariate's distribution in the full panel and are
  identical in all strata (rates and percentages on the natural scale; right-skewed counts, volumes and areas
  log(x + 1); livestock head per 10,000). All pre-specified covariates enter without outcome-based screening;
  multicollinearity is controlled by a variance inflation factor < 5.
- **Covariate effects are estimated in the principal model M4** (BYM + space–time interaction) in every stratum.
- **Estimation settings.** Every INLA call goes through a wrapper that (i) disables the variational-Bayes correction of
  the posterior mean, which in this collinear covariate set displaced some estimates from their maximum-likelihood
  values, and (ii) refits any random-effects model whose DIC exceeds the corresponding non-spatial fit by more than 50
  (a degenerate mode seen occasionally for the space–time models). INLA runs with its default multithreading;
  single-threaded runs were found to degenerate.
- The k-nearest-neighbour graphs use the EPSG:5179 coordinates of the boundary file, and the robustness
  extensions use the credible determinants of the fitted model.
- The paediatric battery (Table 2) is fitted in the script.

## Contents
- `norovirus_spatial_korea.R` — full pipeline: data loading → covariate-year harmonisation and a priori
  transformations → urban/rural stratification (Total / Urban / Rural) → Bayesian negative-binomial models with
  Besag–York–Mollié (BYM) spatial random effects fitted by integrated nested Laplace approximation (INLA), compared
  across six specifications M1–M6 (principal model M4) → Global Moran's I (pre/post) and high-/low-risk districts →
  robustness (BYM2 φ, eight district adjacency graphs, alternative precision priors) → case-versus-outbreak
  counting-unit sensitivity → robustness extensions (CPO, ridge-type regularisation, RW1 temporal effects, unified
  urbanicity-interaction model, k-nearest-neighbour graph, sewerage-coverage adjustment) → paediatric battery →
  result tables.
- `sensitivity_filled_years_influence.R` — sensitivity of the principal estimates to the filled covariate years
  (district mean of observed years, linear interpolation, omission of the filled covariates, observed values only
  for 2020–2023) and to individual years and district-years (each year omitted; the district-years with the largest
  case counts excluded); Supplementary Table S12. Run as `SENS=<specification> Rscript sensitivity_filled_years_influence.R`
  (specifications are listed in the script header).
- `ext_wastewater_kosis.csv` — district industrial wastewater discharge (m³/day), 2020–2023.
- `DATA_DICTIONARY.md` — English glosses for every district-level covariate.

## Data availability
The inputs are aggregated official statistics with **no personal identifiers**. Apart from the wastewater series
above, they are **not redistributed** here. District-year laboratory-confirmed norovirus food-poisoning counts are
from the Korea Ministry of Food and Drug Safety (MFDS) Food Poisoning Statistics System; population denominators
from the Korean Statistical Information Service (KOSIS, Statistics Korea). Covariates are compiled from official
government statistics (Statistics Korea; Ministry of Environment; Korea Water Resources Corporation National
Groundwater Information Center; National Sewerage Information System; Ministry of Land, Infrastructure and
Transport; KDCA Community Health Survey; Health Insurance Review and Assessment Service). District boundaries are
from the Statistics Korea SGIS administrative-boundary service. Edit `BASE_IV` in the script to point to the folder
holding these files, and keep `ext_wastewater_kosis.csv` in the working directory.

## Requirements
R ≥ 4.5 and **R-INLA** (installed automatically from `https://inla.r-inla-download.org/R/stable`).
Other packages auto-installed by the script: MASS, arrow, car, dplyr, openxlsx, sf, spdep, stringr, tidyr, glmmTMB.

## How to run
```
Rscript norovirus_spatial_korea.R
```
The script writes tables and a run log to `output/`. R-INLA fits take roughly 10–30 minutes depending on the machine.

## Reproducibility
The script reproduces the manuscript estimates (principal model M4; incidence rate ratios per 1 SD, 95% CrI):
Total model (1,115 district-years) — total livestock head 1.88 (1.19–2.96), sludge-moisture content 1.61
(1.10–2.36); Urban model (730) — total livestock head 2.62 (1.51–4.52); Rural model (385) — storm-water gullies
0.28 (0.09–0.85). Global Moran's I +0.036 → residual −0.004, with 0 high- and 0 low-risk districts; M1–M6 DIC
2,416.2–2,417.0 (differences within run-to-run variation); BYM2 φ = 0.36 (0.01–0.94); unified
urbanicity-interaction model, livestock × urban 3.98 (1.69–9.41). Two independent runs reproduce every IRR to
within 0.02 and every DIC to within 0.5. An independent maximum-likelihood fit (glmmTMB, written to
`output/mle_crosscheck_glmmTMB.csv`) agrees with the INLA estimates of the non-spatial model to two decimals.
In `sensitivity_filled_years_influence.R`, the sludge-moisture association remains credible under every alternative
handling of the filled covariate years and after excluding the district-years with the largest case counts, whereas the
total livestock head association is not credible when 2024 is omitted or when the district-year with the largest case
count is excluded (Total 1.06, 0.64–1.74).

## Citation
Kim S, Chun BC. Divergent urban and rural environmental drivers of foodborne norovirus infections and their
dissociation from paediatric disease burden in South Korea (manuscript under review). Archived code:
https://doi.org/10.5281/zenodo.20725142 (Zenodo). Citation to be updated on publication.

## License
MIT (code) — see `LICENSE`.
