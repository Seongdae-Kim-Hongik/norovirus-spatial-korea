# Foodborne norovirus — district-level Bayesian spatial analysis, Republic of Korea, 2020–2024

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20725142.svg)](https://doi.org/10.5281/zenodo.20725142)


Reproducible R analysis code for the manuscript:
**"Divergent urban and rural environmental drivers of foodborne norovirus and its dissociation from the paediatric disease burden: a Bayesian spatial analysis across South Korean districts, 2020–2024"** (submitted to *Science of the Total Environment*).

## Contents
- `DATA_DICTIONARY.md` — English glosses for every district-level covariate
  (Korean column key → English label → hypothesised direction), a reading aid for
  the Korean variable names retained in the script.
- `norovirus_spatial_korea.R` — full pipeline:
  data loading → pre-specified covariate construction (transformation + z-standardisation) →
  urban/rural stratification (Total / Urban / Rural) → Bayesian negative-binomial models with
  Besag–York–Mollié (BYM) spatial random effects, a temporal random walk and a space–time
  interaction, fitted by integrated nested Laplace approximation (INLA), compared across six
  specifications M1–M6 (principal model M4) → Global Moran's I (pre/post) and high-/low-risk
  district mapping → paediatric battery (collinearity/VIF + dissociation test) →
  robustness (BYM2 φ, eight neighbourhood graphs, alternative precision priors) →
  case-versus-outbreak counting-unit sensitivity → peer-review robustness extensions (CPO/predictive checks, ridge-type regularization, RW1 temporal effects, a unified urbanicity-interaction model, a k-nearest-neighbour spatial graph, sewerage-coverage adjustment, and per-standard-deviation covariate magnitudes) → result tables (xlsx).

## Data availability
The inputs are aggregated official statistics with **no personal identifiers** and are **not
redistributed** here. District-year laboratory-confirmed norovirus food-poisoning counts are
from the Korea Ministry of Food and Drug Safety (MFDS) Food Poisoning Statistics System;
population denominators from the Korean Statistical Information Service (KOSIS, Statistics Korea).
Pre-specified district covariates are compiled from official government statistics (Statistics
Korea; Korea Water Resources Corporation National Groundwater Information Center; National
Sewerage Information System; Ministry of Land, Infrastructure and Transport; KDCA Community
Health Survey; Health Insurance Review and Assessment Service). District boundaries are from the
Statistics Korea SGIS administrative-boundary service. Edit `BASE_IV` in the script to point to
the folder holding these files.

## Requirements
R ≥ 4.5 and **R-INLA** (installed automatically from `https://inla.r-inla-download.org/R/stable`).
Other packages auto-installed by the script: MASS, arrow, car, dplyr, openxlsx, sf, spdep,
stringr, tidyr.

## How to run
```
Rscript norovirus_spatial_korea.R
```
The script auto-installs packages on first run and writes tables/log to `output/`. R-INLA fits
take roughly 10–30 minutes depending on the machine.

## Reproducibility
The script runs end-to-end and reproduces the manuscript headline numbers (deterministic INLA):
Total-model M4 incidence rate ratios — groundwater reliance 1.52 (1.04–2.23), peri-urban
pastureland 1.57 (1.04–2.37), wastewater-effluent volume 0.47, reservoir area 0.48 (DIC 2,377;
N = 1,112 district-years); Urban-model M4 — sludge-moisture 1.91 (1.15–3.16), pastureland 2.24,
on-site sludge self-treatment 0.56, child share 0–4 y 0.54; Global Moran's I +0.04 → residual
−0.005 with 0 high- and 0 low-risk districts; BYM2 φ = 0.34 (0.01–0.91); case-versus-outbreak
Spearman ρ ≈ 0.95.

## Citation
Kim S, Chun BC. Divergent urban and rural environmental drivers of foodborne norovirus and its
dissociation from the paediatric disease burden: a Bayesian spatial analysis across South Korean
districts, 2020–2024 (manuscript under review). Archived code: https://doi.org/10.5281/zenodo.20725142 (Zenodo). Citation to be updated on publication.

## License
MIT (code) — see `LICENSE`.
