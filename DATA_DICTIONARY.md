# Data dictionary — district-level covariates

English glosses for the pre-specified district-level covariates used in
`norovirus_spatial_korea.R`. The Korean names are the column keys in the
official-statistics input files; they are kept verbatim in the script because
they must match those file headers exactly. This table is a reading aid and is
not used by the code.

Transformations are fixed a priori from each covariate's distribution in the full
2020–2024 panel and are identical in all strata: rates and percentages are entered
on their natural scale ("Continuous"), right-skewed counts, volumes and areas
(skewness > 1) are log(x + 1)-transformed, and livestock head is scaled per 10,000.
All covariates are then z-standardised, so each incidence rate ratio (IRR) is
expressed **per 1 standard deviation** of the transformed covariate. All
pre-specified covariates enter without outcome-based screening, and
multicollinearity is controlled by iterative removal of variables with a variance
inflation factor > 5. Years in which a covariate is zero or missing for all
districts are treated as missing and filled from the nearest available year within
the same district. "Direction" is the *a priori* hypothesised sign, not the fitted
result. Models: T, Total; U, Urban; R, Rural.

| Domain | Manuscript label | script `eng` | Korean (`kr`) | Transform | Models | Hyp. direction |
|---|---|---|---|---|---|---|
| Food source & livestock | Korean beef-cattle farms | `beef_farm` | 한육우농가수 | log(x+1) | T, U | Positive |
| Food source & livestock | Total livestock farms | `farm_total` | 농가수합계 | log(x+1) | T, U | Positive |
| Food source & livestock | Oyster production | `oyster` | 굴 | log(x+1) | T, U | Positive |
| Food source & livestock | Total livestock head | `livestock_total` | 사육두수_합계 | per 10,000 | T, U, R† | Positive |
| Food source & livestock | Korean cattle head | `cattle_beef` | 사육두수_한육우 | per 10,000 | R | Positive |
| Food source & livestock | Poultry head | `cattle_poultry` | 사육두수_가금 | per 10,000 | R | Positive |
| Sludge & waste treatment | On-site sludge self-treatment | `sludge_total` | 자체처리량계 | log(x+1) | T, U | Inverse |
| Sludge & waste treatment | Post-incineration sludge treatment | `sludge_incin` | 소각후처리 | log(x+1) | T, U | Inverse |
| Sludge & waste treatment | Post-drying sludge treatment | `sludge_dry` | 건조후처리 | log(x+1) | T, U | Positive |
| Sludge & waste treatment | Sludge-moisture content | `sludge_moisture` | 함수율_탈수기준 | Continuous | T, U, R | Positive |
| Water & sewerage infrastructure | Wastewater-effluent volume (industrial discharge) | `ww_effluent` | 폐수방류량 | log(x+1) | T, U | Positive |
| Water & sewerage infrastructure | Septic-tank population | `septic_pop` | 정화조인구 | log(x+1) | T, U | Positive |
| Water & sewerage infrastructure | Civil-defence groundwater wells | `gw_civil_count` | 민방위용_개소수 | log(x+1) | T, U, R | Positive |
| Water & sewerage infrastructure | Civil-defence groundwater use | `gw_civil_defense` | 지하수_민방위용 | log(x+1) | R | Positive |
| Water & sewerage infrastructure | School groundwater use | `school_gw_use` | 학교용_지하수이용량 | log(x+1) | R | Positive |
| Water & sewerage infrastructure | Groundwater non-compliance count | `gw_unfit_count` | 지하수_부적합건수 | log(x+1) | R | Positive |
| Water & sewerage infrastructure | Storm-water gullies | `rainwater_gully` | 받이_빗물받이 | log(x+1) | R | Inverse |
| Waterways & land use | Reservoir area | `reservoir` | 유지 | log(x+1) | T, U, R | Positive |
| Waterways & land use | River area | `river` | 하천 | log(x+1) | T, U | Positive |
| Waterways & land use | Peri-urban pastureland | `ranch` | 목장용지 | log(x+1) | T, U | Positive |
| Waterways & land use | Paddy-field area | `paddy` | 답(논) | log(x+1) | T, U | Neutral |
| Waterways & land use | Forest area | `forest` | 임야 | log(x+1) | R | Neutral |
| Waterways & land use | Dry-field area | `dry_field` | 전(밭) | Continuous | R | Positive |
| Hygiene & behaviour | Post-toilet handwashing rate | `handwash_toilet` | 화장실손씻기 | Continuous | T, U, R | Inverse |
| Hygiene & behaviour | Usual handwashing rate | `handwash_usual_adj` | 평소손씻기_조율 | Continuous | R | Inverse |
| Hygiene & behaviour | Healthy-lifestyle practice rate | `health_practice` | 건강생활실천율 | Continuous | T, U | Inverse |
| Hygiene & behaviour | Walking-practice rate | `walking_practice` | 걷기실천율 | Continuous | T, U | Inverse |
| Hygiene & behaviour | Influenza vaccination rate | `flu_vaccination_rate` | 인플루엔자예방접종률 | Continuous | R | Inverse |
| Socioeconomic & vulnerability | Elderly living alone | `elderly_alone` | 독거노인 | Continuous | T, U | Positive |
| Socioeconomic & vulnerability | Rural population | `rural_pop` | 농촌인구수 | log(x+1) | T, U | Positive |
| Socioeconomic & vulnerability | Fiscal autonomy | `fiscal_auto` | 재정자주도 | Continuous | T, U | Inverse |
| Socioeconomic & vulnerability | Fiscal independence | `fiscal_indep` | 재정자립도 | Continuous | T, U | Inverse |
| Socioeconomic & vulnerability | One-person household rate | `alone_rate` | 1인가구율_전체 | Continuous | T, U, R | Neutral |
| Socioeconomic & vulnerability | Basic-livelihood (welfare) recipient rate | `welfare_rate` | 기초생활수급자수율 | Continuous | T, U | Positive |
| Demographic | Households with elderly (65+) | `hh_elderly` | 총가구수65세이상 | Continuous | T, U | Positive |
| Demographic | Sex ratio | `sex_ratio` | 성비 | Continuous | T, U, R | Neutral |
| Demographic | Elderly population share | `elderly_rate` | 고령인구비율 | Continuous | T, U, R | Neutral |
| Healthcare access | In-district outpatient medical cost | `med_out` | 진료비외래 | log(x+1) | T, U | Inverse |
| Healthcare access | In-district total medical cost | `med_total` | 관내진료비전체 | log(x+1) | T, U | Inverse |
| Paediatric | Child population share, 0–4 y | `child_0_4` | 영유아비율_0_4 | Continuous | T, U, R | Positive |
| Paediatric | Child population share, 5–9 y | `child_5_9` | 아동비율_5_9 | Continuous | T, U, R† | Positive |

† removed from the Rural model by variance-inflation pruning.

The paediatric battery adds the two child shares, child-care facility density and
the preventable-hospitalisation rate for paediatric gastroenteritis to the credible
determinants of the principal Total model (total livestock head, sludge-moisture
content and in-district outpatient medical cost) and is fitted in the principal
model M4 in each stratum. The preventable-hospitalisation rate is taken from the
nearest year with published district values.
