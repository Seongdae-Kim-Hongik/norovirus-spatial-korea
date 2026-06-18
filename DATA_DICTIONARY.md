# Data dictionary — district-level covariates

English glosses for the pre-specified district-level covariates used in
`norovirus_spatial_korea.R`. The Korean names are the column keys in the
(non-redistributed) official-statistics input files; they are kept verbatim in
the script because they must match those file headers exactly. This table is a
reading aid for reviewers and is not used by the code.

All covariates are transformed (median-split binary, tertile/quartile category,
per-10,000, or log1p) and then z-standardised, so each incidence rate ratio (IRR)
is expressed **per 1 standard deviation** of the transformed covariate. Covariates
with a bivariate screening p < 0.20 (plus a-priori forced covariates, marked ‡)
are retained, and multicollinearity is controlled by iterative removal of
variables with a variance inflation factor > 10. "Direction" is the *a priori*
hypothesised sign, not the fitted result.

| Domain | Manuscript label | script `eng` | Korean (`kr`) | Hyp. direction |
|---|---|---|---|---|
| Food source & livestock | Korean beef-cattle farms ‡ | `beef_farm` | 한육우농가수 | Positive |
| Food source & livestock | Total livestock farms | `farm_total` | 농가수합계 | Positive |
| Food source & livestock | Oyster production | `oyster` | 굴 | Positive |
| Food source & livestock | Total livestock head | `livestock_total` | 사육두수_합계 | Positive |
| Sludge & waste treatment | On-site sludge self-treatment (total) | `sludge_total` | 자체처리량계 | Inverse |
| Sludge & waste treatment | Sludge, post-incineration treatment | `sludge_incin` | 소각후처리 | Inverse |
| Sludge & waste treatment | Post-drying sludge | `sludge_dry` | 건조후처리 | Positive |
| Sludge & waste treatment | Sludge-moisture content | `sludge_moisture` | 함수율_탈수기준 | Positive |
| Sludge & waste treatment | Outsourced waste (fuel use) | `waste_fuel` | 외부위탁_연료 | Inverse |
| Water & sewerage infrastructure | Wastewater-effluent volume | `ww_effluent` | 폐수방류량 | Positive |
| Water & sewerage infrastructure | Septic-tank population | `septic_pop` | 정화조인구 | Positive |
| Water & sewerage infrastructure | Household groundwater reliance (emergency wells) | `gw_civil_count` | 민방위용_개소수 | Positive |
| Waterways & land use | Reservoir area | `reservoir` | 유지 | Positive |
| Waterways & land use | River area | `river` | 하천 | Positive |
| Waterways & land use | Peri-urban pastureland | `ranch` | 목장용지 | Positive |
| Waterways & land use | Paddy-field area | `paddy` | 답(논) | Neutral |
| Hygiene & behaviour | Post-toilet handwashing rate | `handwash_toilet` | 화장실손씻기 | Inverse |
| Hygiene & behaviour | Healthy-lifestyle practice rate | `health_practice` | 건강생활실천율 | Inverse |
| Hygiene & behaviour | Walking-practice rate | `walking_practice` | 걷기실천율 | Inverse |
| Socioeconomic & vulnerability | Elderly living alone | `elderly_alone` | 독거노인 | Positive |
| Socioeconomic & vulnerability | Rural population | `rural_pop` | 농촌인구수 | Positive |
| Socioeconomic & vulnerability | Fiscal autonomy | `fiscal_auto` | 재정자주도 | Inverse |
| Socioeconomic & vulnerability | Fiscal independence | `fiscal_indep` | 재정자립도 | Inverse |
| Socioeconomic & vulnerability | One-person household rate | `alone_rate` | 1인가구율_전체 | Neutral |
| Socioeconomic & vulnerability | Basic-livelihood (welfare) recipient rate | `welfare_rate` | 기초생활수급자수율 | Positive |
| Demographic | Households with elderly (65+) | `hh_elderly` | 총가구수65세이상 | Positive |
| Demographic | Sex ratio ‡ | `sex_ratio` | 성비 | Neutral |
| Demographic | Elderly population share ‡ | `elderly_rate` | 고령인구비율 | Neutral |
| Healthcare access | Outpatient medical cost | `med_out` | 진료비외래 | Inverse |
| Healthcare access | In-district total medical cost | `med_total` | 관내진료비전체 | Inverse |
| Paediatric | Child population share, 0–4 y | `child_0_4` | 영유아비율_0_4 | Positive |
| Paediatric | Child population share, 5–9 y | `child_5_9` | 아동비율_5_9 | Positive |

‡ a-priori forced covariate (retained regardless of the screening p-value).

The script's `TV_v6` data frame holds, for each covariate, the exact input-file
column key (`code`), the short English name (`eng`), the domain (`cat`) and the
hypothesised direction (`이론방향`: 위험 = Positive, 보호 = Inverse, 중립 = Neutral).
The Urban model drops `farm_total` after urban-stratum screening; the Rural model
uses a reduced set adapted to the sparse rural counts. A separate paediatric
"battery" (the two child shares plus child-care facility density and the
preventable-hospitalisation rate for paediatric gastroenteritis) is added to a
four-variable environmental-anchor model for the dissociation test.
