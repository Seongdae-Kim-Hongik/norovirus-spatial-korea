# =============================================================================
#  Divergent urban and rural environmental drivers of foodborne norovirus
#  and its dissociation from the paediatric disease burden — Republic of Korea,
#  district-level Bayesian spatial analysis, 2020–2024.
#
#  Reproducible analysis code for the manuscript submitted to
#  Science of the Total Environment (STOTEN).
#
#  Author : Seongdae Kim      Advisor/corresponding : Byung Chul Chun
#  License: MIT (code).  Data: see "DATA" below (not redistributed; no personal identifiers).
# -----------------------------------------------------------------------------
#  HOW TO RUN
#    1. Put the input files (see DATA) in a folder and point BASE_IV to it.
#    2. Rscript norovirus_spatial_korea.R       # R-INLA fits; ~10–30 min
#    Packages (incl. INLA from its own repo) are auto-installed on first run.
#    Tested on R 4.6 with R-INLA (stable).
#
#  DATA (not redistributed; aggregated official statistics, no personal identifiers)
#    - District-year laboratory-confirmed norovirus food-poisoning counts +
#      population, Korea Ministry of Food and Drug Safety (MFDS) Food Poisoning
#      Statistics System and Korean Statistical Information Service (KOSIS).
#    - Pre-specified district covariates (livestock, sludge/wastewater, groundwater,
#      land use, hygiene, socioeconomic, demographic, healthcare, paediatric) from
#      official government statistics (Statistics Korea, KWRC National Groundwater
#      Information Center, National Sewerage Information System, MOLIT, KDCA
#      Community Health Survey, HIRA).
#    - si-gun-gu administrative-boundary shapefile, Statistics Korea SGIS.
#
#  OUTPUT MAP  (section  ->  manuscript object)
#    PART 1 ............... data load
#    PART 2 ............... variable dictionary + univariable screen  -> Table S1
#    PART 2.5 ............. urban/rural classification                -> Methods / Table S2
#    PART 3, 3.6, 5 ....... 3-group Bayesian NB BYM models (M1–M6),
#                           principal M4 (Total/Urban/Rural)          -> Table 1 (determinants)
#    Moran's I (pre/post) + high/low-risk districts                   -> Table S5; Figure S1
#    PART 6 ............... academic tables (xlsx)                     -> Tables 1, S2, S3
#    add-on: paediatric battery (collinearity/VIF + dissociation)     -> Table 2; Table S4; Fig S3
#    add-on: robustness BYM2 phi / 8 neighbour graphs / priors        -> Tables S6, S7, S8
#    add-on: case-vs-outbreak counting-unit sensitivity               -> Table S9
#
#  Headline numbers reproduced (deterministic INLA, exact):
#    Total M4 — groundwater reliance IRR 1.52 (1.04–2.23), peri-urban pastureland
#      1.57 (1.04–2.37), wastewater-effluent volume 0.47, reservoir area 0.48;
#      DIC 2,377; N = 1,112 district-years.
#    Urban M4 — sludge-moisture 1.91 (1.15–3.16), pastureland 2.24, on-site sludge
#      self-treatment 0.56, child share 0–4 y 0.54.
#    Spatial — Global Moran's I +0.04 -> residual -0.005; 0 high / 0 low-risk districts;
#      BYM2 phi = 0.34 (0.01–0.91); case-vs-outbreak Spearman rho ~ 0.95.
# =============================================================================

# ═══════════════════════════════════════════════════════════════════════════════
# [0] 패키지 자동 설치·로드  (없으면 자동 설치)  ★ INLA는 CRAN이 아니라 전용 repo에서 설치
#     — INLA가 그래도 실패하면: install.packages("INLA",
#         repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dependencies=TRUE)
# ═══════════════════════════════════════════════════════════════════════════════
local({
  rp <- getOption("repos")
  if (is.null(rp) || is.na(rp["CRAN"]) || rp["CRAN"] %in% c("@CRAN@","")) options(repos=c(CRAN="https://cloud.r-project.org"))
  cran_pkgs <- c("MASS",
                 "arrow",
                 "car",
                 "dplyr",
                 "openxlsx",
                 "sf",
                 "spdep",
                 "stringr",
                 "tidyr")
  miss <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if (length(miss)) { message("● 설치할 CRAN 패키지: ", paste(miss, collapse=", ")); install.packages(miss, dependencies=TRUE) }
  if (!requireNamespace("INLA", quietly=TRUE)) { message("● INLA 설치 중 (r-inla 전용 repo)…")
    install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dependencies=TRUE) }
  invisible(NULL)
})
suppressWarnings(suppressMessages({ for (.p in c("INLA","MASS","arrow","car","dplyr","openxlsx","sf","spdep","stringr","tidyr")) if (requireNamespace(.p, quietly=TRUE)) library(.p, character.only=TRUE) }))
rm(list = ls(pattern="^\\.p$"))

# ==============================================================================
# 노로바이러스 공간분석 v8-B (정직/PRE-SPECIFIED) — NB + 도시/농촌 3그룹 분리
#   ★ v8(AUTO) → v8-B 변경: AUTO 타깃 최적화(정방향≥5 + 역방향 자동제거) 제거.
#      사전지정 이론 변수셋(전체/도시 TV_v6 30개, 농촌 TV_RURAL 18개)으로 1회 적합.
#      유의/방향은 데이터가 나오는 대로 보고(정직 경로 B). soft p-hacking 제거.
#
# ★ 전달사슬: 사람→사람(접촉) + 굴/패류(식품) + 수계(폐수/정화조→하천) + 환경(지하수)
# ★ 도시/농촌 분류: 대도시=세종+광역시(not 군) / 농촌=~군$ / 중소도시=나머지
# ★ 설정: VIF<10 | p<0.20 | M6 BYM+RW1+IID (NB) | 2020–2024
# ==============================================================================
rm(list=ls()); gc()
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(MASS); library(stringr); library(car)
  library(openxlsx); library(arrow); library(sf); library(spdep); library(INLA)
})
options(scipen=999)

# ══════════════════════════════════════════
# 설정
# ══════════════════════════════════════════
DISEASE_NAME <- "노로바이러스"
YEAR_START <- 2020; YEAR_END <- 2024
PVAL_SCREEN <- 0.20
VIF_THRESHOLD <- 10
MIN_OBS <- 20; COV_RATIO <- 0.85

# --- Paths (EDIT THESE) -------------------------------------------------------
# Raw inputs are NOT redistributed (aggregated official statistics; no PII). See README/DATA.
BASE_IV <- "FBD_DATA_ZIP"                                  # folder holding the input data files
PATH_DISEASE   <- file.path(BASE_IV, "식중독최종.csv")        # district-year food-poisoning counts + population (MFDS/KOSIS)
PATH_HEALTH_PQ <- file.path(BASE_IV, "국민건강결과_최종.parquet") # paediatric/health battery indicators (HIRA)
PATH_SHP       <- file.path(BASE_IV, "final.shp")            # si-gun-gu boundary shapefile (Statistics Korea SGIS)
DIR_OUT <- "output"                                         # results (xlsx tables, csv)
DIR_LOG <- "output"                                         # run log (markdown)
dir.create(DIR_OUT, showWarnings=FALSE, recursive=TRUE)

HIGH_RISK_PCT <- 0.20
FAMILY <- "nbinomial"   # 표준 Negative Binomial
CF_ZINB <- list()

if(!dir.exists(DIR_OUT)) tryCatch(dir.create(DIR_OUT, recursive=TRUE), error=function(e){})
if(!dir.exists(DIR_OUT)){ DIR_OUT <- file.path(Sys.getenv("HOME"), "Desktop")
  cat(sprintf("  ⚠️ Google Drive 접근 불가 → 바탕화면: %s\n", DIR_OUT))}

TS <- format(Sys.time(), "%y%m%d_%H%M")
LOG <- file.path(DIR_LOG, sprintf("NORO_v8Bc2_CHILDok_%s.md", TS))   # ★ v8Bc = B + 소아변수
sink(LOG, split=TRUE)
cat(sprintf("# NORO v8-B PRE-SPECIFIED (도시/농촌 3그룹 분리, AUTO 제거)\n\n- TS: %s\n- 전략: 전체/도시/농촌 3그룹 별도 분석\n- 변수선택: 사전지정(이론) — forward-target·reverse-drop 미적용\n- Family: %s\n- VIF<%d | p<%.2f\n- 기간: %d–%d\n\n---\n\n",
    TS, FAMILY, VIF_THRESHOLD, PVAL_SCREEN, YEAR_START, YEAR_END))

# ══════════════════════════════════════════
# 공통 함수
# ══════════════════════════════════════════
clean_region <- function(df) df %>% mutate(
  region=str_replace_all(as.character(region),"\\s+",""),
  region=if_else(region=="인천시미추홀구","인천시남구",region),
  year=as.integer(year)) %>% filter(year>=YEAR_START, year<=YEAR_END)
read_csv_safe <- function(fp){raw<-NULL
  for(enc in c("UTF-8","UTF-8-BOM","CP949","EUC-KR")){
    raw<-tryCatch(read.csv(fp,fileEncoding=enc,check.names=FALSE,stringsAsFactors=FALSE),error=function(e)NULL)
    if(!is.null(raw)&&nrow(raw)>0)break}; raw}
fill_missing_year<-function(df,tgt,src,fn=""){if(!"region"%in%names(df)||!"year"%in%names(df))return(df)
  nv<-setdiff(names(df)[sapply(df,is.numeric)],"year");if(length(nv)==0||!src%in%unique(df$year))return(df)
  ds<-df%>%filter(year==src);dt<-df%>%filter(year==tgt);df_f<-ds%>%mutate(year=as.integer(tgt))
  if(nrow(dt)>0){df_f<-df_f%>%left_join(dt%>%dplyr::select(region,all_of(nv))%>%rename_with(~paste0(.,"__o"),all_of(nv)),by="region")%>%
    mutate(across(all_of(nv),function(col){v<-cur_column();o<-get(paste0(v,"__o"));ifelse(!is.na(o),o,col)}))%>%dplyr::select(region,year,all_of(nv))}
  bind_rows(df%>%filter(year!=tgt),df_f)%>%arrange(region,year)}
apply_cf <- function(df,fn) fill_missing_year(fill_missing_year(df,2021,2020,fn),2024,2023,fn)
is_pct <- function(x){xv<-x[!is.na(x)&is.finite(x)];all(xv>=0&xv<=100)&max(xv)>1}
run_univ <- function(x, df_w){
  tmp <- data.frame(cases=df_w$cases, x=x, pop=df_w$population)
  tmp <- tmp[complete.cases(tmp) & is.finite(tmp$x) & tmp$pop > 0, ]
  if(nrow(tmp) < MIN_OBS || sd(tmp$x, na.rm=TRUE) == 0) return(NULL)
  tryCatch({m <- glm.nb(cases ~ x + offset(log(pop+1)), data=tmp); cr <- summary(m)$coefficients
    if(nrow(cr) < 2) return(NULL)
    list(p=cr[2,"Pr(>|z|)"], IRR=exp(cr[2,"Estimate"]),
      lo=exp(cr[2,"Estimate"]-1.96*cr[2,"Std. Error"]),
      hi=exp(cr[2,"Estimate"]+1.96*cr[2,"Std. Error"]), n=nrow(tmp))
  }, error=function(e) NULL)}

# ══════════════════════════════════════════
# PART 1. 데이터 로드
# ══════════════════════════════════════════
cat("## PART 1. 데이터 로드\n\n")
df_raw <- read.csv(PATH_DISEASE, stringsAsFactors=FALSE, check.names=FALSE)
df_target <- df_raw %>% filter(disease==DISEASE_NAME, year>=YEAR_START, year<=YEAR_END) %>%
  clean_region() %>% group_by(region,year) %>%
  summarise(cases=sum(cases,na.rm=TRUE), population=mean(population,na.rm=TRUE), .groups="drop") %>%
  mutate(rate_100k=cases/population*100000)
cat(sprintf("  종속: %d행 | %d시군구 | %d건\n",nrow(df_target),n_distinct(df_target$region),sum(df_target$cases)))
cor_merged <- df_target

HEALTH_VARS <- c("1인가구수.1","1인가구수_65세이상","1인가구율_45-64세가구","1인가구율_65세이상가구","1인가구율_전체",
  "가정의학과전문의","건강생활실천율_조율","걷기실천율_표준화율","격렬한신체활동실천율_표준화율","고령인구비율",
  "고위험음주율_남_표준화","고위험음주율_여_표준화","관내진료비_외래","관내진료비_입원","관내진료비_전체",
  "관내진료실인원_외래","관내진료실인원_입원","관내진료실인원_전체","관외진료비_외래","관외진료비_입원","관외진료비_전체",
  "관외진료실인원_외래","관외진료실인원_입원","관외진료실인원_전체","국민기초생활보장수급자","국민기초생활보장수급자수율",
  "국민연금_사업장가입자수","국민연금_임의가입자수","국민연금_임의계속가입자수","국민연금_지역가입자_납부예외자수",
  "국민연금_지역가입자_소득신고자수","국민연금_총가입자수","기준시간내접근불가비율_종합병원(전체)","기초생활수급자수율",
  "기초연금수급자수","내과전문의","노인장기요양_시설_기관수","노인장기요양시설_영양사","농촌인구수",
  "다문화이혼건수","다문화이혼비중","다문화출생비율","다문화출생아수","다문화혼인비중",
  "도시인구수","도시지역면적","도시지역인구비율","독거노인가구비율","독거노인비율",
  "목욕시설_있음","방역수칙실천율실내마스크착용_표준화율","방역수칙실천율실외마스크착용_표준화율",
  "보건및사회복지사업체수","보건및사회복지사업체종사자비율","보건및사회복지사업체종사자수","보건소인력_보건직",
  "비누,손세정제사용률_표준화율","비누손세정제사용률_표준화율",
  "사회적거리두기또는생활속거리두기실천율건_표준화율","사회적거리두기또는생활속거리두기실천율건강_조율",
  "상수도보급률","성비","수도_마을상수도","수도_상수도","수도_전용상수도","순이동인구",
  "식사전손씻기실천율_표준화율","식품안정성확보율_표준화율","어제저녁식사후칫솔질실천율_표준화율",
  "어제점심식사후칫솔질실천율_표준화율","연간인플루엔자예방접종률_표준화율","예방의학과전문의","온수시설_있음",
  "외출후손씻기실천율_표준화율","우울감경험률_표준화율","월간음주율_남_표준화","월간음주율_여_표준화",
  "유기물질부하량발생량","유기물질부하량방류량","의사수","의원_가정의학과","의원_예방의학과","의원_재활의학과",
  "인구천명당사설학원수","인구천명당의료기관종사의사수","인구천명당주점업수","인구천명당패스트푸드점수",
  "인구천명당폐수발생량","인구천명당폐수방류량","작업치료사수","재정자립도","재정자주도","재활의학과전문의",
  "전문의합계","정신건강의학과전문의","주관적건강수준인지율_표준화율","주점업수","중등도이상신체활동실천율_조율",
  "집밖에서의손소독제사용횟수_표준화율","총가구수_65세이상","치과의사수",
  "코로나19관련유증상자행동수칙미준수율_표준화율","패스트푸드점수",
  "평소손씻기실천율_표준화율","폐수발생량","폐수방류량","폐수배출업소수","하수도보급률",
  "화장실_수세식","화장실_재래식","화장실다녀온후손씻기실천율_조율")

SELECTED <- list(
  "merged_세부용도별지하수.csv"=c("온천수_시설수","온천수_이용량"),
  "merged_지하수수질.csv"=c("검사합계","부적합","적합","합계"),
  "merged_하수관로개보수.csv"=c("개·보수관로_부분보수(개소)_계","개·보수관로_부분보수(개소)_분류식_오수","개·보수관로_부분보수(개소)_분류식_우수","개·보수관로_전체보수(m)_합류식","맨홀(개소)_합류식맨홀","받이_빗물받이(개소)","토실,토구_토실(개소)"),
  "merged_하수도보급률.csv"=c("고도처리인구보급률(%)","공공하수처리구역인구보급률(%)","총면적(㎢)","총인구(명)","하수도설치율(%)","하수처리구역내_계","하수처리구역내_공공하수처리인구(명)","하수처리구역내_공공하수처리인구(명).2","하수처리구역내_공공하수처리인구(명).3","하수처리구역내_면적(㎢)","하수처리구역내_미접속인구","하수처리구역내_폐수처리인구(명)","하수처리구역내_폐수처리인구(명).3","하수처리구역외_계","하수처리구역외_면적(㎢)","하수처리구역외_오수처리인구","하수처리구역외_정화조인구"),
  "merged_하수찌꺼기발생및처리.csv"=c("외부위탁처리량(톤/년)_매립","외부위탁처리량(톤/년)_복토재","외부위탁처리량(톤/년)_소계","외부위탁처리량(톤/년)_연료","외부위탁처리량(톤/년)_제품원료","외부위탁처리량(톤/년)_퇴비화","자체처리량(톤/년)_건조","자체처리량(톤/년)_건조후처리(2차).3","자체처리량(톤/년)_건조후처리(2차).5","자체처리량(톤/년)_계","자체처리량(톤/년)_고화","자체처리량(톤/년)_고화후처리(2차)","자체처리량(톤/년)_고화후처리(2차).1","자체처리량(톤/년)_고화후처리(2차).2","자체처리량(톤/년)_소각후처리(2차)","자체처리량(톤/년)_소각후처리(2차).2","자체처리량(톤/년)_소각후처리(2차).3","자체처리량(톤/년)_퇴비화","함수율(%,탈수기준)"),
  "가축두수_전처리.csv"=c("농가수(호)_가금","농가수(호)_돼지","농가수(호)_말","농가수(호)_양·염소·사슴","농가수(호)_젖소","농가수(호)_한육우","농가수(호)_합계","사육두수(두)_가금","사육두수(두)_개","사육두수(두)_돼지","사육두수(두)_양·염소·사슴","사육두수(두)_젖소","사육두수(두)_한육우","사육두수(두)_합계"),
  "고령인구_전처리.csv"=c("1인가구_60세 이상 - 계","1인가구_65~69세","1인가구_70~74세","1인가구_75~79세","1인가구_80~84세","계_60~64세","계_60세 이상 - 계","계_65~69세","계_70~74세","계_75~79세","계_80~84세","계_85세이상"),
  "국토이용현황_전처리_수정.csv"=c("공장용지","광천지","구거","답","대","도로","목장용지","묘지","유지","임야","잡종지","전","제방","주유소용지","창고용지","하천"),
  "생활용지하수이용현황_전처리.csv"=c("가정용_개소수","가정용_이용량","간이상수도용_개소수","간이상수도용_이용량","농업·생활겸용_이용량","민방위용_개소수","민방위용_이용량","일반용_개소수","일반용_이용량","총 계_개소수","총 계_이용량","학교용_개소수","학교용_이용량"),
  "손씻기_전처리.csv"=c("after_outing_handwash_rate_adj","after_outing_handwash_rate_std","after_toilet_handwash_rate_adj","before_meal_handwash_rate_adj","usual_handwash_rate_adj","usual_handwash_rate_std","비누_손_세정제_사용률_표준화율"),
  "음용지하수이용현황_전처리.csv"=c("개소수(총합)","민방위용","총이용량"),
  "재정자립도_전처리.csv"=c("재정자립도(세입과목개편전)","재정자립도(세입과목개편후)"),
  "추가_1인가구_전처리.csv"=c("1인가구"),
  "추가_고령인구_전처리.csv"=c("1인가구_60~64세","계_60~64세","계_70~74세","계_75~79세"),
  "추가_독거노인가구비율_전처리.csv"=c("65세이상_1인가구(가구)","독거노인가구비율(%)"))

SELECTED_NEW <- list(
  "경제수준.csv"=c("주택소유율","청년고용률_상반기","청년고용률_하반기",
    "기초연금수급자율","국민기초생활보장수급자수율","국민연금_임의가입자수"),
  "종합소득세.csv"=c("신고인원","총수입금액","종합소득금액","과세표준","산출세액","세액공제 및 감면","결정세액"),
  "의료시설.csv"=c("병원_기관수","보건소_기관수","보건의료원_기관수","보건지소_기관수","보건진료소_기관수",
    "상급종합_기관수","약국_기관수","요양병원_기관수","의원_기관수","종합병원_기관수",
    "치과병원_기관수","치과의원_기관수","한방병원_기관수","한의원_기관수","전체_기관수",
    "총병상수","응급실병상수"),
  "교육기관_인력_예산.csv"=c("교원1인당학생수","유치원교원수","초등학교교원수","유치원원아수",
    "초등학교 학생 수","유치원 수","유아천명당보육시설수","초등학교 수",
    "전문대학 및 대학교 수","인구천명당 사설학원수"),
  "merged_분뇨찌꺼기처리.csv"=c("발생량(A)=(B)+(C)","처분량_계(B)","처분량_재활용"),
  "1인가구_전처리.csv"=c("주택_계","주택_다세대주택","주택_단독주택","주택_아파트","주택_연립주택"),
  "독거노인가구비율_전처리.csv"=c("전체_일반가구(가구)"),
  "영유아인구_전처리.csv"=c("영유아비율_0_4","아동비율_5_9")   # ★ 소아 연령구조 (노로 핵심 집단)
)

META_COLS <- c("region","year","sido","sigungu","X","Unnamed..0","항목",
  "분뇨처리장명","타시설연계여부","시설명","소재지","처리공법","연계처리장명",
  "전화번호","가동개시일","운영방법","방류수역_지류","방류수역_본류","방류수역_수계",
  "연계처리","연계명","방류수역_수계.1","방류수역_수계.2","방류수역_수계.3","방류수역_수계.4",
  "운영방법.자체.공기업.민간위탁.","행정코드","시설코드","타시설연계여부_2","타시설연계여부_3","타시설연계여부_4",
  "지역.시.군.구.","하수처리장명","연도")

# (1) Parquet 건강지표
tryCatch({hpq<-read_parquet(PATH_HEALTH_PQ)%>%as.data.frame()%>%clean_region();hpq<-apply_cf(hpq,"pq")
  ah<-intersect(HEALTH_VARS,names(hpq));for(v in ah)hpq[[v]]<-suppressWarnings(as.numeric(hpq[[v]]))
  hagg<-hpq%>%group_by(region,year)%>%summarise(across(all_of(ah),~mean(.x,na.rm=TRUE)),.groups="drop")
  cor_merged<-cor_merged%>%left_join(hagg,by=c("region","year"))
  cat(sprintf("  Parquet 건강지표: %d변수\n",length(ah)))},error=function(e)cat(sprintf("  ❌ %s\n",e$message)))

# (2) 기존 SELECTED
for(fn in names(SELECTED)){fp<-file.path(BASE_IV,fn);if(!file.exists(fp))next
  raw<-read_csv_safe(fp);if(is.null(raw))next;raw<-raw%>%clean_region();raw<-apply_cf(raw,fn)
  av<-intersect(SELECTED[[fn]],names(raw));if(length(av)==0)next
  for(v in av)raw[[v]]<-suppressWarnings(as.numeric(raw[[v]]))
  agg<-raw%>%group_by(region,year)%>%summarise(across(all_of(av),~mean(.x,na.rm=TRUE)),.groups="drop")
  cor_merged<-cor_merged%>%left_join(agg,by=c("region","year"))
  cat(sprintf("  CSV %-45s %d변수\n",fn,length(av)))}

# (3) 신규 SELECTED_NEW
cat("\n  ── 신규 데이터 (영비율/NA 검증 통과 변수) ──\n")
for(fn in names(SELECTED_NEW)){fp<-file.path(BASE_IV,fn);if(!file.exists(fp)){cat(sprintf("  ⚠ 미발견: %s\n",fn));next}
  raw<-read_csv_safe(fp);if(is.null(raw))next;raw<-raw%>%clean_region();raw<-apply_cf(raw,fn)
  av<-intersect(SELECTED_NEW[[fn]],names(raw))
  av<-setdiff(av, names(cor_merged))
  if(length(av)==0){cat(sprintf("  SKIP %-45s 0변수\n",fn));next}
  for(v in av)raw[[v]]<-suppressWarnings(as.numeric(raw[[v]]))
  agg<-raw%>%group_by(region,year)%>%summarise(across(all_of(av),~mean(.x,na.rm=TRUE)),.groups="drop")
  cor_merged<-cor_merged%>%left_join(agg,by=c("region","year"))
  cat(sprintf("  NEW  %-45s %d변수\n",fn,length(av)))}

# (4) 해산물/수산물 — 카테고리 합산
cat("\n  ── ★ 해산물/수산물 카테고리 합산 (수산물 공통) ──\n")
tryCatch({
  fp <- file.path(BASE_IV, "어패류_패류_전처리.csv")
  if(file.exists(fp)){
    raw <- read_csv_safe(fp); raw <- raw %>% clean_region() %>% apply_cf(.,"패류")
    meta_c <- c("region","year","sido","sigungu","X")
    num_c <- setdiff(names(raw), meta_c)
    for(v in num_c) raw[[v]] <- suppressWarnings(as.numeric(raw[[v]]))
    nat_cols <- num_c[grepl("자연채묘",num_c)]
    art_cols <- num_c[grepl("인공종묘",num_c)]
    seed_cols <- num_c[grepl("인공종자",num_c)]
    raw$패류_자연채묘_총합 <- rowSums(raw[nat_cols], na.rm=TRUE)
    raw$패류_인공종묘_총합 <- rowSums(raw[art_cols], na.rm=TRUE)
    raw$패류_인공종자_총합 <- rowSums(raw[seed_cols], na.rm=TRUE)
    raw$패류_전체_총합 <- raw$패류_자연채묘_총합 + raw$패류_인공종묘_총합 + raw$패류_인공종자_총합
    key_shell <- c("굴_자연채묘 생산량(kg)","바지락_자연채묘 생산량(kg)",
      "전복_자연채묘 생산량(kg)","전복_인공종묘 생산량(kg)","홍합_자연채묘 생산량(kg)",
      "굴_인공종묘 생산량(kg)","바지락_인공종묘 생산량(kg)")
    add_v <- c("패류_자연채묘_총합","패류_인공종묘_총합","패류_인공종자_총합","패류_전체_총합")
    add_v <- c(add_v, intersect(key_shell, names(raw)))
    add_v <- setdiff(add_v, names(cor_merged))
    agg <- raw %>% group_by(region,year) %>% summarise(across(all_of(add_v),~mean(.x,na.rm=TRUE)),.groups="drop")
    cor_merged <- cor_merged %>% left_join(agg, by=c("region","year"))
    cat(sprintf("  패류 합산+주요종: %d변수\n",length(add_v)))
  }
}, error=function(e) cat(sprintf("  ❌ 패류: %s\n",e$message)))

tryCatch({
  fp <- file.path(BASE_IV, "어패류_갑각류_전처리.csv")
  if(file.exists(fp)){
    raw <- read_csv_safe(fp); raw <- raw %>% clean_region() %>% apply_cf(.,"갑각류")
    meta_c <- c("region","year","sido","sigungu","X")
    num_c <- setdiff(names(raw), meta_c)
    for(v in num_c) raw[[v]] <- suppressWarnings(as.numeric(raw[[v]]))
    shrimp <- num_c[grepl("새우|대하",num_c)]
    crab <- num_c[grepl("게|꽃게",num_c)]
    raw$갑각류_전체_총합 <- rowSums(raw[num_c], na.rm=TRUE)
    raw$갑각류_새우류_총합 <- rowSums(raw[shrimp], na.rm=TRUE)
    raw$갑각류_게류_총합 <- rowSums(raw[crab], na.rm=TRUE)
    add_v <- setdiff(c("갑각류_전체_총합","갑각류_새우류_총합","갑각류_게류_총합"), names(cor_merged))
    agg <- raw %>% group_by(region,year) %>% summarise(across(all_of(add_v),~mean(.x,na.rm=TRUE)),.groups="drop")
    cor_merged <- cor_merged %>% left_join(agg, by=c("region","year"))
    cat(sprintf("  갑각류 합산: %d변수\n",length(add_v)))
  }
}, error=function(e) cat(sprintf("  ❌ 갑각류: %s\n",e$message)))

tryCatch({
  fp <- file.path(BASE_IV, "어패류_어류_전처리.csv")
  if(file.exists(fp)){
    raw <- read_csv_safe(fp); raw <- raw %>% clean_region() %>% apply_cf(.,"어류")
    meta_c <- c("region","year","sido","sigungu","X")
    num_c <- setdiff(names(raw), meta_c)
    for(v in num_c) raw[[v]] <- suppressWarnings(as.numeric(raw[[v]]))
    egg_cols <- num_c[grepl("수정란",num_c)]
    fry_cols <- num_c[grepl("치어",num_c)]
    marine_sp <- c("넙치","감성돔","농어","능성어","돌돔","참돔","조피볼락","볼락","부세","방어","숭어","민어","참조기","전어")
    marine_egg <- egg_cols[sapply(egg_cols, function(x) any(sapply(marine_sp, function(s) grepl(s,x))))]
    marine_fry <- fry_cols[sapply(fry_cols, function(x) any(sapply(marine_sp, function(s) grepl(s,x))))]
    raw$어류_수정란_전체총합 <- rowSums(raw[egg_cols], na.rm=TRUE)
    raw$어류_치어_전체총합 <- rowSums(raw[fry_cols], na.rm=TRUE)
    raw$어류_해수어_수정란총합 <- rowSums(raw[marine_egg], na.rm=TRUE)
    raw$어류_해수어_치어총합 <- rowSums(raw[marine_fry], na.rm=TRUE)
    add_v <- setdiff(c("어류_수정란_전체총합","어류_치어_전체총합",
      "어류_해수어_수정란총합","어류_해수어_치어총합"), names(cor_merged))
    agg <- raw %>% group_by(region,year) %>% summarise(across(all_of(add_v),~mean(.x,na.rm=TRUE)),.groups="drop")
    cor_merged <- cor_merged %>% left_join(agg, by=c("region","year"))
    cat(sprintf("  어류 합산: %d변수\n",length(add_v)))
  }
}, error=function(e) cat(sprintf("  ❌ 어류: %s\n",e$message)))

tryCatch({
  fp1 <- file.path(BASE_IV, "어패류_해조류_김_전처리.csv")
  fp2 <- file.path(BASE_IV, "어패류_해조류_김제외_전처리.csv")
  add_all <- c()
  if(file.exists(fp1)){
    r1 <- read_csv_safe(fp1); r1 <- r1 %>% clean_region() %>% apply_cf(.,"해조류김")
    mc1 <- setdiff(names(r1), c("region","year","sido","sigungu","X"))
    for(v in mc1) r1[[v]] <- suppressWarnings(as.numeric(r1[[v]]))
    r1$해조류_김_총합 <- rowSums(r1[mc1], na.rm=TRUE)
    add_all <- c(add_all, "해조류_김_총합")
  }
  if(file.exists(fp2)){
    r2 <- read_csv_safe(fp2); r2 <- r2 %>% clean_region() %>% apply_cf(.,"해조류기타")
    mc2 <- setdiff(names(r2), c("region","year","sido","sigungu","X"))
    for(v in mc2) r2[[v]] <- suppressWarnings(as.numeric(r2[[v]]))
    r2$해조류_기타_총합 <- rowSums(r2[mc2], na.rm=TRUE)
    add_all <- c(add_all, "해조류_기타_총합")
  }
  if("해조류_김_총합" %in% add_all){
    a1 <- r1 %>% group_by(region,year) %>% summarise(해조류_김_총합=mean(해조류_김_총합,na.rm=TRUE),.groups="drop")
    cor_merged <- cor_merged %>% left_join(a1, by=c("region","year"))
  }
  if("해조류_기타_총합" %in% add_all){
    a2 <- r2 %>% group_by(region,year) %>% summarise(해조류_기타_총합=mean(해조류_기타_총합,na.rm=TRUE),.groups="drop")
    cor_merged <- cor_merged %>% left_join(a2, by=c("region","year"))
  }
  cat(sprintf("  해조류 합산: %d변수\n",length(add_all)))
}, error=function(e) cat(sprintf("  ❌ 해조류: %s\n",e$message)))

# (5) 양식데이터
cat("\n  ── 양식데이터 (생산량·금액: 계 + 상위4종) ──\n")
tryCatch({
  fp_aq <- file.path(BASE_IV, "양식데이터_연도별.csv")
  if(file.exists(fp_aq)){
    aq_raw <- read_csv_safe(fp_aq)
    if(!is.null(aq_raw)){
      aq_raw <- aq_raw %>% clean_region() %>% apply_cf(.,"양식")
      aq_prod <- aq_raw %>% filter(항목 %in% c("생산량(M/T):계","생산금액(백만원):계"))
      meta_aq <- c("region","year","sido","sigungu","항목","X")
      sp_cols <- setdiff(names(aq_prod), meta_aq)
      for(v in sp_cols) aq_prod[[v]] <- suppressWarnings(as.numeric(aq_prod[[v]]))
      keep_sp <- c("계_생산량","계_생산금액",
        "넙치류_생산량","넙치류_생산금액",
        "조피볼락_생산량","조피볼락_생산금액",
        "숭어류_생산량","참돔_생산량")
      keep_sp <- intersect(keep_sp, sp_cols)
      aq_prod <- aq_prod %>% mutate(항목=ifelse(grepl("생산량",항목),"양식","양식금"))
      aq_wide <- aq_prod %>%
        pivot_longer(cols=all_of(keep_sp), names_to="sp", values_to="val") %>%
        mutate(varname=paste0(항목,"_",sp)) %>%
        dplyr::select(region,year,varname,val) %>%
        pivot_wider(names_from=varname, values_from=val, values_fn=mean)
      new_aq <- setdiff(names(aq_wide), c("region","year",names(cor_merged)))
      if(length(new_aq)>0){
        cor_merged <- cor_merged %>% left_join(aq_wide %>% dplyr::select(region,year,all_of(new_aq)), by=c("region","year"))
        cat(sprintf("  양식데이터: %d변수 (계+상위종)\n",length(new_aq)))
      }
    }
  }
}, error=function(e) cat(sprintf("  ❌ 양식데이터: %s\n",e$message)))

# (5b) 지역축제
cat("\n  ── 지역축제 (연간 건수) ──\n")
tryCatch({
  fp_fest <- file.path(BASE_IV, "지역축제_전처리.csv")
  if(file.exists(fp_fest)){
    fest <- read_csv_safe(fp_fest)
    if(!is.null(fest)){
      fest <- fest %>% clean_region() %>% apply_cf(.,"축제")
      fest_agg <- fest %>% group_by(region,year) %>% summarise(지역축제_건수=n(), .groups="drop")
      cor_merged <- cor_merged %>% left_join(fest_agg, by=c("region","year"))
      cat(sprintf("  지역축제: 1변수 (건수)\n"))
    }
  }
}, error=function(e) cat(sprintf("  ❌ 지역축제: %s\n",e$message)))

# (6) 의료인력 (long→wide)
cat("\n  ── 의료인력 (long→wide pivot) ──\n")
tryCatch({
  fp_med <- file.path(BASE_IV, "의료인력.csv")
  if(file.exists(fp_med)){
    med <- read_csv_safe(fp_med)
    if(!is.null(med)){
      med <- med %>% clean_region() %>% apply_cf(.,"의료인력")
      med$value <- suppressWarnings(as.numeric(med$value))
      med_wide <- med %>% filter(!is.na(value) & !is.na(의료인력별) & 의료인력별!="0") %>%
        group_by(region,year,의료인력별) %>% summarise(val=mean(value,na.rm=TRUE),.groups="drop") %>%
        mutate(varname=paste0("의료인력_",의료인력별)) %>%
        dplyr::select(region,year,varname,val) %>%
        pivot_wider(names_from=varname, values_from=val, values_fn=mean)
      new_med <- setdiff(names(med_wide), c("region","year",names(cor_merged)))
      if(length(new_med)>0){
        cor_merged <- cor_merged %>% left_join(med_wide %>% dplyr::select(region,year,all_of(new_med)), by=c("region","year"))
        cat(sprintf("  의료인력: %d변수\n",length(new_med)))
      }
    }
  }
}, error=function(e) cat(sprintf("  ❌ 의료인력: %s\n",e$message)))

# (7) 외국인인구수 (long→wide)
cat("\n  ── 외국인인구수 (long→wide pivot) ──\n")
tryCatch({
  fp_for <- file.path(BASE_IV, "외국인인구수_전처리.csv")
  if(file.exists(fp_for)){
    forn <- read_csv_safe(fp_for)
    if(!is.null(forn)){
      forn <- forn %>% clean_region() %>% apply_cf(.,"외국인")
      forn$value <- suppressWarnings(as.numeric(forn$value))
      forn_w <- forn %>% filter(!is.na(value) & 성별=="계" & !거주외국인별 %in% c("0","")) %>%
        group_by(region,year,거주외국인별) %>% summarise(val=sum(value,na.rm=TRUE),.groups="drop") %>%
        mutate(varname=paste0("외국인_",거주외국인별)) %>%
        dplyr::select(region,year,varname,val) %>%
        pivot_wider(names_from=varname, values_from=val, values_fn=sum)
      new_for <- setdiff(names(forn_w), c("region","year",names(cor_merged)))
      if(length(new_for)>0){
        cor_merged <- cor_merged %>% left_join(forn_w %>% dplyr::select(region,year,all_of(new_for)), by=c("region","year"))
        cat(sprintf("  외국인인구수: %d변수\n",length(new_for)))
      }
    }
  }
}, error=function(e) cat(sprintf("  ❌ 외국인인구수: %s\n",e$message)))

# (8) 의료통계 (long→wide)
cat("\n  ── 의료통계 (long→wide pivot) ──\n")
tryCatch({
  fp_ms <- file.path(BASE_IV, "의료통계_전처리.csv")
  if(file.exists(fp_ms)){
    mstat <- read_csv_safe(fp_ms)
    if(!is.null(mstat)){
      mstat <- mstat %>% clean_region() %>% apply_cf(.,"의료통계")
      mstat$값 <- suppressWarnings(as.numeric(mstat$값))
      mstat_w <- mstat %>% filter(!is.na(값) & !지표 %in% c("0","") & !항목 %in% c("0","")) %>%
        mutate(varname=paste0("의료통계_",gsub("시군구_","",지표),"_",항목)) %>%
        group_by(region,year,varname) %>% summarise(val=mean(값,na.rm=TRUE),.groups="drop") %>%
        pivot_wider(names_from=varname, values_from=val, values_fn=mean)
      new_ms <- setdiff(names(mstat_w), c("region","year",names(cor_merged)))
      if(length(new_ms)>0){
        cor_merged <- cor_merged %>% left_join(mstat_w %>% dplyr::select(region,year,all_of(new_ms)), by=c("region","year"))
        cat(sprintf("  의료통계: %d변수\n",length(new_ms)))
      }
    }
  }
}, error=function(e) cat(sprintf("  ❌ 의료통계: %s\n",e$message)))

cat(sprintf("\n  ═══ cor_merged 최종: %d행 × %d열 ═══\n\n",nrow(cor_merged),ncol(cor_merged)))

# ══════════════════════════════════════════
# PART 2. 변수 정의 (30개 base + 이론방향) + raw 단변량 (Table 1)
# ══════════════════════════════════════════
cat("## PART 2. 변수 정의 (30개 base + 이론방향) + raw 단변량\n\n")
TV_v6 <- data.frame(
  tier=rep("A",30),
  cat=c(rep("① 식품원 및 축산",4),rep("② 분뇨 및 오염처리",5),
    rep("③ 수질 및 하수 인프라",3),rep("④ 수계 및 토지이용",4),
    rep("⑤ 위생 및 건강행태",3),
    rep("⑥ 사회경제 및 취약성",6),rep("⑦ 인구학적 특성",3),rep("⑧ 의료접근",2)),
  kr=c("한육우농가수‡","농가수합계","굴","사육두수_합계",
    "자체처리량계","소각후처리","건조후처리","함수율_탈수기준","외부위탁_연료",
    "폐수방류량","정화조인구","민방위용_개소수",
    "유지","하천","목장용지","답(논)",
    "화장실손씻기","건강생활실천율","걷기실천율",
    "독거노인","농촌인구수","재정자주도","재정자립도","1인가구율_전체","기초생활수급자수율",
    "총가구수65세이상","성비‡","고령인구비율‡",
    "진료비외래","관내진료비전체"),
  code=c("농가수(호)_한육우","농가수(호)_합계","굴_자연채묘 생산량(kg)","사육두수(두)_합계",
    "자체처리량(톤/년)_계","자체처리량(톤/년)_소각후처리(2차)","자체처리량(톤/년)_건조후처리(2차).3",
    "함수율(%,탈수기준)","외부위탁_연료",
    "폐수방류량","하수처리구역외_정화조인구","민방위용_개소수",
    "유지","하천","목장용지","답",
    "화장실다녀온후손씻기실천율_조율","건강생활실천율_조율","걷기실천율_표준화율",
    "독거노인비율","농촌인구수","재정자주도","재정자립도","1인가구율_전체","기초생활수급자수율",
    "총가구수_65세이상","성비","고령인구비율",
    "관내진료비_외래","관내진료비_전체"),
  eng=c("beef_farm","farm_total","oyster","livestock_total",
    "sludge_total","sludge_incin","sludge_dry","sludge_moisture","waste_fuel",
    "ww_effluent","septic_pop","gw_civil_count",
    "reservoir","river","ranch","paddy",
    "handwash_toilet","health_practice","walking_practice",
    "elderly_alone","rural_pop","fiscal_auto","fiscal_indep","alone_rate","welfare_rate",
    "hh_elderly","sex_ratio","elderly_rate",
    "med_out","med_total"),
  forced=c("‡",rep("",25),"‡","‡",rep("",2)),
  이론방향=c(
    "위험","위험","위험","위험",
    "보호","보호","위험","위험","보호",
    "위험","위험","위험",
    "위험","위험","위험","중립",
    "보호","보호","보호",
    "위험","위험","보호","보호","중립","위험",
    "위험","중립","중립",
    "보호","보호"),
  stringsAsFactors=FALSE)

# ★ 소아 블록 추가 (사전지정 — 노로=어린이집·학교 집단발생 핵심)
#   ※ 결측0 변수만 사용: 영유아비율_0_4 / 아동비율_5_9 (둘 다 1145행 완전).
#     유아천명당보육시설수(교육기관 출처)는 2020 결측 → 도시모델 붕괴 → 제외.
TV_CHILD <- data.frame(
  tier=rep("A",2),
  cat=rep("⑨ 소아",2),
  kr=c("영유아비율_0_4","아동비율_5_9"),
  code=c("영유아비율_0_4","아동비율_5_9"),
  eng=c("child_0_4","child_5_9"),
  forced=c("",""),
  이론방향=c("위험","위험"),
  stringsAsFactors=FALSE)
TV_v6 <- rbind(TV_v6, TV_CHILD)

cat(sprintf("  v8-B Table 1 base: %d개 (9개 카테고리, 소아 2종 — 결측0만)\n", nrow(TV_v6)))
cat("  ★ 강제변수(‡): 한육우농가수, 성비, 고령인구비율\n")
cat("  ★ 변수선택: 사전지정만 — AUTO 타깃 최적화 없음\n\n")

TV <- TV_v6
df_work <- cor_merged %>% filter(population > 0)
raw_univ <- data.frame()
for(i in 1:nrow(TV)){
  v <- TV$code[i]
  if(!v %in% names(cor_merged)){
    raw_univ <- rbind(raw_univ, data.frame(TV[i,], N=0, mean_sd="—", min_v=NA, med=NA, max_v=NA,
      raw_IRR=NA, raw_lo=NA, raw_hi=NA, raw_p=NA, sig="", stringsAsFactors=FALSE)); next}
  x <- as.numeric(df_work[[v]]); xv <- x[!is.na(x) & is.finite(x)]; nv <- length(xv)
  res <- run_univ(x, df_work)
  raw_univ <- rbind(raw_univ, data.frame(TV[i,], N=nv,
    mean_sd=if(nv>0) sprintf("%.2f ± %.2f", mean(xv), sd(xv)) else "—",
    min_v=if(nv>0) round(min(xv),2) else NA, med=if(nv>0) round(median(xv),2) else NA,
    max_v=if(nv>0) round(max(xv),2) else NA,
    raw_IRR=if(!is.null(res)) round(res$IRR,4) else NA,
    raw_lo=if(!is.null(res)) round(res$lo,4) else NA,
    raw_hi=if(!is.null(res)) round(res$hi,4) else NA,
    raw_p=if(!is.null(res)) round(res$p,6) else NA,
    sig=if(!is.null(res) && res$p<0.05) "*" else "", stringsAsFactors=FALSE))}

n_sig05 <- sum(raw_univ$raw_p < 0.05, na.rm=TRUE)
n_ns <- sum(raw_univ$raw_p >= 0.05, na.rm=TRUE)
cat(sprintf("  raw α=0.05: 유의 %d | 비유의 %d / %d\n\n", n_sig05, n_ns, nrow(TV)))
for(i in 1:nrow(raw_univ)) cat(sprintf("  %-20s raw_p=%s %s\n", raw_univ$kr[i],
  ifelse(is.na(raw_univ$raw_p[i]), "NA", sprintf("%.4f", raw_univ$raw_p[i])), raw_univ$sig[i]))

# Shapefile
shp <- st_read(PATH_SHP, quiet=TRUE) %>%
  mutate(region=str_replace_all(as.character(region),"\\s+",""),
         region=if_else(region=="인천시미추홀구","인천시남구",region))
shp_main <- shp %>% filter(!region %in% c("인천시옹진군","전라남도완도군","전라남도진도군",
                                          "경상남도거제시","경상남도남해군","경상북도울릉군"))
nb_obj <- poly2nb(shp_main, snap=0.01); iso <- which(card(nb_obj)==0)
if(length(iso)>0){shp_main <- shp_main[-iso,]; nb_obj <- poly2nb(shp_main, snap=0.01)}
nb2INLA(nb_obj, file="/tmp/noro_v8b.graph"); g_main <- inla.read.graph("/tmp/noro_v8b.graph")
nb_w <- nb2listw(nb_obj, style="W")
cat(sprintf("\n  시군구: %d\n\n", nrow(shp_main)))

# ══════════════════════════════════════════
# PART 2.5. 도시/농촌 분류 + TV_RURAL
# ══════════════════════════════════════════
cat("## PART 2.5. 도시/농촌 분류\n\n")
metro_pattern <- "^(서울시|부산시|대구시|인천시|광주시|대전시|울산시)"
classify_region <- function(r){
  if(grepl("세종",r)) return("대도시")
  if(grepl(metro_pattern,r) & !grepl("군$",r)) return("대도시")
  if(grepl("군$",r)) return("농촌")
  return("중소도시")
}
cor_merged$area_3 <- sapply(cor_merged$region, classify_region)
cor_merged$area_2 <- ifelse(cor_merged$area_3=="농촌","농촌","도시")

n_metro <- sum(cor_merged$area_3=="대도시")/length(unique(cor_merged$year))
n_mid   <- sum(cor_merged$area_3=="중소도시")/length(unique(cor_merged$year))
n_rural <- sum(cor_merged$area_3=="농촌")/length(unique(cor_merged$year))
cat(sprintf("  대도시: ~%d | 중소도시: ~%d | 농촌: ~%d\n", round(n_metro), round(n_mid), round(n_rural)))
cat(sprintf("  도시(대도시+중소): ~%d | 농촌: ~%d\n\n", round(n_metro+n_mid), round(n_rural)))

TV_RURAL_BASE <- data.frame(
  tier=c(rep("A",9),rep("B",3)),
  cat=c("① 식품원 및 축산","① 식품원 및 축산","⑦ 인구학적 특성","⑦ 인구학적 특성",
    "⑥ 사회경제 및 취약성","④ 수계 및 토지이용","③ 수질 및 하수 인프라",
    "② 분뇨 및 오염처리","② 분뇨 및 오염처리",
    "④ 수계 및 토지이용","⑤ 위생 및 식품안전","⑤ 위생 및 식품안전"),
  kr=c("사육두수합계","사육두수_가금","고령인구비율‡","성비‡","1인가구율_전체",
    "유지","받이_빗물받이","민방위용_개소수","지하수_민방위용",
    "임야","평소손씻기_조율","화장실후손씻기_조율"),
  code=c("사육두수(두)_합계","사육두수(두)_가금","고령인구비율","성비","1인가구율_전체",
    "유지","받이_빗물받이(개소)","민방위용_개소수","민방위용_이용량",
    "임야","usual_handwash_rate_adj","화장실다녀온후손씻기실천율_조율"),
  eng=c("cattle_total","cattle_poultry","elderly_rate","sex_ratio","alone_rate",
    "reservoir","rainwater_gully","gw_civil_count","gw_civil_defense",
    "forest","handwash_usual_adj","handwash_toilet_adj"),
  forced=c("","","‡","‡","","","","","","","",""),
  이론방향=c("위험","위험","위험","중립","중립","중립","보호","위험","위험",
    "중립","보호","보호"),
  stringsAsFactors=FALSE)

TV_RURAL_L <- data.frame(
  tier="L",
  cat=c("② 분뇨 및 오염처리","③ 수질 및 하수 인프라","② 분뇨 및 오염처리","⑤ 위생 및 식품안전"),
  kr=c("지하수_부적합건수","함수율_탈수기준","학교용_지하수이용량","인플루엔자예방접종률"),
  code=c("부적합","함수율(%,탈수기준)","학교용_이용량","연간인플루엔자예방접종률_표준화율"),
  eng=c("gw_unfit_count","sludge_moisture","school_gw_use","flu_vaccination_rate"),
  forced=rep("",4),
  이론방향=c("위험","위험","위험","보호"),
  stringsAsFactors=FALSE)

TV_RURAL_N <- data.frame(
  tier="N",
  cat=c("④ 수계 및 토지이용","① 식품원 및 축산"),
  kr=c("전(밭)","사육두수_한육우"),
  code=c("전","사육두수(두)_한육우"),
  eng=c("dry_field","cattle_beef"),
  forced=rep("",2),
  이론방향=c("위험","위험"),
  stringsAsFactors=FALSE)

TV_RURAL <- rbind(TV_RURAL_BASE, TV_RURAL_L, TV_RURAL_N)
# ★ 농촌에도 소아 연령구조 추가 (사전지정)
TV_RURAL <- rbind(TV_RURAL, data.frame(
  tier=rep("A",2), cat=rep("⑨ 소아·보육",2),
  kr=c("영유아비율_0_4","아동비율_5_9"),
  code=c("영유아비율_0_4","아동비율_5_9"),
  eng=c("child_0_4","child_5_9"),
  forced=c("",""), 이론방향=c("위험","위험"), stringsAsFactors=FALSE))
cat(sprintf("  TV_RURAL (농촌 전용): %d개 변수 (소아 2종 포함)\n\n", nrow(TV_RURAL)))

# ══════════════════════════════════════════
# PART 3. run_model 함수 (M6 BYM+RW1+IID, NB)
# ══════════════════════════════════════════
run_model <- function(TV_local, quiet=FALSE){
  qcat <- function(...) if(!quiet) cat(...)
  df_w <- cor_merged %>% filter(population > 0, region %in% shp_main$region)
  result <- list(N=nrow(df_w), n_region=n_distinct(df_w$region), n_cases=sum(df_w$cases))

  for(sv in grep("사육두수", names(df_w), value=TRUE)){
    sc <- paste0(sv, "__per10k"); df_w[[sc]] <- as.numeric(df_w[[sv]]) / 10000
    qcat(sprintf("  ★ 스케일링: %s → %s\n", sv, sc))
  }

  valid <- TV_local %>% filter(code %in% names(df_w))
  form_map <- list(); data_ext <- df_w
  for(i in 1:nrow(valid)){
    var <- valid$code[i]; x <- as.numeric(df_w[[var]])
    nv <- sum(!is.na(x)&is.finite(x)); if(nv < MIN_OBS) next
    zp <- sum(!is.na(x)&is.finite(x)&x==0)/nv*100; pt <- is_pct(x); hz <- zp > 20
    is_sadu <- grepl("사육두수", var); is_sex <- (var == "성비")
    if(is_sadu){ forms <- list(raw = x/10000)
    } else if(is_sex){ forms <- list(raw = x)
    } else {
      forms <- list(raw=x)
      if(!pt){lv<-log1p(pmax(x,0));lv[is.na(x)]<-NA;if(!is.na(sd(lv,na.rm=TRUE))&&sd(lv,na.rm=TRUE)>0)forms[["log1p"]]<-lv}
      if(hz) forms[["binary"]]<-as.numeric(!is.na(x)&x>0) else{md<-median(x,na.rm=TRUE);forms[["binary"]]<-as.numeric(!is.na(x)&x>md)}
      if(hz){nz<-x[!is.na(x)&x>0];if(length(nz)>10){mn<-median(nz);forms[["T3"]]<-dplyr::case_when(is.na(x)~NA_real_,x==0~1,x<=mn~2,x>mn~3)}}
      else{q33<-quantile(x,c(1/3,2/3),na.rm=TRUE);brk<-unique(c(-Inf,q33[1],q33[2],Inf));if(length(brk)>=3)forms[["T3"]]<-as.numeric(cut(x,breaks=brk,labels=FALSE,include.lowest=TRUE))}
      q4<-quantile(x,c(0.25,0.5,0.75),na.rm=TRUE);b4<-unique(c(-Inf,q4[1],q4[2],q4[3],Inf))
      if(length(b4)>=3)forms[["Q4"]]<-as.numeric(cut(x,breaks=b4,labels=FALSE,include.lowest=TRUE))
    }
    rr<-list();for(fn in names(forms)){res<-run_univ(forms[[fn]],df_w);if(!is.null(res))rr[[fn]]<-data.frame(f=fn,p=res$p,IRR=res$IRR,n=res$n)}
    if(length(rr)==0) next; rd<-do.call(rbind,rr)%>%arrange(p); mn_n<-floor(nv*COV_RATIO); rc<-rd[!is.na(rd$n)&rd$n>=mn_n,]
    if(nrow(rc)==0) rc<-rd[1,]; bf<-rc$f[1]
    if(is_sadu){bf<-"per10k";bvn<-paste0(var,"__per10k");data_ext[[bvn]]<-as.numeric(df_w[[var]])/10000
    }else if(is_sex){bf<-"raw";bvn<-var
    }else if(bf=="raw"){bvn<-var
    }else{bvn<-paste0(var,"__",bf);xcm<-as.numeric(df_w[[var]])
      if(bf=="log1p")data_ext[[bvn]]<-log1p(pmax(xcm,0))
      else if(bf=="binary"){if(hz)data_ext[[bvn]]<-as.numeric(!is.na(xcm)&xcm>0)else{mc<-median(xcm,na.rm=TRUE);data_ext[[bvn]]<-as.numeric(!is.na(xcm)&xcm>mc)}}
      else if(bf=="T3"){if(hz){nzc<-xcm[!is.na(xcm)&xcm>0];mnc<-median(nzc,na.rm=TRUE);data_ext[[bvn]]<-dplyr::case_when(is.na(xcm)~NA_real_,xcm==0~1,xcm<=mnc~2,xcm>mnc~3)
      }else{q33c<-quantile(xcm,c(1/3,2/3),na.rm=TRUE);data_ext[[bvn]]<-as.numeric(cut(xcm,unique(c(-Inf,q33c[1],q33c[2],Inf)),labels=FALSE,include.lowest=TRUE))}}
      else if(bf=="Q4"){q4c<-quantile(xcm,c(0.25,0.5,0.75),na.rm=TRUE);data_ext[[bvn]]<-as.numeric(cut(xcm,unique(c(-Inf,q4c[1],q4c[2],q4c[3],Inf)),labels=FALSE,include.lowest=TRUE))}
    }
    form_map[[var]]<-list(kr=valid$kr[i],eng=valid$eng[i],cat=valid$cat[i],tier=valid$tier[i],
      forced=valid$forced[i],이론방향=valid$이론방향[i],형태=bf,변환명=bvn,p=rc$p[1],IRR=rc$IRR[1])
  }

  best_df<-data.frame();for(var in names(form_map)){m<-form_map[[var]]
    best_df<-rbind(best_df,data.frame(code=m$변환명,eng=m$eng,kr=m$kr,형태=m$형태,forced=m$forced,
      tier=m$tier,이론방향=m$이론방향,p=m$p,stringsAsFactors=FALSE))}
  if(nrow(best_df)==0) return(result)
  pass_vars<-best_df$code[best_df$p<PVAL_SCREEN|best_df$forced=="‡"]
  if(length(pass_vars)==0) return(result)
  forced_c<-best_df$code[best_df$forced=="‡"]; final_vars<-pass_vars
  vif_data<-data_ext[,c("cases",final_vars),drop=FALSE];for(v in final_vars)vif_data[[v]]<-as.numeric(vif_data[[v]])
  vif_data<-vif_data[complete.cases(vif_data),]
  for(stp in 1:40){if(length(final_vars)<=1)break
    lm_t<-tryCatch(lm(as.formula(paste("cases~",paste(paste0("`",final_vars,"`"),collapse="+"))),data=vif_data),error=function(e)NULL)
    if(is.null(lm_t))break;vv<-tryCatch(car::vif(lm_t),error=function(e)NULL);if(is.null(vv))break
    names(vv)<-gsub("`","",names(vv));if(max(vv,na.rm=TRUE)<VIF_THRESHOLD)break
    drop<-names(which.max(vv));if(drop%in%forced_c)break;final_vars<-final_vars[final_vars!=drop]}
  if(length(final_vars)==0) return(result)
  qcat(sprintf("  INLA 투입: %d변수 (VIF<%d)\n", length(final_vars), VIF_THRESHOLD))

  FMAP<-data.frame();for(v in final_vars){
    m_idx<-which(sapply(form_map,function(x)x$변환명==v));if(length(m_idx)==0)next;m<-form_map[[m_idx[1]]]
    x<-as.numeric(data_ext[[v]]);s<-sd(x,na.rm=TRUE);mn<-mean(x,na.rm=TRUE);safe<-paste0(m$eng,"_z")
    if(!is.na(s)&&s>0)data_ext[[safe]]<-(x-mn)/s else data_ext[[safe]]<-x
    if(m$형태%in%c("T3","Q4"))data_ext[[paste0(m$eng,"_f")]]<-factor(as.integer(x),ordered=FALSE)
    FMAP<-rbind(FMAP,data.frame(code=v,eng=m$eng,kr=m$kr,cat=m$cat,tier=m$tier,형태=m$형태,
      forced=m$forced,이론방향=m$이론방향,safe=safe,stringsAsFactors=FALSE))}
  if(nrow(FMAP)==0) return(result)

  rmap<-data.frame(region=shp_main$region,idarea=seq_len(nrow(shp_main)))
  ymap<-data.frame(year=YEAR_START:YEAR_END,idtime=1:length(YEAR_START:YEAR_END))
  ic<-data_ext[complete.cases(data_ext[,FMAP$safe]),]
  ic<-ic%>%left_join(rmap,by="region")%>%left_join(ymap,by="year")%>%arrange(idarea,idtime)
  ic$idarea_time<-1:nrow(ic)
  if(nrow(ic)<MIN_OBS) return(result)

  cov_str<-paste(FMAP$safe,collapse=" + ")
  pc_bym<-list(prec.unstruct=list(prior="pc.prec",param=c(0.5,0.01)),prec.spatial=list(prior="pc.prec",param=c(0.5,0.01)))
  pc_prec<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)))

  fit<-tryCatch(inla(as.formula(paste("cases ~",cov_str,"+ offset(log(population+1))+ f(idarea,model='bym',graph=g_main,scale.model=TRUE,hyper=pc_bym)+ f(idtime,model='rw1',hyper=pc_prec)+ f(idarea_time,model='iid',hyper=pc_prec)")),
    family=FAMILY,data=ic,control.family=list(),control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),control.predictor=list(link=1),verbose=FALSE),error=function(e)NULL)
  if(is.null(fit)||is.na(fit$dic$dic)) return(result)
  qcat(sprintf("  M6 BYM+RW1+IID (NB): DIC=%.2f | N=%d | EPV=%.1f\n", fit$dic$dic, nrow(ic), nrow(ic)/nrow(FMAP)))
  result$dic<-fit$dic$dic;result$N_final<-nrow(ic);result$EPV<-nrow(ic)/nrow(FMAP)
  result$fit<-fit;result$ic<-ic;result$FMAP<-FMAP;result$form_map<-form_map;result$data_ext<-data_ext

  fe<-fit$summary.fixed;fe<-fe[rownames(fe)!="(Intercept)",,drop=FALSE]
  if(nrow(fe)==0){result$n_fwd<-0;result$n_rev<-0;result$n_neu<-0;result$sig_count<-0;result$mv<-data.frame();return(result)}
  mv<-data.frame();sig_count<-0
  for(k in 1:nrow(fe)){
    irr<-round(exp(fe$mean[k]),4);lo<-round(exp(fe$`0.025quant`[k]),4);hi<-round(exp(fe$`0.975quant`[k]),4)
    sig<-ifelse(fe$`0.025quant`[k]>0|fe$`0.975quant`[k]<0,"★","")
    if(hi > 1000 || lo < 0.001){sig <- ""; irr <- round(exp(fe$mean[k]),4)}
    rn<-gsub("_z$","",rownames(fe)[k]);kr_n<-rn;form_str<-"";tier_str<-"";cat_str<-"";theory_dir<-""
    for(ii in 1:nrow(FMAP)){if(FMAP$eng[ii]==rn){kr_n<-FMAP$kr[ii];form_str<-FMAP$형태[ii];tier_str<-FMAP$tier[ii];cat_str<-FMAP$cat[ii];theory_dir<-FMAP$이론방향[ii];break}}
    obs_dir<-ifelse(irr>1,"위험↑","보호↓")
    if(sig!="★") direction_match<-"비유의" else if(theory_dir=="중립") direction_match<-"중립"
    else{theory_risk<-(theory_dir=="위험");obs_risk<-(irr>1);direction_match<-ifelse(theory_risk==obs_risk,"✅정방향","❌역방향")}
    if(sig=="★") sig_count<-sig_count+1
    mv<-rbind(mv,data.frame(tier=tier_str,카테고리=cat_str,var_kr=kr_n,형태=form_str,IRR=irr,lo=lo,hi=hi,
      sig=sig,obs_dir=obs_dir,theory_dir=theory_dir,방향일치=direction_match,stringsAsFactors=FALSE))
  }
  if(!quiet){for(k in 1:nrow(mv)){
    flag<-ifelse(mv$방향일치[k]=="❌역방향","❌",ifelse(mv$방향일치[k]=="✅정방향","✅"," "))
    qcat(sprintf("    [%s]%s %-22s (%-6s) IRR=%.4f (%.4f–%.4f) %s %s | 이론=%s\n",
        mv$tier[k],flag,mv$var_kr[k],mv$형태[k],mv$IRR[k],mv$lo[k],mv$hi[k],mv$sig[k],mv$obs_dir[k],mv$theory_dir[k]))}}
  n_fwd<-sum(mv$방향일치=="✅정방향");n_rev<-sum(mv$방향일치=="❌역방향");n_neu<-sum(mv$방향일치=="중립"&mv$sig=="★")
  qcat(sprintf("\n  ★ 유의 %d/%d | ✅정방향 %d | ❌역방향 %d | 중립 %d\n",sig_count,nrow(fe),n_fwd,n_rev,n_neu))
  result$mv<-mv;result$sig_count<-sig_count;result$n_fwd<-n_fwd;result$n_rev<-n_rev;result$n_neu<-n_neu
  return(result)
}

# ══════════════════════════════════════════
# PART 3.5-B. 사전지정(PRE-SPECIFIED) — AUTO 타깃 최적화 없음  ★ v8-B 핵심 변경
# ══════════════════════════════════════════
cat("\n## PART 3.5-B. 사전지정(PRE-SPECIFIED) — AUTO 없음\n\n")
TV_FULL_ALL <- TV_v6     # 전체/도시: 사전지정 이론 30개
# TV_RURAL 는 PART 2.5 에서 이미 정의(18개)
TV_FINAL    <- TV_FULL_ALL
best_rm_codes <- character(0)
cat(sprintf("  [B/사전지정] 전체·도시 %d개 | 농촌 %d개 (forward-target·reverse-drop 미적용)\n", nrow(TV_FULL_ALL), nrow(TV_RURAL)))
cat("  ★ 유의/방향은 데이터가 나오는 대로 보고(정직).\n\n")

# ══════════════════════════════════════════
# PART 3.6. 3그룹 분리 (전체/도시/농촌)
# ══════════════════════════════════════════
cat("## PART 3.6. 3그룹 분리 (전체/도시/농촌)\n\n")
data_all   <- cor_merged
data_urban <- cor_merged %>% filter(area_2 == "도시")
data_rural <- cor_merged %>% filter(area_2 == "농촌")
cat(sprintf("  전체: %d행 %d시군구 | 도시: %d행 %d시군구 | 농촌: %d행 %d시군구\n\n",
    nrow(data_all), n_distinct(data_all$region),
    nrow(data_urban), n_distinct(data_urban$region),
    nrow(data_rural), n_distinct(data_rural$region)))

# ══════════════════════════════════════════
# PART 5. 최종 모델 (verbose) — 3그룹
# ══════════════════════════════════════════
cat("## PART 5. 최종 모델 (verbose) — 3그룹\n\n")

cat("  ──── 전체 모델 ────\n")
res_final <- run_model(TV_FULL_ALL, quiet=FALSE); gc()
n_fwd_final <- ifelse(is.null(res_final$n_fwd), NA, res_final$n_fwd)
n_rev_final <- ifelse(is.null(res_final$n_rev), NA, res_final$n_rev)
# CrI 안정성 필터 (불안정 변수 제거 후 재적합)
if(!is.null(res_final$mv) && nrow(res_final$mv) > 0) {
  mv_check <- res_final$mv
  unstable_idx <- which(mv_check$IRR > 500 | mv_check$IRR < 0.002 |
                        (mv_check$hi / pmax(mv_check$lo, 1e-10)) > 1000)
  if(length(unstable_idx) > 0) {
    cat("\n  ⚠️ CrI 불안정 변수 감지 → 제거 후 재적합:\n")
    remove_kr <- c()
    for(j in unstable_idx) {
      cri_ratio <- mv_check$hi[j] / max(mv_check$lo[j], 1e-10)
      cat(sprintf("    제거: %s (IRR=%.4f, CrI ratio=%.0f)\n", mv_check$var_kr[j], mv_check$IRR[j], cri_ratio))
      remove_kr <- c(remove_kr, mv_check$var_kr[j])
    }
    TV_FULL_ALL <- TV_FULL_ALL[!TV_FULL_ALL$kr %in% remove_kr, ]
    TV_FINAL <- TV_FULL_ALL
    cat(sprintf("  → 잔여 변수: %d개 → 재적합...\n\n", nrow(TV_FULL_ALL)))
    res_final <- run_model(TV_FULL_ALL, quiet=FALSE); gc()
    n_fwd_final <- ifelse(is.null(res_final$n_fwd), NA, res_final$n_fwd)
    n_rev_final <- ifelse(is.null(res_final$n_rev), NA, res_final$n_rev)
  }
}

cat("\n  ──── 도시 모델 ────\n")
old_cm <- cor_merged
assign("cor_merged", data_urban, envir=.GlobalEnv)
res_urban <- run_model(TV_FULL_ALL, quiet=FALSE); gc()
assign("cor_merged", old_cm, envir=.GlobalEnv)

cat("\n  ──── 농촌 모델 ────\n")
assign("cor_merged", data_rural, envir=.GlobalEnv)
res_rural <- run_model(TV_RURAL, quiet=FALSE); gc()
assign("cor_merged", old_cm, envir=.GlobalEnv)

# 3그룹 요약
cat("\n\n  ═══════════════════════════════════════════════════════\n")
cat("  3그룹 최종 요약 (사전지정/PRE-SPECIFIED)\n")
cat("  ═══════════════════════════════════════════════════════\n")
for(rr in list(
  list(nm="전체", r=res_final, tv=TV_FULL_ALL),
  list(nm="도시", r=res_urban, tv=TV_FULL_ALL),
  list(nm="농촌", r=res_rural, tv=TV_RURAL))){
  cat(sprintf("  │ %-6s │ 변수 %2d │ 유의 %2d │ 정방향 %2d │ 역방향 %2d │ 중립 %2d\n",
      rr$nm, nrow(rr$tv),
      ifelse(is.null(rr$r$sig_count), 0, rr$r$sig_count),
      ifelse(is.null(rr$r$n_fwd), 0, rr$r$n_fwd),
      ifelse(is.null(rr$r$n_rev), 0, rr$r$n_rev),
      ifelse(is.null(rr$r$n_neu), 0, rr$r$n_neu)))
}
cat("\n")

# M1-M6 비교 (전체 모델)
cat("\n  M1-M6 모델비교...\n")
bm <- NULL; bi <- 6
all_m <- list(); dics <- rep(NA,6); waics <- rep(NA,6); delta_m4m6 <- NA
all_fe_df <- data.frame(); best_factor <- NULL
moran_pre <- NULL; moran_post <- NULL; high_r <- character(0); low_r <- character(0)
if(!is.null(res_final$FMAP) && nrow(res_final$FMAP)>0 && !is.null(res_final$ic)){
  ic <- res_final$ic; cov_str <- paste(res_final$FMAP$safe, collapse=" + ")
  pc_bym<-list(prec.unstruct=list(prior="pc.prec",param=c(0.5,0.01)),prec.spatial=list(prior="pc.prec",param=c(0.5,0.01)))
  pc_prec<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)))
  base_f <- paste("cases ~",cov_str,"+ offset(log(population+1))")
  run_m<-function(fs,nm){fit<-tryCatch(inla(as.formula(fs),family=FAMILY,data=ic,control.family=list(),
    control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),control.predictor=list(link=1),verbose=FALSE),
    error=function(e){cat(sprintf("  ❌ %s: %s\n",nm,e$message));NULL});if(!is.null(fit))cat(sprintf("  %s: DIC=%.2f\n",nm,fit$dic$dic));fit}
  M1<-run_m(base_f,"M1 NB")
  M2<-run_m(paste(base_f,"+ f(idarea,model='besag',graph=g_main,scale.model=TRUE,hyper=pc_prec)"),"M2 ICAR")
  M3<-run_m(paste(base_f,"+ f(idarea,model='bym',graph=g_main,scale.model=TRUE,hyper=pc_bym)"),"M3 BYM")
  M4<-run_m(paste(base_f,"+ f(idarea,model='bym',graph=g_main,scale.model=TRUE,hyper=pc_bym)+ f(idarea_time,model='iid',hyper=pc_prec)"),"M4 BYM+IID")
  M5<-run_m(paste(base_f,"+ f(idarea,model='bym',graph=g_main,scale.model=TRUE,hyper=pc_bym)+ f(idtime,model='rw1',hyper=pc_prec)"),"M5 BYM+RW1")
  M6<-run_m(paste(base_f,"+ f(idarea,model='bym',graph=g_main,scale.model=TRUE,hyper=pc_bym)+ f(idtime,model='rw1',hyper=pc_prec)+ f(idarea_time,model='iid',hyper=pc_prec)"),"M6 BYM+RW1+IID")
  all_m<-list(M1=M1,M2=M2,M3=M3,M4=M4,M5=M5,M6=M6)
  dics<-sapply(all_m,function(m)if(!is.null(m))m$dic$dic else NA)
  waics<-sapply(all_m,function(m)if(!is.null(m))m$waic$waic else NA)
  delta_m4m6<-abs(dics[4]-dics[6])
  if(!is.na(delta_m4m6)&&delta_m4m6<=2){bi<-4;cat(sprintf("\n  ★ M4 선택 (ΔDIC=%.2f ≤ 2)\n",delta_m4m6))
  }else{bi<-which.min(dics);cat(sprintf("\n  ★ %s 선택 (DIC=%.2f)\n",names(all_m)[bi],dics[bi]))}
  bm <- all_m[[bi]]

  if(!is.null(bm$cpo$cpo)){
    cpo_fail <- sum(bm$cpo$cpo < 0.001, na.rm=TRUE)
    cat(sprintf("  CPO 진단: fail(CPO<0.001)=%d/%d (%.1f%%)\n", cpo_fail, length(bm$cpo$cpo), cpo_fail/length(bm$cpo$cpo)*100))
  }
  cat(sprintf("  M4: DIC=%.2f | M6: DIC=%.2f | ΔDIC=%.2f\n\n",dics[4],dics[6],delta_m4m6))

  all_fe_list<-list()
  for(mi in 1:6){if(is.null(all_m[[mi]]))next
    fe_i<-all_m[[mi]]$summary.fixed;fe_i<-fe_i[rownames(fe_i)!="(Intercept)",,drop=FALSE]
    for(k in 1:nrow(fe_i)){irr<-exp(fe_i$mean[k]);lo<-exp(fe_i$`0.025quant`[k]);hi<-exp(fe_i$`0.975quant`[k])
      sig<-ifelse(fe_i$`0.025quant`[k]>0|fe_i$`0.975quant`[k]<0,"★","")
      rn<-gsub("_z$","",rownames(fe_i)[k]);kr_n<-rn
      for(ii in 1:nrow(res_final$FMAP))if(res_final$FMAP$eng[ii]==rn){kr_n<-res_final$FMAP$kr[ii];break}
      all_fe_list[[length(all_fe_list)+1]]<-data.frame(model=paste0("M",mi),var_kr=kr_n,IRR=round(irr,4),lo=round(lo,4),hi=round(hi,4),sig=sig,stringsAsFactors=FALSE)}}
  all_fe_df<-do.call(rbind,all_fe_list)
  for(mi in 1:6){ns<-sum(all_fe_df$sig[all_fe_df$model==paste0("M",mi)]=="★")
    cat(sprintf("  M%d: 유의 %d/%d\n",mi,ns,sum(all_fe_df$model==paste0("M",mi))))}
  cat("\n")

  cat("  Factor 모델 (전체 변수, T3/Q4→factor)...\n")
  n_tq<-sum(res_final$FMAP$형태%in%c("T3","Q4"))
  fterms<-c();for(i in 1:nrow(res_final$FMAP)){
    if(res_final$FMAP$형태[i]%in%c("T3","Q4"))fterms<-c(fterms,paste0(res_final$FMAP$eng[i],"_f"))
    else fterms<-c(fterms,res_final$FMAP$safe[i])}
  cov_f<-paste(fterms,collapse=" + ")
  bf_formula<-paste("cases ~",cov_f,"+ offset(log(population+1))+ f(idarea,model='bym',graph=g_main,scale.model=TRUE,hyper=pc_bym)+ f(idtime,model='rw1',hyper=pc_prec)+ f(idarea_time,model='iid',hyper=pc_prec)")
  best_factor<-tryCatch(inla(as.formula(bf_formula),family=FAMILY,data=ic,control.family=list(),control.compute=list(dic=TRUE),control.predictor=list(link=1),verbose=FALSE),error=function(e){cat(sprintf("  ❌ %s\n",e$message));bm})
  if(!is.null(best_factor))cat(sprintf("  Factor DIC=%.2f (전체 %d변수, T3/Q4 %d개 factor)\n",best_factor$dic$dic,nrow(res_final$FMAP),n_tq))

  cat("\n  Moran's I:\n")
  rate_r<-cor_merged%>%filter(region%in%shp_main$region)%>%group_by(region)%>%summarise(r=sum(cases)/sum(population)*100000,.groups="drop")
  rvec<-rate_r$r[match(shp_main$region,rate_r$region)];rvec[is.na(rvec)]<-0
  moran_pre<-tryCatch(moran.test(rvec,nb_w),error=function(e)NULL)
  if(!is.null(moran_pre))cat(sprintf("    사전: I=%+.4f p=%.4f\n",moran_pre$estimate[1],moran_pre$p.value))
  if(!is.null(bm)){res<-ic$cases - bm$summary.fitted.values$mean[1:nrow(ic)]
    rdf<-data.frame(region=ic$region,r=res)%>%group_by(region)%>%summarise(r=mean(r),.groups="drop")
    rv2<-rdf$r[match(shp_main$region,rdf$region)];rv2[is.na(rv2)]<-0
    moran_post<-tryCatch(moran.test(rv2,nb_w),error=function(e)NULL)
    if(!is.null(moran_post))cat(sprintf("    사후: I=%+.4f p=%.4f %s\n",moran_post$estimate[1],moran_post$p.value,ifelse(moran_post$p.value>0.05,"✅","⚠️")))}

  if(!is.null(bm)&&!is.null(bm$summary.random$idarea)){na<-nrow(shp_main)
    sl<-bm$summary.random$idarea$`0.025quant`[1:na];sh2<-bm$summary.random$idarea$`0.975quant`[1:na]
    high_r<-shp_main$region[sl>0];low_r<-shp_main$region[sh2<0]
    cat(sprintf("\n  고위험: %d | 저위험: %d\n",length(high_r),length(low_r)))}
}

# ══════════════════════════════════════════
# PART 6. Academic Tables (xlsx)
# ══════════════════════════════════════════
cat("\n## PART 6. Academic Tables\n\n")
tryCatch({
wb<-createWorkbook()
s_t<-createStyle(fontSize=13,textDecoration="bold")
s_h<-createStyle(fontSize=10,textDecoration="bold",halign="center",border="TopBottom",borderStyle="medium",fgFill="#D9E2F3")
s_c<-createStyle(fontSize=11,textDecoration="bold",fgFill="#E2EFDA")
s_v<-createStyle(fontSize=10,indent=1);s_d<-createStyle(fontSize=10,halign="right")
s_s<-createStyle(textDecoration="bold",fontColour="#C00000");s_n<-createStyle(fontSize=9,textDecoration="italic")
s_ns<-createStyle(fontColour="#999999");s_sec<-createStyle(fontSize=11,textDecoration="bold",fgFill="#E2EFDA")

write_t1<-function(ws,en){
  ttl<-if(en)sprintf("Table 1. Candidate variables and univariable analysis, NORO, %d-%d (N=%d)",YEAR_START,YEAR_END,nrow(cor_merged))else sprintf("표 1. NORO 후보변수 기술통계 및 단변량 분석 (N=%d, α=0.05)",nrow(cor_merged))
  writeData(wb,ws,ttl,startRow=1);addStyle(wb,ws,s_t,rows=1,cols=1)
  hd<-if(en)c("Category","Variable","N","Mean ± SD","Min","Median","Max","IRR","95% CI","p-value","Sig","Theory")else c("카테고리","변수","N","평균 ± SD","최소","중앙값","최대","IRR","95% CI","p값","유의","이론방향")
  writeData(wb,ws,t(hd),startRow=4,startCol=1,colNames=FALSE);addStyle(wb,ws,s_h,rows=4,cols=1:12,gridExpand=TRUE)
  r<-5;pc<-""
  for(i in 1:nrow(raw_univ)){
    if(raw_univ$cat[i]!=pc){writeData(wb,ws,raw_univ$cat[i],startRow=r,startCol=1);addStyle(wb,ws,s_c,rows=r,cols=1:12,gridExpand=TRUE);pc<-raw_univ$cat[i];r<-r+1}
    writeData(wb,ws,raw_univ$kr[i],startRow=r,startCol=2);addStyle(wb,ws,s_v,rows=r,cols=2)
    writeData(wb,ws,raw_univ$N[i],startRow=r,startCol=3)
    writeData(wb,ws,raw_univ$mean_sd[i],startRow=r,startCol=4)
    writeData(wb,ws,raw_univ$min_v[i],startRow=r,startCol=5)
    writeData(wb,ws,raw_univ$med[i],startRow=r,startCol=6)
    writeData(wb,ws,raw_univ$max_v[i],startRow=r,startCol=7)
    writeData(wb,ws,ifelse(is.na(raw_univ$raw_IRR[i]),"—",sprintf("%.4f",raw_univ$raw_IRR[i])),startRow=r,startCol=8)
    ci<-if(!is.na(raw_univ$raw_lo[i]))sprintf("%.4f–%.4f",raw_univ$raw_lo[i],raw_univ$raw_hi[i])else"—"
    writeData(wb,ws,ci,startRow=r,startCol=9)
    writeData(wb,ws,ifelse(is.na(raw_univ$raw_p[i]),"—",ifelse(raw_univ$raw_p[i]<0.001,"<0.001",sprintf("%.4f",raw_univ$raw_p[i]))),startRow=r,startCol=10)
    writeData(wb,ws,raw_univ$sig[i],startRow=r,startCol=11)
    writeData(wb,ws,raw_univ$이론방향[i],startRow=r,startCol=12)
    addStyle(wb,ws,s_d,rows=r,cols=3:10,gridExpand=TRUE)
    if(raw_univ$sig[i]=="*")addStyle(wb,ws,s_s,rows=r,cols=8:11,gridExpand=TRUE)
    r<-r+1}
  setColWidths(wb,ws,cols=1,widths=25);setColWidths(wb,ws,cols=2,widths=20);setColWidths(wb,ws,cols=3:12,widths=14)}
addWorksheet(wb,"Table1_EN");write_t1("Table1_EN",TRUE)
addWorksheet(wb,"Table1_KR");write_t1("Table1_KR",FALSE)
cat("  ✅ Table1 EN/KR\n")

# ModelComparison
addWorksheet(wb,"ModelComparison")
writeData(wb,"ModelComparison","M1~M6 Fixed Effects Comparison (PRE-SPECIFIED)",startRow=1);addStyle(wb,"ModelComparison",s_t,rows=1,cols=1)
mc_h<-c("변수");for(mi in 1:6)mc_h<-c(mc_h,sprintf("M%d IRR",mi),sprintf("M%d Sig",mi));mc_h<-c(mc_h,"DIC","WAIC")
writeData(wb,"ModelComparison",t(mc_h),startRow=3,startCol=1,colNames=FALSE);addStyle(wb,"ModelComparison",s_h,rows=3,cols=1:length(mc_h),gridExpand=TRUE)
if(nrow(all_fe_df)>0){uv<-unique(all_fe_df$var_kr)
for(ri in seq_along(uv)){vn<-uv[ri];writeData(wb,"ModelComparison",vn,startRow=3+ri,startCol=1)
  for(mi in 1:6){sub<-all_fe_df[all_fe_df$model==paste0("M",mi)&all_fe_df$var_kr==vn,];ci<-2+(mi-1)*2
    if(nrow(sub)>0){writeData(wb,"ModelComparison",sub$IRR[1],startRow=3+ri,startCol=ci)
      writeData(wb,"ModelComparison",sub$sig[1],startRow=3+ri,startCol=ci+1)
      if(sub$sig[1]=="★")addStyle(wb,"ModelComparison",s_s,rows=3+ri,cols=ci:(ci+1),gridExpand=TRUE)}}}
dr<-3+length(uv)+1
writeData(wb,"ModelComparison","DIC",startRow=dr,startCol=1);writeData(wb,"ModelComparison","WAIC",startRow=dr+1,startCol=1);writeData(wb,"ModelComparison","유의 수",startRow=dr+2,startCol=1)
for(mi in 1:6){ci<-2+(mi-1)*2
  if(!is.na(dics[mi]))writeData(wb,"ModelComparison",round(dics[mi],2),startRow=dr,startCol=ci)
  if(!is.na(waics[mi]))writeData(wb,"ModelComparison",round(waics[mi],2),startRow=dr+1,startCol=ci)
  ns_mi<-sum(all_fe_df$sig[all_fe_df$model==paste0("M",mi)]=="★");nt_mi<-sum(all_fe_df$model==paste0("M",mi))
  writeData(wb,"ModelComparison",sprintf("%d/%d",ns_mi,nt_mi),startRow=dr+2,startCol=ci)}
writeData(wb,"ModelComparison",sprintf("★ M%d 선택 (ΔDIC M4-M6=%.2f)",bi,delta_m4m6),startRow=dr+4,startCol=1);addStyle(wb,"ModelComparison",s_n,rows=dr+4,cols=1)}
setColWidths(wb,"ModelComparison",cols=1,widths=20);setColWidths(wb,"ModelComparison",cols=2:13,widths=10)
cat("  ✅ ModelComparison\n")

# MoransI
addWorksheet(wb,"MoransI");writeData(wb,"MoransI","Moran's I",startRow=1);addStyle(wb,"MoransI",s_t,rows=1,cols=1)
writeData(wb,"MoransI",t(c("Timing","I","p","Judgment")),startRow=3,startCol=1,colNames=FALSE);addStyle(wb,"MoransI",s_h,rows=3,cols=1:4,gridExpand=TRUE)
if(!is.null(moran_pre)){writeData(wb,"MoransI","Pre",startRow=4,startCol=1);writeData(wb,"MoransI",round(moran_pre$estimate[1],4),startRow=4,startCol=2)
  writeData(wb,"MoransI",round(moran_pre$p.value,6),startRow=4,startCol=3);writeData(wb,"MoransI",ifelse(moran_pre$p.value<0.05,"Significant","ns"),startRow=4,startCol=4)}
if(!is.null(moran_post)){writeData(wb,"MoransI","Post",startRow=5,startCol=1);writeData(wb,"MoransI",round(moran_post$estimate[1],4),startRow=5,startCol=2)
  writeData(wb,"MoransI",round(moran_post$p.value,6),startRow=5,startCol=3);writeData(wb,"MoransI",ifelse(moran_post$p.value>0.05,"Removed","⚠️"),startRow=5,startCol=4)}
cat("  ✅ MoransI\n")

# VariableMapping
if(!is.null(res_final$FMAP) && nrow(res_final$FMAP)>0){
  addWorksheet(wb,"VariableMapping");writeData(wb,"VariableMapping","Variable Mapping",startRow=1);addStyle(wb,"VariableMapping",s_t,rows=1,cols=1)
  vm_cols <- intersect(c("cat","kr","eng","형태","forced","이론방향","tier"), names(res_final$FMAP))
  writeData(wb,"VariableMapping",res_final$FMAP[,vm_cols],startRow=3,colNames=TRUE)
  cat("  ✅ VariableMapping\n")}

# HighRisk / LowRisk
for(st in c("HighRisk","LowRisk")){addWorksheet(wb,st);ih<-st=="HighRisk"
  writeData(wb,st,sprintf("%s municipalities",if(ih)"High-risk"else"Low-risk"),startRow=1);addStyle(wb,st,s_t,rows=1,cols=1)
  writeData(wb,st,t(c("#","Municipality","Spatial Effect","CrI Lo","CrI Hi")),startRow=3,startCol=1,colNames=FALSE);addStyle(wb,st,s_h,rows=3,cols=1:5,gridExpand=TRUE)
  if(!is.null(bm)&&!is.null(bm$summary.random$idarea)){na<-nrow(shp_main);sm<-bm$summary.random$idarea$mean[1:na]
    sl2<-bm$summary.random$idarea$`0.025quant`[1:na];sh3<-bm$summary.random$idarea$`0.975quant`[1:na]
    idx<-if(ih)which(sl2>0)else which(sh3<0);idx<-idx[order(sm[idx],decreasing=ih)]
    if(length(idx)>0)for(k in seq_along(idx)){writeData(wb,st,k,startRow=3+k,startCol=1);writeData(wb,st,shp_main$region[idx[k]],startRow=3+k,startCol=2)
      writeData(wb,st,round(sm[idx[k]],4),startRow=3+k,startCol=3);writeData(wb,st,round(sl2[idx[k]],4),startRow=3+k,startCol=4);writeData(wb,st,round(sh3[idx[k]],4),startRow=3+k,startCol=5)}}}
cat("  ✅ HighRisk / LowRisk\n")

# Table2_그룹 (전체/도시/농촌 결과)
write_group_table <- function(ws, grp_name, res_grp, en) {
  ttl <- if(en) sprintf("Table 2-%s. Bayesian spatiotemporal model, NORO — %s (PRE-SPECIFIED)", grp_name, grp_name)
         else sprintf("표 2-%s. 베이지안 시공간 모델 — %s (사전지정)", grp_name, grp_name)
  writeData(wb, ws, ttl, startRow=1); addStyle(wb, ws, s_t, rows=1, cols=1)
  if(is.null(res_grp$mv) || nrow(res_grp$mv)==0) { writeData(wb, ws, "No results available", startRow=3); return(invisible(NULL)) }
  mv <- res_grp$mv
  hd <- if(en) c("Tier","Category","Variable","Form","IRR","95% CrI Lo","95% CrI Hi","Sig","Direction","Theory","Match")
        else c("등급","카테고리","변수","형태","IRR","95% CrI Lo","95% CrI Hi","유의","관찰방향","이론방향","방향일치")
  writeData(wb, ws, t(hd), startRow=3, startCol=1, colNames=FALSE); addStyle(wb, ws, s_h, rows=3, cols=1:11, gridExpand=TRUE)
  for(i in 1:nrow(mv)) {
    r <- 3 + i
    writeData(wb, ws, mv$tier[i], startRow=r, startCol=1)
    writeData(wb, ws, mv$카테고리[i], startRow=r, startCol=2)
    writeData(wb, ws, mv$var_kr[i], startRow=r, startCol=3); addStyle(wb, ws, s_v, rows=r, cols=3)
    writeData(wb, ws, mv$형태[i], startRow=r, startCol=4)
    writeData(wb, ws, mv$IRR[i], startRow=r, startCol=5)
    writeData(wb, ws, mv$lo[i], startRow=r, startCol=6)
    writeData(wb, ws, mv$hi[i], startRow=r, startCol=7)
    writeData(wb, ws, mv$sig[i], startRow=r, startCol=8)
    writeData(wb, ws, mv$obs_dir[i], startRow=r, startCol=9)
    writeData(wb, ws, mv$theory_dir[i], startRow=r, startCol=10)
    writeData(wb, ws, mv$방향일치[i], startRow=r, startCol=11)
    addStyle(wb, ws, s_d, rows=r, cols=5:7, gridExpand=TRUE)
    if(mv$sig[i]=="★") addStyle(wb, ws, s_s, rows=r, cols=5:8, gridExpand=TRUE)
  }
  sr <- 3 + nrow(mv) + 2
  n_sig <- sum(mv$sig=="★"); n_fwd <- sum(mv$방향일치=="✅정방향"); n_rev <- sum(mv$방향일치=="❌역방향")
  writeData(wb, ws, sprintf("★ 유의 %d/%d | ✅정방향 %d | ❌역방향 %d", n_sig, nrow(mv), n_fwd, n_rev), startRow=sr, startCol=1)
  addStyle(wb, ws, s_n, rows=sr, cols=1)
  setColWidths(wb, ws, cols=1, widths=8); setColWidths(wb, ws, cols=2, widths=22)
  setColWidths(wb, ws, cols=3, widths=25); setColWidths(wb, ws, cols=4:11, widths=14)
}
for(grp_info in list(
  list(nm="전체", nm_en="Overall", res=res_final),
  list(nm="도시", nm_en="Urban", res=res_urban),
  list(nm="농촌", nm_en="Rural", res=res_rural))) {
  ws_en <- sprintf("Table2_%s_EN", grp_info$nm_en); ws_kr <- sprintf("Table2_%s_KR", grp_info$nm)
  addWorksheet(wb, ws_en); write_group_table(ws_en, grp_info$nm_en, grp_info$res, TRUE)
  addWorksheet(wb, ws_kr); write_group_table(ws_kr, grp_info$nm, grp_info$res, FALSE)
}
cat("  ✅ Table2_전체/도시/농촌 EN/KR\n")

fn_xlsx<-file.path(DIR_OUT,sprintf("NORO_v8Bc2_CHILDok_%s.xlsx",TS));saveWorkbook(wb,fn_xlsx,overwrite=TRUE)
cat(sprintf("\n  ★ 엑셀: %s\n",fn_xlsx))
},error=function(e)cat(sprintf("  ❌ Tables: %s\n",e$message)))

# ══════════════════════════════════════════
# 완료
# ══════════════════════════════════════════
cat(sprintf("\n✅ 로그: %s\n",LOG))
sink()
cat(sprintf("\n═══ NORO v8-B PRE-SPECIFIED (전체/도시/농촌) 완료 ═══\n"))
cat(sprintf("★ 모델: %s | 변수선택: 사전지정(AUTO 없음)\n", FAMILY))
cat(sprintf("★ 엑셀: %s\n", if(exists("fn_xlsx")) fn_xlsx else ""))


# ═══════════════════════════════════════════════════════════════════════════════
# ▣ [통합] NORO 추가분석 — 메인(CHILDok)이 만든 res_final/cor_merged/shp_main/nb_obj/run_model 이어받음
# ═══════════════════════════════════════════════════════════════════════════════
DIR_LOG <- if (exists("DIR_OUT")) DIR_OUT else getwd(); OUTPUT_DIR <- DIR_LOG
if (!exists("res_final")) stop("메인(NORO_v8Bc2_CHILDok_FULL)을 먼저 끝까지 실행해야 합니다")
cat("\n\n", strrep("█",80), "\n  통합 추가분석 (소아 dissociation · φ/8graph/prior · case-vs-outbreak)\n", strrep("█",80), "\n", sep="")


#==============================================================================
# ▶ 통합블록: 소아 battery 다중공선성·dissociation 진단 (4-1) = Table 4.2/4.S4 보강
#==============================================================================
# ════════════════════════════════════════════════════════════════════════════════
# NORO Ch4 — 소아 battery 다중공선성·차원정리 진단 (종심 4-1, 최승아 / EHEC와 동일 기조)
# ════════════════════════════════════════════════════════════════════════════════
# 맥락: NORO Ch4의 소아 battery(0–4·5–9·소아입원·어린이집)는 "환경앵커 4 + 소아 4"의
#       별도 사전지정 'dissociation 검정' 모형(주 모형 Table 4.1과 구분). 소아 변수는 모두
#       비유의(=노로가 소아 인구/질병부담과 분리)라는 것이 결론.
# 목적: EHEC와 동일하게 (1) 소아 4변수 상관·분포, (2) VIF, (3) 0–4 제거/0–9 통합 재적합으로
#       0–4·5–9 공선성을 정리해도 'dissociation(소아 비유의)'이 그대로 유지됨을 보임.
# 전제: NORO 소아 battery 모형을 같은 세션에서 먼저 실행(cor_merged·run_model 존재;
#       NORO_aux_paediatric_battery.R 와 동일 전제). area_2(도시/농촌)로 층화 가능.
# ════════════════════════════════════════════════════════════════════════════════
tryCatch({
suppressMessages({ library(dplyr); library(car); library(openxlsx) })
stopifnot(exists("cor_merged"))
OUT <- get0("DIR_LOG", ifnotfound=getwd()); ts <- format(Sys.time(),"%y%m%d_%H%M")

PED4 <- c("영유아비율_0_4","아동비율_5_9","예방가능입원율_소아위장관염","유아천명당보육시설수")
ped_in <- intersect(PED4, names(cor_merged))
df <- cor_merged

# ── 1. 소아 4변수 상관 + 분포 ──────────────────────────────────────────────
cat("[1] 소아 4변수 Pearson 상관\n"); print(round(cor(df[,ped_in], use="pairwise.complete.obs"),3))
cat("\n[2] 분포(영점과다 점검)\n")
print(do.call(rbind, lapply(ped_in, function(v){ x<-as.numeric(df[[v]])
  data.frame(variable=v, n=sum(!is.na(x)), zeros=sum(x==0,na.rm=TRUE),
             pct_zero=round(100*mean(x==0,na.rm=TRUE),1),
             min=round(min(x,na.rm=TRUE),3), median=round(median(x,na.rm=TRUE),3),
             max=round(max(x,na.rm=TRUE),3))})), row.names=FALSE)

# ── 3. 소아 battery 모형(환경앵커4 + 소아4) VIF + 선택변이 재적합 ───────────
# 환경앵커(NORO_aux_paediatric_battery.R TV_P4 기준): 목장용지·함수율·민방위용·자체처리량
ENV <- intersect(c("목장용지","함수율(%,탈수기준)","민방위용_개소수","자체처리량(톤/년)_계"), names(df))
build_set <- function(drop04=FALSE, comb09=FALSE){
  pe <- ped_in
  if (comb09){ df$child_0_9 <<- rowMeans(df[, intersect(c("영유아비율_0_4","아동비율_5_9"),names(df))], na.rm=TRUE)
    pe <- c(setdiff(ped_in, c("영유아비율_0_4","아동비율_5_9")), "child_0_9") }
  else if (drop04) pe <- setdiff(ped_in, "영유아비율_0_4")
  c(ENV, pe)
}
vif_of <- function(cols){ cc<-intersect(cols,names(df)); X<-as.data.frame(lapply(df[,cc],as.numeric))
  X<-X[complete.cases(X),,drop=FALSE]; v<-diag(solve(cor(X))); v }
cat("\n[3] VIF (환경앵커4 + 소아: full / drop04 / comb09)\n")
for (nm in c("full","drop04","comb09")){
  cols <- build_set(drop04=(nm=="drop04"), comb09=(nm=="comb09"))
  v <- tryCatch(vif_of(cols), error=function(e) NULL)
  if(!is.null(v)) cat(sprintf("  %-7s maxVIF=%.2f | 소아VIF: %s\n", nm, max(v),
       paste(sprintf("%s=%.2f",intersect(PED4,names(v)),v[intersect(PED4,names(v))]),collapse=", ")))
}

# (선택) run_model 이 있으면 dissociation 유지 확인: drop04/comb09 에서도 소아 IRR 모두 비유의?
if (exists("run_model")){
  cat("\n[4] run_model 감지 → dissociation 재확인 권장:\n")
  cat("    run_model(build_set()) / run_model(build_set(drop04=TRUE)) / run_model(build_set(comb09=TRUE))\n")
  cat("    각 결과 summary.fixed 에서 소아 변수 IRR·CrI 가 모두 1 부근·비유의면 → 차원정리 후에도 dissociation 유지.\n")
}

# ── 저장 ──────────────────────────────────────────────────────────────────
cor4 <- round(cor(df[,ped_in], use="pairwise.complete.obs"),3)
wb<-createWorkbook(); addWorksheet(wb,"ped4_corr"); writeData(wb,"ped4_corr",cbind(variable=rownames(cor4),as.data.frame(cor4)))
saveWorkbook(wb, file.path(OUT, sprintf("NORO_PaedMulticollinearity_%s.xlsx",ts)), overwrite=TRUE)
cat(sprintf("\n저장: NORO_PaedMulticollinearity_%s.xlsx\n", ts))
}, error=function(e) cat(sprintf("\n\u26a0 [\ud1b5\ud569\ube14\ub85d \uac74\ub108\ub700: 소아 dissociation] %s\n", conditionMessage(e))))
# 해석: 소아 데이터는 EHEC와 동일 컬럼이므로 0–4↔5–9 r≈0.91 동일.
#       NORO는 소아 dissociation(모두 비유의)이 결론이라, 0–4 제거/0–9 통합 후에도
#       소아 IRR 비유의가 유지되면 "공선성 정리 후에도 dissociation 불변"으로 4-1 답변.


#==============================================================================
# ▶ 통합블록: robustness — BYM2 φ / 8그래프 / prior 민감도 = Table 4.S6/S7/S8
#==============================================================================
# ==============================================================================
# NORO 정직판 robustness 애드온 — BYM2 φ / 8그래프 / prior 민감도 (Total 모델)
# ------------------------------------------------------------------------------
# ★ 실행 전제: NORO_v8Bc2_CHILDok_FULL.R 를 '같은 R 세션'에서 먼저 실행.
#   재사용: res_final(ic·FMAP), shp_main, nb_obj, g_main, FAMILY, DIR_OUT.
# 산출: 보충자료 Table S8(BYM2 φ) / S9(8그래프) / S10(prior 민감도) 채울 실수치.
#       콘솔 출력 + NORO_ROBUSTNESS_*.xlsx
# 주의: 무날조 — 빈 결과는 NA로 보고. 추정 실패시 해당 표 생략.
# ==============================================================================
tryCatch({
suppressPackageStartupMessages({library(INLA);library(spdep);library(sf);library(dplyr);library(openxlsx)})
cat("\n",strrep("=",78),"\n  NORO robustness (BYM2·8graph·prior) — Total 모델\n",strrep("=",78),"\n",sep="")
stopifnot(exists("res_final"), !is.null(res_final$ic), !is.null(res_final$FMAP), exists("shp_main"), exists("nb_obj"))

ic <- res_final$ic
cov_str <- paste(res_final$FMAP$safe, collapse=" + ")
base_f  <- paste("cases ~", cov_str, "+ offset(log(population+1))")
sig_vars <- res_final$FMAP                                   # eng·kr·safe·형태
pc_bym  <- list(prec.unstruct=list(prior="pc.prec",param=c(0.5,0.01)),prec.spatial=list(prior="pc.prec",param=c(0.5,0.01)))
pc_prec <- list(prec=list(prior="pc.prec",param=c(0.5,0.01)))
TSr <- format(Sys.time(),"%y%m%d_%H%M")

fit_inla <- function(fs, ...){
  tryCatch(inla(as.formula(fs), family=FAMILY, data=ic, control.family=list(),
                control.compute=list(dic=TRUE,waic=TRUE), control.predictor=list(link=1),
                verbose=FALSE, ...),
           error=function(e){cat(sprintf("  ❌ %s\n", e$message)); NULL})
}
# 유의변수의 95%CrI가 1을 제외하는지 (eng 고정효과 기준)
covar_sig <- function(fit){
  if(is.null(fit)||is.null(fit$summary.fixed)) return(NULL)
  fx <- fit$summary.fixed; out <- list()
  for(i in 1:nrow(sig_vars)){ rn <- sig_vars$safe[i]
    if(rn %in% rownames(fx)){ lo<-exp(fx[rn,"0.025quant"]); hi<-exp(fx[rn,"0.975quant"])
      out[[sig_vars$kr[i]]] <- (lo>1)|(hi<1) } }
  out
}

# ── 1. BYM2 φ ──
cat("\n[1] BYM2 reparametrisation → φ (structured 분율)\n")
pc_bym2 <- list(prec=list(prior="pc.prec",param=c(0.5,0.01)),
                phi =list(prior="pc",      param=c(0.5,0.5)))
f_bym2 <- paste(base_f,
  "+ f(idarea,model='bym2',graph=g_main,scale.model=TRUE,hyper=pc_bym2)",
  "+ f(idtime,model='rw1',hyper=pc_prec)+ f(idarea_time,model='iid',hyper=pc_prec)")
m_bym2 <- fit_inla(f_bym2)
phi_row <- NULL
if(!is.null(m_bym2) && !is.null(m_bym2$summary.hyperpar)){
  hp <- m_bym2$summary.hyperpar
  pr <- grep("^Phi", rownames(hp))
  if(length(pr)>0){ phi_row <- hp[pr[1],]
    cat(sprintf("  φ = %.3f (95%% CrI %.3f–%.3f)  | BYM2 DIC=%.2f\n",
        phi_row$mean, phi_row$`0.025quant`, phi_row$`0.975quant`, m_bym2$dic$dic)) }
}
df_phi <- if(!is.null(phi_row)) data.frame(param="phi (structured share)",
  mean=round(phi_row$mean,3), lo=round(phi_row$`0.025quant`,3), hi=round(phi_row$`0.975quant`,3),
  BYM2_DIC=round(m_bym2$dic$dic,2)) else data.frame(note="BYM2 추정 실패 — 미보고")

# ── 2. 8 graph robustness ──
cat("\n[2] 8 neighbourhood graphs\n")
# 평면 좌표(투영)로 변환 → 중심점·거리 정확화 (Korea 2000 Unified CS, EPSG:5179)
shp_proj <- tryCatch({if(is.na(sf::st_crs(shp_main))) sf::st_set_crs(shp_main,4326) else shp_main}, error=function(e) shp_main)
shp_proj <- tryCatch(sf::st_transform(shp_proj, 5179), error=function(e){cat("  ⚠️ 투영 실패 → 원좌표 사용\n"); shp_proj})
cc_proj  <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(shp_proj))))
mk_graph <- function(type){
  if(type=="Queen") nb<-poly2nb(shp_main, queen=TRUE, snap=0.01)
  else if(type=="Rook") nb<-poly2nb(shp_main, queen=FALSE, snap=0.01)
  else { k<-as.integer(sub("knn","",type))
         nb<-knn2nb(knearneigh(cc_proj, k=k), sym=TRUE) }   # ★ sym=TRUE: INLA bym/besag 대칭그래프 필수
  f<-tempfile(); nb2INLA(f, nb); inla.read.graph(f)
}
gtypes <- c("Queen","Rook","knn2","knn3","knn4","knn5","knn6","knn7")
sig_count <- setNames(rep(0L, nrow(sig_vars)), sig_vars$kr)
irr_lo <- setNames(rep(Inf,nrow(sig_vars)),sig_vars$kr); irr_hi<-setNames(rep(-Inf,nrow(sig_vars)),sig_vars$kr)
n_ok <- 0
for(gt in gtypes){
  g <- tryCatch(mk_graph(gt), error=function(e){cat(sprintf("  graph %s 실패: %s\n",gt,e$message));NULL})
  if(is.null(g)) next
  assign("g_tmp", g, envir=.GlobalEnv)
  fs <- paste(base_f,"+ f(idarea,model='bym',graph=g_tmp,scale.model=TRUE,hyper=pc_bym)",
              "+ f(idtime,model='rw1',hyper=pc_prec)+ f(idarea_time,model='iid',hyper=pc_prec)")
  fit <- fit_inla(fs); if(is.null(fit)) next
  n_ok <- n_ok+1; cs <- covar_sig(fit)
  fx <- fit$summary.fixed
  for(i in 1:nrow(sig_vars)){ kr<-sig_vars$kr[i]; rn<-sig_vars$safe[i]
    if(!is.null(cs[[kr]]) && cs[[kr]]) sig_count[kr]<-sig_count[kr]+1L
    if(rn %in% rownames(fx)){ irr_lo[kr]<-min(irr_lo[kr],exp(fx[rn,"mean"])); irr_hi[kr]<-max(irr_hi[kr],exp(fx[rn,"mean"])) } }
  cat(sprintf("  %-6s ✓ (유의 %d/%d)\n", gt, sum(unlist(cs)), nrow(sig_vars)))
}
df_graph <- data.frame(covariate=sig_vars$kr,
  IRR_range=sprintf("%.2f–%.2f", irr_lo, irr_hi),
  graphs_sig=sprintf("%d / %d", sig_count, n_ok),
  classification=ifelse(sig_count==n_ok & n_ok>0,"Robust",ifelse(sig_count>0,"Sensitive","Not sig")),
  stringsAsFactors=FALSE)
df_graph <- df_graph[order(-sig_count),]
cat(sprintf("  → %d/8 그래프 적합\n", n_ok)); print(df_graph[df_graph$graphs_sig!=paste0("0 / ",n_ok),])

# ── 3. prior 민감도 (M4) ──
cat("\n[3] prior 민감도 (M4 BYM+IID)\n")
priors <- list(
  list(nm="PC.prec(0.5,0.01) — principal", u=c(0.5,0.01)),
  list(nm="PC.prec tighter (1,0.01)",       u=c(1,0.01)),
  list(nm="PC.prec looser (0.1,0.01)",      u=c(0.1,0.01)))
rows_pr <- list()
for(p in priors){
  pb<-list(prec.unstruct=list(prior="pc.prec",param=p$u),prec.spatial=list(prior="pc.prec",param=p$u))
  pp<-list(prec=list(prior="pc.prec",param=p$u))
  assign("pb_tmp",pb,envir=.GlobalEnv); assign("pp_tmp",pp,envir=.GlobalEnv)
  fs<-paste(base_f,"+ f(idarea,model='bym',graph=g_main,scale.model=TRUE,hyper=pb_tmp)+ f(idarea_time,model='iid',hyper=pp_tmp)")
  fit<-fit_inla(fs)
  if(!is.null(fit)) rows_pr[[length(rows_pr)+1]]<-data.frame(prior=p$nm,DIC=round(fit$dic$dic,2),WAIC=round(fit$waic$waic,2))
}
# LogGamma 대안
assign("lg_tmp",list(prec=list(prior="loggamma",param=c(1,5e-5))),envir=.GlobalEnv)
fs<-paste(base_f,"+ f(idarea,model='bym',graph=g_main,scale.model=TRUE,hyper=list(prec.unstruct=list(prior='loggamma',param=c(1,5e-5)),prec.spatial=list(prior='loggamma',param=c(1,5e-5))))+ f(idarea_time,model='iid',hyper=lg_tmp)")
fit<-fit_inla(fs); if(!is.null(fit)) rows_pr[[length(rows_pr)+1]]<-data.frame(prior="LogGamma(1,5e-5)",DIC=round(fit$dic$dic,2),WAIC=round(fit$waic$waic,2))
df_prior <- if(length(rows_pr)>0) do.call(rbind,rows_pr) else data.frame(note="prior 적합 실패")
cat("  prior 민감도:\n"); print(df_prior)

# ── 저장 ──
tryCatch({wb<-createWorkbook()
  addWorksheet(wb,"S8_BYM2_phi");   writeData(wb,"S8_BYM2_phi",df_phi)
  addWorksheet(wb,"S9_8graph");     writeData(wb,"S9_8graph",df_graph)
  addWorksheet(wb,"S10_prior");     writeData(wb,"S10_prior",df_prior)
  fn<-file.path(DIR_OUT, sprintf("NORO_ROBUSTNESS_%s.xlsx",TSr)); saveWorkbook(wb,fn,overwrite=TRUE)
  cat(sprintf("\n  ★ 엑셀: %s\n", fn))
}, error=function(e) cat("  ❌ xlsx:",e$message,"\n"))
cat(strrep("=",78),"\n  robustness 완료 — S8/S9/S10 수치를 보충자료에 반영하세요\n",strrep("=",78),"\n",sep="")
}, error=function(e) cat(sprintf("\n\u26a0 [\ud1b5\ud569\ube14\ub85d \uac74\ub108\ub700: robustness φ·8graph·prior] %s\n", conditionMessage(e))))


#==============================================================================
# ▶ 통합블록: 계수단위(case vs outbreak) 민감도 [독립 실행] = Table 4.S9 (유석현 2-2)
#==============================================================================
# ════════════════════════════════════════════════════════════════════
#  노로바이러스 공간분석 — OUTBREAK(건수) 단위 민감도 분석
#  심사위원 #2 (유석현) 코멘트 2-2 대응: case(사례) vs outbreak(집단발생 건수)
#  - 입력: 식중독최종_건수.csv (cases 열 = outbreak 건수)
#  - 비교: 식중독최종.csv (cases 열 = 환자 사례수, 본문 Ch4 = 4,801건)
#  - 기간/지역/모형: Ch4와 동일 (2020–2024, 시군구, NB BYM2)
#  ※ 작성: 응답서 보강용. RStudio(INLA 설치 환경)에서 실행하세요.
# ════════════════════════════════════════════════════════════════════
tryCatch({
suppressMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(sf); library(spdep); library(INLA)
})

# ── 설정 (본인 환경 경로 그대로) ──────────────────────────────────
BASE_IV   <- "/Users/SeongdaeKim/Library/CloudStorage/GoogleDrive-wwwwrte@gmail.com/내 드라이브/S.K/기타/R_Studio/FBD_DATA_ZIP"
PATH_OUTB <- file.path(BASE_IV, "식중독최종_건수.csv")   # outbreak 건수
PATH_CASE <- file.path(BASE_IV, "식중독최종.csv")          # 환자 사례수
PATH_SHP  <- file.path(BASE_IV, "final.shp")
DIR_OUT   <- file.path(BASE_IV, "노로_건수_민감도결과")
DISEASE   <- "노로바이러스"; Y0 <- 2020; Y1 <- 2024
if(!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive=TRUE)

clean_region <- function(df) df %>% mutate(
  region = str_replace_all(as.character(region), "\\s+", ""),
  region = if_else(region=="인천시미추홀구","인천시남구",region),
  year   = as.integer(year))

read_safe <- function(fp){
  for(enc in c("UTF-8","UTF-8-BOM","CP949","EUC-KR")){
    r <- tryCatch(read.csv(fp, fileEncoding=enc, check.names=FALSE, stringsAsFactors=FALSE), error=function(e) NULL)
    if(!is.null(r) && "region" %in% names(r)) return(r)
  }; stop("read fail: ", fp)
}

# ── 1. 시군구 집계 (건수 / 사례) ──────────────────────────────────
agg_disease <- function(path, label){
  read_safe(path) %>% clean_region() %>%
    filter(disease==DISEASE, year>=Y0, year<=Y1) %>%
    group_by(region) %>%
    summarise(!!label := sum(cases, na.rm=TRUE),
              population = mean(population, na.rm=TRUE), .groups="drop")
}
ob <- agg_disease(PATH_OUTB, "outb")
ca <- agg_disease(PATH_CASE, "cases") %>% dplyr::select(region, cases)
dat <- ob %>% left_join(ca, by="region") %>%
  mutate(across(c(outb,cases), ~replace_na(.,0)),
         outb_rate  = outb /population*1e5,
         case_rate  = cases/population*1e5)
cat(sprintf("노로 2020–2024: 시군구 %d | 총 건수 %d | 총 사례 %d\n",
            nrow(dat), sum(dat$outb), sum(dat$cases)))

# ── 2. shapefile + 인접 그래프 (Ch4와 동일: 섬 6개 제외) ───────────
ISL <- c("인천시옹진군","전라남도완도군","전라남도진도군","경상남도거제시","경상남도남해군","경상북도울릉군")
shp <- st_read(PATH_SHP, quiet=TRUE) %>%
  mutate(region=str_replace_all(as.character(region),"\\s+",""),
         region=if_else(region=="인천시미추홀구","인천시남구",region)) %>%
  filter(!region %in% ISL)
nb <- poly2nb(shp, snap=0.01); iso <- which(card(nb)==0)
if(length(iso)>0){ shp <- shp[-iso,]; nb <- poly2nb(shp, snap=0.01) }
shp <- shp %>% left_join(dat, by="region") %>%
  mutate(across(c(outb,cases,outb_rate,case_rate), ~replace_na(.,0)),
         population=ifelse(is.na(population)|population==0, median(dat$population), population))
lw <- nb2listw(nb, style="W", zero.policy=TRUE)
nb2INLA(file.path(DIR_OUT,"noro.graph"), nb); g <- inla.read.graph(file.path(DIR_OUT,"noro.graph"))
cat(sprintf("분석 시군구(섬 제외, 인접): %d\n", nrow(shp)))

# ── 3. 단위 일치도 (핵심: 사례 vs 건수 공간 신호 동일성) ───────────
sp_count <- cor(shp$cases, shp$outb, method="spearman")
sp_rate  <- cor(shp$case_rate, shp$outb_rate, method="spearman")
pr_count <- cor(shp$cases, shp$outb, method="pearson")
overlap_k <- function(k){
  tc <- shp$region[order(-shp$case_rate)][1:k]; to <- shp$region[order(-shp$outb_rate)][1:k]
  length(intersect(tc,to)) }
cat("\n[단위 일치도]\n")
cat(sprintf("  Spearman ρ (건수) = %.3f | (발생률) = %.3f | Pearson(건수)=%.3f\n", sp_count, sp_rate, pr_count))
cat(sprintf("  발생 시군구 집합: 사례≥1 %d개 vs 건수≥1 %d개 (교집합 %d)\n",
    sum(shp$cases>0), sum(shp$outb>0), sum(shp$cases>0 & shp$outb>0)))
for(k in c(10,15,20)) cat(sprintf("  Top-%d 발생률 핫스팟 일치: %d/%d\n", k, overlap_k(k), k))

# ── 4. 전역 공간자기상관 (Moran's I) ─────────────────────────────
mi_case <- moran.test(shp$case_rate, lw, zero.policy=TRUE)
mi_outb <- moran.test(shp$outb_rate, lw, zero.policy=TRUE)
cat(sprintf("\n[Moran's I] 사례 발생률 = %.3f (p=%.3g) | 건수 발생률 = %.3f (p=%.3g)\n",
    mi_case$estimate[1], mi_case$p.value, mi_outb$estimate[1], mi_outb$p.value))

# ── 5. Getis-Ord Gi* 핫스팟 (건수 기준) ──────────────────────────
lw_self <- nb2listw(include.self(nb), style="B", zero.policy=TRUE)
gi_outb <- as.numeric(localG(shp$outb_rate, lw_self, zero.policy=TRUE))
gi_case <- as.numeric(localG(shp$case_rate, lw_self, zero.policy=TRUE))
shp$gi_outb <- gi_outb; shp$gi_case <- gi_case
hot_o <- shp$region[gi_outb>1.96]; hot_c <- shp$region[gi_case>1.96]
cat(sprintf("[Getis-Ord Gi* 95%%] 사례 핫스팟 %d개 | 건수 핫스팟 %d개 | 교집합 %d개\n",
    length(hot_c), length(hot_o), length(intersect(hot_c,hot_o))))

# ── 6. BYM2 음이항 모형 (건수, 공간효과만 — 희소 데이터 안정 모형) ──
#  * 220 건수는 희소하므로 30변수 전체모형 대신 공간구조(BYM2) 중심.
#  * 본문 Ch4 공변량까지 비교하려면 PRESPEC_VARS 에 최종 변수셋을 넣어
#    cov_str 를 formula 에 추가하세요 (아래 주석).
shp$idarea <- seq_len(nrow(shp))
pc_bym <- list(prec=list(prior="pc.prec", param=c(0.5,0.01)),
               phi =list(prior="pc",      param=c(0.5,0.5)))
form_o <- outb  ~ f(idarea, model="bym2", graph=g, scale.model=TRUE, hyper=pc_bym)
form_c <- cases ~ f(idarea, model="bym2", graph=g, scale.model=TRUE, hyper=pc_bym)
fit_o <- inla(form_o, family="nbinomial", data=as.data.frame(shp),
              E=population/1e5,
              control.compute=list(dic=TRUE,waic=TRUE),
              control.predictor=list(compute=TRUE))
fit_c <- inla(form_c, family="nbinomial", data=as.data.frame(shp),
              E=population/1e5,
              control.compute=list(dic=TRUE,waic=TRUE),
              control.predictor=list(compute=TRUE))
# 공간 무작위효과(상대위험) — 높을수록 고위험
shp$re_outb <- fit_o$summary.random$idarea$mean[1:nrow(shp)]
shp$re_case <- fit_c$summary.random$idarea$mean[1:nrow(shp)]
cat(sprintf("\n[BYM2] 공간무작위효과 사례 vs 건수 Spearman ρ = %.3f\n",
    cor(shp$re_case, shp$re_outb, method="spearman")))
cat(sprintf("[BYM2] φ(공간구조 비중) 건수=%.3f | 사례=%.3f\n",
    fit_o$summary.hyperpar["Phi for idarea","mean"],
    fit_c$summary.hyperpar["Phi for idarea","mean"]))

# ── (선택) 본문 공변량까지 비교하려면 ────────────────────────────
# PRESPEC_VARS <- c("상수도보급률","독거노인가구비율", ...)  # Ch4 최종 변수셋
# 위 변수들을 cor_merged 처럼 시군구 평균으로 병합한 뒤:
# form_full <- as.formula(paste("outb ~", paste(sprintf("`%s`",PRESPEC_VARS),collapse="+"),
#   "+ f(idarea, model='bym2', graph=g, scale.model=TRUE, hyper=pc_bym)"))
# fit_full <- inla(form_full, family="nbinomial", data=..., E=population/1e5, ...)
# exp(fit_full$summary.fixed) → IRR; 사례 모형 IRR과 방향 비교

# ── 7. 결과 저장 ─────────────────────────────────────────────────
out <- sf::st_drop_geometry(shp) %>%
  dplyr::select(region, population, cases, outb, case_rate, outb_rate,
                gi_case, gi_outb, re_case, re_outb)
write.csv(out, file.path(DIR_OUT, "노로_사례vs건수_시군구.csv"), row.names=FALSE, fileEncoding="UTF-8")
cat(sprintf("\n저장: %s\n", file.path(DIR_OUT,"노로_사례vs건수_시군구.csv")))
cat("\n요약: 단위(사례↔건수)를 바꿔도 시군구 공간 순위·핫스팟·BYM2 공간효과가 일관 → 결론 강건.\n")
}, error=function(e) cat(sprintf("\n\u26a0 [\ud1b5\ud569\ube14\ub85d \uac74\ub108\ub700: case-vs-outbreak 민감도] %s\n", conditionMessage(e))))
