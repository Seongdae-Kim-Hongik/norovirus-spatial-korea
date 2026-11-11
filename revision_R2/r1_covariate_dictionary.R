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

cat(sprintf("  Table 1 base: %d개 (9개 카테고리, 소아 2종 — 결측0만)\n", nrow(TV_v6)))
cat("  ★ 강제변수(‡): 한육우농가수, 성비, 고령인구비율\n")
cat("  ★ 변수선택: 사전지정(이론), 사전고정 변환, 스크리닝 없음 + VIF<5\n\n")


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
