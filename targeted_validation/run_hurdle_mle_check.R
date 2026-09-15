# Post hoc targeted analyses of sludge-moisture content (Supplementary Table S13). Requires analysis_input.rds from prepare_analysis_input.R.
suppressPackageStartupMessages({ library(dplyr); library(stringr); library(glmmTMB) })
x <- readRDS("analysis_input.rds"); d <- as.data.frame(x$cor_merged)
clean_region <- function(z) str_replace_all(as.character(z), "\\s+", "")
d$region <- clean_region(d$region); d <- d[d$region %in% x$shp_main$region & d$population > 0, ]
d$urban <- as.integer(d$area_2 == "도시"); d$year_f <- factor(d$year)
base_data <- "../FBD_DATA_ZIP"
ob <- read.csv(file.path(base_data, "식중독최종_건수.csv"), check.names = FALSE)
ob$region <- clean_region(ob$region)
ob <- ob %>% filter(disease == "노로바이러스", year %in% 2020:2024) %>% group_by(region, year) %>%
  summarise(outbreaks = sum(cases, na.rm = TRUE), .groups = "drop")
d <- left_join(d, ob, by = c("region", "year")); d$outbreaks[is.na(d$outbreaks)] <- 0
sl <- read.csv(file.path(base_data, "merged_하수찌꺼기발생및처리.csv"), check.names = FALSE)
sl$region <- clean_region(sl$region); sl$year <- as.integer(sl$year)
sl$gen <- suppressWarnings(as.numeric(sl[["찌꺼기발생량(톤/년,탈수기준)_합계"]]))
sl$moist <- suppressWarnings(as.numeric(sl[["함수율(%,탈수기준)"]]))
sg <- sl %>% filter(year %in% 2020:2024) %>% group_by(region, year) %>% summarise(
  gen = sum(gen, na.rm = TRUE), moist = if(any(is.finite(moist) & moist > 0)) mean(moist[is.finite(moist) & moist > 0]) else NA_real_, .groups = "drop")
d <- left_join(d, sg, by = c("region", "year")); d$facility <- as.integer(d$gen > 0)
mu <- mean(d$moist[d$facility == 1], na.rm = TRUE); ss <- sd(d$moist[d$facility == 1], na.rm = TRUE)
d$moisture_within_z <- ifelse(d$facility == 1 & is.finite(d$moist), (d$moist - mu)/ss, ifelse(d$facility == 0, 0, NA_real_))
zg <- function(v) { a <- as.numeric(d[[v]]); (a-mean(a,na.rm=TRUE))/sd(a,na.rm=TRUE) }
d$fiscal_z <- zg("재정자립도"); d$elderly_z <- zg("고령인구비율"); d$child_z <- zg("영유아비율_0_4")
rows <- list(); k <- 0
for (outcome in c("cases", "outbreaks")) for (sn in c("principal", "drop_namwon2024", "drop_namwon_anseong2024")) {
  dd <- d
  if(sn == "drop_namwon2024") dd <- dd[!(dd$region == "전라북도남원시" & dd$year == 2024),]
  if(sn == "drop_namwon_anseong2024") dd <- dd[!((dd$region %in% c("전라북도남원시","경기도안성시")) & dd$year == 2024),]
  fm <- as.formula(paste(outcome, "~ facility + moisture_within_z + urban + year_f + fiscal_z + elderly_z + child_z + offset(log(population+1)) + (1|region)"))
  fit <- glmmTMB(fm, family = nbinom2, data = dd)
  cs <- summary(fit)$coefficients$cond
  for(tm in c("facility","moisture_within_z")) { k <- k+1; rows[[k]] <- data.frame(outcome=outcome,specification=sn,term=tm,N=nobs(fit),IRR=exp(cs[tm,"Estimate"]),lo=exp(cs[tm,"Estimate"]-1.96*cs[tm,"Std. Error"]),hi=exp(cs[tm,"Estimate"]+1.96*cs[tm,"Std. Error"])) }
}
out <- bind_rows(rows); out$credible <- out$lo > 1 | out$hi < 1
write.csv(out, "moisture_hurdle_mle_results.csv", row.names = FALSE)
print(out, row.names = FALSE)
