# Post hoc targeted analyses of sludge-moisture content (Supplementary Table S13). Requires analysis_input.rds from prepare_analysis_input.R.
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(INLA); library(glmmTMB)
})


# Same estimation settings as the public analysis script: VB correction off + degeneracy guard
inla <- function(formula, ..., control.inla = list()) {
  control.inla$control.vb <- list(enable = FALSE)
  fit <- INLA::inla(formula, ..., control.inla = control.inla)
  fs <- paste(deparse(formula, width.cutoff = 500L), collapse = " ")
  if (grepl("f\\(", fs) && !is.null(fit$dic) && is.finite(fit$dic$dic)) {
    ns <- gsub("\\+\\s*f\\((?:[^()]|\\([^()]*\\))*\\)", "", fs, perl = TRUE)
    ref <- tryCatch(INLA::inla(as.formula(ns), ..., control.inla = control.inla), error = function(e) NULL)
    tries <- 0
    while (!is.null(ref) && is.finite(ref$dic$dic) && fit$dic$dic > ref$dic$dic + 50 && tries < 5) {
      tries <- tries + 1; cat(sprintf("  [guard] degenerate fit (DIC %.1f vs non-spatial %.1f): refit %d\n", fit$dic$dic, ref$dic$dic, tries))
      fit <- INLA::inla(formula, ..., control.inla = control.inla) }
  }
  fit }

x <- readRDS("analysis_input.rds")
d <- as.data.frame(x$cor_merged)
shp_main <- x$shp_main
g_main <- x$g_main

clean_region <- function(z) str_replace_all(as.character(z), "\\s+", "")
d$region <- clean_region(d$region)
d <- d[d$region %in% shp_main$region & d$population > 0, ]
d$urban <- as.integer(d$area_2 == "도시")
d$year_f <- factor(d$year)

# Outbreak-event counts, aggregated to the same district-year panel.
base_data <- "../FBD_DATA_ZIP"
ob <- read.csv(file.path(base_data, "식중독최종_건수.csv"), check.names = FALSE, stringsAsFactors = FALSE)
ob$region <- clean_region(ob$region)
ob <- ob %>% filter(disease == "노로바이러스", year %in% 2020:2024) %>%
  group_by(region, year) %>% summarise(outbreaks = sum(cases, na.rm = TRUE), .groups = "drop")
d <- left_join(d, ob, by = c("region", "year"))
d$outbreaks[is.na(d$outbreaks)] <- 0

# Facility indicator from positive sludge generation. A recorded zero moisture value
# at a facility is treated as unavailable for the within-facility moisture slope.
sl <- read.csv(file.path(base_data, "merged_하수찌꺼기발생및처리.csv"), check.names = FALSE, stringsAsFactors = FALSE)
sl$region <- clean_region(sl$region); sl$year <- as.integer(sl$year)
sl$gen <- suppressWarnings(as.numeric(sl[["찌꺼기발생량(톤/년,탈수기준)_합계"]]))
sl$moist_raw <- suppressWarnings(as.numeric(sl[["함수율(%,탈수기준)"]]))
sg <- sl %>% filter(year %in% 2020:2024) %>% group_by(region, year) %>%
  summarise(sludge_generated = sum(gen, na.rm = TRUE),
            moisture_measured = if (any(is.finite(moist_raw) & moist_raw > 0)) mean(moist_raw[is.finite(moist_raw) & moist_raw > 0]) else NA_real_,
            .groups = "drop")
d <- left_join(d, sg, by = c("region", "year"))
d$facility <- as.integer(is.finite(d$sludge_generated) & d$sludge_generated > 0)

z_global <- function(v, log_form = FALSE) {
  a <- suppressWarnings(as.numeric(d[[v]])); if (log_form) a <- log1p(pmax(a, 0))
  (a - mean(a, na.rm = TRUE)) / sd(a, na.rm = TRUE)
}
d$ranch_z <- z_global("목장용지", TRUE)
d$moisture_z <- z_global("함수율(%,탈수기준)", FALSE)
d$groundwater_z <- z_global("민방위용_개소수", TRUE)
d$onsite_z <- z_global("자체처리량(톤/년)_계", TRUE)
d$fiscal_z <- z_global("재정자립도", FALSE)
d$elderly_z <- z_global("고령인구비율", FALSE)
d$child_z <- z_global("영유아비율_0_4", FALSE)

mm <- d$moisture_measured
mm_mean <- mean(mm[d$facility == 1], na.rm = TRUE)
mm_sd <- sd(mm[d$facility == 1], na.rm = TRUE)
d$moisture_within_z <- ifelse(d$facility == 1 & is.finite(mm), (mm - mm_mean) / mm_sd,
                             ifelse(d$facility == 0, 0, NA_real_))

rmap <- data.frame(region = shp_main$region, idarea = seq_len(nrow(shp_main)))
d <- left_join(d, rmap, by = "region")
d <- d[order(d$idarea, d$year), ]
d$idarea_time <- seq_len(nrow(d))

pc_bym <- list(prec.unstruct = list(prior = "pc.prec", param = c(0.5, 0.01)),
               prec.spatial = list(prior = "pc.prec", param = c(0.5, 0.01)))
pc_prec <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))

exposures <- c(ranch = "ranch_z", moisture = "moisture_z",
               groundwater = "groundwater_z", onsite = "onsite_z")
specs <- list(
  minimal = character(),
  sociodemographic = c("fiscal_z", "elderly_z", "child_z")
)

term_row <- function(fit, term, engine, outcome, exposure, spec, n) {
  s <- fit$summary.fixed
  if (!term %in% rownames(s)) return(NULL)
  data.frame(engine = engine, outcome = outcome, exposure = exposure, specification = spec,
             term = term, N = n, IRR = exp(s[term, "mean"]),
             lo = exp(s[term, "0.025quant"]), hi = exp(s[term, "0.975quant"]),
             stringsAsFactors = FALSE)
}

fit_one <- function(outcome, exposure_name, exposure_var, spec_name, adjust) {
  needed <- unique(c(outcome, exposure_var, "urban", "year_f", adjust, "population", "idarea", "idarea_time", "region"))
  dat <- d[complete.cases(d[, needed]), needed]
  dat$exposure_rural <- dat[[exposure_var]] * (1 - dat$urban)
  dat$exposure_urban <- dat[[exposure_var]] * dat$urban
  fixed <- c("exposure_rural", "exposure_urban", "urban", "year_f", adjust)
  spatial <- "f(idarea, model='bym', graph=g_main, scale.model=TRUE, hyper=pc_bym) + f(idarea_time, model='iid', hyper=pc_prec)"
  fm <- as.formula(paste(outcome, "~", paste(fixed, collapse = " + "), "+ offset(log(population + 1)) +", spatial))
  fit <- inla(fm, family = "nbinomial", data = dat,
              control.inla = list(control.vb = list(enable = FALSE)),
              control.compute = list(dic = TRUE, waic = TRUE),
              control.predictor = list(link = 1), verbose = FALSE)
  rows <- rbind(term_row(fit, "exposure_rural", "INLA_M4", outcome, exposure_name, spec_name, nrow(dat)),
                term_row(fit, "exposure_urban", "INLA_M4", outcome, exposure_name, spec_name, nrow(dat)))

  # Non-spatial district-random-intercept check with the same fixed effects.
  fm_mle <- as.formula(paste(outcome, "~", paste(fixed, collapse = " + "),
                             "+ offset(log(population + 1)) + (1|region)"))
  mle <- tryCatch(glmmTMB(fm_mle, family = nbinom2, data = dat), error = function(e) NULL)
  if (!is.null(mle)) {
    cs <- summary(mle)$coefficients$cond
    for (tm in c("exposure_rural", "exposure_urban")) if (tm %in% rownames(cs)) {
      rows <- rbind(rows, data.frame(engine = "glmmTMB_RE", outcome = outcome, exposure = exposure_name,
        specification = spec_name, term = tm, N = nrow(dat), IRR = exp(cs[tm, "Estimate"]),
        lo = exp(cs[tm, "Estimate"] - 1.96 * cs[tm, "Std. Error"]),
        hi = exp(cs[tm, "Estimate"] + 1.96 * cs[tm, "Std. Error"])))
    }
  }
  rows
}

results <- list(); k <- 0
for (outcome in c("cases", "outbreaks")) for (nm in names(exposures)) for (sp in names(specs)) {
  cat("fit", outcome, nm, sp, "\n"); flush.console()
  k <- k + 1; results[[k]] <- fit_one(outcome, nm, exposures[[nm]], sp, specs[[sp]])
}

# Moisture hurdle model: facility status and within-facility moisture are separate.
fit_hurdle <- function(outcome, sensitivity = "principal") {
  dd <- d
  if (sensitivity == "drop_namwon2024") dd <- dd[!(dd$region == "전라북도남원시" & dd$year == 2024), ]
  if (sensitivity == "drop_namwon_anseong2024") dd <- dd[!((dd$region == "전라북도남원시" | dd$region == "경기도안성시") & dd$year == 2024), ]
  adjust <- c("fiscal_z", "elderly_z", "child_z")
  needed <- c(outcome, "facility", "moisture_within_z", "urban", "year_f", adjust,
              "population", "idarea", "idarea_time", "region")
  dat <- dd[complete.cases(dd[, needed]), needed]
  fixed <- c("facility", "moisture_within_z", "urban", "year_f", adjust)
  fm <- as.formula(paste(outcome, "~", paste(fixed, collapse = " + "),
    "+ offset(log(population + 1)) + f(idarea, model='bym', graph=g_main, scale.model=TRUE, hyper=pc_bym) + f(idarea_time, model='iid', hyper=pc_prec)"))
  fit <- inla(fm, family = "nbinomial", data = dat,
              control.inla = list(control.vb = list(enable = FALSE)),
              control.compute = list(dic = TRUE, waic = TRUE), control.predictor = list(link = 1), verbose = FALSE)
  out <- rbind(term_row(fit, "facility", "INLA_M4", outcome, "moisture_hurdle", sensitivity, nrow(dat)),
               term_row(fit, "moisture_within_z", "INLA_M4", outcome, "moisture_hurdle", sensitivity, nrow(dat)))
  out
}
for (outcome in c("cases", "outbreaks")) { cat("fit", outcome, "moisture_hurdle\n"); flush.console(); k <- k + 1; results[[k]] <- fit_hurdle(outcome, "principal") }

res <- bind_rows(results) %>% mutate(credible = lo > 1 | hi < 1)
write.csv(res, "targeted_model_results.csv", row.names = FALSE)

audit <- data.frame(
  item = c("analytic district-years", "case total", "outbreak total", "facility district-years",
           "nonfacility district-years", "facility rows with positive measured moisture",
           "facility rows with zero/unavailable moisture", "within-facility moisture mean", "within-facility moisture SD"),
  value = c(nrow(d), sum(d$cases), sum(d$outbreaks), sum(d$facility == 1), sum(d$facility == 0),
            sum(d$facility == 1 & is.finite(d$moisture_measured)),
            sum(d$facility == 1 & !is.finite(d$moisture_measured)), mm_mean, mm_sd)
)
write.csv(audit, "targeted_model_audit.csv", row.names = FALSE)

influence <- list(); j <- 0
for (outcome in c("cases", "outbreaks")) for (sn in c("principal", "drop_namwon2024", "drop_namwon_anseong2024")) {
  cat("influence", outcome, sn, "\n"); flush.console()
  j <- j + 1; influence[[j]] <- fit_hurdle(outcome, sn)
}
influence <- bind_rows(influence) %>% mutate(credible = lo > 1 | hi < 1)
write.csv(influence, "moisture_influence_results.csv", row.names = FALSE)
print(res, row.names = FALSE)
