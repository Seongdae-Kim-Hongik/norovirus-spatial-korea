# =============================================================================
#  Sensitivity of the principal estimates to the filled covariate years and to
#  individual years and district-years (Supplementary Table S12).
#
#  Usage (from the repository folder, with FBD_DATA_ZIP/ and ext_wastewater_kosis.csv
#  in place, as for norovirus_spatial_korea.R):
#      SENS=<specification> Rscript sensitivity_filled_years_influence.R
#
#  <specification>
#    principal         nearest-year filling (main analysis)
#    timeinv           filled covariates replaced by the district mean of observed years
#                      (livestock head / farm counts and oyster production included)
#    linear            linear interpolation between observed years
#    excl              covariates requiring filling omitted
#    obs2023           2020-2023 only; covariates with any filled year in 2020-2023 omitted
#    yrdrop_<year>     all district-years of <year> excluded (2020 ... 2024)
#    drop_namwon2024   district-year with the largest case count excluded
#    drop_anseong2024  district-year with the second-largest case count excluded
#    drop_bothall      both districts excluded in all years
#
#  The main script is run up to (not including) PART 6, after modifying the covariate
#  panel immediately after year harmonisation. Estimates of the Total, Urban and Rural
#  principal models (M4) are written to output/sens_<specification>_fixed_effects.csv.
#  Run specifications in separate folders when running them in parallel.
# =============================================================================
SENS <- Sys.getenv("SENS", "principal")
L <- readLines("norovirus_spatial_korea.R", encoding = "UTF-8")

if (SENS == "obs2023") {
  i <- grep("^YEAR_START <- 2020; YEAR_END <- 2024$", L); stopifnot(length(i) == 1)
  L[i] <- "YEAR_START <- 2020; YEAR_END <- 2023"
}
# per-run adjacency-graph file
L <- gsub('"/tmp/noro.graph"', "g_file", L, fixed = TRUE)
i <- grep("nb2INLA(nb_obj, file=g_file)", L, fixed = TRUE); stopifnot(length(i) == 1)
L[i] <- paste0("g_file <- tempfile(fileext = \".graph\"); ", L[i])

# record the filled years of each harmonised covariate and its unfilled values
i1 <- grep("x[cm$year %in% bad] <- NA; cm[[v]] <- x", L, fixed = TRUE); stopifnot(length(i1) == 1)
L[i1] <- sub("x[cm$year %in% bad] <- NA; cm[[v]] <- x",
             "x[cm$year %in% bad] <- NA; cm[[v]] <- x; IMPUTED[[v]] <- bad; RAWCM[[v]] <- cm[, c(\"region\",\"year\",v)]", L[i1], fixed = TRUE)
i0 <- grep("^HARMONISE <- intersect\\(HARMONISE, names\\(cm\\)\\)", L); stopifnot(length(i0) == 1)
L[i0] <- paste0(L[i0], "\nIMPUTED <- list(); RAWCM <- list()")

# modify the harmonised panel
i2 <- grep("^cor_merged <- cm$", L); stopifnot(length(i2) == 1)
L[i2] <- paste0(L[i2], '
SENS <- Sys.getenv("SENS", "principal")
OYSTER <- "굴_자연채묘 생산량(kg)"
if (SENS == "obs2023") {
  cor_merged <- cor_merged[cor_merged$year <= 2023, ]
  for (v in c(names(IMPUTED)[sapply(IMPUTED, function(b) any(b <= 2023))], OYSTER)) cor_merged[[v]] <- NA_real_
} else if (SENS == "timeinv") {
  for (v in names(IMPUTED)) {
    raw <- RAWCM[[v]]; raw <- raw[raw$year <= 2023 & !(raw$year %in% IMPUTED[[v]]), ]
    mu <- tapply(suppressWarnings(as.numeric(raw[[v]])), raw$region, mean, na.rm = TRUE)
    cor_merged[[v]] <- unname(mu[as.character(cor_merged$region)])
  }
  for (v in intersect(c("사육두수(두)_합계","사육두수(두)_한육우","사육두수(두)_가금","농가수(호)_합계","농가수(호)_한육우", OYSTER), names(cor_merged))) {
    s <- cor_merged[cor_merged$year <= 2023 & !(v == OYSTER & cor_merged$year == 2021), ]
    mu <- tapply(suppressWarnings(as.numeric(s[[v]])), s$region, mean, na.rm = TRUE)
    cor_merged[[v]] <- unname(mu[as.character(cor_merged$region)])
  }
} else if (SENS == "linear") {
  for (v in names(IMPUTED)) {
    raw <- RAWCM[[v]]; x <- suppressWarnings(as.numeric(raw[[v]])); out <- x
    for (rg in unique(raw$region)) { k <- which(raw$region == rg); ok <- !is.na(x[k])
      if (sum(ok) >= 2) out[k] <- approx(raw$year[k][ok], x[k][ok], xout = raw$year[k], rule = 2)$y else if (sum(ok) == 1) out[k] <- x[k][ok] }
    lk <- setNames(out, paste(raw$region, raw$year)); cor_merged[[v]] <- unname(lk[paste(cor_merged$region, cor_merged$year)])
  }
} else if (SENS == "excl") {
  for (v in names(IMPUTED)) cor_merged[[v]] <- NA_real_
} else if (grepl("^yrdrop_", SENS)) {
  cor_merged <- cor_merged[cor_merged$year != as.integer(sub("yrdrop_", "", SENS)), ]
} else if (grepl("^drop_", SENS)) {
  dy <- list(drop_namwon2024 = list(r = "전라북도남원시", y = 2024), drop_anseong2024 = list(r = "경기도안성시", y = 2024),
             drop_bothall = list(r = c("전라북도남원시", "경기도안성시"), y = 2020:2024))[[SENS]]
  cor_merged <- cor_merged[!(cor_merged$region %in% dy$r & cor_merged$year %in% dy$y), ]
}
cat("  sensitivity specification:", SENS, "| district-year rows:", nrow(cor_merged), "\\n")')

p6 <- grep("^# PART 6\\. ", L)[1] - 2
code <- c(L[1:p6], '
fe_out <- function(r, nm) { f <- r$fit$summary.fixed; f <- f[rownames(f) != "(Intercept)", ]
  data.frame(model = nm, N = r$N_final, DIC = round(r$dic, 2), term = rownames(f), IRR = round(exp(f$mean), 2),
             lo = round(exp(f[, "0.025quant"]), 2), hi = round(exp(f[, "0.975quant"]), 2)) }
dir.create("output", showWarnings = FALSE)
write.csv(rbind(fe_out(res_final, "Total"), fe_out(res_urban, "Urban"), fe_out(res_rural, "Rural")),
          file.path("output", paste0("sens_", SENS, "_fixed_effects.csv")), row.names = FALSE)
cat("done:", SENS, "\\n")')
eval(parse(text = code, encoding = "UTF-8"))
