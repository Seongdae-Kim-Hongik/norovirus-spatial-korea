# Evidence grading fits (Supplementary Table S14). Run from a folder containing norovirus_spatial_korea.R, ext_wastewater_kosis.csv and FBD_DATA_ZIP/:
#   MODE=principal Rscript evidence_fit.R ; MODE=ridge31 Rscript evidence_fit.R
# The sensitivity table combines output/sens_<spec>_fixed_effects.csv from sensitivity_filled_years_influence.R.
# MODE=principal : principal specification (VIF<5, default fixed-effect prior) -> exact posterior direction probabilities
# MODE=ridge31   : all pre-specified covariates retained (no VIF pruning) with a ridge-type N(0,1) prior on fixed effects
MODE <- Sys.getenv("MODE", "principal")
L <- readLines("norovirus_spatial_korea.R", encoding = "UTF-8")
L <- gsub('"/tmp/noro.graph"', "g_file", L, fixed = TRUE)
i <- grep("nb2INLA(nb_obj, file=g_file)", L, fixed = TRUE); stopifnot(length(i) == 1)
L[i] <- paste0("g_file <- tempfile(fileext = \".graph\"); ", L[i])
if (MODE == "ridge31") {
  i <- grep("^VIF_THRESHOLD <- 5$", L); stopifnot(length(i) == 1); L[i] <- "VIF_THRESHOLD <- 1e9"
  j <- grep('fit<-tryCatch(inla(as.formula(paste("cases ~",cov_str,', L, fixed = TRUE); stopifnot(length(j) == 1)
  stopifnot(grepl("control.family=list(),", L[j + 1], fixed = TRUE))
  L[j + 1] <- sub("control.family=list(),", "control.family=list(),control.fixed=list(prec=1,prec.intercept=1e-3),", L[j + 1], fixed = TRUE)
}
p6 <- grep("^# PART 6\\. ", L)[1] - 2
code <- c(L[1:p6], '
MODE <- Sys.getenv("MODE", "principal")
out <- do.call(rbind, lapply(list(Total = res_final, Urban = res_urban, Rural = res_rural), function(r) NULL))
rows <- list()
for (nm in c("Total","Urban","Rural")) { r <- get(c(Total="res_final",Urban="res_urban",Rural="res_rural")[[nm]])
  f <- r$fit$summary.fixed; mf <- r$fit$marginals.fixed
  for (t in setdiff(rownames(f), "(Intercept)")) {
    p_neg <- tryCatch(INLA::inla.pmarginal(0, mf[[t]]), error = function(e) NA)
    rows[[length(rows) + 1]] <- data.frame(mode = MODE, model = nm, N = r$N_final, term = t, IRR = exp(f[t, "mean"]), lo = exp(f[t, "0.025quant"]), hi = exp(f[t, "0.975quant"]),
                                           P_pos = 1 - p_neg, P_neg = p_neg) } }
res <- do.call(rbind, rows); dir.create("output", showWarnings = FALSE)
write.csv(res, file.path("output", paste0("evidence_", MODE, ".csv")), row.names = FALSE)
cat("done", MODE, nrow(res), "\\n")')
eval(parse(text = code, encoding = "UTF-8"))
