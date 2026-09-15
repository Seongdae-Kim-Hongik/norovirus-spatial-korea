# Rebuilds the harmonised district-year panel and spatial graph from the main script (stops before model fitting).
# Run from this folder with the main script and ext_wastewater_kosis.csv one level up and FBD_DATA_ZIP/ available there.
L <- readLines("../norovirus_spatial_korea.R", encoding = "UTF-8")
L <- gsub('nb2INLA(nb_obj, file="/tmp/noro.graph"); g_main <- inla.read.graph("/tmp/noro.graph")',
          'g_file <- tempfile(fileext=".graph"); nb2INLA(nb_obj, file=g_file); g_main <- inla.read.graph(g_file)', L, fixed = TRUE)
owd <- setwd(".."); on.exit(setwd(owd))
cut <- grep("^# PART 3\\.5\\. ", L)[1] - 2
eval(parse(text = L[1:cut], encoding = "UTF-8"))
while (sink.number() > 0) sink()
setwd(owd)
saveRDS(mget(c("cor_merged", "shp_main", "g_main", "nb_obj", "FIX_FORMS", "FAMILY")), "analysis_input.rds")
cat("saved analysis_input.rds:", nrow(cor_merged), "rows\n")
