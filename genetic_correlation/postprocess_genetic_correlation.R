res.dt <- data.table::data.table()
for (i in c("T2D.eur", "T2D.multiancestry", "SCZ.eur", "SCZ.multiancestry")){
  dt <- data.table::fread(paste0("/project_data/GeneticCorrelation/", i, ".correlation.result"))
  res.dt <- rbind(res.dt, dt)
}

res.dt[, p1:=sub(".ldsc.sumstats.gz", "", p1)]
res.dt[, p2:=sub(".ldsc.sumstats.gz", "", p2)]

data.table::fwrite(res.dt, "/project_data/GeneticCorrelation/results_ldsc_analysis.csv")
data.table::fwrite(res.dt[, .(p1,p2,rg,se,p)], "/project_data/GeneticCorrelation/short_results_ldsc_analysis.csv")
