project_folder <- "C:/Users/ana.arruda/OneDrive - Helmholtz Zentrum München/Projects/SCZ_T2D/"
setwd(project_folder)

#----------- Druggable genome -----------
druggable.genome <- data.table::as.data.table(readxl::read_xlsx(paste0(project_folder, "druggable_genome.xlsx"), sheet=1, skip=0))
hc.genes <- c(read.table(paste0(project_folder, "hc.txt")))$V1
likely.genes <- c(read.table(paste0(project_folder, "likely.txt")))$V1

likely.dt <- data.table::data.table(genes=likely.genes)
likely.dt <- merge(likely.dt, druggable.genome[, .(hgnc_names, druggability_tier)], by.x="genes", by.y="hgnc_names", all.x=TRUE)
likely.dt[is.na(druggability_tier), druggability_tier:="not included"]

hc.dt <- druggable.genome[hgnc_names %in% hc.genes]

#----------- Pharos -----------
pharos.dt <- data.table::fread("pharos_query_results.csv")
pharos.dt[Symbol %in% likely.genes]
pharos.dt <- pharos.dt[Symbol %in% likely.genes]
likely.dt <- merge(likely.dt, pharos.dt[, .(Symbol, T_level=`Target Development Level`)], by.x="genes", by.y="Symbol")

#----------- Manual lookup of Tier 1 or Tclin genes
likely.dt[druggability_tier=="Tier 1" | T_level=="Tclin"]

data.table::fwrite(likely.dt, "lookup_druggability.csv")


