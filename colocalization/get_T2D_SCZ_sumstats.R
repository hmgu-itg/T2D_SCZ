library(dplyr)

project_folder <- "/lustre/groups/itg/teams/zeggini/projects/SCZ_T2D/"
source("/project_data/scripts/read_files_config.R")

region <- data.table::fread(paste0(output.path, "GWAS_regions.csv")) %>% .[CHR!=23]
credset <- data.table::fread(paste0(output.path, "credible_set.csv"))
credset[, SNP:=ID]
result <- data.table::fread(paste0(output.path, "final_coloc_result_pp4.csv"))

read.data <- function(data.path, gwas.coloc.result, chr.col="CHR", pos.col="POS", range=1e+6){
  dt <- data.table::fread(data.path)
  res.dt <- data.table::data.table()
  
  for(i in 1:nrow(gwas.coloc.result)){
    reg <- gwas.coloc.result[i]
    res.dt <- rbind(res.dt, dt[get(chr.col)==reg$CHR
                               & data.table::inrange(get(pos.col), reg$POS-range, reg$POS+range)])
  }
  return(res.dt)
}

gwas <- lapply(traits, function(i){read.data(data.path=paste0(output.path, "GWAS_", i, "_precoloc.txt"), gwas.coloc.result=result, chr.col="chr", pos.col="pos")})
names(gwas) <- traits

data.table::fwrite(gwas[[1]], paste0(data.path, "T2D_gtex_coloc.txt"))
data.table::fwrite(gwas[[2]], paste0(data.path, "SCZ_gtex_coloc.txt"))

# Get regions for colocalization
region <- data.table::fread(paste0(output.path, "GWAS_regions.csv")) %>% .[CHR!=23]
result <- data.table::fread(paste0(output.path, "final_coloc_result_pp4.csv"))
range=1e+6
coloc.regions <- result[, .(region, CHR, POS, lead.SNP, rsID)]
coloc.regions[, `:=`(start=POS-range, end=POS+range)]

data.table::fwrite(coloc.regions, paste0(output.path, "coloc_regions.csv"))

# Mapping data
mapping.file <- data.table::fread("data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz")
mapping.file <- data.table::data.table()

for (i in 1:nrow(coloc.regions)){
  pos <- coloc.regions[i, POS]
  chr <- coloc.regions[i, CHR]
  mapping.dt <- rbind(mapping.dt, mapping.file[sub("chr", "", CHR)==chr & data.table::inrange(POS, pos-range, pos+range)])
}

data.table::fwrite(mapping.dt, "data/GTEX_mapping_file.txt")


