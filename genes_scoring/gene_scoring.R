# compare Andrei and mine mQTL fetal rsID
# dt1 <- data.table::fread("/project_data/data/bulk/mqtl/fetal/with.rs.txt.gz")
# dt1[, rsID:=unlist(stringr::str_split(V10, " "))[2], by=seq_len(nrow(dt1))]
# dt2 <- data.table::fread("/project_data/data/bulk/mqtl/mQTL_fetal.csv")

setwd("C:/Users/ana.arruda/OneDrive - Helmholtz Zentrum München/Projects/SCZ_T2D/")
project_folder <- "C:/Users/ana.arruda/OneDrive - Helmholtz Zentrum München/Projects/SCZ_T2D/"
gwas_coloc_folder <- paste0(project_folder, "gwas_coloc/")  # paste0(project_folder, "GWASColoc/")
lookup_folder <- paste0(project_folder, "lookups/")  # paste0(project_folder, "GWASColoc/lookups/")

res.qtl.gtex.coloc <- data.table::fread(paste0(gwas_coloc_folder, "qtl_gtex_hyprcoloc_results.csv"))
res.qtl.coloc <- data.table::fread(paste0(gwas_coloc_folder, "qtl_hyprcoloc_results.csv"))
res.qtl.coloc <- rbind(res.qtl.coloc, res.qtl.gtex.coloc, fill=TRUE)

# Get true hyprcoloc regions
gwas.coloc.credset <- data.table::fread(paste0(gwas_coloc_folder, "final_coloc_credible_set.csv"))
credset.qtl.gtex.coloc <- data.table::fread(paste0(gwas_coloc_folder, "qtl_gtex_hyprcoloc_credible_set.csv"))
credset.qtl.coloc <- data.table::fread(paste0(gwas_coloc_folder, "qtl_hyprcoloc_credible_set.csv"))
credset.qtl.coloc <- rbind(credset.qtl.coloc, credset.qtl.gtex.coloc, fill=TRUE)

for (reg in unique(gwas.coloc.credset$region)){
  credset.qtl.coloc[region==reg, VarInGWASColoc:=ifelse(ID %in% gwas.coloc.credset[region==reg, snp], 1, 0), by=c("qtl.type", "tissue", "geneID")]
  credset.qtl.coloc[region==reg, NVarInGWASColoc:=sum(VarInGWASColoc), by=c("qtl.type", "tissue", "geneID")]
  }

res.qtl.coloc <- merge(res.qtl.coloc, unique(credset.qtl.coloc[, .(qtl.type, tissue, geneID, NVarInGWASColoc)]), by=c("qtl.type", "tissue", "geneID"))
res.qtl.coloc[, TrueColoc:=ifelse(NVarInGWASColoc>0, 1, 0)]
res.qtl.coloc <- res.qtl.coloc[TrueColoc==1]

# Manual lookup of gene IDs to convert from b38 to b37
res.qtl.coloc[geneID=="ENSG00000273611", geneID:="ENSG00000108278"]
res.qtl.coloc[geneID=="ENSG00000277161", geneID:="ENSG00000184886"]
res.qtl.coloc[geneID=="ENSG00000278259", geneID:="ENSG00000141140"]
res.qtl.coloc[geneID=="ENSG00000278311", geneID:="ENSG00000260508"]

# for mqtl --> find genes related to CpG
# BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
# BiocManager::install("missMethyl")
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(org.Hs.eg.db)
library(limma)
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

entrez.id <- data.table::data.table()
for (cpg in res.qtl.coloc[qtl.type=="mqtl", geneID]){
  tryCatch({
    entrez.id <- rbind(entrez.id, data.table::data.table(cpg=cpg, entrez.id=missMethyl::getMappedEntrezIDs(cpg, rownames(ann), array.type="EPIC")$sig.eg))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# Get Ensembl ID from Entrez ID
entrez.id$ensembl.id <- mapIds(org.Hs.eg.db, 
                               keys=entrez.id$entrez.id, 
                               keytype="ENTREZID", 
                               column="ENSEMBL")

res.qtl.coloc <- merge(res.qtl.coloc, entrez.id[, .(ensembl.id, cpg)], by.x="geneID", by.y="cpg", all.x=TRUE)

res.qtl.coloc[geneID=="Q9HD26", ensembl.id:="ENSG00000047932"]
res.qtl.coloc[qtl.type %in% c("eqtl","sc_eqtl", "sqtl", "tqtl"), ensembl.id:=geneID]

data.table::fwrite(res.qtl.coloc, paste0(gwas_coloc_folder, "qtl_hyprcoloc_results_geneid.csv"))

deg <- data.table::fread(paste0(lookup_folder, "DEG.csv"))
omim <- data.table::fread(paste0(lookup_folder, "OMIM.csv"))
ko.mice <- data.table::fread(paste0(lookup_folder, "KOmice.csv"))
hc.genes <- data.table::fread(paste0(lookup_folder, "HC.csv"))

all.genes <- data.table::fread(paste0(gwas_coloc_folder, "all_genes.csv"))

# Merge lookup results to one table
lookup <- merge(deg, hc.genes, by=c("gene", "ID"), suffixes=c(".deg", ".hc"))
lookup <- merge(lookup, ko.mice[, .(gene, T2D, SCZ)], by="gene", all=TRUE)
omim <- unique(omim, by="id")
lookup <- merge(lookup, omim[, .(name, id,T2D=t1, SCZ=t2)], by.x=c("gene", "ID"), by.y=c("name", "id"), suffixes=c(".mice", ".omim"), all.x=T)

# Add qtl colocalization information
library(dplyr)
scz.tissues <- c("brain", "Excitatory.neurons", "Oligodendrocytes", "cortex_eur", "cerebellum", "Brain_Cerebellum",
                 "Brain_Frontal_Cortex", "Brain_Hypothalamus", "fetal")
t2d.tissues <- c("islets", "adipose", "Adipose_Visceral", "Adipose_Subcutaneous")
qtl <- all.genes[, .(gene=name, ID=id)] %>% .[, `:=` (T2D.qtl=ifelse(ID %in% res.qtl.coloc[tissue %in% t2d.tissues, ensembl.id], 1, 0), 
                                                      SCZ.qtl=ifelse(ID %in% res.qtl.coloc[tissue %in% scz.tissues, ensembl.id], 1, 0))]
lookup <- merge(lookup, qtl, by=c("ID", "gene"))

# Add missense variants information
library(dplyr) 
missense <- data.table::fread(paste0(lookup_folder, "missense.csv"))
missense <- all.genes[, .(gene=name, ID=id)] %>% .[, missense:=ifelse(ID %in% missense$ensembl_gene_name, 1, 0)]

lookup <- merge(lookup, missense, by=c("ID", "gene"))

lookup[, `:=` (T2D.score=sum(T2D.deg, T2D.mice, T2D.omim, T2D.qtl, na.rm=TRUE),
               SCZ.score=sum(SCZ.deg, SCZ.mice, SCZ.omim, SCZ.qtl, na.rm=TRUE)), by=seq_len(nrow(lookup))]
lookup[, Total.score:=sum(T2D.score, SCZ.score), by=seq_len(nrow(lookup))]
lookup[, `:=` (Final.T2D.score=ifelse(T2D.score==0, T2D.hc, T2D.score),
               Final.SCZ.score=ifelse(SCZ.score==0, SCZ.hc, SCZ.score)), by=seq_len(nrow(lookup))]
lookup[, Final.score:=sum(Final.T2D.score, Final.SCZ.score, missense), by=seq_len(nrow(lookup))]

data.table::setcolorder(lookup, c("ID", "gene", "T2D.deg", "SCZ.deg", "T2D.mice", "SCZ.mice", "T2D.omim", "SCZ.omim", "T2D.qtl", "SCZ.qtl", "T2D.score", "SCZ.score", "Total.score",
                      "T2D.hc", "SCZ.hc", "Final.T2D.score", "Final.SCZ.score", "missense", "Final.score"))

data.table::fwrite(lookup, paste0(project_folder, "lookup_results.csv"))

lookup[Final.T2D.score>=1 & Final.SCZ.score>=1, .N]
lookup[Final.T2D.score>=1 & Final.SCZ.score>=1 & Final.score>=3, .N]
lookup[Final.T2D.score>=1 & Final.SCZ.score>=1 & Final.score>=4, .N]

data.table::fwrite(list(lookup[Final.T2D.score>=1 & Final.SCZ.score>=1, gene]), paste0(project_folder, "potential.txt"))
data.table::fwrite(list(lookup[Final.T2D.score>=1 & Final.SCZ.score>=1 & Final.score>=3, gene]), paste0(project_folder, "likely.txt"))
data.table::fwrite(list(lookup[Final.T2D.score>=1 & Final.SCZ.score>=1 & Final.score>=4, gene]), paste0(project_folder, "hc.txt"))

# sensitivity analysis for KO mice
with.motoric <- data.table::fread(paste0(lookup_folder, "KO_lookup_withmotoric_020823.csv"))
without.motoric <- data.table::fread(paste0(lookup_folder, "KO_lookup_nomotoric_020823.csv"))
ko.motoric.lookup <- merge(with.motoric[, .(gene, SCZ)], without.motoric[, .(gene, SCZ)], by="gene", suffixes=c(".with", ".without"))
ko.motoric.lookup[, SCZ.sum:=sum(SCZ.without, SCZ.with, na.rm=TRUE), by=seq_len(nrow(ko.motoric.lookup))]
ko.motoric.lookup[SCZ.sum==1]
lookup[gene %in% ko.motoric.lookup[SCZ.sum==1, gene]]
print("Including motoric traits, KMT2E would be a HC gene and GDAP2 a likely one")
