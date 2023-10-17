create_bed <- function(input_file, output_file){
  if(is.character(input_file)) out.dt <- data.table::fread(input_file)
  else out.dt <- input_file
  out.dt$start <- out.dt$POS-1
  data.table::fwrite(unique(out.dt[, .(CHR=as.integer(CHR), start=as.integer(start), POS=as.integer(POS), EA, beta)], by=c("CHR", "POS")), 
                     output_file, col.names = F, row.names = F, quote = F, sep = "\t")
}

process_liftdown <- function(input_file, bed_file_hg19, output_file){
  if(is.character(input.file)) out.dt <- data.table::fread(input_file)
  else out.dt <- input_file
  liftedover = data.table::fread(bed_file_hg19)
  colnames(liftedover) <- c("CHR", "start", "POS", "EA", "beta", "ID")
  liftedover[, start:=NULL]
  dim(liftedover)
  
  # Remove the positions not mapping to chromosomes 1 to 22
  liftedover <- subset(liftedover, CHR %in% 1:22)
  dim(liftedover)
  # Do variants map to several positions ?
  nb.positions <- table(liftedover$variant_description)
  head(nb.positions)
  table(nb.positions)
  
  # Create the final file
  liftedover[, `:=` (beta=NULL, EA=NULL)]
  liftedover <- unique(liftedover)
  res <- merge(out.dt, liftedover, by="ID")
  rm(out.dt)
  rm(liftedover)
  
  res$CHR=res$CHR.y
  res$POS=res$POS.y
  res[, `:=` (CHR.x=NULL, CHR.y=NULL, POS.x=NULL, POS.y=NULL)]
  
  data.table::fwrite(res, output_file)
}

preprocess.gtex <- function(tissue, mqtl.type="eQTL", data.path=data.path, tmp.path=tmp.path){
  # https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/README_eQTL_v8.txt
  dt <- data.table::fread(paste0("/reference_data/GTEx/GTEx_Analysis_v8_", mqtl.type, "/", tissue, ".allpairs.txt.gz"), tmpdir=tmp.path)
  sig.QTL <- data.table::fread(paste0("/reference_data/GTEx/GTEx_Analysis_v8_", mqtl.type, "/", tissue, ".v8.signif_variant_gene_pairs.txt.gz"), tmpdir=tmp.path)[, sig.QTL:=1]
  
  out.dt <- merge(dt, sig.QTL, by=c("variant_id", "gene_id", "tss_distance", "ma_samples", "ma_count",
                                     "pval_nominal", "slope", "slope_se"), all.x=TRUE)
  rm(dt)
  rm(sig.QTL)
  
  out.dt <- out.dt[, .(CHR=as.integer(sub("chr", "", sub("_.*", "", variant_id))), 
                       POS=as.integer(stringr::str_match(variant_id, ".*_(.*?)_.*_.*_.*")[,2]), 
                       geneID=gsub("\\..*", "", gene_id), pval=pval_nominal, beta=slope, se=slope_se, MAF=maf.x, 
                       EA=stringr::str_match(variant_id, ".*_.*_.*_(.*?)_.*")[,2], 
                       NEA=stringr::str_match(variant_id, ".*_.*_(.*?)_.*_.*")[,2], indep.QTL=ifelse(is.na(sig.QTL), 0, 1))]
  out.dt <- out.dt[!is.na(CHR)]

  
  data.table::fwrite(out.dt, paste0(data.path, "bulk/eqtl/eQTL_", tissue, "_GTEx_hg38.csv"))
  # rm(out.dt)
  
  #### Liftover ####
  # Create .bed file
  # create_bed(paste0(data.path, "bulk/eqtl/eQTL_", tissue, "_GTEx_hg38.csv"), paste0(data.path, "eqtl_", tissue, "_GTEx.hg38.bed"))
  create_bed(input_file=out.dt, output_file=paste0(data.path, "bulk/", mqtl.type, "/", "_", tissue, ".hg38.bed"))
  
  ### ON WORKER: Run CrossMap to convert the positions on 
  if(!file.exists(paste0(data.path, "bulk/", mqtl.type, "/", mqtl.type, "_", tissue, ".hg19.bed"))){
    # system(paste0("CrossMap.py bed ", data.path, "hg38ToHg19.over.chain.gz ", data.path, "bulk/", mqtl.type, "/", mqtl.type, "_", tissue, "_GTEx.hg38.bed ", data.path, "bulk/", mqtl.type, "/", mqtl.type, "_", tissue, "_GTEx.hg19.bed 2>&1"))
    system(paste0("CrossMap.py bed ", data.path, "hg38ToHg19.over.chain.gz ", data.path, "bulk/", mqtl.type, "/", mqtl.type, "_", tissue, "_GTEx.hg38.bed ", data.path, "bulk/", mqtl.type, "/", mqtl.type, "_", tissue, "_GTEx.hg19.bed"))
  }
  
  # process_liftdown(paste0(data.path, "bulk/", mqtl.type, "/", mqtl.type, "_", tissue, "_hg38.csv"), 
  process_liftdown(out.dt, 
                   paste0(data.path, "bulk/", mqtl.type, "/", mqtl.type, "_", tissue, ".hg19.bed"),
                   paste0(data.path, "bulk/", mqtl.type, "/", mqtl.type, "_", tissue, "hg19.csv"))
  
  # Import the lifted positions and keep only autosomes
  liftedover = data.table::fread(paste0(data.path, mqtl.type, "_", tissue, "_GTEx.hg19.bed"))
  colnames(liftedover) <- c("CHR", "start", "POS", "EA", "beta", "ID")
  liftedover[, start:=NULL]
  dim(liftedover)
  #Remove the positions not mapping to chromosomes 1 to 22
  liftedover <- subset(liftedover, CHR %in% 1:22)
  dim(liftedover)
  # Do variants map to several positions ?
  nb.positions <- table(liftedover$variant_description)
  head(nb.positions)
  table(nb.positions)
  
  # Create the final file
  liftedover[, `:=` (beta=NULL, EA=NULL)]
  liftedover <- unique(liftedover)
  res <- merge(out.dt, liftedover, by="ID")
  rm(out.dt)
  rm(liftedover)
  
  res$CHR=res$CHR.y
  res$POS=res$POS.y
  res[, `:=` (CHR.x=NULL, CHR.y=NULL, POS.x=NULL, POS.y=NULL)]
  
  data.table::fwrite(res[, ID:=NULL], paste0(data.path, "bulk/", mqtl.type, "/", mqtl.type, "_", tissue, "_GTEx_hg19_noid.csv"))
  
  # add ID
  system()
}

project_folder <- "/lustre/groups/itg/teams/zeggini/projects/SCZ_T2D/"
source(paste0(project_folder, "scripts/read_files_config.R"))

#######################################################################################
#--------------------------------- sQTL and eQTL GTEx --------------------------------#
#######################################################################################
tissues.qtl <- list(c("eQTL", "Adipose_Visceral_Omentum"), c("sQTL", "Adipose_Visceral_Omentum"), 
                    c("sQTL", "Adipose_Subcutaneous"), c("sQTL", "Liver"))  #  c("eQTL", "Adipose_Subcutaneous"), c("eQTL", "Liver")

tissue="Adipose_Visceral_Omentum"
mqtl.type="eQTL"

for (t in tissues.qtl){
  if(!file.exists(paste0(data.path, "bulk/", t[1], "/", t[1], "_", t[2], "_GTEx_hg19_noid.csv")) ||
     !file.exists(paste0(data.path, "bulk/", t[1], "/", t[1], "_", t[2], "_GTEx_hg19.csv"))){
    preprocess.gtex(tissue=t[2], mqtl.type=t[1])
  }
}

#######################################################################################
#--------------------------- sQTL from pancreatic islets -----------------------------#
#######################################################################################
# dt <- data.table::fread("/project_data/data/bulk/sqtl/adipose/spliceQTL_METSIM_n426_adipose_summaryStats_250kb_all.txt.gz", tmpdir="/project_data/tempdata")
dt <- data.table::fread("data/NominallySignificant_sQTLs.txt.gz")
out.dt <- dt[, .(CHR=sub(":.*", "", Variant_id), POS=stringr::str_match(Variant_id, ".*:(.*?):.*:.*")[,2], 
                 geneName=Gene, pval=Nominal_pvalue, beta=Slope, EA=sub(".*:", "", Variant_id), NEA=stringr::str_match(Variant_id, ".*:.*:(.*?):.*")[,2], 
                 indep.QTL=ifelse(Variant_id==Lead_sSNP, 1, 0), Splice=Junction)]
out.dt[, se:=abs(beta/qnorm(pval/2))]
head(out.dt)
rm(dt)

# data.table::fwrite(out.dt, "/project_data/data/bulk/sqtl/sQTL_islets.csv")
data.table::fwrite(out.dt, "data/sQTL_islets.csv")
rm(out.dt)

# add ID

#######################################################################################
#------------------------- sQTL from METSIM adipose tissue ---------------------------#
#######################################################################################
dt <- data.table::fread("/project_data/data/bulk/sqtl/adipose/spliceQTL_METSIM_n426_adipose_summaryStats_250kb_all.txt.gz", tmpdir="/project_data/tempdata")
# dt[pval<9.6e-6, .(count=.N), by=ENSG]

out.dt <- dt[, .(CHR=sub(":.*", "", Variant), POS=sub(".*:", "", sub("_.*", "", Variant)), rsID, 
                 geneID=gsub("\\..*", "", Gene), pval, beta, se=NA, MAF=ifelse(EAF<0.5, EAF, 1-EAF), 
                 EAF, EA=sub(".*/", "", Variant), NEA=sub(".*_", "", sub("/.*", "", Variant)), 
                 indep.eQTL=NA, Splice)]
out.dt[, se:=abs(beta/qnorm(pval/2))]
head(out.dt)
rm(dt)

# check indep.sQTLs

data.table::fwrite(out.dt, "/project_data/data/bulk/sqtl/sQTL_adipose.csv")
rm(out.dt)

# add ID

#######################################################################################
#------------------------- eQTL from METSIM adipose tissue ---------------------------#
#######################################################################################
dt <- data.table::fread("/project_data/data/bulk/eqtl/METSIM_adipose_ciseQTLs_summaryStats_all.txt.gz", tmpdir="/project_data/tempdata")
# dt[pval<9.6e-6, .(count=.N), by=ENSG]

out.dt <- dt[, .(CHR=sub(":.*", "", Variant), POS=sub(".*:", "", sub("_.*", "", Variant)), rsID, 
                 geneID=gsub("\\..*", "", ENSG), pval, beta, se=SE, MAF=ifelse(EAF<0.5, EAF, 1-EAF), 
                 EAF, EA=sub(".*/", "", Variant), NEA=sub(".*_", "", sub("/.*", "", Variant)), 
                 indep.eQTL=NA)]
                 # indep.eQTL=ifelse(pval<9.6e-6, 1, 0))]
rm(dt)

# check indep.eQTLs

data.table::fwrite(out.dt, "/project_data/data/bulk/eqtl/eQTL_adipose.csv")
rm(out.dt)

# add ID

#################################################################################
#------------------------- eQTL from brain MetaBrain ---------------------------#
#################################################################################
for (tissue in c("basalganglia-EUR-30", "hippocampus-EUR-30", "cerebellum-EUR-60", "cortex-AFR-40", "cortex-EAS-30", "cortex-EUR-80")){
  # tissue="hippocampus-EUR-30"
  dt <- data.table::fread(paste0("/project_data/data/bulk/eqtl/MetaBrain/", tissue, "Pcs.txt.gz"), tmpdir="/project_data/tempdata")
  indep.eQTL <- data.table::as.data.table(readxl::read_xlsx(paste0("/project_data/data/bulk/eqtl/MetaBrain/indep_eqtls_", tissue, ".xlsx"))) %>% .[, .(Gene, SNP, indep.eQTL=1)]
  
  out.dt <- merge(dt, indep.eQTL, by=c("Gene", "SNP"), all.x=TRUE)
  out.dt[is.na(indep.eQTL), indep.eQTL:=0]
  out.dt[indep.eQTL==1, .N]
  out.dt[indep.eQTL==0, .N]
  
  out.dt <- out.dt[, .(CHR=as.integer(SNPChr), POS=as.integer(SNPPos), geneID=gsub("\\..*", "", Gene), pval=MetaP, beta=MetaBeta, se=MetaSE,
                       MAF=ifelse(SNPEffectAlleleFreq<0.5, as.numeric(SNPEffectAlleleFreq), 1-as.numeric(SNPEffectAlleleFreq)), EAF=SNPEffectAlleleFreq, 
                       EA=SNPEffectAllele, NEA=ifelse(SNPEffectAllele==gsub(".*_", "", SNP), gsub(".*:", "", gsub("_.*", "", SNP)), gsub(".*_", "", SNP)), 
                       ID=stringr::str_match(SNP, ".*:.*:(.*?):.*")[,2], indep.eQTL)]
  rm(dt)
  rm(indep.eQTL)
  
  #### Lift-down from b38 to b37 ####
  # Create .bed 
  out.dt <- out.dt[!is.na(POS)]
  out.dt$start <- out.dt$POS-1 
  data.table::fwrite(unique(out.dt[, .(CHR=as.integer(CHR), start=as.integer(start), POS=as.integer(POS), ALT, beta, ID)], by=c("CHR", "POS", "ID")), paste0("/project_data/data/eqtl_", tissue, ".hg38.bed"), col.names = F, row.names = F, quote = F, sep = "\t")
  data.table::fwrite(out.dt[, start:=NULL], paste0("/project_data/data/bulk/eqtl/eQTL_", tissue, "_hg38.csv"))
  rm(out.dt)
}

### ON WORKER: Run CrossMap to convert the positions on 
# for i in cerebellum-EUR-60 cortex-AFR-40 cortex-EAS-30 cortex-EUR-80 # basalganglia-EUR-30 hippocampus-EUR-30
# do
# CrossMap.py bed /lustre/groups/itg/teams/zeggini/projects/SCZ_T2D/data/hg38ToHg19.over.chain.gz /lustre/groups/itg/teams/zeggini/projects/SCZ_T2D/data/eqtl_${i}.hg38.bed /lustre/groups/itg/teams/zeggini/projects/SCZ_T2D/data/eqtl_${i}.hg19.bed 2>&1
# done

for (tissue in c("basalganglia-EUR-30", "hippocampus-EUR-30", "cerebellum-EUR-60", "cortex-AFR-40", "cortex-EAS-30", "cortex-EUR-80")){
  # tissue = "hippocampus-EUR-30"
  out.dt <- data.table::fread(paste0("/project_data/data/bulk/eqtl/eQTL_", tissue, "_hg38.csv"))
  # Import the lifted positions and keep only autosomes
  liftedover = data.table::fread(paste0("/project_data/data/eqtl_", tissue, ".hg19.bed"))
  colnames(liftedover) <- c("CHR", "start", "POS", "ALT", "beta", "ID")
  liftedover[, start:=NULL]
  dim(liftedover)
  # Do variants map to several positions ?
  nb.positions <- table(liftedover$variant_description)
  head(nb.positions)
  table(nb.positions)
  
  # Create the final file
  liftedover[, `:=` (beta=NULL, ALT=NULL)]
  liftedover <- unique(liftedover)
  liftedover <- liftedover[ID!="nors"]
  res <- merge(out.dt, liftedover, by="ID")
  rm(out.dt)
  rm(liftedover)
  
  res$CHR=res$CHR.y
  res$POS=res$POS.y
  res[, `:=` (CHR.x=NULL, CHR.y=NULL, POS.x=NULL, POS.y=NULL)]
  
  colnames(res) <- c("rsID", "geneID", "pval", "beta", "se", "MAF", "EAF", "EA", "NEA",
                     "indep.eQTL", "CHR", "POS", "ID")
  res[, CHR:=as.integer(CHR)]
  res[, POS:=as.integer(POS)]
  res[, rsID:=ID]
  res[, ID:=NULL]
  
  data.table::fwrite(res, paste0("/project_data/data/bulk/eqtl/eQTL_", tissue, "_hg19.csv"))
}



################################################################################
#------------------------ eQTL from brain PsychENCODE -------------------------#
################################################################################
dt <- data.table::fread("/project_data/data/bulk/eqtl/Full_hg19_cis-eQTL.txt.gz")
colnames(dt) <- c("gene_id", "gene_chr", "gene_start", "gene_end", "strand", "number_of_SNPs_tested", 
                  "SNP_distance_to_TSS", "SNP_id", "SNP_chr", "SNP_start", "SNP_end", "nominal_pval", 
                  "regression_slope", "top_SNP")

snps.dt <- data.table::fread("/project_data/data/bulk/eqtl/SNP_Information_Table_with_Alleles.txt")

out.dt <- merge(dt, snps.dt, by.x="SNP_id", by.y="PEC_id")
rm(dt)
rm(snps.dt)

out.dt <- out.dt[, .(CHR=strsplit(gene_chr, "r")[[1]][2], POS=SNP_start, rsID=Rsid, geneID=gene_id, pval=nominal_pval, beta=regression_slope, 
                     se=NA, MAF=NA, EAF=NA, EA=ALT, NEA=REF, ID=NA, indep.eQTL=top_SNP)]
out.dt <- out.dt[, se:=abs(regression_slope/qnorm(nominal_pval/2)), by=seq_len(nrow(out.dt))]

# add ID
out.dt <- out.dt[, ID:=paste(paste(CHR, POS, sep = ":"), sort_alleles(EA, NEA), sep = "_")]

data.table::fwrite(out.dt, "/project_data/data/bulk/eqtl/eQTL_brain.csv")
rm(out.dt)

################################################################################
#-------------------------- eQTL from liver from GTEx -------------------------#
################################################################################
dt <- data.table::fread("/reference_data/GTEx/GTEx_Analysis_v8_eQTL/Liver.allpairs.txt.gz")
sig.pairs <- data.table::fread("/reference_data/GTEx/GTEx_Analysis_v8_eQTL/Liver.v8.signif_variant_gene_pairs.txt.gz") %>% .[, .(variant_id, gene_id, sig.eQTL=1)]

out.dt <- merge(dt, sig.pairs, by=c("variant_id", "gene_id"), all.x=TRUE)
out.dt[is.na(sig.eQTL), sig.eQTL:=0]
out.dt[sig.eQTL==1, .N]
out.dt[sig.eQTL==0, .N]

out.dt <- out.dt[, .(CHR=as.integer(gsub("chr(\\d+|X+)_.*", "\\1", variant_id)), POS=as.integer(gsub(".*_(\\d+)_.*", "\\1", variant_id)), 
                     geneID=gene_id, pval=pval_nominal, beta=slope, se=slope_se, MAF=maf, EAF=NA, ALT=gsub("chr(\\d+|X+)_(\\d+)_([A-Z]*)_([A-Z]*)_b38", "\\4", variant_id), 
                     REF=gsub("chr(\\d+|X+)_(\\d+)_([A-Z]*)_([A-Z]*)_b38", "\\3", variant_id), ID=variant_id, sig.eQTL)]
rm(dt)
rm(sig.pairs)

# delete indels and X chromosome
out.dt <- out.dt[CHR!="X" & nchar(ALT)==1 & nchar(REF)==1] 

data.table::fwrite(out.dt, "/project_data/data/bulk/eqtl/eQTL_Liver.csv")
rm(out.dt)

#### Liftover ####
# Create .bed file
out.dt$start <- out.dt$POS-1
out.dt <- out.dt[, start:=as.integer(start)]
write.table(out.dt[, .(CHR, start, POS, ALT, beta, ID)], "/project_data/data/eqtl_liver.bed", col.names = F, row.names = F, quote = F, sep = "\t")

### ON WORKER: Run CrossMap to convert the positions on 
# CrossMap.py bed /lustre/groups/itg/teams/zeggini/projects/SCZ_T2D/data/hg38ToHg19.over.chain.gz /lustre/groups/itg/teams/zeggini/projects/SCZ_T2D/data/eqtl_liver.bed /lustre/groups/itg/teams/zeggini/projects/SCZ_T2D/data/eqtl_liver.hg19.bed 2>&1

# Import the lifted positions and keep only autosomes
liftedover = data.table::fread("/project_data/data/eqtl_liver.hg19.bed")
colnames(liftedover) <- c("CHR", "start", "POS", "ALT", "beta", "ID")
liftedover[, start:=NULL]
dim(liftedover)
#Remove the positions not mapping to chromosomes 1 to 22
liftedover <- subset(liftedover, CHR %in% 1:22)
dim(liftedover)
# Do variants map to several positions ?
nb.positions <- table(liftedover$variant_description)
head(nb.positions)
table(nb.positions)

# Create the final file
liftedover[, `:=` (beta=NULL, ALT=NULL)]
liftedover <- unique(liftedover)
res <- merge(out.dt, liftedover, by="ID")
rm(out.dt)
rm(liftedover)

res$CHR=res$CHR.y
res$POS=res$POS.y
res[, `:=` (CHR.x=NULL, CHR.y=NULL, POS.x=NULL, POS.y=NULL)]

# add ID (in AWK!)
# res <- res[, ID:=paste(paste(CHR, POS, sep = ":"), sort_alleles(EA, NEA), sep = "_")]

data.table::fwrite(res, "/project_data/data/bulk/eqtl/eQTL_Liver_hg19.csv")

################################################################################
#-------------------------- eQTL from fetal brain -----------------------------#
################################################################################
dt <- data.table::fread("/project_data/data/bulk/eqtl/fetal_brain_eqtl.txt.gz")
snps.dt <- data.table::fread("/project_data/data/bulk/eqtl/fetal_brain_snp_info.bed.gz") %>% .[, .(V1,V2,V4,V7,V8,V9)]
colnames(snps.dt) <- c("chr", "pos", "rsid", "REF", "ALT", "info")

# slope = from alternative allele

out.dt <- merge(dt, snps.dt, by.x="variant_id", by.y="rsid")
rm(dt)
rm(snps.dt)

out.dt <- out.dt[, `:=` (CHR=strsplit(chr, "r")[[1]][2], EAF=strsplit(strsplit(info, ";")[[1]][1], "=")[[1]][2]), by=seq_len(nrow(out.dt))]
out.dt <- out.dt[, .(CHR=as.integer(CHR), POS=pos, rsID=variant_id, geneID=gene_id, pval=pval_nominal, beta=slope, se=slope_se, MAF=maf, EAF=as.integer(EAF), 
                     EA=ALT, NEA=REF, ID=NA)]

data.table::fwrite(out.dt, "/project_data/data/bulk/eqtl/eQTL_fetal_brain_hg38.csv")
rm(out.dt)

#### Liftover ####
dt <- data.table::fread("/project_data/data/bulk/eqtl/old/eQTL_fetal_brain_hg38_id.csv")

# Create .bed file
bed.dt <- unique(dt[, .(CHR, POS, EA, beta, ID=newID)], by="ID")
head(bed.dt)
nrow(bed.dt)
bed.dt$start <- bed.dt$POS-1
bed.dt <- bed.dt[, start:=as.integer(start)]
write.table(bed.dt[, .(CHR, start, POS, EA, beta, ID)], "/project_data/data/eqtl_fetal_brain.bed", col.names = F, row.names = F, quote = F, sep = "\t")

### ON WORKER: Run CrossMap to convert the positions on 
# CrossMap.py bed /lustre/groups/itg/teams/zeggini/projects/SCZ_T2D/data/hg38ToHg19.over.chain.gz /lustre/groups/itg/teams/zeggini/projects/SCZ_T2D/data/eqtl_fetal_brain.bed /lustre/groups/itg/teams/zeggini/projects/SCZ_T2D/data/eqtl_fetal_brain.hg19.bed 2>&1

# Import the lifted positions and keep only autosomes
liftedover = read.table("/project_data/data/eqtl_fetal_brain.hg19.bed")
colnames(liftedover) <- c("CHR", "start", "POS", "ALT", "beta", "ID")
liftedover = liftedover[,-which(colnames(liftedover) == "start")]
dim(liftedover)
#Remove the positions not mapping to chromosomes 1 to 22
liftedover <- subset(liftedover, CHR %in% 1:22)
dim(liftedover)

# Do variants map to several positions ?
nb.positions <- table(liftedover$variant_description)
head(nb.positions)
table(nb.positions)

# Create the final file 
liftedover <- data.table::as.data.table(liftedover)
res <- merge(dt, liftedover[, .(CHR, POS, ID)], by.x="newID", by.y="ID")
rm(dt)
rm(liftedover)

res$CHR=res$CHR.y
res$POS=res$POS.y
res[, `:=` (CHR.x=NULL, CHR.y=NULL, POS.x=NULL, POS.y=NULL, newID=NULL)]

data.table::fwrite(res, "/project_data/data/bulk/eqtl/eQTL_fetal_brain.csv")


#------------------ Add indep.eQTL information ---------------------#
dt <- data.table::fread("/project_data/data/bulk/eqtl/eQTL_fetal_brain.csv")
top.eqtl <- data.table::as.data.table(readxl::read_xlsx("/project_data/data/bulk/eqtl/13059_2018_1567_MOESM1_ESM.xlsx", sheet=2, skip=2))

dt[, indep.eQTL:=ifelse(rsID %in% top.eqtl$rsID, 1, 0), by=geneID]

data.table::fwrite(dt, "/project_data/data/bulk/eqtl/eQTL_fetal_brain_indep_eQTL.csv")

################################################################################
#------------------------------ pQTL from brain -------------------------------#
################################################################################
dt <- data.table::fread("/project_data/data/bulk/pqtl/ROSMAP_DLPFC_pQTLs.csv")
out.dt <- dt[, .(CHR=CHR, POS=POS, rsID=NA, proteinID=UNIPROT, pval=P, adj.pval=FDR,
                 beta=BETA, se=SE, MAF=NA, EAF=NA, EA=ALT, NEA=REF, N=N, ID=NA)]
rm(dt)
data.table::fwrite(out.dt, "/project_data/data/bulk/pqtl/pQTL_brain.csv")
rm(out.dt)

# add ID
out.dt <- data.table::fread("/project_data/data/bulk/pqtl/pQTL_brain.csv")
out.dt <- out.dt[, ID:=paste(paste(CHR, POS, sep = ":"), sort_alleles(EA, NEA), sep = "_")]
data.table::fwrite(out.dt, "/project_data/data/bulk/pqtl/pQTL_brain.csv")
rm(out.dt)

################################################################################
#------------------------------ pQTL from liver -------------------------------#
################################################################################
dt <- as.data.table(readxl::read_xlsx("/project_data/data/bulk/pqtl/12915_2020_830_MOESM4_ESM.xlsx"))

# beta for Minor allele
out.dt <- dt[, .(CHR=as.integer(Chromosome), POS=as.integer(Position), rsID=Rs.number, geneID=Protein, pval=as.numeric(`P-value`), 
                 beta=as.numeric(Beta), se=as.numeric(NA), MAF=as.numeric(Minor.allele.frequency), EAF=as.numeric(Minor.allele.frequency), 
                 EA=Minor.allele, NEA=NA, indep.pQTL=ifelse(Independent.loci.marker=="Yes", 1, 0))]
out.dt <- out.dt[, se:=abs(beta/qnorm(pval/2)), by=seq_len(nrow(out.dt))]
out.dt <- out.dt[EA!="D" & !is.na(CHR)]
rm(dt)

data.table::fwrite(out.dt, "/project_data/data/bulk/pqtl/pQTL_liver_raw.csv")

# NEA missing!!!!
# https://support.bioconductor.org/p/9144015/
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(BSgenome.Hsapiens.UCSC.hg19)
my_rsids <- unique(out.dt$rsID)
gpos <- snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, my_rsids[!is.na(my_rsids)],ifnotfound="drop")
seqlevelsStyle(gpos) <- "UCSC"
z <- inferRefAndAltAlleles(gpos, BSgenome.Hsapiens.UCSC.hg19)
mcols(gpos) <- cbind(mcols(gpos), z)
gpos.dt <- as.data.table(gpos)
rm(gpos)
rm(z)

merged.dt <- merge(out.dt, gpos.dt, by.x="rsID", by.y="RefSNP_id") # , all.x=TRUE)
merged.dt <- merged.dt[nchar(alt_alleles)==1]
merged.dt <- merged.dt[, alt_allele:=unlist(alt_alleles)]
merged.dt <- merged.dt[, .(CHR, POS, rsID, geneID=proteinID, pval, beta, se, MAF, EAF, EA, NEA=ifelse(EA==ref_allele, alt_allele, ref_allele), indep.pQTL)]


merged.dt <- merged.dt[, ID:=paste(paste(CHR, POS, sep = ":"), sort_alleles(EA, NEA), sep = "_")]
rm(out.dt)

data.table::fwrite(merged.dt, "/project_data/data/bulk/pqtl/pQTL_liver.csv")

################################################################################
#----------------------- caQTL from brain PsychENCODE -------------------------#
################################################################################
dt <- data.table::fread("/project_data/data/bulk/caqtl/DER-09_hg19_cQTL.significant.txt")
dt <- dt[, `:=` (CHR=strsplit(Peak_chr, "r")[[1]][2], se=abs(regression_slope/qnorm(nominal_pval/2))), by=seq_len(nrow(dt))]
dt <- dt[, .(CHR=CHR, POS=SNP_start, SNP_id, rsID=NA, peakID=Peak_id, peakCenter=Peak_center, pval=nominal_pval, 
             beta=regression_slope, se=se, MAF=NA, EAF=NA, EA=NA, NEA=NA, ID=NA, indep.caQTL=top_SNP)]
snps.dt <- data.table::fread("/project_data/data/bulk/eqtl/SNP_Information_Table_with_Alleles.txt")

out.dt <- merge(dt, snps.dt, by.x="SNP_id", by.y="PEC_id")
rm(dt)
rm(snps.dt)

out.dt <- out.dt[, .(CHR, POS, rsID=Rsid, peakID, peakCenter, pval, beta, se, MAF=NA, EAF=NA, EA=ALT, NEA=REF, ID=NA, indep.caQTL)]

data.table::fwrite(out.dt, "/project_data/data/caQTL_brain_significant.csv")
rm(out.dt)

# add ID

################################################################################
#------------------------------ caQTL from liver ------------------------------#
################################################################################
dt <- data.table::fread("/project_data/data/bulk/caqtl/combined_RASQUAL-results-2PCs_100kbFromPeakCenters-MAF0.1_liver-secondPass-20samples.txt.gz")
dt <- dt[var_ID!="SKIPPED"]


# dt <- dt[, ID:=paste(var_ID, paste(ref_allele, alt_allele, sep="/"), sep="_"), by=seq_len(nrow(dt))]
# pval.dt <- data.table::fread("/project_data/data/bulk/caqtl/eigenMT-results_combined_RASQUAL-100kbFromPeakCenters_liver-secondPass-20samples.txt.gz")
# dt[var_ID=="1:730087" & peak=="peak11", "chi_square_statistic_(2xlogLikelihood_ratio)"]
# pval.dt[var=="1:730087_T/C" & gene=="peak11"]
# out.dt <- merge(dt, pval.dt, by.x=c("ID", "peak"), by.y=c("var", "gene"))
# rm(dt)
# rm(pval.dt)`

out.dt <- dt[, .(CHR=chr, POS=var_position_hg19, rsID=var_ID, peakID=peak, pval=pchisq(`chi_square_statistic_(2xlogLikelihood_ratio)`,df=1, lower.tail=FALSE), 
                 `effect_size_(Pi)`, `sequencing_mapping_error_rate_(Delta)`, 
                 `chi_square_statistic_(2xlogLikelihood_ratio)`, `reference_allele_mapping_bias_(Phi)`, overdispersion,
                 MAF=ifelse(alt_allele_frequency<=0.5, alt_allele_frequency, 1-alt_allele_frequency), EAF=alt_allele_frequency, EA=alt_allele, 
                 NEA=ref_allele, ID=NA)]

rm(dt)
data.table::fwrite(out.dt, "/project_data/data/bulk/caqtl/caQTL_liver_RASQUAL.csv")
rm(out.dt)

################################################################################
#----------------------------- caQTL from islets ------------------------------#
################################################################################
dt <- data.table::as.data.table(readxl::read_xlsx("/project_data/data/bulk/caqtl/db180393supplementarytable8.xlsx", skip=1))
out.dt <- dt[, .(CHR=Chromosome, POS=SNP_Position, rsID=rsID, peakID=Feature_Coordinates, pval=PValue, EffectSize, ChiSq_Statistic, Delta, Phi, Overdispersion, 
                 MAF=ifelse(AlelleFreq<=0.5, AlelleFreq, 1-AlelleFreq), EAF=AlelleFreq, EA=Alt, NEA=Ref, ID=NA)]
rm(dt)

data.table::fwrite(out.dt, "/project_data/data/bulk/caqtl/caQTL_islets_RASQUAL.csv")
rm(out.dt)

################################################################################
#--------------------------- single-cell from brain ---------------------------#
################################################################################
snp.pos <- data.table::fread("/project_data/data/single_cell/brain/snp_pos.txt.gz")

for (cell.type in c("Astrocytes", "Excitatory.neurons", "Inhibitory.neurons", "Microglia", "Oligodendrocytes", "OPCs...COPs", "Pericytes")){  # "Endothelial.cells", 
  dt <- data.table::fread(paste0("/project_data/data/single_cell/brain/", cell.type, ".gz"))
  colnames(dt) <- c("geneID", "rsID", "Distance_to_TSS", "pval", "beta")
  
  out.dt <- merge(dt, snp.pos, by.x="rsID", by.y="SNP")
  out.dt[, se:=abs(beta/qnorm(pval/2)), by=seq_len(nrow(out.dt))]
  
  out.dt <- out.dt[, .(CHR=as.integer(gsub("chr(\\d+|X+):.*", "\\1", SNP_id_hg19)), POS=as.integer(gsub(".*:(\\d+)", "\\1", SNP_id_hg19)), 
                       rsID, geneID=sub(".*_", "", geneID), gene_name=sub("_.*", "", geneID), pval, beta, se, MAF, EA=effect_allele, NEA=other_allele)]
  
  data.table::fwrite(out.dt, paste0("/project_data/data/single_cell/brain/sc_eqtl_", cell.type, ".csv"))
}

indep.qtl <- data.table::as.data.table(readxl::read_xlsx("/project_data/data/single_cell/brain/sc_eqtl_indep_info.xlsx", skip=3))
indep.qtl[cell_type=="OPCs / COPs", cell_type:="OPCs...COPs"]
indep.qtl[, cell_type:=sub(" ", "\\.", cell_type)]
for (c in c("Endothelial.cells", "Excitatory.neurons", "Inhibitory.neurons", "Microglia", "Oligodendrocytes", "OPCs...COPs", "Pericytes")){ # "Astrocytes", 
  dt <- data.table::fread(paste0("/project_data/data/single_cell/brain/sc_eqtl_", c, ".csv"))
  dt[, indep.eQTL:=ifelse(rsID %in% indep.qtl[cell_type==c, SNP], TRUE, FALSE)]
}

################################################################################
#-------------------------- mQTL from fetal brain -----------------------------#
################################################################################
dt <- data.table::fread("/project_data/data/bulk/mqtl/fetal/All_Imputed_BonfSignificant_mQTLs.csv.gz")

# t <- MungeSumstats::format_sumstats(path=dt, ref_genome="GRCh37")
out.dt <- dt[, .(CHR=SNP_Chr, POS=SNP_BP, rsID=NA, geneID=ProbeID, CpGpos=DNAm_BP, pval=p.value, beta,
                 MAF=NA, EAF=NA, EA=SNP_Allele, NEA=NA)]
out.dt[, se:=abs(beta/qnorm(pval/2)), by=seq_len(nrow(out.dt))]
rm(dt)

# add ID on shell

data.table::fwrite(out.dt, "/project_data/data/bulk/mqtl/mQTL_fetal.csv")
rm(out.dt)


################################################################################
#----------------------------- mQTL from brain --------------------------------#
################################################################################
dt <- data.table::fread("/project_data/data/bulk/mqtl/brain/mQTL_brain_50Kb.csv")
out.dt <- dt[, .(CHR=chr, POS=pos, rsID=rsID, CpG, CpGpos, pval=p, beta, se,
                 EAF=NA, EA=A1, NEA=A2)]
rm(dt)

snp.maf <- data.table::fread("/project_data/data/bulk/mqtl/brain/snp_maf.txt")
colnames(snp.maf) <- c("rsID", "MAF")
snp.maf[, rsID:=sub("chr", "", sub("\\.", ":", rsID))]
out.dt <- merge(out.dt, snp.maf, by="rsID", all.x=TRUE)
rm(snp.maf)

data.table::fwrite(out.dt, "/project_data/data/bulk/mqtl/mQTL_brain.csv")
rm(out.dt)