#####################################################################################
#---------------------------------- Functions --------------------------------------#
#####################################################################################
plot.hyprcoloc <- function(traits, tissue, qtl.type, gene, data, reg, res, plot.path, genes.info, ld.path=NULL, range=1e+6){
  sel.snp <- res$results$candidate_snp[1]
  sel.chr <- as.integer(gsub(":.*", "", sel.snp))
  sel.pos <- as.integer(gsub("_.*", "", gsub(".*:", "", sel.snp)))
  sel.rsid <- data[ID==sel.snp, snp.t1]
  
  # Get data
  trait1.region <- data[chr==sel.chr & dplyr::between(pos, sel.pos-range, sel.pos+range) & !is.null(snp.t1),
                        .(CHR=chr, SNP=snp.t1, P=pval.t1, BP=pos, logP=-log10(as.numeric(pval.t1)))]
  trait2.region <- data[chr==sel.chr & dplyr::between(pos, sel.pos-range, sel.pos+range),
                        .(CHR=chr, SNP=snp.t2, P=pval.t2, BP=pos, logP=-log10(as.numeric(pval.t2)))]
  qtl.region <- data[chr==sel.chr & dplyr::between(pos, sel.pos-range, sel.pos+range),
                     .(CHR=chr, SNP=snp.t1, P=pvalue, BP=pos, logP=-log10(as.numeric(pvalue)))]
  
  # Get LD matrix
  # ld.file <- data.table::fread(paste0(ld.path, "region", reg , "_LD.ld")) %>%
  #   .[SNP_A==sel.rsid | SNP_B==sel.rsid] %>% .[SNP_B==sel.rsid, c("SNP_B", "SNP_A") := .(SNP_A, SNP_B)]
  ld.file <- NA
  
  # Regional association plot
  data.lst <- list(trait1.region, trait2.region, qtl.region)
  names(data.lst) <-  c(traits[[1]], traits[[2]], paste(qtl.type, tissue, gene, sep="_"))
  locus.zoom(data = data.lst,
             offset_bp = range,
             genes.data = genes.info,
             file.name = paste0(plot.path, paste(traits[1], traits[2], reg, qtl.type, tissue, gene, sep="_"), ".png"),
             #secondary.snp = ifelse(credible.set<2, NA, credible.set),
             snp=sel.rsid, 
             ignore.lead=TRUE,
             #ld.file=ld.file,
             pp="PP4",
             pp.value=round(res$results$posterior_prob[1], digits=3),
             nplots=TRUE)
}

request_datasets_from_api <- function(study_id = "",
                                      quant_method = "",
                                      sample_group = "",
                                      tissue_id = "",
                                      study_label = "",
                                      tissue_label = "",
                                      condition_label = "") {
  size = 1000 #Page size
  start = 0 #Page start
  
  parameter_values = c(study_id,quant_method,sample_group,tissue_id,study_label, 
                       tissue_label,condition_label)
  parameter_names = c('study_id','quant_method','sample_group','tissue_id',
                      'study_label','tissue_label','condition_label')
  
  while (T) {
    URL = glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/?size={size}&start={start}")
    
    #Adding defined parameters to the request
    for (i in 1:length(parameter_values)) {
      par = parameter_values[i]
      par_name = parameter_names[i]
      if (par != "")
        URL = glue("{URL}&{par_name}={par}")
    }
    
    r <- GET(URL, accept_json())
    cont <- content(r, "text", encoding = "UTF-8")
    
    # If the request was unsuccessful
    if (status_code(r) != 200) {
      #If we get no results at all, print error
      if (start == 0) {
        print(glue("Error {status_code(r)}"))
        print(cont)
        return ()
      }
      #else just break
      break
    }
    
    cont_df <- fromJSON(cont)
    
    if (start == 0) {
      responses <- cont_df
    }
    else{
      responses <- rbind(responses, cont_df)
    }
    start <- start + size
  }
  return(responses)
}

request_associations_around_position <- function(dataset_id, position, chromosome_id, offset = 500000){
  size = 1000
  start = 0
  range_start = position - offset
  range_end = position + offset
  
  
  while (TRUE){
    URL = glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/{dataset_id}/associations?size={size}&start={start}&pos={chromosome_id}:{range_start}-{range_end}")
    
    r <- GET(URL, accept_json())
    cont <- content(r, "text", encoding = "UTF-8")
    
    if (status_code(r) != 200) {
      # Loop will break if the request was unsuccessful
      if(start==0) {
        print(glue("Error {status_code(r)}"))
        print(cont)
        return()}
      break
    }
    
    
    cont_df <- fromJSON(cont)
    
    if (start == 0){
      responses <- cont_df
    }
    else{
      responses <- rbind(responses, cont_df)
    }
    start <- start + size
  }
  return(responses)
}

get.credset <- function(snpscores, data, reg, gene, tissue, qtl.type, value=0.95){
  sorted <- sort(snpscores, decreasing=TRUE)
  cs <- cumsum(sorted)
  w <- which(cs > value)[1]
  credset <- data.table::data.table(region=reg, 
                                    geneID=gene, 
                                    qtl.type=qtl.type, 
                                    tissue=tissue,
                                    ID=names(sorted)[1:w], 
                                    SNP.PP4=sorted[1:w],
                                    rsID=data[ID %in% names(sorted)[1:w], snp.t1],
                                    chr=data[ID %in% names(sorted)[1:w], chr],
                                    pos=data[ID %in% names(sorted)[1:w], pos])
  return(credset)
}

library(httr)
library(jsonlite)
library(glue)
library(rtracklayer)
library(dplyr)
library(hyprcoloc)

setwd("C:/Users/ana.arruda/OneDrive - Helmholtz Zentrum München/Projects/SCZ_T2D/")

# source("scripts/read_files_config.R")
source("scripts_bckup/plot_functions.R")
GRCh37_Genes <- read.delim("data/UCSC_GRCh37_Genes_UniqueList.txt", stringsAsFactors = FALSE, header = TRUE)

t2d <- unique(data.table::fread("data/T2D_gtex_coloc.txt"))
scz <- unique(data.table::fread("data/SCZ_gtex_coloc.txt"))

tissues.qtl <- data.table::fread("tissue_qtl.csv")
tissues.gtex <- tissues.qtl[source=="GTEx"]

# Get dataset ID from eQTL catalogue
tissues.gtex[qtl.type=="eqtl", dataset.id:=request_datasets_from_api(quant_method="ge", sample_group=tolower(tissue))$dataset_id, by=tissue]
# tissues.gtex[qtl.type=="sqtl", dataset.id:=request_datasets_from_api(quant_method="txrev", sample_group=tolower(tissue))$dataset_id, by=tissue]
tissues.gtex[qtl.type=="sqtl", dataset.id:=request_datasets_from_api(quant_method="leafcutter", sample_group=tolower(tissue))$dataset_id, by=tissue]
tissues.gtex[qtl.type=="tqtl", dataset.id:=request_datasets_from_api(quant_method="txrev", sample_group=tolower(tissue))$dataset_id, by=tissue]

# Define the regions of interest
coloc.regions <- data.table::fread("gwas_coloc/final_coloc_result_pp4.csv")

# Lift over regional data from build 37 to build 38
print("lifting over genomic risk coordinates for coloc region")
system("gzip -d data/hg19ToHg38.over.chain.gz")
ch <- import.chain("data/hg19ToHg38.over.chain")
for (i in 1:nrow(coloc.regions)){
  regions.b37 <- GRanges(coloc.regions[i, .(CHR, start=POS, end=POS)])
  seqlevelsStyle(regions.b37) = "UCSC" # necessary
  # regions.b38 <- data.table::as.data.table(unlist(liftOver(regions.b37, ch)))
  coloc.regions[i, POS.b38:=data.table::as.data.table(unlist(liftOver(regions.b37, ch)))$start]
}

# mapping.file <- data.table::fread("data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz")

# Perform colocalization, filtering, or any other analysis as needed
res.dt <- data.table::data.table()
credset.dt <- data.table::data.table()

for (i in 15:nrow(tissues.gtex)){
  t <- tissues.gtex[i]
  
  for (j in 1:nrow(coloc.regions)){
    reg <- coloc.regions[j]
    
    # Get data
    data <- merge(t2d[chr==reg$CHR & data.table::inrange(pos, reg$POS-1e6, reg$POS+1e6)], 
                  scz[chr==reg$CHR & data.table::inrange(pos, reg$POS-1e6, reg$POS+1e6)], 
                  by=c("ID", "chr", "pos"), allow.cartesian=T, suffixes=c(".t1", ".t2"))
    data[snp.t1=="", snp.t1:=snp.t2]
    data[snp.t2=="", snp.t2:=snp.t1]
  
    # Fetch the tissue-specific eQTL data from the API for all genes in the region (b38)
    qtl.data <- rbind(request_associations_around_position(t$dataset.id, reg$POS.b38-500000, reg$CHR),
                      request_associations_around_position(t$dataset.id, reg$POS.b38+500000, reg$CHR))
    qtl.data <- unique(data.table::as.data.table(qtl.data))
    
    if(nrow(qtl.data)==0) {
      print("No SNPs in qtl data")
      next
    }
    
    qtl.data <- qtl.data[, chromosome:=as.integer(chromosome)]
    
    # Lift down data from build 38 to build 37
    print("Lifting down genomic risk coordinates of qtl data")
    path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
    ch = import.chain(path)
    qtl.data.b38 <- GRanges(qtl.data[, `:=` (start=position, end=position)])
    seqlevelsStyle(qtl.data.b38) = "UCSC" # necessary
    qtl.data <- data.table::as.data.table(unlist(liftOver(qtl.data.b38, ch)))
    qtl.data[, `:=`(POS.b37=start, POS.b38=position, chromosome=as.integer(sub("chr", "", seqnames)))]
    qtl.data[, `:=`(start=NULL, end=NULL, seqnames=NULL, width=NULL, strand=NULL, position=NULL)]
    qtl.data <- qtl.data[nchar(ref)==1 & nchar(alt)==1]
   
    # for (gene in unique(data$gene_id)){
    for (gene in unique(qtl.data$gene_id)){
      tmp.data <- merge(data, qtl.data[gene_id==gene], by.x=c("chr", "pos"), by.y=c("chromosome", "POS.b37"))
      
      if(nrow(tmp.data)==0) {
        print("No SNPs in common between GWAS and qtl data")
        next
      }
      
      # Flip alleles based on reference GWAS
      tmp.data[, `:=` (ea=ifelse(alt==ea.t1, alt, ref), 
                       nea=ifelse(alt==ea.t1, ref, alt), 
                       beta=ifelse(alt==ea.t1, beta, -beta)), by=seq_len(nrow(tmp.data))]
      tmp.data[, `:=`(alt=NULL, ref=NULL)]
      
      # ----------- Prepare beta and ses matrices -----------
      print("Preparing input matrices")
      ses <- tmp.data[gene_id==gene, .(ID, rsID=snp.t1, se.t1, se.t2, se)]
      colnames(ses) <- c("ID", "rsID", "t2d", "scz", paste(t$qtl.type, t$tissue, gene, sep="_"))
      betas <- tmp.data[gene_id==gene, .(ID, rsID=snp.t1, beta.t1, beta.t2, beta)]
      colnames(betas) <- c("ID", "rsID", "t2d", "scz", paste(t$qtl.type, t$tissue, gene, sep="_"))
      
      # ----------- Run HyPrColoc -----------
      print("Colocalization analysis using HyPrColoc")
      id <- betas$ID
      rsid <- betas$rsID
      betas_mat <- as.matrix(betas[, c('ID', 'rsID'):=NULL])
      rownames(betas_mat) <- id
      ses_mat <- as.matrix(ses[, c('ID', 'rsID'):=NULL])
      rownames(ses_mat) <- id
      binary.traits = c(1,1,0)
      betas_mat <- na.omit(betas_mat)
      ses_mat <- na.omit(ses_mat)
      
      if (nrow(betas_mat)>1){
        res <- hyprcoloc::hyprcoloc(betas_mat, ses_mat, trait.names=colnames(betas_mat), snp.id=id, bb.alg=FALSE,
                                    binary.outcomes=binary.traits, prior.1=1e-4, prior.c=0.05, snpscores = TRUE)
        
        
        if(res$results$posterior_prob>0.8) {
          print(paste("Plotting for", t$qtl.type, "in", t$tissue, "and", gene, sep=" "))
          plot.hyprcoloc(traits=c("t2d", "scz"), t$tissue, t$qtl.type, gene, unique(tmp.data[gene_id==gene], by="snp.t1"), 
                         reg$region, res, "gwas_coloc/qtl_plots/", genes.info=GRCh37_Genes)
          
          # get 95% credible set
          credset.dt <- rbind(credset.dt, get.credset(snpscores=res$snpscores, 
                                                      data=tmp.data[gene_id==gene], 
                                                      reg=reg$region, 
                                                      gene=gene, 
                                                      tissue=t$tissue, 
                                                      qtl.type=t$qtl.type))
                
          res.dt <- rbind(res.dt, data.table::data.table(region=reg$region, 
                                                         geneID=gene, 
                                                         qtl.type=t$qtl.type, 
                                                         tissue=t$tissue, 
                                                         posterior_prob=res$results$posterior_prob,
                                                         candidate_snp=res$results$candidate_snp, 
                                                         rsID=tmp.data[ID==res$results$candidate_snp & gene_id==gene]$snp.t1,
                                                         posterior_explained_by_snp=res$results$posterior_explained_by_snp,
                                                         qtl.beta=tmp.data[ID==res$results$candidate_snp & gene_id==gene]$beta, 
                                                         qtl.se=tmp.data[ID==res$results$candidate_snp & gene_id==gene]$se, 
                                                         qtl.pval=tmp.data[ID==res$results$candidate_snp & gene_id==gene]$pval), fill=TRUE)
        }
      }
    }
  }
}

res.dt <- res.dt[order(region)]
res.dt[, geneID:=gsub("(.*)\\..*", "\\1", gsub(".*_(.*)", "\\1", geneID))]
data.table::fwrite(res.dt, "gwas_coloc/qtl_gtex_hyprcoloc_results.csv")
data.table::fwrite(credset.dt, "gwas_coloc/qtl_gtex_hyprcoloc_credible_set.csv")




