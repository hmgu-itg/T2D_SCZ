library(dplyr)
#####################################################################################
#---------------------------------- Functions --------------------------------------#
#####################################################################################
plot.hyprcoloc <- function(traits, tissue, qtl.type, gene, data, reg, res, plot.path, ld.path, genes.info, range=1e+6){
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
                     .(CHR=chr, SNP=snp.t1, P=pval, BP=pos, logP=-log10(as.numeric(pval)))]
  
  # Get LD matrix
  #ld.file <- NA
  ld.file <- data.table::fread(paste0(ld.path, "region", reg , "_LD.ld")) %>%
    .[SNP_A==sel.rsid | SNP_B==sel.rsid] %>% .[SNP_B==sel.rsid, c("SNP_B", "SNP_A") := .(SNP_A, SNP_B)]
  
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
             ld.file=ld.file,
             pp="PP4",
             pp.value=round(res$results$posterior_prob[1], digits=3),
             nplots=TRUE)
}

read.data <- function(data.path, gwas.coloc.result, chr.col="CHR", pos.col="POS", range=1e+6){
  dt <- data.table::fread(data.path)
  res.dt <- data.table::data.table()

  for(i in 1:nrow(gwas.coloc.result)){
    reg <- gwas.coloc.result[i]
    res.dt <- rbind(res.dt, dt[get(chr.col)==reg$CHR
                               & data.table::inrange(get(pos.col), reg$POS-range, reg$POS+range)])
  }
  return(unique(res.dt))
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

gwas.qtl.hyprcoloc <- function(traits, tissues.qtl, out.path, plot.path, genes.info, range=1e+6){
  gwas.coloc.result <- data.table::fread(paste0(output.path, "final_coloc_result_pp4.csv"))
  
  # Read GWAS data sets
  gwas <- lapply(traits, function(i){read.data(data.path=paste0(out.path, "GWAS_", i, "_precoloc.txt"), 
                                               gwas.coloc.result=gwas.coloc.result, chr.col="chr", 
                                               pos.col="pos")})
  names(gwas) <- traits
  
  res.dt <- data.table::data.table()
  credset.dt <- data.table::data.table()
  
  for (i in 1:nrow(tissues.qtl)){
    tissue <- tissues.qtl[i]$tissue
    qtl.type <- tissues.qtl[i]$qtl.type
    
    #for (tissue in tissues){
    # Read QTL data sets
    if(!file.exists(get(paste(tissue, qtl.type, "path", sep=".")))) {
      print(paste0("We do not have any data for ", qtl.type, " in ", tissue))
    }
    
    else{
      qtl <- read.data(data.path=get(paste(tissue, qtl.type, "path", sep=".")), gwas.coloc.result)
      qtl <- qtl[!is.na(se)]
      # colnames(qtl) <- c("rsID", "geneID", "pval", "beta", "se", "MAF", "EAF", "EA", "NEA",
      #                    "indep.eQTL", "CHR", "POS", "ID")
      # qtl[, CHR:=as.integer(CHR)]
      # qtl[, POS:=as.integer(POS)]
      # qtl[, SNP:=rsID]
      
      for (reg in gwas.coloc.result$region){
        sel.snp <- gwas.coloc.result[region==reg]$lead.SNP
        sel.chr <- gwas.coloc.result[region==reg]$CHR
        lead.pos <- gwas.coloc.result[region==reg]$POS
        
        data <- merge(gwas[[traits[1]]][chr==sel.chr & data.table::inrange(pos, lead.pos-range, lead.pos+range)], 
                      gwas[[traits[2]]][chr==sel.chr & data.table::inrange(pos, lead.pos-range, lead.pos+range)], 
                      by=c("ID", "chr", "pos"), allow.cartesian=T, suffixes=c(".t1", ".t2"))
        data[snp.t1=="", snp.t1:=snp.t2]
        data[snp.t2=="", snp.t2:=snp.t1]
        
        qtl.tmp <- qtl[CHR==sel.chr & data.table::inrange(POS, lead.pos-range, lead.pos+range)]
        for (gene in unique(qtl.tmp$geneID)){
          # print(gene)
          if(!("ID" %in% colnames(qtl)) || sum(is.na(qtl$ID))!=0) tmp.data <- merge(data, qtl.tmp[geneID==gene], by.x=c("chr", "pos"), by.y=c("CHR", "POS"))
          else tmp.data <- merge(data, qtl.tmp[geneID==gene], by.x=c("ID", "chr", "pos"), by.y=c("ID", "CHR", "POS"))
          
          if(nrow(tmp.data)==0) {
            print("No SNPs in common between GWAS and qtl data")
            next
          }
          
          # Flip alleles based on reference GWAS
          tmp.data[, `:=` (ea=ifelse(EA!=ea.t1, NEA, EA), 
                           nea=ifelse(EA!=ea.t1, EA, NEA), 
                           beta=ifelse(EA!=ea.t1, -beta, beta)), by=seq_len(nrow(tmp.data))]
          tmp.data[, `:=`(EA=NULL, NEA=NULL)]
          
          # ----------- Prepare beta and ses matrices -----------
          print("Preparing input matrices")
          ses <- tmp.data[geneID==gene, .(ID, rsID=snp.t1, se.t1, se.t2, se)]
          colnames(ses) <- c("ID", "rsID", traits[1], traits[2], paste(qtl.type, tissue, gene, sep="_"))
          betas <- tmp.data[geneID==gene, .(ID, rsID=snp.t1, beta.t1, beta.t2, beta)]
          colnames(betas) <- c("ID", "rsID", traits[1], traits[2], paste(qtl.type, tissue, gene, sep="_"))
          
          # ----------- Run HyPrColoc -----------
          print("Colocalization analysis using HyPrColoc")
          id <- betas$ID
          rsid <- betas$rsID
          betas_mat <- as.matrix(betas[, c('ID', 'rsID'):=NULL])
          rownames(betas_mat) <- id
          ses_mat <- as.matrix(ses[, c('ID', 'rsID'):=NULL])
          rownames(ses_mat) <- id
          binary.traits = c(rep(1,length(traits)), rep(0,ncol(betas_mat)-length(traits)))
          betas_mat <- na.omit(betas_mat)
          ses_mat <- na.omit(ses_mat)
          
          if (nrow(betas_mat)>1){
            res <- hyprcoloc::hyprcoloc(betas_mat, ses_mat, trait.names=colnames(betas_mat), snp.id=id, bb.alg=FALSE,
                                        binary.outcomes=binary.traits, prior.1=1e-4, prior.2=0.95, snpscores=TRUE)  #prior.c=0.05
            
            if(res$results$posterior_prob>0.8) {
              # get 95% credible set
              print(paste0("Plotting for ", gene))
              plot.hyprcoloc(traits, tissue, qtl.type, gene, unique(tmp.data[geneID==gene], by="snp.t1"), reg, res, 
                             plot.path, ld.path, genes.info)
              
              credset.dt <- rbind(credset.dt, get.credset(snpscores=res$snpscores, 
                                                          data=tmp.data[geneID==gene], 
                                                          reg=reg, 
                                                          gene=gene, 
                                                          tissue=tissue, 
                                                          qtl.type=qtl.type))
              
              res.dt <- rbind(res.dt, data.table::data.table(region=reg, 
                                                             geneID=gene, 
                                                             qtl.type=qtl.type, 
                                                             tissue=tissue, 
                                                             posterior_prob=res$results$posterior_prob,
                                                             candidate_snp=res$results$candidate_snp, 
                                                             rsID=tmp.data[ID==res$results$candidate_snp & geneID==gene]$snp.t1,
                                                             posterior_explained_by_snp=res$results$posterior_explained_by_snp,
                                                             # qtl.beta=qtl[ID==res$results$candidate_snp & geneID==gene]$beta, 
                                                             # qtl.se=qtl[ID==res$results$candidate_snp & geneID==gene]$se, 
                                                             # qtl.pval=qtl[ID==res$results$candidate_snp & geneID==gene]$pval))
                                                             qtl.beta=tmp.data[ID==res$results$candidate_snp & geneID==gene]$beta, 
                                                             qtl.se=tmp.data[ID==res$results$candidate_snp & geneID==gene]$se, 
                                                             qtl.pval=tmp.data[ID==res$results$candidate_snp & geneID==gene]$pval), fill=TRUE)
            }
          }
        }
      }
    }
  }
  return(list(res.dt, credset.dt))
}

#####################################################################################
#------------------------------------ Main -----------------------------------------#
#####################################################################################
project_folder="/project_data/"
source(paste0(project_folder, "scripts/read_config_file.R"))
source("/project_data/scripts/plot_functions.R")
GRCh37_Genes <- read.delim(genes.data, stringsAsFactors = FALSE, header = TRUE)

tissues.qtl <- data.table::fread("/project_data/GWASColoc/tissue_qtl.csv")
tissues.qtl <- tissues.qtl[source!="GTEx" & source!="no" & tissue!=""]

all.dt <- gwas.qtl.hyprcoloc(traits, tissues.qtl, output.path, qtl.coloc.plots.path, GRCh37_Genes)

result.dt <- all.dt[[1]][order(region)]
result.dt[, geneID:=gsub("(.*)\\..*", "\\1", gsub(".*_(.*)", "\\1", geneID))]
data.table::fwrite(result.dt, paste0(output.path, "qtl_hyprcoloc_results.csv"))
data.table::fwrite(all.dt[[2]], paste0(output.path, "qtl_hyprcoloc_credible_set.csv"))
