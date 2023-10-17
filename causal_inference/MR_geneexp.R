################################################
#------------------- Functions -----------------
################################################
library(dplyr)

sort_alleles <- Vectorize(function(x,y) {
  paste(sort(c(x, y)), collapse = "_")
}) 

local.clump <- function(data, pval=5e-08){
  # https://mrcieu.github.io/ieugwasr/articles/local_ld.html
  clumped.dt <- ieugwasr::ld_clump(dplyr::tibble(rsid=data$SNP, pval=data$pval.exposure), #, id=data$id.exposure),
                         plink_bin=genetics.binaRies::get_plink_binary(),
                         bfile="/reference_data/1kG/EUR_1kg_v3/EUR",
                         clump_p=pval)
  return(data[data$SNP %in% clumped.dt$rsid, ])
}

proxy.lookup <- function(snp, api.token, population="EUR"){
  my_proxies <- LDlinkR::LDproxy(snp = snp, 
                        pop = population, 
                        r2d = "r2", 
                        token = api.token,
                        genome_build = "grch37")
  return(my_proxies$RS_Number)
}

get.filename <- function(trait, geneexp.exp=TRUE, hc.genes.id=NA) {  #, genes=NA){
  if(trait %in% c("T2D_european", "SCZ_european", "T2D_multiancestry", "SCZ_multiancestry")){
    dt <- data.table::fread(paste0("/project_data/data/", trait, ".txt"))
    return(dt)
  }
  else {
    dt <- data.table::fread(get(paste(trait, "path", sep="."))) %>% .[geneID %in% hc.genes.id]
    # if (geneexp.exp) {
    #   # dt <- tryCatch(
    #   #   dt[indep.eQTL == TRUE],
    #   #   error=function(e) e
    #   # )
    #   dt <- dt[indep.eQTL == TRUE | indep.eQTL==1]
    # }
    # if (!is.na(genes)) dt <- dt[geneID %in% genes]
    return(dt)
  }
  print("ERROR: Trait not known")
}

get.data <- function(file, data.type, data.dir){
  if(data.type=="quant"){
    data <- TwoSampleMR::format_data(type = data.dir,
                                     dat = file,
                                     snp_col = "rsID",
                                     beta_col = "beta",
                                     se_col = "se",
                                     effect_allele_col = "EA",
                                     other_allele_col = "NEA",
                                     eaf_col = "EAF",
                                     pval_col = "pval",
                                     chr_col="CHR",
                                     pos_col="POS",
                                     id_col="ID",
                                     samplesize_col="N",
                                     gene_col = "geneID",
                                     info_col="indep.eQTL")
  }
  else if(data.type=="binary"){
    data <- TwoSampleMR::format_data(type = data.dir,
                                     dat = file,
                                     snp_col = "rsID",
                                     beta_col = "beta",
                                     se_col = "se",
                                     effect_allele_col = "EA",
                                     other_allele_col = "NEA",
                                     eaf_col = "EAF",
                                     pval_col = "pval",
                                     chr_col="CHR",
                                     pos_col="POS",
                                     id_col="ID",
                                     ncase_col="Ncases",
                                     samplesize_col="N")
  }
  else {
    print("Data type not recognized (options: quant or binary).")
  }
  
  return(data)
}

run.TwoSampleMR <- function(data, filename, folder="geneexp"){
  if(nrow(data)==0) {
    print("No SNP in common between exposure and outcome :(")
    tmp.dt <- data.table::data.table()
  }
  else if (nrow(data)==1) {
    print("Only one SNP in common between exposure and outcome, running Wald ratio method")
    res <- TwoSampleMR::mr(data, method_list="mr_wald_ratio")
    res <- TwoSampleMR::generate_odds_ratios(res)
    
    fstat <- (data$beta.exposure)^2/(data$se.exposure)^2
    snp <- data$SNP
    tmp.dt <- data.table::as.data.table(res) 
    tmp.dt <- tmp.dt[, `:=`(SNP=snp,Fstat=fstat,Q=as.numeric(NA),Q_df=as.numeric(NA),Q_pval=as.numeric(NA),plt.egger_intercept=as.numeric(NA),plt.pval=as.numeric(NA),plt.se=as.numeric(NA))]
  }
  else {
    print("Multiple SNPs in common between exposure and outcome, running MR")
    # Assumption 1: check strength of IVs with F-statistics nd remove the ones with Fstat<10
    data$fstat_per_snp <- (data$beta.exposure)^2/(data$se.exposure)^2
    data <- data[data$fstat_per_snp>10,]
    
    if(nrow(data)==0)  {
      print("No SNP in common between exposure and outcome :(")
      tmp.dt <- data.table::data.table()
      return(tmp.dt)
    }
    else if (nrow(data)==1) {
      print("Only one SNP in common between exposure and outcome, running Wald ratio method")
      res <- TwoSampleMR::mr(data, method_list="mr_wald_ratio")
      res <- TwoSampleMR::generate_odds_ratios(res)
      
      fstat <- (data$beta.exposure)^2/(data$se.exposure)^2
      tmp.dt <- data.table::as.data.table(res) %>% .[, `:=`(Fstat=fstat,Q=as.numeric(NA),Q_df=as.numeric(NA),Q_pval=as.numeric(NA),plt.egger_intercept=as.numeric(NA),plt.pval=as.numeric(NA),plt.se=as.numeric(NA))]
    }
    else {
      fstat_overall <- mean((data$beta.exposure)^2/(data$se.exposure)^2)
      
      res <- TwoSampleMR::mr(data, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
      res <- TwoSampleMR::generate_odds_ratios(res)
      
      if(nrow(res)>0){
        #------------ Sensitivity analysis
        TwoSampleMR::mr_scatter_plot(res, data)   #scatter plot
        ggplot2::ggsave(paste0("/project_data/CausalInference/", folder, "/TwoSampleMR_", filename, ".jpg"))
        
        het <- TwoSampleMR::mr_heterogeneity(data)   #heterogeneity statistics
        plt <- TwoSampleMR::mr_pleiotropy_test(data)  #MR Egger intercept for directional pleiotropy
        tmp.dt <- data.table::as.data.table(merge(res, het, by=c("id.exposure", "id.outcome", "outcome", "exposure", "method"), all.x=T))
        tmp.dt <- tmp.dt[, `:=`(plt.egger_intercept=as.numeric(NA),plt.pval=as.numeric(NA),plt.se=as.numeric(NA))]
        tmp.dt <- tmp.dt[method=="MR Egger", `:=`(plt.egger_intercept=plt$egger_intercept,plt.pval=plt$pval,plt.se=plt$se)]
        
        TwoSampleMR::mr_forest_plot(TwoSampleMR::mr_leaveoneout(data))
        ggplot2::ggsave(paste0("/project_data/CausalInference/", folder, "/LeaveOneOut_", filename, ".jpg"))
        TwoSampleMR::mr_forest_plot(TwoSampleMR::mr_singlesnp(data))
        ggplot2::ggsave(paste0("/project_data/CausalInference/", folder, "/SingleSNP_", filename, ".jpg"))
      }
      else{
        tmp.dt <- data.table::data.table()
      }
    }
  }
  return(tmp.dt)
}

output.wide.result <- function(dt){
  dt[, `:=` (id.exposure=NULL, id.outcome=NULL)]
  dt <- data.table::dcast(dt, exposure + outcome + nsnp + Fstat ~ method, 
                          value.var=c("b", "se", "pval", "lo_ci", "up_ci", "or", "or_lci95", "or_uci95", "Q", "Q_df", "Q_pval", "plt.egger_intercept", "plt.pval", "plt.se", "p.adj.fdr", "SNP")) 

  # remove plt from IVW
  dt[ ,`:=`(`plt.egger_intercept_Inverse variance weighted`=NULL, `plt.pval_Inverse variance weighted`=NULL, `plt.se_Inverse variance weighted`=NULL)]
  # remove Q and plt from WM
  dt[ ,`:=`(`plt.egger_intercept_Weighted median`=NULL, `plt.pval_Weighted median`=NULL, `plt.se_Weighted median`=NULL,
            `Q_Weighted median`=NULL, `Q_df_Weighted median`=NULL, `Q_pval_Weighted median`=NULL)]
  # remove Q and plt from Wald ratio
  dt[ ,`:=`(`plt.egger_intercept_Wald ratio`=NULL, `plt.pval_Wald ratio`=NULL, `plt.se_Wald ratio`=NULL,
            `Q_Wald ratio`=NULL, `Q_df_Wald ratio`=NULL, `Q_pval_Wald ratio`=NULL)]
  
  data.table::setcolorder(dt, c(colnames(dt)[1:3], 
                                stringr::str_subset(colnames(dt), "Inverse variance weighted"),
                                stringr::str_subset(colnames(dt), "Wald ratio"),
                                stringr::str_subset(colnames(dt), "Weighted median"),
                                stringr::str_subset(colnames(dt), "MR Egger")))
  return(dt)
}

wrap.MR <- function(exposure.lst, outcome.lst, exp.phen, out.phen, folder, output.file, dt, genes, p.thres, geneexp.exp, api.token, dt.steiger=NULL, ivs=NULL){
  if (geneexp.exp==TRUE){
    for (out in outcome.lst) {
      #----------- Get outcome -----------
      outcome <- get.data(get.filename(trait=out[1]), data.dir="outcome", data.type=out[2])
      outcome$outcome <- out[1]
      
      for (exp in exposure.lst){
        #----------- Read and clump exposure -----------
        exposure.raw <- get.filename(trait=exp[1], geneexp.exp=geneexp.exp, hc.genes.id=genes$id)  
        
        if(nrow(exposure.raw[indep.eQTL==TRUE | indep.eQTL==1])==0){
          print(paste0("No eQTL for any HC gene in ", exp[1]))
          next
        }
  
        for (cur.gene in genes$id){
          cur.gene.name <- genes[id==cur.gene, name]
          
          if (nrow(exposure.raw[geneID==cur.gene & (indep.eQTL==TRUE | indep.eQTL==1)])==0){
            print(paste0("No eQTL for gene ", cur.gene.name, " in ", exp[1]))
            next
          }
          
          exposure <- get.data(exposure.raw[geneID==cur.gene & (indep.eQTL==TRUE | indep.eQTL==1)], data.dir="exposure", data.type=exp[2])
          exposure$exposure <- paste(exp[1], cur.gene.name, sep="_")
          
          if(nrow(exposure)==0) {
            print(paste0("No eQTL for gene ", cur.gene.name, " in ", exp[1]))
            next
          }
          
          if(nrow(exposure)>1){
            if(!is.null(ivs)) clump.exposure <- exposure[exposure$SNP %in% ivs,]
            else clump.exposure <- local.clump(data=exposure, pval=p.thres)
          }
          else clump.exposure <- exposure
    
          #----------- Harmonize data -----------
          if(nrow(outcome[outcome$SNP %in% clump.exposure$SNP,])==0) {
            print("No valid SNP in common between exposure and outcome, looking for proxies")
            # Find proxies
            if(exp[1]=="cortex_eas.eqtl") my.pop="EAS"
            else if(exp[1]=="cortex_afr.eqtl") my.pop="AFR"
            else my.pop="EUR"
            
            proxy.snps <- proxy.lookup(clump.exposure$SNP, api.token, pop=my.pop)
            
            for(i in 2:length(proxy.snps)){
              if(nrow(outcome[outcome$SNP %in% proxy.snps[i],])>0 & nrow(exposure.raw[geneID==cur.gene & rsID==proxy.snps[i]])>0){
            #     print(paste0("i is ", i))
            #   }
            # }
                clump.exposure <- get.data(exposure.raw[geneID==cur.gene & rsID==proxy.snps[i]], data.dir="exposure", data.type=exp[2])
                clump.exposure$exposure <- paste(exp[1], cur.gene.name, sep="_")
                break
              }
            }
          }
          
          data <- tryCatch(
            TwoSampleMR::harmonise_data(exposure_dat=clump.exposure, outcome_dat=outcome),
            error=function(e) e
          )
          
          if(inherits(data, "error")) {
            print("Harmonization failed.")
            next 
          }
        
          if(nrow(subset(data, mr_keep==TRUE))==0) {
            print("No SNP left after data harmonization, not able to run MR :(")
            next
          }
          
          data <- subset(data, mr_keep==TRUE)
         
          #----------- Run TwoSampleMR for unfiltered data -----------
          dt <- rbind(dt, run.TwoSampleMR(data, filename=paste(exp[1], cur.gene.name, out[1], sep="_"), folder=folder), fill=TRUE)
          
          if(nrow(data)>1 & !is.null(dt.steiger)) {
            #----------- Apply steiger filtering on data -----------
            filtered.data <- TwoSampleMR::steiger_filtering(data)
            filtered.data <- filtered.data[filtered.data$steiger_dir==TRUE,]  # & filtered.data$steiger_pval<0.05,]
            
            #----------- Run TwoSampleMR for Steiger-filtered data -----------
            dt.steiger <- rbind(dt.steiger, run.TwoSampleMR(filtered.data, filename=paste(exp[1], cur.gene.name, out[1], "Steiger", sep="_"), folder=folder), fill=TRUE)
          }
        }
      }
    }
  }
  else{
    for (exp in exposure.lst) {
      #----------- Read and clump exposure -----------
      exposure <- get.data(get.filename(trait=exp[1]), data.dir="exposure", data.type=exp[2])
      exposure$exposure <- exp[1]
      if(!is.null(ivs)) clump.exposure <- exposure[exposure$SNP %in% ivs,]
      else clump.exposure <- local.clump(data=exposure, pval=p.thres)
      
      for (out in outcome.lst){
        #----------- Get outcome -----------
        outcome.raw <- get.filename(trait=out[1], geneexp.exp=FALSE, hc.genes.id=genes$id)
        
        for (cur.gene in genes$id){
          cur.gene.name <- genes[id==cur.gene, name]
          
          if (nrow(outcome.raw[geneID==cur.gene])==0){
            print("No SNP for this gene in the qtl data.")
            next
          }
          
          outcome <- get.data(outcome.raw[geneID==cur.gene], data.dir="outcome", data.type=out[2])
          outcome$outcome <- paste(out[1], cur.gene.name, sep="_")
          
          if (nrow(outcome[outcome$id.outcome %in% clump.exposure$id.exposure,])==0) {
            print("No valid SNP in common between exposure and outcome, not able to run MR :(")
            next
          }
          
          #----------- Harmonize data -----------
          data <- tryCatch(
            TwoSampleMR::harmonise_data(exposure_dat=clump.exposure, outcome_dat=outcome),
            error=function(e) e
          )
          
          if(inherits(data, "error")) {
            print("Harmonization failed.")
            next 
          }
          
          if(nrow(subset(data, mr_keep==TRUE))==0) {
            print("No SNP left after data harmonization, not able to run MR :(")
            next
          }
          
          data <- subset(data, mr_keep==TRUE)
          
          #----------- Run TwoSampleMR for unfiltered data -----------
          dt <- rbind(dt, run.TwoSampleMR(data, filename=paste(exp[1], cur.gene.name, out[1], sep="_"), folder=folder), fill=TRUE)
          
          
          if(nrow(data)>1 & !is.null(dt.steiger)) {
            #----------- Apply steiger filtering on data -----------
            filtered.data <- TwoSampleMR::steiger_filtering(data)
            filtered.data <- filtered.data[filtered.data$steiger_dir==TRUE,]  # & filtered.data$steiger_pval<0.05,]
            
            #----------- Run TwoSampleMR for Steiger-filtered data -----------
            dt.steiger <- rbind(dt.steiger,run.TwoSampleMR(filtered.data, filename=paste(exp[1], cur.gene.name, out[1], "Steiger", sep="_"), folder=folder), fill=TRUE)
          }
        }
      }
    }
  }
  
  if(nrow(dt)>0){
    dt[, p.adj.fdr:=p.adjust(pval, method = "fdr"), by=method]
    data.table::fwrite(dt, paste0("/project_data/CausalInference/", folder, "/", output.file, "_analysis_MR.csv"))
    data.table::fwrite(output.wide.result(dt), paste0("/project_data/CausalInference/", folder, "/", output.file, "_analysis_MR_wide.csv"))
  }

  if(!is.null(dt.steiger)){
    if(nrow(dt.steiger)>0){
      dt.steiger[, p.adj.fdr:=p.adjust(pval, method = "fdr"), by=method]
      data.table::fwrite(dt.steiger, paste0("/project_data/CausalInference/", folder, "/Steiger_", output.file, "_analysis_MR.csv"))
      data.table::fwrite(output.wide.result(dt.steiger), paste0("/project_data/CausalInference/", folder, "/Steiger_", output.file, "_analysis_MR_wide.csv"))
    }
  }

  return(dt)
}

project_folder="/project_data/"
source("/project_data/scripts/read_files_config.R")
hc.genes <- data.table::data.table(name=c(read.table("/project_data/hc_symbol.txt"))$V1, 
                                   id=c(read.table("/project_data/hc_ensembl.txt"))$V1)

tissues.qtl <- data.table::fread("/project_data/GWASColoc/tissue_qtl.csv")
tissues.qtl <- tissues.qtl[source!="GTEx" & source!="METSIM"]
tissues.qtl <- tissues.qtl[-c(8,9,23)]
tissues.lst <- lapply(tissues.qtl[, paste(tissue, qtl.type, sep="."), by=seq_len(nrow(tissues.qtl))]$V1, function(x) c(x, "quant"))

api.token <- ""

###########################################################################
#--------------- Exposure: molecular qtl, outcome: disease ---------------#
###########################################################################
dt <- wrap.MR(exposure.lst=tissues.lst,  # [-c(2,3)] 
              outcome.lst=list(c("SCZ_multiancestry", "binary")), # c("T2D_multiancestry", "binary")),
              folder="geneexp_new",
              dt=data.table::data.table(),
              dt.steiger=NULL, #data.table::data.table(),
              output.file="scz_geneexp_disease_multiancestry",
              p.thres=1,
              genes=hc.genes,
              geneexp.exp=TRUE,
              api.token=api.token)

dt <- wrap.MR(exposure.lst=tissues.lst[-c(2,3)], 
              outcome.lst=list(c("SCZ_european", "binary"), c("T2D_european", "binary")),
              folder="geneexp_new",
              dt=data.table::data.table(),
              dt.steiger=NULL, #data.table::data.table(),
              output.file="geneexp_disease_european",
              p.thres=1,
              genes=hc.genes,
              geneexp.exp=TRUE,
              api.token=api.token)

# adipose tissue METSIM
dt <- wrap.MR(exposure.lst=list(c("adipose.eqtl", "quant")),
              outcome.lst=list(c("SCZ_multiancestry", "binary"), c("T2D_multiancestry", "binary")),
              folder="geneexp_new",
              dt=data.table::data.table(),
              dt.steiger=NULL, #data.table::data.table(),
              output.file="adipose_geneexp_disease_multiancestry",
              p.thres=5e-5,
              genes=hc.genes,
              geneexp.exp=TRUE,
              api.token=api.token)

dt <- wrap.MR(exposure.lst=list(c("adipose.eqtl", "quant")),
              outcome.lst=list(c("SCZ_european", "binary"), c("T2D_european", "binary")),
              folder="geneexp_new",
              dt=data.table::data.table(),
              dt.steiger=NULL, #data.table::data.table(),
              output.file="adipose_geneexp_disease_european",
              p.thres=5e-5,
              genes=hc.genes,
              geneexp.exp=TRUE,
              api.token=api.token)

#------------------- Concatenate all results -----------------------
scz.dt <- data.table::fread("/project_data/CausalInference/geneexp_new/scz_geneexp_disease_multiancestry_analysis_MR.csv")
t2d.dt <- data.table::fread("/project_data/CausalInference/geneexp_new/t2d_geneexp_disease_multiancestry_analysis_MR.csv")
t2d.dt2 <- data.table::fread("/project_data/CausalInference/geneexp_new/t2d_cortex_afr_geneexp_disease_multiancestry_analysis_MR.csv")
gtex.dt <- data.table::fread("/project_data/CausalInference/geneexp_new/gtex_geneexp_disease_multiancestry_analysis_MR.csv")
dt <- rbind(scz.dt, t2d.dt, t2d.dt2, gtex.dt)  # adipose.dt
dt[, p.adj.fdr:=p.adjust(pval, method = "fdr"), by=method]
dt[, `:=` (id.exposure=NULL, id.outcome=NULL)]
data.table::fwrite(dt, "/project_data/CausalInference/geneexp_new/geneexp_disease_multiancestry_analysis_MR.csv")

eur.dt <- data.table::fread("/project_data/CausalInference/geneexp_new/geneexp_disease_european_analysis_MR.csv")
gtex.dt <- data.table::fread("/project_data/CausalInference/geneexp_new/gtex_geneexp_disease_european_analysis_MR.csv")
dt <- rbind(eur.dt, gtex.dt)  # adipose.dt
dt[, p.adj.fdr:=p.adjust(pval, method = "fdr"), by=method]
dt[, `:=` (id.exposure=NULL, id.outcome=NULL)]
data.table::fwrite(dt, "/project_data/CausalInference/geneexp_new/all_geneexp_disease_european_analysis_MR.csv")

############## Add indep.eQTL info to sc data ################
indep.qtl <- data.table::as.data.table(readxl::read_xlsx("/project_data/data/single_cell/brain/sc_eqtl_indep_info.xlsx"))
indep.qtl[cell_type=="OPCs / COPs", cell_type:="OPCs...COPs"]
indep.qtl[, cell_type:=sub(" ", "\\.", cell_type)]

for (c in c("Astrocytes", "Endothelial.cells", "Excitatory.neurons", "Inhibitory.neurons", "Microglia", "Oligodendrocytes", "OPCs...COPs", "Pericytes")){
  dt <- data.table::fread(paste0("/project_data/data/single_cell/brain/sc_eqtl_", c, ".csv"))
  dt[, indep.eQTL:=ifelse(rsID %in% indep.qtl[cell_type==c, SNP], TRUE, FALSE)]
  data.table::fwrite(dt, paste0("/project_data/data/single_cell/brain/sc_eqtl_", c, ".csv"))
}



