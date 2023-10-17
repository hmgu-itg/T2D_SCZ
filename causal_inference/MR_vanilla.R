library(dplyr)

################################################
#------------------- Functions -----------------
################################################
local.clump <- function(data, pval=5e-08){
  # https://mrcieu.github.io/ieugwasr/articles/local_ld.html
  clumped.dt <- ieugwasr::ld_clump(dplyr::tibble(rsid=data$SNP, pval=data$pval.exposure), #id=data$id.exposure),
                                   plink_bin=genetics.binaRies::get_plink_binary(),
                                   bfile="/reference_data/1kG/EUR_1kg_v3/EUR",
                                   clump_p=pval)
  return(data[data$SNP %in% clumped.dt$rsid, ])
}

read.pgc <- function(trait){
  if(trait=="ADHD") {
    dt <- data.table::fread(paste0("/project_data/data/PGC_", trait, ".txt.gz"))
    return(dt[, .(CHR, POS=BP, rsID=SNP, EA=A1, NEA=A2, beta=log(OR), se=SE, pval=P, Ncases=Nca, N=Nca+Nco)])  
  }
  else if(trait=="AnorexiaNervosa") {
    dt <- data.table::fread(paste0("/project_data/data/PGC_", trait, ".vcf.tsv.gz"))
    return(dt[, .(CHR=CHROM, POS, rsID=ID, EA=REF, NEA=ALT, beta=BETA, se=SE, pval=PVAL, Ncases=NCAS, N=NCAS+NCON)]) 
  }
  else if(trait=="Autism") {
    dt <- data.table::fread(paste0("/project_data/data/PGC_", trait, ".txt.gz"))
    return(dt[, .(CHR, POS=BP, rsID=SNP, EA=A1, NEA=A2, beta=log(OR), se=SE, pval=P, Ncases=18381, N=46350)]) 
  }
  else if(trait=="BipolarDisorder") {
    dt <- data.table::fread(paste0("/project_data/data/PGC_", trait, ".vcf.tsv.gz"))
    return(dt[, .(CHR=`#CHROM`, POS, rsID=ID, EA=A1, NEA=A2, beta=BETA, se=SE, pval=PVAL, Ncases=NCAS, N=NCAS+NCON)]) 
  }
  else if(trait=="Depression") {
    dt <- data.table::fread(paste0("/project_data/data/PGC_", trait, ".txt"))
    return(dt[, .(CHR=NA, POS=NA, rsID=MarkerName, EA=toupper(A1), NEA=toupper(A2), beta=LogOR, se=StdErrLogOR, pval=P, Ncases=170756, N=500199)]) 
  }
  else if(trait=="PanicDisorder"){
    dt <- data.table::fread(paste0("/project_data/data/PGC_", trait, ".vcf.tsv.gz"))
    return(dt[, .(CHR=`#CHROM`, POS, rsID=ID, EA=A1, NEA=A2, beta=BETA, se=SE, pval=PVAL, Ncases=NCAS, N=NCAS+NCON)])
  }
  else if(trait=="PTSD_eur") {
    dt <- data.table::fread(paste0("/project_data/data/PGC_", trait, ".txt.gz"))
    return(dt[, .(CHR, POS=BP, rsID=SNP, EA=A1, NEA=A2, beta=log(OR), se=SE, pval=P, Ncases=Nca, N=Nca+Nco)]) 
  }
  else {
    return(print("Trait not recognized.")) 
  }
}

read.magic <- function(trait, subset){
  if(subset=="european") pop="EUR"
  else if(subset=="multiancestry") pop="TA"
  
  dt <- data.table::fread(paste0("/project_data/data/MAGIC1000G_", trait, "_", pop, ".tsv.gz"))
  return(dt[, .(CHR=chromosome, POS=base_pair_location, rsID=variant, EA=effect_allele, NEA=other_allele, 
                EAF=effect_allele_frequency, beta=beta, se=standard_error, pval=p_value, N=sample_size)])
}

read.adiposity <- function(trait, subset){
  if (trait %in% c("bmi", "whr")){
    dt <- data.table::fread(paste0("/project_data/data/", trait, ".giant-ukbb.meta-analysis.combined.23May2018.txt.gz"))
    return(dt[, .(CHR, POS, rsID=sub(":.*", "", SNP), EA=Tested_Allele, NEA=Other_Allele, EAF=Freq_Tested_Allele, beta=BETA, 
                  se=SE, pval=P, N)])  
  }
  else if (trait %in% c("body_fat_percentage", "whole_body_fat_mass")){
    dt <- data.table::fread(paste0("/project_data/data/", trait, "_neale_rsid.txt.gz"))
    return(dt[, .(CHR, POS, rsID=SNP, EA=minor_allele, EAF=minor_AF, NEA=major_allele, beta, se, pval, N=n_complete_samples)])  
  }
  else if (subset=="european" & trait %in% c("ldl", "tgs", "hdl", "apoa", "apob")){
    dt <- data.table::fread(paste0("/project_data/data/", trait, "_neale_rsid.txt.gz"))
    return(dt[, .(CHR, POS, rsID, EA=minor_allele, EAF=minor_AF, NEA=major_allele, beta, se, pval, N=n_complete_samples, ID)])  
  }
  # else if (subset=="multiancestry" & trait %in% c("body_fat_percentage", "whole_body_fat_mass")){
  #   dt <- data.table::fread(paste0("/project_data/data/", trait, "_multiancestry_neale_rsid.txt.gz"))
  #   return(dt[, .(CHR, POS, rsID=SNP, EA=minor_allele, EAF=minor_AF, NEA=major_allele, beta, se, pval, N=n_complete_samples)])  
  # }
  else if (subset=="multiancestry" & trait %in% c("ldl", "tgs", "hdl", "apoa", "apob")){
    dt <- data.table::fread(paste0("/project_data/data/", trait, "_multiancestry_neale_rsid.txt.gz"))
    return(dt[, .(CHR, POS, rsID, EA=minor_allele, EAF=minor_AF, NEA=major_allele, beta, se, pval, N=n_complete_samples, ID)])  
  }
  else {
    return(print("Trait not recognized.")) 
  }
}

get.filename <- function(trait, subset="european"){
  if(trait %in% c("T2D", "SCZ")){
    dt <- data.table::fread(paste0("/project_data/data/", trait, "_", subset, ".txt"))
  }
  else if(trait %in% c("2hGlu", "FI", "FG", "HbA1c")){
    dt <- read.magic(trait, subset)
  }
  else if (trait %in% c("ADHD", "AnorexiaNervosa", "Autism", "BipolarDisorder", "Depression", "PanicDisorder", "PTSD_eur")){
    dt <- read.pgc(trait)
  }
  else if (trait %in% c("bmi", "whr", "whole_body_fat_mass", "body_fat_percentage", "ldl", "tgs", "hdl", "apoa", "apob")){
    dt <- read.adiposity(trait, subset)
  }
  else print("ERROR: Subset not known")
  
  return(dt)
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
                                     samplesize_col="N")
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
                                     ncase_col="Ncases",
                                     samplesize_col="N")
  }
  else {
    print("Data type not recognized (options: quant or binary).")
  }
  
  return(data)
}

run.TwoSampleMR <- function(data, filename, folder="vanilla"){
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
    print("Multiple SNPs in common between exposure and outcome, running MR")
    
    # Assumption 1: check strength of IVs with F-statistica nd remove the ones with Fstat<10
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
      
      #------------ Sensitivity analysis
      #TwoSampleMR::mr_scatter_plot(res, data)   #scatter plot
      #ggplot2::ggsave(paste0("/project_data/CausalInference/", folder, "/TwoSampleMR_", filename, ".jpg"))
      
      # check heterogeneity with Cochran's Q-statistic
      het <- TwoSampleMR::mr_heterogeneity(data)   #heterogeneity statistics
      
      # Assumption 2: check pleiotropy with MR-Egger intercept
      plt <- TwoSampleMR::mr_pleiotropy_test(data)  #MR Egger intercept for directional pleiotropy
  
      tmp.dt <- data.table::as.data.table(merge(res, het, by=c("id.exposure", "id.outcome", "outcome", "exposure", "method"), all.x=T)) %>% 
        .[, `:=`(plt.egger_intercept=as.numeric(NA),plt.pval=as.numeric(NA),plt.se=as.numeric(NA))] %>%
        .[method=="MR Egger", `:=`(plt.egger_intercept=plt$egger_intercept,plt.pval=plt$pval,plt.se=plt$se)]
      
      tmp.dt[, Fstat:=fstat_overall]
      
      # Sensitivity analyses: leave-one-out and forest plot
      #TwoSampleMR::mr_forest_plot(TwoSampleMR::mr_leaveoneout(data))
      #ggplot2::ggsave(paste0("/project_data/CausalInference/", folder, "/LeaveOneOut_", filename, ".jpg"))
      #TwoSampleMR::mr_forest_plot(TwoSampleMR::mr_singlesnp(data))
      #ggplot2::ggsave(paste0("/project_data/CausalInference/", folder, "/SingleSNP_", filename, ".jpg"))
    }
  }
  return(tmp.dt)
}

output.wide.result <- function(dt){
  dt[, `:=` (id.exposure=NULL, id.outcome=NULL)]
  dt <- data.table::dcast(dt, exposure + outcome + nsnp + Fstat ~ method, 
                          value.var=c("b", "se", "pval", "lo_ci", "up_ci", "or", "or_lci95", "or_uci95", "Q", "Q_df", "Q_pval", "plt.egger_intercept", "plt.pval", "plt.se", "p.adj.fdr")) 
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

wrap.MR <- function(exposure.lst, outcome.lst, exp.phen, out.phen, folder, output.file, dt, dt.steiger=NULL, ivs=NULL, custom.clump.p=5e-08){
  for (exp in exposure.lst){
    #----------- Read and clump exposure -----------
    exposure <- get.data(get.filename(trait=exp[1], subset=exp[2]), data.dir="exposure", data.type=exp[3])
    exposure$exposure <- paste(exp[1], exp[2], sep=".")
    if(!is.null(ivs)) clump.exposure <- exposure[exposure$SNP %in% ivs,]
    else clump.exposure <- local.clump(data=exposure, pval=custom.clump.p)
    
    for (out in outcome.lst) {
      #----------- Get outcome -----------
      outcome <- get.data(get.filename(trait=out[1], subset=out[2]), data.dir="outcome", data.type=out[3])
      outcome$outcome <- paste(out[1], out[2], sep=".")
      
      #----------- Harmonize data -----------
      data <- subset(TwoSampleMR::harmonise_data(exposure_dat=clump.exposure, outcome_dat=outcome), mr_keep==TRUE)
      
      #----------- Run TwoSampleMR for unfiltered data -----------
      dt <- rbind(dt, run.TwoSampleMR(data, filename=paste0(paste(exp, collapse="."), "_", paste(out, collapse=".")), folder=folder), fill=TRUE)
      
      if(!is.null(dt.steiger)) {
        #----------- Apply steiger filtering on data -----------
        # data$prevalence.exposure <- 0.1
        # data$prevalence.outcome <- 0.01
        # data$ncontrol.exposure <- data$samplesize.exposure - data$ncase.exposure
        # data$ncontrol.outcome <- data$samplesize.outcome - data$ncase.outcome
        # data$units.exposure <- "log odds"
        # data$units.outcome <- "log odds"
        filtered.data <- TwoSampleMR::steiger_filtering(data)
        filtered.data <- filtered.data[filtered.data$steiger_dir==TRUE & filtered.data$steiger_pval<0.05,]
        
        #----------- Run TwoSampleMR for Steiger-filtered data -----------
        dt.steiger <- rbind(dt.steiger,run.TwoSampleMR(filtered.data, filename=paste0(paste(exp, collapse="."), "_", paste(out, collapse="."), "_Steiger"), folder=folder), fill=TRUE)
      }
    }
  }
  
  dt[, p.adj.fdr:=p.adjust(pval, method = "fdr"), by=method]
  data.table::fwrite(output.wide.result(dt), paste0("/project_data/CausalInference/", folder, "/", output.file, "_analysis_MR.csv"))
  
  if(!is.null(dt.steiger)){
    dt.steiger[, p.adj.fdr:=p.adjust(pval, method = "fdr"), by=method]
    data.table::fwrite(output.wide.result(dt.steiger), paste0("/project_data/CausalInference/", folder, "/Steiger_", output.file, "_analysis_MR.csv"))
  }   
}

magic.traits <- c("2hGlu", "FG", "FI", "HbA1c")
pgc.traits <- c("ADHD", "AnorexiaNervosa", "Autism", "BipolarDisorder", "Depression")  # , "PanicDisorder", "PTSD_eur")
adiposity.traits <- c("bmi", "whr", "whole_body_fat_mass", "body_fat_percentage", "ldl", "tgs", "hdl", "apoa", "apob")

########################################################################################
#------------------ Exposure: T2D.multiancestry; Outcome: all SCZ ---------------------#
########################################################################################
dt <- data.table::data.table()
dt.steiger <- data.table::data.table()

wrap.MR(exposure.lst=list(c("T2D", "multiancestry", "binary")),  
        outcome.lst=c(list(c("SCZ", "multiancestry", "binary")),
                     lapply(pgc.traits, function(i){i <- c(i, "european", "binary")})),
        # ivs=T2D.ivs, 
        # custom.clump.p=5e-09,
        folder="test",
        dt=dt,
        dt.steiger=dt.steiger,
        output.file="T2D_all_multiancestry")

####################################################################################
#--------------- Exposure: SCZ.multiancestry; Outcome: all T2D --------------------#
####################################################################################
dt <- data.table::data.table()
dt.steiger <- data.table::data.table()

wrap.MR(exposure.lst=list(c("SCZ", "multiancestry", "binary")),
        outcome.lst=c(list(c("T2D", "multiancestry", "binary")), 
                      lapply(magic.traits, function(i){i <- c(i, "european", "quant")}),
                      lapply(adiposity.traits, function(i){i <- c(i, "european", "quant")})), 
        folder="test",
        dt=dt,
        dt.steiger=dt.steiger,
        output.file="SCZ_all_multiancestry")

####################################################################################
#----------------- Exposure: T2D all; Outcome: SCZ multiancestry ------------------#
####################################################################################
dt <- data.table::data.table()
dt.steiger <- data.table::data.table()

wrap.MR(exposure.lst=c(list(c("T2D", "multiancestry", "binary")),
                       lapply(magic.traits, function(i){i <- c(i, "european", "binary")}),
                       lapply(adiposity.traits, function(i){i <- c(i, "european", "quant")})), 
        outcome.lst=list(c("SCZ", "multiancestry", "binary")), 
        folder="test",
        dt=dt,
        dt.steiger=dt.steiger,
        output.file="all_SCZ_multiancestry")

####################################################################################
#----------------- Exposure: SCZ all; Outcome: T2D multiancestry ------------------#
####################################################################################
dt <- data.table::data.table()
dt.steiger <- data.table::data.table()

wrap.MR(exposure.lst=c(list(c("SCZ", "multiancestry", "binary")),
                       lapply(pgc.traits, function(i){i <- c(i, "european", "binary")})), 
        outcome.lst=list(c("T2D", "multiancestry", "binary")), 
        folder="test",
        dt=dt,
        dt.steiger=dt.steiger,
        output.file="all_T2D_multiancestry")

#################### Sensitivity analysis using European data only #################

########################################################################################
#------------------ Exposure: T2D European; Outcome: all SCZ ---------------------#
########################################################################################
dt <- data.table::data.table()
dt.steiger <- data.table::data.table()

wrap.MR(exposure.lst=list(c("T2D", "european", "binary")),  
        outcome.lst=c(list(c("SCZ", "european", "binary")),
                      lapply(pgc.traits, function(i){i <- c(i, "european", "binary")})),
        # ivs=T2D.ivs, 
        folder="test",
        # custom.clump.p=5e-09,
        dt=dt,
        dt.steiger=dt.steiger,
        output.file="T2D_all_european")

####################################################################################
#--------------- Exposure: SCZ European; Outcome: all T2D --------------------#
####################################################################################
dt <- data.table::data.table()
dt.steiger <- data.table::data.table()

wrap.MR(exposure.lst=list(c("SCZ", "european", "binary")),
        outcome.lst=c(list(c("T2D", "european", "binary")), 
                      lapply(magic.traits, function(i){i <- c(i, "european", "quant")}),
                      lapply(adiposity.traits, function(i){i <- c(i, "european", "quant")})), 
        folder="test",
        dt=dt,
        dt.steiger=dt.steiger,
        output.file="SCZ_all_european") 

####################################################################################
#----------------- Exposure: T2D all; Outcome: SCZ European ------------------#
####################################################################################
dt <- data.table::data.table()
dt.steiger <- data.table::data.table()

wrap.MR(exposure.lst=c(list(c("T2D", "european", "binary")),
                       lapply(magic.traits, function(i){i <- c(i, "european", "binary")}),
                       lapply(adiposity.traits, function(i){i <- c(i, "european", "quant")})), 
        outcome.lst=list(c("SCZ", "european", "binary")), 
        folder="test",
        dt=dt,
        dt.steiger=dt.steiger,
        output.file="all_SCZ_european")

####################################################################################
#----------------- Exposure: SCZ all; Outcome: T2D European ------------------#
####################################################################################
dt <- data.table::data.table()
dt.steiger <- data.table::data.table()

wrap.MR(exposure.lst=c(list(c("SCZ", "european", "binary")),
                       lapply(pgc.traits, function(i){i <- c(i, "european", "binary")})), 
        outcome.lst=list(c("T2D", "european", "binary")), 
        folder="test",
        dt=dt,
        dt.steiger=dt.steiger,
        output.file="all_T2D_european")

########################################################################################
#------------------ Exposure: adiposity; Outcome: SCZ/T2D ---------------------#
########################################################################################
dt <- data.table::data.table()
dt.steiger <- data.table::data.table()

wrap.MR(exposure.lst=lapply(adiposity.traits[5:9], function(i){i <- c(i, "european", "quant")}), 
        outcome.lst=c(list(c("T2D", "european", "binary")), 
                      list(c("SCZ", "european", "binary"))), 
        folder="test",
        dt=dt,
        dt.steiger=dt.steiger,
        ivs=data.table::fread("/project_data/CausalInference/mvmr/adiposity_ivs.csv")$snp, 
        output.file="adiposity_disease_european")

############## Concatenate all results #####################
dt <- data.table::data.table()
dt.steiger <- data.table::data.table()
folder = "test"
for(file in c("T2D_all_multiancestry", "SCZ_all_multiancestry", "all_SCZ_multiancestry", "all_T2D_multiancestry")){
  dt <- rbind(dt, data.table::fread(paste0("/project_data/CausalInference/", folder, "/", file, "_analysis_MR.csv")), fill=TRUE)
  dt.steiger <- rbind(dt.steiger, data.table::fread(paste0("/project_data/CausalInference/", folder, "/Steiger_", file, "_analysis_MR.csv")), fill=TRUE)
}
data.table::fwrite(dt, paste0("/project_data/CausalInference/", folder, "/MR_vanilla_multiancestry.csv"))
data.table::fwrite(dt.steiger, paste0("/project_data/CausalInference/", folder, "/MR_vanilla_multiancestry_Steiger.csv"))

dt <- data.table::data.table()
dt.steiger <- data.table::data.table()
for(file in c("T2D_all_european", "SCZ_all_european", "all_SCZ_european", "all_T2D_european")){
  dt <- rbind(dt, data.table::fread(paste0("/project_data/CausalInference/", folder, "/", file, "_analysis_MR.csv")), fill=TRUE)
  dt.steiger <- rbind(dt.steiger, data.table::fread(paste0("/project_data/CausalInference/", folder, "/Steiger_", file, "_analysis_MR.csv")), fill=TRUE)
  
}
data.table::fwrite(dt, paste0("/project_data/CausalInference/", folder, "/MR_vanilla_european.csv"))
data.table::fwrite(dt.steiger, paste0("/project_data/CausalInference/", folder, "/MR_vanilla_european_Steiger.csv"))



