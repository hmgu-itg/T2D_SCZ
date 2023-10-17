library(data.table)
library(dplyr)
library(Hmisc)

setwd("C:/Users/ana.arruda/OneDrive - Helmholtz Zentrum München/Projects/SCZ_T2D")

################################################
#------------------- Functions -----------------
################################################
sort_alleles <- Vectorize(function(x, y) {
  paste(sort(c(x,y)), collapse = "_")
})

local.clump <- function(data, pval=5e-08){
  # https://mrcieu.github.io/ieugwasr/articles/local_ld.html
  clumped.dt <- ieugwasr::ld_clump(dplyr::tibble(rsid=data$SNP, pval=data$pval.exposure), #id=data$id.exposure),
                                   plink_bin=genetics.binaRies::get_plink_binary(),
                                   bfile="/reference_data/1kG/EUR_1kg_v3/EUR",
                                   clump_p=pval)
  return(data[data$SNP %in% clumped.dt$rsid, ])
}

get.filename <- function(trait, subset="european", type="univariate"){
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
  else if (trait %in% c("child_bmi", "adult_bmi")){
    dt <- read.child.adult.bmi(trait, type=type)
  }
  else print("ERROR: Subset not known")
  
  return(dt)
}

read.child.adult.bmi <- function(trait, type="univariate"){
  if(trait=="child_bmi" & type=="univariate"){
    dt <- data.table::fread("/project_data/CausalInference/child_bmi_ivs.csv")
    return(dt[, .(CHR=Chromosome, POS=`Base position`, rsID=SNP, EA=`Effect allele`, 
                  EAF=`Effect allele frequency`, NEA=`Other allele`, beta=Beta, 
                  se=`Standard Error`, pval=P, N=463005)])  
  }
  else if(trait=="adult_bmi" & type=="univariate"){
    dt <- data.table::fread("/project_data/CausalInference/adult_bmi_ivs.csv")
    return(dt[, .(CHR=Chromosome, POS=`Base position`, rsID=SNP, EA=`Effect allele`, 
                  EAF=`Effect allele frequency`, NEA=`Other allele`, beta=Beta, 
                  se=`Standard Error`, pval=P, N=463005)])
  }
  else if(trait=="child_bmi" & type=="multivariate"){
    dt <- data.table::fread("/project_data/CausalInference/comparison_bmi_ivs.csv")
    return(dt[, .(CHR=Chromosome, POS=`Base position`, rsID=SNP, EA=`Effect allele`, 
                  NEA=`Other allele`, beta=`Beta (Age 10)`, se=`SE (Age 10)`, pval=`P (Age 10)`,
                  N=463005, EAF)])
  }
  else if(trait=="adult_bmi" & type=="multivariate"){
    dt <- data.table::fread("/project_data/CausalInference/comparison_bmi_ivs.csv")
    return(dt[, .(CHR=Chromosome, POS=`Base position`, rsID=SNP, EA=`Effect allele`, 
                  NEA=`Other allele`, beta=`Beta (Adult)`, se=`SE (Adult)`, pval=`P (Adult)`, 
                  N=463005), EAF])
  }
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

get_odds_ratios_mvmr <- function(mr_res){
  mr_res$or <- exp(mr_res$b)
  mr_res$or_lci95 <- exp(mr_res$b - 1.96 * mr_res$se)
  mr_res$or_uci95 <- exp(mr_res$b + 1.96 * mr_res$se)
  return(mr_res)
}

run.mvmregger <- function(data, phen){
  #----------- Run MV MR-Egger from MVMR package -----------
  mv_exp_beta <- as.data.frame(data$exposure_beta)
  mv_exp_se <- as.data.frame(data$exposure_se)
  for(i in 1:length(colnames(mv_exp_beta))){
    assign(paste0("exposure", i), data$expname[data$expname$id.exposure==colnames(mv_exp_beta)[i],]$exposure)
    colnames(mv_exp_beta)[1] <- paste0(data$expname[data$expname$id.exposure==colnames(mv_exp_beta)[i],]$exposure, "_beta")
    colnames(mv_exp_se)[1] <- paste0(data$expname[data$expname$id.exposure==colnames(mv_exp_beta)[i],]$exposure, "_se")
  }
  # colnames(mv_exp_beta) <- c(paste0(exposure1, "_beta"), paste0(exposure2, "_beta"))
  # colnames(mv_exp_se) <- c(paste0(exposure1, "_se"), paste0(exposure2, "_se"))
  
  mv_out_beta <- as.data.frame(data$outcome_beta)
  mv_out_se <- as.data.frame(data$outcome_se)
  colnames(mv_out_beta) <- c(paste(phen, "beta", sep="_"))
  colnames(mv_out_se) <- c(paste(phen, "se", sep="_"))
  
  mv_exp_beta$SNP <- row.names(mv_exp_beta)
  mv_exp_se$SNP <- row.names(mv_exp_se)
  
  mv_exp <- merge(mv_exp_beta, mv_exp_se, by="SNP")
  mv_exp <- cbind(mv_exp, mv_out_beta)
  mv_exp <- cbind(mv_exp, mv_out_se)
  
  F.data <- MVMR::format_mvmr(BXGs = mv_exp[,c(2,3)],
                              BYG = mv_exp[,6],
                              seBXGs = mv_exp[,c(4,5)],
                              seBYG = mv_exp[,7],
                              RSID = mv_exp[,1])
  
  sres <- MVMR::strength_mvmr(r_input = F.data, gencov = 0)
  pres <- MVMR::pleiotropy_mvmr(r_input = F.data, gencov = 0)
  res_mvmr <- MVMR::ivw_mvmr(r_input = F.data)
  res_mvmr2 <- MVMR::qhet_mvmr(F.data, matrix(c(1,0,0,1), nrow=2, ncol=2), CI = F, iterations = 100)
  
  dt.egger <- merge(data.table::as.data.table(res_mvmr) %>% .[, .(exposure=c(exposure1, exposure2), effect=Estimate, se=`Std. Error`, tvalue=`t value`, `Pr(>|t|)`)],
                    data.table::as.data.table(res_mvmr2) %>% .[, .(exposure=c(exposure1, exposure2), effect.qhet=`Effect Estimates`)],
                    by="exposure") %>% 
    .[, `:=`(outcome=phen,Fstat=paste(sres$exposure1, sres$exposure2, sep="_"), Qstat=pres$Qstat, Qpval=pres$Qpval)]
  return(dt.egger)
}

run.mvmr <- function(data, phen, filename, folder="mvmr"){
  if(nrow(data$exposure_beta)<2) {
    print("After harmonizing data, not enough SNPs left :(")
    dt <- data.table::data.table()
    return(dt)
  }
  else {
    print("Multiple SNPs in common between exposure and outcome, running MVMR")
    
    # Assumption 1: check strength of IVs with F-statistica nd remove the ones with Fstat<10
    # data$fstat_per_snp <- (data$beta.exposure)^2/(data$se.exposure)^2
    # data <- data[data$fstat_per_snp>10,]
    
    #----------- Run MV TwoSampleMR ivw -----------
    res <- TwoSampleMR::mv_multiple(data, pval_threshold = 5e-01, plots=TRUE)
    
    #------------ Save output and plots
    ggplot2::ggsave(paste0("/project_data/CausalInference/", folder, "/mvmr_adiposity_", phen, ".jpg"), gridExtra::arrangeGrob(grobs = res$plots))
    dt <- data.table::as.data.table(res$result[,-c(1,3)])
  }
  return(dt)
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

wrap.MVMR <- function(exposure.lst, outcome.lst, exp.phen, out.phen, folder, output.file, dt, dt.egger=NULL, ivs=NULL, custom.clump.p=5e-08, egger=TRUE, type="multivariate"){
  # exposure.lst = list of lists of vectors
  clump.exposure <- data.frame()
  
  for (exp in exposure.lst){
    #----------- Read and clump exposure -----------
    exposure <- get.data(get.filename(trait=exp[1], subset=exp[2], type=type), data.dir="exposure", data.type=exp[3])
    exposure$exposure <- paste(exp[1], exp[2], sep=".")
    if(!is.null(ivs)) clump.exp <- exposure[exposure$SNP %in% ivs,]
    else clump.exp <- local.clump(data=exposure, pval=custom.clump.p)
    
    clump.exposure <- rbind(clump.exposure, clump.exp)
  }
  
  for (out in outcome.lst) {
    #----------- Get outcome -----------
    outcome <- get.data(get.filename(trait=out[1], subset=out[2]), data.dir="outcome", data.type=out[3])
    outcome$outcome <- paste(out[1], out[2], sep=".")
    
    #----------- Harmonize data -----------
    data <- TwoSampleMR::mv_harmonise_data(exposure_dat=clump.exposure, outcome_dat=outcome)
    
    #----------- Run TwoSampleMR for unfiltered data -----------
    dt <- rbind(dt, run.mvmr(data, filename=paste0(paste(exp, collapse="."), "_", paste(out, collapse=".")), folder=folder, phen=out[1]), fill=TRUE)
    if(egger) dt.egger <- rbind(dt.egger, run.mvmregger(data, phen=out[1]), fill=TRUE)
  }
  
  data.table::fwrite(dt, paste0("/project_data/CausalInference/", folder, "/", output.file, "_analysis_MR.csv"))
  if(egger) data.table::fwrite(dt.egger, paste0("/project_data/CausalInference/", folder, "/egger_", output.file, "_analysis_MR.csv"))
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
      TwoSampleMR::mr_scatter_plot(res, data)   #scatter plot
      ggplot2::ggsave(paste0("/project_data/CausalInference/", folder, "/TwoSampleMR_", filename, ".jpg"))
      
      # check heterogeneity with Cochran's Q-statistic
      het <- TwoSampleMR::mr_heterogeneity(data)   #heterogeneity statistics
      
      # Assumption 2: check pleiotropy with MR-Egger intercept
      plt <- TwoSampleMR::mr_pleiotropy_test(data)  #MR Egger intercept for directional pleiotropy
      
      tmp.dt <- data.table::as.data.table(merge(res, het, by=c("id.exposure", "id.outcome", "outcome", "exposure", "method"), all.x=T)) %>% 
        .[, `:=`(plt.egger_intercept=as.numeric(NA),plt.pval=as.numeric(NA),plt.se=as.numeric(NA))] %>%
        .[method=="MR Egger", `:=`(plt.egger_intercept=plt$egger_intercept,plt.pval=plt$pval,plt.se=plt$se)]
      
      tmp.dt[, Fstat:=fstat_overall]
      
      # Sensitivity analyses: leave-one-out and forest plot
      TwoSampleMR::mr_forest_plot(TwoSampleMR::mr_leaveoneout(data))
      ggplot2::ggsave(paste0("/project_data/CausalInference/", folder, "/LeaveOneOut_", filename, ".jpg"))
      TwoSampleMR::mr_forest_plot(TwoSampleMR::mr_singlesnp(data))
      ggplot2::ggsave(paste0("/project_data/CausalInference/", folder, "/SingleSNP_", filename, ".jpg"))
    }
  }
  return(tmp.dt)
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
        filtered.data <- filtered.data[filtered.data$steiger_dir==TRUE,]  # & filtered.data$steiger_pval<0.05,]
        
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

############################################################################
#------------------------ Local data preprocessing ------------------------#
############################################################################
child <- as.data.table(readxl::read_xlsx("data/rict052821.wt.xlsx", sheet=3, skip=3))
adult <- as.data.table(readxl::read_xlsx("data/rict052821.wt.xlsx", sheet=4, skip=3))
comparison <- as.data.table(readxl::read_xlsx("data/rict052821.wt.xlsx", sheet=5, skip=3))

child[SNP %in% adult$SNP, .N]
adult[SNP %in% child$SNP, .N]

child[SNP %nin% comparison$SNP, .N]
adult[SNP %nin% comparison$SNP, .N]

comparison.eaf <- merge(comparison, child[, .(SNP, EAF=`Effect allele frequency`)], by="SNP", all.x=TRUE)
comparison.eaf <- merge(comparison.eaf, adult[, .(SNP, EAF=`Effect allele frequency`)], by="SNP", all.x=TRUE)
comparison.eaf[, EAF:=ifelse(is.na(EAF.x), EAF.y, EAF.x)]
comparison.eaf[, `:=` (EAF.x=NULL, EAF.y=NULL)]

fwrite(child, "data/child_bmi_ivs.csv")
fwrite(adult, "data/adult_bmi_ivs.csv")
fwrite(comparison.eaf[SNP %in% c(adult$SNP, child$SNP)], "data/comparison_bmi_ivs.csv")

#######################################################################
#----------------------------- Read data -----------------------------#
#######################################################################
child.bmi <- fread("/project_data/CausalInference/child_bmi_ivs.csv")
adult.bmi <- fread("/project_data/CausalInference/adult_bmi_ivs.csv")
comparison.bmi <- fread("/project_data/CausalInference/comparison_bmi_ivs.csv")

############################################################################
#----------------------------- Univariable MR -----------------------------#
############################################################################
dt <- data.table::data.table()
dt.steiger <- data.table::data.table()

wrap.MR(exposure.lst=list(c("child_bmi", "european", "quant")),  
        outcome.lst=c(list(c("SCZ", "multiancestry", "binary")),
                      list(c("T2D", "multiancestry", "binary"))),
        ivs=child.bmi$SNP, 
        folder="child_adult_bmi",
        dt=dt,
        dt.steiger=dt.steiger,
        output.file="univariate_child_bmi")


dt <- data.table::data.table()
dt.steiger <- data.table::data.table()

wrap.MR(exposure.lst=list(c("adult_bmi", "european", "quant")),  
        outcome.lst=c(list(c("SCZ", "multiancestry", "binary")),
                      list(c("T2D", "multiancestry", "binary"))),
        ivs=adult.bmi$SNP, 
        folder="child_adult_bmi",
        dt=dt,
        dt.steiger=dt.steiger,
        output.file="univariate_adult_bmi")


#---------------- Sensitivity analysis with European subset only
dt <- data.table::data.table()
dt.steiger <- data.table::data.table()

wrap.MR(exposure.lst=list(c("child_bmi", "european", "quant")),  
        outcome.lst=c(list(c("SCZ", "european", "binary")),
                      list(c("T2D", "european", "binary"))),
        ivs=child.bmi$SNP, 
        folder="child_adult_bmi",
        dt=dt,
        dt.steiger=dt.steiger,
        output.file="univariate_child_bmi_european")


dt <- data.table::data.table()
dt.steiger <- data.table::data.table()

wrap.MR(exposure.lst=list(c("adult_bmi", "european", "quant")),  
        outcome.lst=c(list(c("SCZ", "european", "binary")),
                      list(c("T2D", "european", "binary"))),
        ivs=adult.bmi$SNP, 
        folder="child_adult_bmi",
        dt=dt,
        dt.steiger=dt.steiger,
        output.file="univariate_adult_bmi_european")

############################################################################
#---------------------------- Multivariable MR ----------------------------#
############################################################################
dt <- data.table::data.table()
dt.egger <- data.table::data.table()

wrap.MVMR(exposure.lst=list(c("child_bmi", "european", "quant"),
                            c("adult_bmi", "european", "quant")),
          outcome.lst=c(list(c("SCZ", "multiancestry", "binary")),
                        list(c("T2D", "multiancestry", "binary"))),
          folder="child_adult_bmi",
          dt=dt,
          dt.egger=dt.egger,
          ivs=comparison.bmi$SNP,
          output.file="multivariate_adult_child_bmi")

# add OR and adjusted pvalue
project_folder <- "C:/Users/ana.arruda/OneDrive - Helmholtz Zentrum München/Projects/SCZ_T2D/"
setwd(project_folder)
dt <- data.table::fread("causal_inference/child_adult_bmi/multivariate_adult_child_bmi_analysis_MR.csv")
dt <- get_odds_ratios_mvmr(dt)
dt[, p.adj := p.adjust(pval, method = "fdr")]
data.table::fwrite(dt, paste0(project_data, "causal_inference/child_adult_bmi/multivariate_adult_child_bmi_analysis_MR.csv"))

#---------------- Sensitivity analysis with European subset only
dt <- data.table::data.table()
dt.egger <- data.table::data.table()

wrap.MVMR(exposure.lst=list(c("child_bmi", "european", "quant"),
                            c("adult_bmi", "european", "quant")),
          outcome.lst=c(list(c("SCZ", "european", "binary")),
                        list(c("T2D", "european", "binary"))),
          folder="child_adult_bmi",
          dt=dt,
          dt.egger=dt.egger,
          ivs=comparison.bmi$SNP,
          output.file="multivariate_adult_child_bmi_european")

# add OR and adjusted pvalue
project_folder <- "C:/Users/ana.arruda/OneDrive - Helmholtz Zentrum München/Projects/SCZ_T2D/"
setwd(project_folder)
dt <- data.table::fread("causal_inference/child_adult_bmi/multivariate_adult_child_bmi_european_analysis_MR.csv")
dt <- get_odds_ratios_mvmr(dt)
dt[, p.adj := p.adjust(pval, method = "fdr")]
data.table::fwrite(dt, paste0(project_data, "causal_inference/child_adult_bmi/multivariate_adult_child_bmi_european_analysis_MR.csv"))
