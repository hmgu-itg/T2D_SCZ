#####################################################################################
#---------------------------------- Functions --------------------------------------#
#####################################################################################
sort_alleles <- Vectorize(function(x, y) {
  paste(sort(c(x,y)), collapse = "_")
})

# Add ID column and get list of IDs with corresponding traits
add_ID <- function(dict, trait) {
  dict[[trait]] <- dict[[trait]][, ID := paste(paste(CHR, POS, sep = ":"), sort_alleles(EA, NEA), sep = "_")]
}

#' Define regions aroung distinct association signals from all traits of interest
#'
#' @param file.lst list of paths for excel file containing independent signals for the traits of interest (columns: CHR, POS, rsID, ID)
#' @param out.path filepath to save regions for trait
#' @param wd window around signals, defines regions' size
#' @return regions to perform colocalization
#' @export
get.regions <- function(file.lst, out.path, out.name="GWAS_regions.csv", wd=1e+6) {
  signals <- data.table::data.table()
  for (file in file.lst) {
    signals <- rbind(signals, data.table::fread(file))
  }
  regions <- unique(na.omit(signals))[, ":="(start = ifelse(POS-wd>0, POS-wd, 0), end = POS+wd)]
  regions <- regions[order(CHR, POS)]
  data.table::fwrite(regions, paste0(out.path, out.name))
  # GWAS_overlap_regions <- regions[, data.table::as.data.table(IRanges::reduce(IRanges::IRanges(start, end), min.gapwidth = 0L)), CHR]
  
  return(regions)
}

#' Define regions around distinct association signals from all traits of interest
#'
#' @param gwas.file excluding indels (columns: snp, chr, pos, ea, nea, eaf, maf, beta, se, pval, n, ncases)
#' @param wd window around signals, defines regions' size
#' @return regions to perform colocalization
#' @export
read.gwas <- function(trait, gwas.file, dict, snp="rsID", chr="CHR", pos="POS", ea="EA", nea="NEA", eaf="EAF", maf="MAF", beta="beta", se="se", pval="pval", n="N", ncases="Ncases", ID="ID") {
  dict[[trait]] <- data.table::fread(gwas.file, select=c(snp, chr, pos, ea, nea, eaf, maf, beta, se, pval, n, ncases, ID), verbose = FALSE)
  data.table::setcolorder(dict[[trait]], c(snp, chr, pos, ea, nea, eaf, maf, beta, se, pval, n, ncases, ID))
  data.table::setnames(dict[[trait]], c("snp", "chr", "pos", "ea", "nea", "eaf", "maf", "beta", "se", "pval", "n", "ncases", "ID"))
  
  # TO DO: add rsID if missing
  
  # Add ID: 
  if (is.na(dict[[trait]]$ID)) add_ID(dict, trait)
}

# Extract GWAS variants within 1Mb from merged independent signals
select.gwas <- function(dict, regions, trait) {
  clumped.gwas <- data.table::data.table()
  for (i in 1:nrow(regions)) {
    clumped.gwas <- rbind(clumped.gwas, dict[[trait]][chr %in% regions$CHR[i] &  dplyr::between(pos, regions$start[i], regions$end[i])])
    dict[[paste0(trait, "_coloc")]] <- unique(clumped.gwas)
  }
}

# Flip alleles based on reference GWAS
flip.alleles <- function(dict, ref_trait, alt_trait) {
  tmp.dt <- merge(dict[[ref_trait]][, .(chr, pos, ID, ea, nea)], dict[[paste0(alt_trait, "_coloc")]], by=c("chr", "pos", "ID"), suffixes=c(".ref", ""), all.y=TRUE)
  tmp.dt[, `:=` (ea.new=ifelse(ea!=ea.ref, nea, ea), nea.new=ifelse(ea!=ea.ref, ea, nea), beta=ifelse(ea!=ea.ref, -beta, beta)), by=seq_len(nrow(tmp.dt))]
  tmp.dt[, `:=` (ea=ea.new, nea=nea.new), by=seq_len(nrow(tmp.dt))]
  dict[[paste0(alt_trait, "_coloc")]] <- tmp.dt[, `:=`(ea.new=NULL, nea.new=NULL, ea.ref=NULL, nea.ref=NULL)]
}

preprocess.gwas <- function(signals.lst, gwas.lst, traits, out.path, out.file.suffix="_precoloc.txt", region.file.name="GWAS_regions.csv", wd=1e+6) {
  GWAS <- hash::hash()
  
  print("Loading independent GWAS signals")
  regions <- get.regions(signals.lst, out.path, region.file.name, wd)
  
  # for (trait in traits){
  idx=1
  for (trait in c("t1", "t2")){  # ES FUNKTIONIERT :O
    print(paste0("Loading GWAS data for ", traits[idx]))
    read.gwas(trait, get(paste0(trait, ".filepath")), GWAS)
    
    print("Selecting GWAS signals around merged independent signals")
    select.gwas(GWAS, regions, trait)
    
    if(trait!="t1") {
      # rm(dict[[trait]])
      print("Flipping alleles if needed")
      flip.alleles(dict=GWAS, ref_trait="t1", alt_trait=trait)
    }
    
    print("Outputing GWAS files for colocalization")
    data.table::fwrite(GWAS[[paste0(trait, "_coloc")]], paste0(out.path, "GWAS_", traits[idx], out.file.suffix))
    idx=idx+1
  }
  # del(traits[1], dict)
}

#' Output a concise data table for the 95% credible set of colocalized regions
#'
#' @param res list of colocalization result for each region
#' @param credset list of 95% credible set for each region
#' @return data table for 95% credible set with PP4 for each variant
#' @export
output.coloc.result <- function(data, res, region, all.regions, pp4=TRUE) {
  if(pp4) {
    # output credible set as data.table for GWAS eQTL colocalization
    o <- order(res$results$SNP.PP.H4,decreasing=TRUE)
    cs <- cumsum(res$results$SNP.PP.H4[o])
    w <- which(cs > 0.95)[1]
    credset <- res$results[o,][1:w,]$snp
    lead.snp <- credset[1]
    credset.dt <- data.table::data.table(region=region, 
                                         snp=credset,
                                         ID=data[, ID[match(credset, snp)]],
                                         CHR=data[, chr[match(credset, snp)]],
                                         POS=data[, pos[match(credset, snp)]],
                                         PP4=rep(round(unname(res$summary[paste0("PP.H4.abf")]), digits=3), length(credset)),
                                         SNP.PP4=res$results$SNP.PP.H4[o][1:length(credset)],
                                         N_credset=length(credset))
    credset.dt <- merge(credset.dt, data[, .(ID, snp, ea, nea, get(paste0("beta.", traits[1])), get(paste0("se.", traits[1])), get(paste0("pval.", traits[1])), get(paste0("beta.", traits[2])), get(paste0("se.", traits[2])), get(paste0("pval.", traits[2])))], by=c("snp", "ID"))
    colnames(credset.dt) <- c(colnames(credset.dt)[1:10], eval(paste0("beta.", traits[1])), eval(paste0("se.", traits[1])), eval(paste0("pval.", traits[1])), eval(paste0("beta.", traits[2])), eval(paste0("se.", traits[2])), eval(paste0("pval.", traits[2])))
    
    res.dt <- data.table::data.table(region=region,
                                     N_credset=length(credset),
                                     lead.SNP=lead.snp,
                                     ID=data[, snp[match(lead.snp, snp)]],
                                     CHR=data[, chr[match(lead.snp, snp)]],
                                     POS=data[, pos[match(lead.snp, snp)]],
                                     PP4=round(unname(res$summary[paste0("PP.H4.abf")]), digits=3),
                                     SNP.PP4=res$results$SNP.PP.H4[o][1])
    return(list(result=res.dt, credset=credset.dt))
  }
  else{
    res.dt <- data.table::data.table(region=region,
                                     CHR=all.regions[region, CHR],
                                     start=all.regions[region, start],
                                     end=all.regions[region, end],
                                     PP3=round(unname(res$summary[paste0("PP.H3.abf")]), digits=3))
    return(res.dt)
  }
}

get.ld <- function(credset, data, out.path, ld.path, range=5e5){
  for (reg in as.integer(unique(credset$region))){
    region.data <- data[chr==credset[region==reg]$CHR[1] & dplyr::between(pos, credset[region==reg]$POS[1] - range, credset[region==reg]$POS[1] + range) & snp.x!=".", snp.x]
    chr <- data[snp.x %in% region.data, unique(chr)]
    fwrite(as.list(region.data), paste0(ld.path, "region", reg, ".chr", chr), sep="\n")
  }

  ################### PLINK ####################
  # REGION=(6 108 161 171 248 264 308)
  # CHR=(1 6 9 10 12 14 17)
  # REGION=(388 554 588 635)
  # CHR=(7 12 17 6)
  # 
  # for index in {0..3}
  # do
  # plink2
  #    --threads 10
  #    --bgen /lustre/groups/itg/shared/referenceData/ukbiobank/chip/imputed/ukb_imp_chr"${CHR[${index}]}"_v3.bgen ref-first
  #    --sample /lustre/groups/itg/shared/referenceData/ukbiobank/chip/imputed/ukb1020_imp_chr1_v2_s487406.sample
  #    --extract /lustre/groups/itg/teams/zeggini/projects/SCZ_T2D/GWASColoc/LDvariants/region${REGION[${index}]}.chr${CHR[${index}]}
  #    --make-bed
  #    --out /lustre/groups/itg/teams/zeggini/projects/SCZ_T2D/GWASColoc/LDvariants/region${REGION[${index}]}.chr${CHR[${index}]}_plink
  # done
  # 
  # for i in {0..3}
  # do
  # plink
  #    --threads 10
  #    --bfile region${REGION[$i]}.chr${CHR[$i]}_plink
  #    --r2 inter-chr
  #    --ld-window-r2 0
  #    --out region${REGION[$i]}_LD &
  # done &
}

plot.coloc <- function(data, traits, res.coloc, credset, region, genes.data, coloc.path, plot.path, ld.path, out.file.suffix="", pp4=TRUE, range=5e+5){
  source(paste0(project_folder, "scripts/plot_functions.R"))
  GRCh37_Genes <- read.delim(genes.data, stringsAsFactors = FALSE, header = TRUE)
  
  # Get data
  data=unique(data, by="snp")
  trait1.region <- data[, .(CHR=chr, SNP=snp, P=get(paste0("pval.", traits[1])), BP=pos, logP=-log10(as.numeric(get(paste0("pval.", traits[1])))))]
  trait2.region <- data[, .(CHR=chr, SNP=snp, P=get(paste0("pval.", traits[2])), BP=pos, logP=-log10(as.numeric(get(paste0("pval.", traits[2])))))]
  data.lst <- list(trait1.region, trait2.region)
  names(data.lst) <-  c(traits[1], traits[2])
  
  if(pp4){
    sel.snp <- res.coloc$lead.SNP  # rsid
    sel.chr <- res.coloc$CHR
    sel.pos <- res.coloc$POS
    sel.PP4 <- res.coloc$PP4
    
    # Get LD matrix
    if(!file.exists(paste0(ld.path, "region", region, "_LD.ld"))){
      print("Calculating LD with plink")
      ld.file <- NULL
      system(paste0("plink --bim /lustre/groups/itg/teams/data-informatics/projects/ukbb/genotypes/imputed/EGAD00010001474/GO2_var_ids/ukbb.imputed.v3.chr", sel.chr, ".bim.orig --bed /lustre/groups/itg/teams/data-informatics/projects/ukbb/genotypes/imputed/EGAD00010001474/GO2_var_ids/ukbb.imputed.v3.chr", sel.chr, ".bed --fam /lustre/groups/itg/teams/data-informatics/projects/ukbb/genotypes/imputed/EGAD00010001474/GO2_var_ids/ukbb.imputed.v3.chr", sel.chr, ".fam --from ",  data[pos==min(data$pos), snp], " --to ", data[pos==max(data$pos), snp], " --threads 10 --r2 inter-chr --ld-snp ", sel.snp, "  --ld-window-r2 0 --out ", ld.path, "region", region, "_LD"))
    }
    ld.file <- data.table::fread(paste0(ld.path, "region", region, "_LD.ld"))[SNP_A==sel.snp | SNP_B==sel.snp]
    # ld.file <- data.table::fread(paste0(ld.path, "region108_LD.ld"))[SNP_A==sel.snp | SNP_B==sel.snp]
    ld.file[SNP_B==sel.snp, c("SNP_B", "SNP_A") := .(SNP_A, SNP_B)]
  
    # Regional association plot
    locus.zoom(data = data.lst,
               offset_bp = range,
               genes.data = GRCh37_Genes,
               file.name = paste0(plot.path, paste(traits[1], traits[2], region, sep="_"), out.file.suffix, ".png"),
               # secondary.snp = ifelse(credible.set<2, NA, credible.set),
               snp=sel.snp,
               ignore.lead=TRUE,
               ld.file=ld.file,
               pp="PP4",
               pp.value=round(unname(sel.PP4), digits=3),
               nplots=TRUE)
    
    # Plot PP4 of credible set
    locus.zoom(data = credset[, .(CHR, SNP=rsID, PP4=SNP.PP4, BP=POS, logP=SNP.PP4)],
               region = c(sel.chr, credset[, min(POS)]-100,credset[, max(POS)]+100),
               genes.data = GRCh37_Genes,
               file.name = paste0(plot.path, paste(traits[1], traits[2], "PP4", region, sep="_"), out.file.suffix, ".png"),
               snp=sel.snp,
               ignore.lead=TRUE,
               ld.file=ld.file,
               sig.type="PP4")
  }
  else{
    # Regional association plot without ld colouring
    locus.zoom(data = data.lst,
               region = c(res.coloc$CHR, res.coloc$start, res.coloc$end),
               offset_bp = 0,
               genes.data = GRCh37_Genes,
               file.name = paste0(plot.path, paste(traits[1], traits[2], region, "PP3", out.file.suffix, ".png", sep="_")),
               pp="PP3",
               pp.value=round(unname(res.coloc$PP3), digits=3),
               nplots=TRUE)
  }
}

#' Perform pairwise colocalization analysis on multiple loci using coloc.abf
#'
#' @param regions regions to perform colocalization (columns: CHR, start end)
#' @param trait1 summary stats of trait1 as a data.table (essential columns: )
#' @param trait2 summary stats of trait2 as a data.table (essential columns: )
#' @param cc_ratio1 case vs. control ratio for trait1
#' @param cc_ratio2 case vs. control ratio for trait2
#' @param pp4.thres PP4 threshold for defining colocalized regions (default=0.8)
#' @return result of colocalization for each region along with 95% credible set for regions that reach PP4 threshold
#' @export
gwas.coloc <- function(traits, N_trait1, N_trait2, Ncases_trait1, Ncases_trait2, out.path, plot.path, ld.path, genes.data, pp4.thres=0.8, out.file.suffix=""){
  # regions <- data.table::fread(paste0(out.path, "GWAS_regions", out.file.suffix, ".csv"))[CHR!=23]
  regions <- data.table::fread(paste0(project_folder, "GWASColoc/GWAS_regions", out.file.suffix, ".csv"))[CHR!=23]
  data <- merge(data.table::fread(paste0(out.path, "GWAS_", traits[1], "_precoloc", out.file.suffix, ".txt")),
                data.table::fread(paste0(out.path, "GWAS_", traits[2], "_precoloc", out.file.suffix, ".txt")), 
                by=c("ID", "chr", "pos", "ea", "nea"), allow.cartesian = T, suffixes=paste0(".", traits))
  data[, snp:=ifelse(eval(paste0("snp.", traits[1]))=="", get(paste0("snp.", traits[2])), get(paste0("snp.", traits[1]))), by=seq_len(nrow(data))]
  data[, eval(paste0("snp.", traits[1])):=NULL]
  data[, eval(paste0("snp.", traits[2])):=NULL]
  
  print("Colocalization analysis using coloc.abf")
  credset.dt <- data.table::data.table()
  res.pp4.dt <- data.table::data.table()
  # res.pp3.dt <- data.table::data.table()
  
  coloc.regions <- c(6,108, 161,171,248,264,308,388,554,588,635)
  coloc.regions <- c(6,171,308,388)
  
  # for (i in 1:nrow(regions)){
  for (i in coloc.regions){
    data_tmp <- data[chr==regions[i]$CHR & dplyr::between(pos, regions[i]$start, regions[i]$end)]
    data_tmp <- data_tmp[!is.na(get(paste0("beta.", traits[1])))]
    data_tmp <- data_tmp[!is.na(get(paste0("beta.", traits[2])))]
    data_tmp <- unique(data_tmp, by="snp")
    
    res <- coloc::coloc.abf(dataset1=list(beta=data_tmp[, get(paste0("beta.", traits[1]))], varbeta=(data_tmp[, get(paste0("se.", traits[1]))])^2, type="cc", s=Ncases_trait1/N_trait1, N=N_trait1, snp=data_tmp$snp, position=data_tmp$pos),
                            dataset2=list(beta=data_tmp[, get(paste0("beta.", traits[2]))], varbeta=(data_tmp[, get(paste0("se.", traits[2]))])^2, type="cc", s=Ncases_trait2/N_trait2, N=N_trait2, snp=data_tmp$snp, position=data_tmp$pos))
    
   if (!is.na(res$summary[6]) & res$summary[6]>pp4.thres){
      output.coloc <- output.coloc.result(data=data_tmp, res=res, region=i, all.regions=regions, pp4=TRUE)
      # Get 95% credible set for regions where PP4 > 0.8
      credset.dt <- rbind(credset.dt, output.coloc$credset)
      res.pp4.dt <- rbind(res.pp4.dt, output.coloc$result)
      
      # plot with LD colouring
      print("Plotting colocalization result")
      plot.coloc(data=data_tmp, traits=traits, res.coloc=output.coloc$result, credset=output.coloc$credset, region=i, genes.data, coloc.path,
                 plot.path, ld.path, out.file.suffix)
    }
    # else if (!is.na(res$summary[5]) & res$summary[5]>pp4.thres){
    #   output.coloc <- output.coloc.result(data=data_tmp, res=res, region=i, all.regions=regions, pp4=FALSE)
    #   res.pp3.dt <- rbind(res.pp3.dt, output.coloc)
    #   
    #   # plot without LD colouring
    #   plot.coloc(data=data_tmp, traits=traits, res.coloc=output.coloc, credset=NA, region=i, genes.data, coloc.path, 
    #              plot.path=paste0(plot.path, "PP3/"), ld.path, out.file.suffix, pp4=FALSE)
    # }
  }
  
  res.pp4.dt <- unique(res.pp4.dt, by=colnames(res.pp4.dt)[2:6])
  data.table::fwrite(credset.dt[region %in% res.pp4.dt$region], paste0(out.path, "final_coloc_credible_set", out.file.suffix, ".csv"))
  data.table::fwrite(res.pp4.dt, paste0(out.path, "final_coloc_result_pp4", out.file.suffix, ".csv"))
  # data.table::fwrite(res.pp3.dt[, data.table::as.data.table(IRanges::reduce(IRanges::IRanges(start, end), min.gapwidth = 0L)), CHR], paste0(out.path, "final_coloc_result_pp3", out.file.suffix, ".csv"))
}

overlap_ensembl <- function(chr, start, end){
  server <- "https://grch37.rest.ensembl.org"
  ext <- paste0("/overlap/region/human/", chr, ":", start, "-", end, "?feature=gene")
  r <- httr::GET(paste(server, ext, sep = ""), httr::content_type("application/json"))
  httr::stop_for_status(r)
  
  gene.dt <- data.table::as.data.table(jsonlite::fromJSON(jsonlite::toJSON(httr::content(r))))[, .(name=as.character(external_name), id=as.character(id))]
  
  return(gene.dt)
}

gene.overlap <- function(out.path, genes.path, out.file.suffix="", wd=1e+06){
  genes.dt <- data.table::data.table(name=character(), id=character())
  dt <- data.table::fread(paste0(out.path, "final_coloc_result_pp4", out.file.suffix, ".csv"))
  dt[, `:=` (start=ifelse(POS-wd<0, 0, POS-wd), end=POS+wd)]
  
  tmp <- dt[, overlap_ensembl(CHR, start, end), by=seq_len(nrow(dt))]
  genes.dt <- rbind(genes.dt, tmp, fill=TRUE)
  genes.dt <- unique(genes.dt, by=c('name', 'id'))
  
  data.table::fwrite(genes.dt, genes.path)
}

################################################################################
#---------------------------------- Main --------------------------------------#
################################################################################
project_folder <- "/lustre/groups/itg/teams/zeggini/projects/SCZ_T2D/"
source(paste0(project_folder, "scripts/read_files_config.R"))
preprocess.gwas(signals.lst=list(t1.indep.signals, t2.indep.signals), out.path=output.path, traits=traits)
gwas.coloc(traits, t1.samplesize, t2.samplesize, t1.ncases, t2.ncases, output.path, plot.path=coloc.plots.path, ld.path, genes.data)
gene.overlap(out.path=output.path, genes.path=all.genes)

credset.dt <- data.table::fread(paste0(output.path, "credible_set.csv"))
credset.dt <- credset.dt[order(region, -SNP.PP4)]
data.table::fwrite(credset.dt, paste0(output.path, "credible_set.csv"))

# Sensitivity analysis: Europeans only
t1.filepath <- paste0(data.path, traits[1], "_european.txt")
t2.filepath <- paste0(data.path, traits[2], "_european.txt")
preprocess.gwas(signals.lst=list(t1.indep.signals, t2.indep.signals), 
                out.path=paste0(project_folder, "GWASColoc/european/"), traits=traits)
gwas.coloc(traits, t1.samplesize, t2.samplesize, t1.ncases, t2.ncases, 
           out.path=paste0(project_folder, "GWASColoc/european/"), 
           plot.path=paste0(output.path, "plots/european/"), ld.path, genes.data)
N_trait1=t1.samplesize
N_trait2=t2.samplesize
Ncases_trait1=t1.ncases
Ncases_trait2=t2.ncases
out.path=paste0(project_folder, "GWASColoc/european/")
plot.path=paste0(output.path, "plots/european/")
pp4.thres=0.8
out.file.suffix=""
