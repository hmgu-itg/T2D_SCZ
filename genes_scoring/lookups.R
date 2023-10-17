# https://github.com/nalu357/ComorbidityAnalysis/tree/master/R 

#' @title Look up genes in DEGs
#' @description Check if any gene in your list is a DEG for any trait
#'
#' @param DEG named list of data tables containing DEGs for a trait (columns: gene, ID, logFC, pval.adj)
#' @param genes data table of genes of interest (columns: name, ID)
#' @param file.path path to save result table, default="lookups/DEG.csv"
#' @param logFC.theta threshold for abs(logFC), default=log(1.5, base=2)
#' @param alpha significance threshold for p-value, default=0.05
#' @return a data table of genes with 1 or 0 in each DEG column
#' @export
#'
#' @importFrom data.table fwrite
#' @examples

#' @title Search OMIM
#' @description Get OMIM's entry for a gene
#'
#' @param gene gene symbol
#' @return data table of OMIM's entry for the gene
#' @export
#'
#' @importFrom httr GET content stop_for_status content_type
#' @examples
search.OMIM <- function(gene) {
  stopifnot(is.character(gene), length(gene)==1)
  gene <- tolower(gene)
  server <- "https://api.omim.org"
  r <- httr::GET(paste(server, "/api/geneMap/search?search=", gene, "&include=geneMap&apiKey=cBg3cWK0TxibU0wvtkwtMA&format=json&start=0&limit=100", sep = ""), httr::content_type("application/json"))
  httr::stop_for_status(r)
  if (length(jsonlite::fromJSON(jsonlite::toJSON(httr::content(r)))$omim$searchResponse$geneMapList$geneMap$phenotypeMapList) != 0) {
    all.phen <- jsonlite::fromJSON(jsonlite::toJSON(httr::content(r)))$omim$searchResponse$geneMapList$geneMap$phenotypeMapList[[1]]$phenotypeMap$phenotype
    phen <- all.phen[sapply(all.phen, function(i) {
      !grepl("\\{|\\[|\\?", i)
    })]
    if (length(phen) == 0) {
      return(NA)
    } else {
      return(unlist(phen))
    }
  }
  else {
    return(NA)
  }
}

#' @title Look up genes in OMIM
#' @description if any gene in your list has an entry in OMIM showing relevant phenotypes for a trait of interest
#'
#' @param terms named list of terms related to each trait of interest
#' @param genes data table of genes of interest (columns: name, ID)
#' @param file.path path to save result table, default="lookups/OMIM.csv"
#' @return a data table of genes with 1 or 0 for each trait column
#' @export
#'
#' @importFrom data.table fwrite
#' @examples
lookup.OMIM <- function(terms.lst, genes.path, file.path) {
  genes <- fread(genes.path)
  genes <- unique(genes, by="id")
  
  stopifnot(length(terms.lst)>=1, is.data.frame(genes), nrow(genes)>0)
  genes$OMIM <- lapply(genes$name, search.OMIM)
  genes2 <- genes[!is.na(OMIM)]
  genes2 <- genes2[,.(OMIM = unlist(OMIM)), by = setdiff(names(genes), 'OMIM')]
  
  for (trait in names(terms.lst)) {
    pattern <- paste(terms.lst[[trait]], collapse = "|")
    genes <- genes[, eval(trait) := as.integer(grepl(pattern, OMIM, ignore.case = TRUE)), by = name]
  }
  
  # data.table::fwrite(genes[, OMIM := NULL], file.path)
  
  # Check terms manually for phenotypes related to T2D or SCZ
  data.table::fwrite(list(unique(genes2$OMIM)), "/project_data/omim_terms.txt")
  
  genes2 <- genes[!is.na(OMIM)]
  genes2 <- genes2[,.(OMIM = unlist(OMIM)), by = setdiff(names(genes), 'OMIM')]
  data.table::fwrite(genes2, file.path)
  return(genes2)
  
}

# https://github.com/nalu357/ComorbidityAnalysis/tree/master/R 

#' @title Look up genes in DEGs
#' @description Check if any gene in your list is a DEG for any trait
#'
#' @param DEG named list of data tables containing DEGs for a trait (columns: gene, ID, logFC, pval.adj)
#' @param genes data table of genes of interest (columns: name, ID)
#' @param file.path path to save result table, default="lookups/DEG.csv"
#' @param logFC.theta threshold for abs(logFC), default=log(1.5, base=2)
#' @param alpha significance threshold for p-value, default=0.05
#' @return a data table of genes with 1 or 0 in each DEG column
#' @export
#'
#' @importFrom data.table fwrite
#' @examples

#' @title Search OMIM
#' @description Get OMIM's entry for a gene
#'
#' @param gene gene symbol
#' @return data table of OMIM's entry for the gene
#' @export
#'
#' @importFrom httr GET content stop_for_status content_type
#' @examples

#' @title Look up genes in OMIM
#' @description if any gene in your list has an entry in OMIM showing relevant phenotypes for a trait of interest
#'
#' @param terms named list of terms related to each trait of interest
#' @param genes data table of genes of interest (columns: name, ID)
#' @param file.path path to save result table, default="lookups/OMIM.csv"
#' @return a data table of genes with 1 or 0 for each trait column
#' @export
#'
#' @importFrom data.table fwrite
#' @examples

# library(xml2)
# library(dplyr)
# library(openxlsx)

.simpleCap <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

search.IMPC <- function(gene){
  server <- "https://www.ebi.ac.uk"
  ext <- paste0("/mi/impc/solr/genotype-phenotype/select?q=marker_symbol:", gene, "&rows=500")
  r <- httr::GET(paste(server, ext, sep = ""), httr::content_type("application/json"))
  httr::stop_for_status(r)
  
  IMPC.dt <- data.table::as.data.table(jsonlite::fromJSON(jsonlite::toJSON(httr::content(r)))$response$docs) 
  if(length(IMPC.dt)!=0){
    # IMPC.dt <- IMPC.dt[, .(phenotype=mp_term_name, effect=effect_size, pval=p_value, class=parameter_name, source="IMPC")]
    IMPC.dt <- IMPC.dt[, .(phenotype=as.character(mp_term_name), source="IMPC")]
    return(IMPC.dt)
  }
  else{
    return(IMPC.dt)
  }
}

mgi.ret <- function(gene, mgi.path){
  # mgi <- as.data.table(read.xlsx(mgi.path)) %>% .[, .(Input, Term)]
  mgi <- data.table::fread(mgi.path)[, .(Input, Term)]
  if (is.na(mgi[Input==gene][1]$Term)) {
    return(data.table::data.table())
  }
  else{
    return(mgi[Input==gene, .(phenotype=Term, source="MGI")])
  }
}

search.RGD <- function(gene, specie=2, src="MP"){
  # species: 1=human, 2=mouse, 3=rat
  # source: mouse phenotype=MP, human pheotype=HP, diseases=DOID
  
  gene <- tolower(gene)
  
  server <- "https://rest.rgd.mcw.edu"
  r <- httr::GET(paste(server, "/rgdws/genes/", gene, "/", specie, sep = ""), httr::content_type("application/json"))
  httr::stop_for_status(r)
  if (length(jsonlite::fromJSON(jsonlite::toJSON(httr::content(r))))>1){
    rgdID <- jsonlite::fromJSON(jsonlite::toJSON(httr::content(r)))$rgdId
    
    ext <- paste0("/rgdws/annotations/rgdId/", rgdID, "/", src)
    r <- httr::GET(paste(server, ext, sep = ""), httr::content_type("application/json"))
    httr::stop_for_status(r)
    
    if(length(jsonlite::fromJSON(jsonlite::toJSON(httr::content(r))))!=0){
      phen <- jsonlite::fromJSON(jsonlite::toJSON(httr::content(r)))$term
      RGD.dt <- data.table::data.table(phenotype=as.character(phen), source="RGD")
      return(RGD.dt)
    }
    else return(data.table::data.table())
  }
  else return(data.table::data.table())
}

knockout.lookup <- function(data, sources, terms, trait){
  result <- data.table::data.table(t(rep(NA, length(sources))))
  colnames(result) <- sources
  if(nrow(data)!=0){
    pattern <- paste(terms, collapse="|")
    tmp <- data[, grepl(pattern, phenotype, ignore.case=TRUE), by=.(source, seq_len(nrow(data)))]
    for (src in sources) {
      if(nrow(data[source==src])>0) result[, eval(src):=0]
      if(sum(tmp[source==src]$V1)>0) result <- result[, eval(src):=1]
    }
    result[, N:=rowSums(.SD, na.rm=TRUE), .SDcols=1:3]
  }
  colnames(result) <- paste(colnames(result), trait, sep="_")
  return(result)
}

lookup.KOmice <- function(genes.path, traits, terms.lst, ko.path, mgi.result){
  genes.dt <- data.table::fread(genes.path)
  genes.dt <- unique(genes.dt, by="id")
  
  #------------------- add IMPC ----------------------#
  impc.lst <- lapply(.simpleCap(genes.dt$name), search.IMPC)
  names(impc.lst) <- genes.dt$name
  
  #------------------- add MGI ----------------------#
  # http://www.informatics.jax.org/batch/summary
  mgi.lst <- lapply(genes.dt$name, mgi.ret, mgi.path=mgi.result)
  names(mgi.lst) <- genes.dt$name
  
  #------------------- add RGD ----------------------#
  rgd.lst <- lapply(genes.dt$name, search.RGD)
  names(rgd.lst) <- genes.dt$name
  
  #------------------ Combine all ------------------ #
  knockout.lst <- Map(rbind, Map(rbind, impc.lst, mgi.lst), rgd.lst)
  names(knockout.lst) <- genes.dt$name
  
  #------------- check entries to expand terms ----------------#
  # t <- unique(unlist(lapply(knockout.lst, function(i){terms <- list(); if(nrow(i)>1) terms <- c(terms, i$phenotype); return(unlist(terms))})))
  # pattern <- paste(unlist(terms.lst), collapse="|")
  # t[grepl(pattern, t, ignore.case=TRUE)]
  # t2 <- t[!grepl(pattern, t, ignore.case=TRUE)]
  
  #------------------- check for OA and T2D terms ----------------------#
  t1.result <- lapply(knockout.lst, knockout.lookup, sources=c("IMPC", "MGI", "RGD"), terms=terms.lst[[1]], trait=traits[1])
  t2.result <- lapply(knockout.lst, knockout.lookup, sources=c("IMPC", "MGI", "RGD"), terms=terms.lst[[2]], trait=traits[2])
  combined.result <- Map(cbind, t1.result, t2.result)
  
  #---------------- output result -------------------#
  knockout.table <- cbind(data.table::data.table(gene=names(combined.result)), data.table::rbindlist(combined.result, fill=TRUE)) %>% 
    .[, T2D:=ifelse(N_T2D>0,1,N_T2D)] %>%
    .[, SCZ:=ifelse(N_SCZ>0,1,N_SCZ)]
    # .[, T2D:=ifelse(rowSums(.SD, na.rm=TRUE)>0,1,rowSums(.SD)), .SDcols = 2:4] %>% 
    # .[, SCZ:=ifelse(rowSums(.SD, na.rm=TRUE)>0,1,rowSums(.SD)), .SDcols = 6:8]
  data.table::fwrite(knockout.table, ko.path)
  return(knockout.table)
}

lookup.missense <- function(snp.list, missense.path){
  snpMart = biomaRt::useEnsembl(biomart = "snps", 
                                dataset = "hsapiens_snp",
                                GRCh=37)
  missense.dt <- biomaRt::getBM(attributes = c('refsnp_id', 'chr_name', 'ensembl_gene_name', 'consequence_type_tv', 'consequence_allele_string', 'sift_prediction', 'polyphen_prediction'),
                                filters = 'snp_filter', 
                                values = snp.list, 
                                mart = snpMart) %>% as.data.table() %>% .[consequence_type_tv=="missense_variant"]
  data.table::fwrite(missense.dt, missense.path)
  return(missense.dt)
}

lookup.DEG <- function(deg.list, deg.path, genes.path, traits, gene.id.col, gene.name.col){
  genes.dt <- data.table::fread(genes.path)
  genes.dt <- unique(genes.dt, by="id")
  res.dt <- data.table::data.table(gene=genes.dt$name, ID=genes.dt$id)
  
  for (i in 1:length(deg.list)){
    deg <- data.table::fread(deg.list[i])
    tmp <- deg[get(gene.id.col[i]) %in% genes.dt$id] %>% .[, .(gene=get(gene.name.col[i]), ID=get(gene.id.col[i]))]
    res.dt[, eval(traits[i]):=ifelse(ID %in% tmp$ID, 1, 0)]
  }
  
  data.table::fwrite(res.dt, deg.path)
  return(res.dt)
}

preprocess.hc.genes <- function(trait){
  if (trait=="T2D"){
    trait="T2D"
    dt <- data.table::fread("/project_data/data/t2d_hc.txt") 
    ensembl = biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
    dt <- data.table::as.data.table(biomaRt::getBM(attributes=c('external_gene_name', 'ensembl_gene_id'), 
                                                   filters = 'external_gene_name', 
                                                   values = dt$gene.name, 
                                                   mart = ensembl))
    colnames(dt) <- c("gene.name", "gene.id")
    fwrite(unique(dt, by = "gene.name"), paste0("/project_data/data/", trait, "_hc.txt"))
  }
  else if (trait=="SCZ"){
    trait="SCZ"
    dt <- data.table::fread("/project_data/data/INT-17_SCZ_High_Confidence_Gene_List.csv")
    dt <- dt[, .(gene.name=sczgenenames, gene.id=ensembl_names)]
    fwrite(unique(dt, by = "gene.name"), paste0("/project_data/data/", trait, "_hc.txt"))
  }
  else print(paste0("Trait ", trait, " not supported"))
}

lookup.HC <- function(hc.list, hc.path, genes.path, traits, gene.id.col, gene.name.col){
  genes.dt <- fread(genes.path)
  genes.dt <- unique(genes.dt, by="id")
  res.dt <- data.table::data.table(gene=genes.dt$name, ID=genes.dt$id)
  
  for (i in 1:length(hc.list)){
    hc <- data.table::fread(hc.list[i])
    tmp <- hc[get(gene.id.col[i]) %in% genes.dt$id] %>% .[, .(gene=get(gene.name.col[i]), ID=get(gene.id.col[i]))]
    res.dt[, eval(traits[i]):=ifelse(ID %in% tmp$ID, 1, 0)]
  }
  
  data.table::fwrite(res.dt, hc.path)
  return(res.dt)
}


################################################################################
#---------------------------------- Main --------------------------------------#
################################################################################
source("/project_data/scripts/read_files_config.R")
credset <- data.table::fread(paste0(output.path, "credible_set.csv"))

omim <- lookup.OMIM(terms.lst=list(t1=t1.terms, t2=t2.terms), genes.path=all.genes, omim.path)
ko.mice <- lookup.KOmice(genes.path=all.genes, traits, terms.lst=list(t1.terms, t2.terms), ko.path, mgi.result)
missense <- lookup.missense(snp.list=credset$snp, missense.path)
deg <- lookup.DEG(deg.list=deg.list, genes.path=all.genes, deg.path=deg.path, traits=traits, 
                  gene.id.col=c("Ensembl gene", "Ensembl_Name"), gene.name.col=c("Gene name", "Gene_Name"))
hc <- lookup.HC(hc.list=hc.list, hc.path=hc.path, genes.path=all.genes, traits=traits, 
                gene.id.col=c("gene.id", "gene.id"), gene.name.col=c("gene.name", "gene.name"))


