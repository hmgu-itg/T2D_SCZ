#############################################################################################
#---------------------------------- Manual preprocess --------------------------------------#
#############################################################################################
scz <- data.table::fread("/project_data/data/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv.gz")
scz <- scz[, `:=`(N=NCAS+NCON, EAF=((FCAS*NCAS)+(FCON*NCON))/(NCAS+NCON), rsID=ID)]
scz <- scz[, MAF:=ifelse(EAF<0.5, EAF, (1-EAF))]
scz <- scz[, ID:=paste(paste(CHROM, POS, sep=":"), sort_alleles(A1, A2), sep="_")]
scz <- scz[, .(CHR=CHROM, POS, rsID, pval=PVAL, beta=BETA, se=SE, MAF, EAF, N, Ncases=NCAS, EA=A1, NEA=A2, ID)]
scz <- scz[, CHR:=as.integer(CHR)]
scz <- scz[nchar(CHR)<3]
data.table::fwrite(scz, "/project_data/data/SCZ_primary.txt")

t2d <- data.table::fread('unzip -p /project_data/data/Mahajan.NatGen2022.DIAMANTE-TA.sumstat.zip')
t2d <- t2d[, .(CHR=`chromosome(b37)`, POS=`position(b37)`, rsID, pval=`MR-MEGA_p-value_association`, beta=`Fixed-effects_beta`,
               se=`Fixed-effects_SE`, MAF=NA, EAF=NA, N=1339889, Ncases=180834, EA=toupper(effect_allele), NEA=toupper(other_allele))]
t2d[, ID:=paste(paste(CHR, POS, sep=":"), sort_alleles(EA, NEA), sep="_")]
data.table::fwrite(t2d, "/project_data/data/T2D_multiancestry.txt")

t2d.indep.signals <- data.table::as.data.table(readxl::read_xlsx("/project_data/data/T2D_indep_signals.xlsx", skip=3)) %>%
  .[, .(CHR=`...3`, POS=`...4`, rsID=`...2`, ID=paste(paste(`...3`, `...4`, sep=":"), sort_alleles(c(Risk, Other)), sep="_"))]
data.table::fwrite(t2d.indep.signals, "/project_data/data/T2D_indep_signals.txt")

t2d.indep.signals.5e8 <- data.table::as.data.table(readxl::read_xlsx("data/41588_2022_1058_MOESM3_ESM.xlsx", skip=2, sheet=6)) %>%
  .[, .(CHR=Chr, POS=`Position (bp, b37)`, rsID=`Lead SNV`)]
data.table::fwrite(t2d.indep.signals.5e8, "data/T2D_indep_signals_5e-8.txt")

# merge with T2D multiancestry summary stats to get allele info for ID
indep.signals <- data.table::fread("/project_data/data/T2D_indep_signals_5e-8.txt")
t2d <- data.table::fread("/project_data/data/T2D_multiancestry.txt")
dt <- merge(indep.signals, t2d[, .(rsID, ID)], all.x=TRUE, by="rsID")
data.table::fwrite(dt, "/project_data/data/T2D_indep_signals_5e-8.txt")

scz.indep.signals <- data.table::as.data.table(readxl::read_xls("/project_data/data/SCZ_indep_signals.xls")) %>%
  .[, c("A1", "A2") := tstrsplit(A1A2, "/", fixed=TRUE)] %>%
  .[P<5e-8, .(CHR, POS=BP, rsID=SNP, ID=paste(paste(CHR, BP, sep=":"), sort_alleles(A1, A2), sep="_"))]
data.table::fwrite(scz.indep.signals, "/project_data/data/SCZ_indep_signals.txt")

t2d.eur <- data.table::fread('/project_data/data/DIAMANTE-EUR.sumstat.txt')
t2d.eur <- t2d.eur[, .(CHR=`chromosome(b37)`, POS=`position(b37)`, rsID, pval=`Fixed-effects_p-value`, beta=`Fixed-effects_beta`,
                       se=`Fixed-effects_SE`, MAF=ifelse(`effect_allele_frequency`<0.5, `effect_allele_frequency`, 1-`effect_allele_frequency`),
                       EAF=`effect_allele_frequency`, N=933970, Ncases=80154, EA=toupper(effect_allele), NEA=toupper(other_allele))]
# t2d.eur[, ID:=paste(paste(CHR, POS, sep=":"), sort_alleles(EA, NEA), sep="_")]
data.table::fwrite(t2d.eur, "/project_data/data/T2D_european_noID.txt")

scz.eur <- data.table::fread("/project_data/data/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz")
scz.eur <- scz.eur[, `:=`(N=NCAS+NCON, EAF=((FCAS*NCAS)+(FCON*NCON))/(NCAS+NCON), rsID=ID)]
scz.eur <- scz.eur[, MAF:=ifelse(EAF<0.5, EAF, (1-EAF))]
# scz <- scz[, ID:=paste(paste(CHROM, POS, sep=":"), sort_alleles(A1, A2), sep="_")]
scz.eur <- scz.eur[, .(CHR=CHROM, POS, rsID, pval=PVAL, beta=BETA, se=SE, MAF, EAF, N, Ncases=NCAS, EA=A1, NEA=A2)]
scz.eur <- scz.eur[, CHR:=as.integer(CHR)]
scz.eur <- scz.eur[nchar(CHR)<3]
data.table::fwrite(scz.eur, "/project_data/data/SCZ_european_noid.txt")
