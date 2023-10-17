DATA_DIR="/lustre/groups/itg/teams/zeggini/projects/SCZ_T2D"

#####################################################################################
#--------------------------------- Format data -------------------------------------#
#####################################################################################
# ------------------- SCZ ------------------- 
for i in primary european asian
do
    case $i in
        primary)
            j="multiancestry"
            n_cases=74776
            n=175799
            ;;
        european)
            j="eur"
            n_cases=53386
            n=130644
            ;;
        asian)
            j="asian"
            n_cases=14004
            n=30761
            ;;
    esac
    echo "snpid chr bp a1 a2 beta pval" > ${DATA_DIR}/GeneticCorrelation/SCZ.$j.for-ldscr
    zcat ${DATA_DIR}/data/PGC3_SCZ_wave3.$i.autosome.public.v3.vcf.tsv.gz | grep -v "^#" | grep -v "POS" | awk 'OFS="\t"{print $2,$1,$3,$4,$5,$9,$11}' >> ${DATA_DIR}/GeneticCorrelation/SCZ.$j.for-ldscr
    /lustre/groups/itg/shared/software/ldsc/munge_sumstats.py --sumstats ${DATA_DIR}/GeneticCorrelation/SCZ.$j.for-ldscr --out ${DATA_DIR}/GeneticCorrelation/SCZ.$j.ldsc --merge-alleles  ${DATA_DIR}/GeneticCorrelation/w_hm3.snplist --N $n --N-cas ${n_cases}
done

# ------------------- T2D ------------------- 
for i in EUR EAS TA
do
    case $i in
        TA)
            j="multiancestry"
            n_cases=180834
            n=1339889
            ;;
        EUR)
            j="eur"
            n_cases=80154
            n=933970
            ;;
        EAS)
            j="asian"
            n_cases=72808
            n=332915
            ;;
    esac
    
    echo "snpid chr bp a1 a2 beta pval" > ${DATA_DIR}/GeneticCorrelation/T2D.$j.for-ldscr
    
    if [ $j == "multiancestry" ]
    then
        cat ${DATA_DIR}/data/DIAMANTE-$i.sumstat.txt | grep -v "rsID" | grep -v "nane-nan" | awk 'OFS="\t"{print $4,$1,$2,$5,$6,$10,$7}' >> ${DATA_DIR}/GeneticCorrelation/T2D.$j.for-ldscr
    else
        cat ${DATA_DIR}/data/DIAMANTE-$i.sumstat.txt | grep -v "rsID" | grep -v "nane-nan" | awk 'OFS="\t"{print $4,$1,$2,$5,$6,$8,$10}' >> ${DATA_DIR}/GeneticCorrelation/T2D.$j.for-ldscr
    fi
    
    /lustre/groups/itg/shared/software/ldsc/munge_sumstats.py --sumstats ${DATA_DIR}/GeneticCorrelation/T2D.$j.for-ldscr --out ${DATA_DIR}/GeneticCorrelation/T2D.$j.ldsc --merge-alleles  ${DATA_DIR}/GeneticCorrelation/w_hm3.snplist --N $n --N-cas ${n_cases}
done

# ------------------- MAGIC ------------------- 
for i in 2hGlu FG FI HbA1c
do
    for j in EUR  # EUR EAS 
    do
        case $j in
            TA)
                k="multiancestry"
                ;;
            EUR)
                k="eur"
                ;;
            EAS)
                k="asian"
                ;;
        esac
        
        echo "snpid chr bp a1 a2 beta pval N" > ${DATA_DIR}/GeneticCorrelation/${i}.${k}.for-ldscr
        
        if [ $k == "multiancestry" ]
        then
            # zcat ${DATA_DIR}/data/MAGIC1000G_${i}_${j}.tsv.gz | grep -v "variant" | awk 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$8,$7}' >> ${DATA_DIR}/GeneticCorrelation/${i}.${k}.for-ldscr
            echo "makes no sense beacuse there is no beta available for multiancestry"
        else
            zcat ${DATA_DIR}/data/MAGIC1000G_${i}_${j}.tsv.gz | grep -v "variant" | awk 'OFS="\t"{print $1,$2,$3,$4,$5,$7,$9,$10}' >> ${DATA_DIR}/GeneticCorrelation/${i}.${k}.for-ldscr
        fi
        
        /lustre/groups/itg/shared/software/ldsc/munge_sumstats.py --sumstats ${DATA_DIR}/GeneticCorrelation/${i}.${k}.for-ldscr --out ${DATA_DIR}/GeneticCorrelation/${i}.${k}.ldsc --merge-alleles  ${DATA_DIR}/GeneticCorrelation/w_hm3.snplist
    done
done

# ------------------- PGC ------------------- 
# for i in ADHD AnorexiaNervosa Austism BipolarDisorder Depression PanicDisorder PTSD_eur

echo "snpid chr bp a1 a2 or pval n_cases n_controls" > ${DATA_DIR}/GeneticCorrelation/ADHD.for-ldscr
zcat ${DATA_DIR}/data/PGC_ADHD.txt.gz | grep -v "A1" | awk 'OFS="\t"{print $2,$1,$3,$4,$5,$9,$11,$17,$18}' >> ${DATA_DIR}/GeneticCorrelation/ADHD.for-ldscr
/lustre/groups/itg/shared/software/ldsc/munge_sumstats.py --sumstats ${DATA_DIR}/GeneticCorrelation/ADHD.for-ldscr --out ${DATA_DIR}/GeneticCorrelation/ADHD.ldsc --merge-alleles ${DATA_DIR}/GeneticCorrelation/w_hm3.snplist

echo "snpid chr bp a1 a2 beta pval n_cases n_controls" > ${DATA_DIR}/GeneticCorrelation/AnorexiaNervosa.for-ldscr
zcat ${DATA_DIR}/data/PGC_AnorexiaNervosa.vcf.tsv.gz| grep -v "^#" | grep -v "POS" | awk 'OFS="\t"{print $3,$1,$2,$4,$5,$6,$8,$12,$13}' >> ${DATA_DIR}/GeneticCorrelation/AnorexiaNervosa.for-ldscr
/lustre/groups/itg/shared/software/ldsc/munge_sumstats.py --sumstats ${DATA_DIR}/GeneticCorrelation/AnorexiaNervosa.for-ldscr --out ${DATA_DIR}/GeneticCorrelation/AnorexiaNervosa.ldsc --merge-alleles ${DATA_DIR}/GeneticCorrelation/w_hm3.snplist

echo "snpid chr bp a1 a2 or pval" > ${DATA_DIR}/GeneticCorrelation/Autism.for-ldscr
zcat ${DATA_DIR}/data/PGC_Autism.txt.gz| grep -v "A1" | awk 'OFS="\t"{print $2,$1,$3,$4,$5,$7,$9}' >> ${DATA_DIR}/GeneticCorrelation/Autism.for-ldscr
/lustre/groups/itg/shared/software/ldsc/munge_sumstats.py --sumstats ${DATA_DIR}/GeneticCorrelation/Autism.for-ldscr --out ${DATA_DIR}/GeneticCorrelation/Autism.ldsc --merge-alleles ${DATA_DIR}/GeneticCorrelation/w_hm3.snplist --N-con 27969 --N-cas 18381

echo "snpid chr bp a1 a2 beta pval n_cases n_controls" > ${DATA_DIR}/GeneticCorrelation/BipolarDisorder.for-ldscr
zcat ${DATA_DIR}/data/PGC_BipolarDisorder.vcf.tsv.gz| grep -v "^#" | awk 'OFS="\t"{print $3,$1,$2,$4,$5,$6,$8,$14,$15}' >> ${DATA_DIR}/GeneticCorrelation/BipolarDisorder.for-ldscr
/lustre/groups/itg/shared/software/ldsc/munge_sumstats.py --sumstats ${DATA_DIR}/GeneticCorrelation/BipolarDisorder.for-ldscr --out ${DATA_DIR}/GeneticCorrelation/BipolarDisorder.ldsc --merge-alleles ${DATA_DIR}/GeneticCorrelation/w_hm3.snplist

echo "snpid a1 a2 beta pval" > ${DATA_DIR}/GeneticCorrelation/Depression.for-ldscr
cat ${DATA_DIR}/data/PGC_Depression.txt | grep -v "A1" | awk 'OFS="\t"{print $1,$2,$3,$5,$7}' >> ${DATA_DIR}/GeneticCorrelation/Depression.for-ldscr
/lustre/groups/itg/shared/software/ldsc/munge_sumstats.py --sumstats ${DATA_DIR}/GeneticCorrelation/Depression.for-ldscr --out ${DATA_DIR}/GeneticCorrelation/Depression.ldsc --merge-alleles ${DATA_DIR}/GeneticCorrelation/w_hm3.snplist --N-con 329443 --N-cas 170756

echo "snpid chr bp a1 a2 beta pval n_cases n_controls" > ${DATA_DIR}/GeneticCorrelation/PanicDisorder.for-ldscr
zcat ${DATA_DIR}/data/PGC_PanicDisorder.vcf.tsv.gz| grep -v "^#" | awk 'OFS="\t"{print $3,$1,$2,$4,$5,$6,$8,$14,$15}' >> ${DATA_DIR}/GeneticCorrelation/PanicDisorder.for-ldscr
/lustre/groups/itg/shared/software/ldsc/munge_sumstats.py --sumstats ${DATA_DIR}/GeneticCorrelation/PanicDisorder.for-ldscr --out ${DATA_DIR}/GeneticCorrelation/PanicDisorder.ldsc --merge-alleles ${DATA_DIR}/GeneticCorrelation/w_hm3.snplist

echo "snpid chr bp a1 a2 or pval n_cases n_controls" > ${DATA_DIR}/GeneticCorrelation/PTSD_eur.for-ldscr
zcat ${DATA_DIR}/data/PGC_PTSD_eur.txt.gz| grep -v "A1" | awk 'OFS="\t"{print $2,$1,$3,$4,$5,$9,$11,$17,$18}' >> ${DATA_DIR}/GeneticCorrelation/PTSD_eur.for-ldscr
/lustre/groups/itg/shared/software/ldsc/munge_sumstats.py --sumstats ${DATA_DIR}/GeneticCorrelation/PTSD_eur.for-ldscr --out ${DATA_DIR}/GeneticCorrelation/PTSD_eur.ldsc --merge-alleles ${DATA_DIR}/GeneticCorrelation/w_hm3.snplist

# ------------------- adiposity ------------------- 
for i in bmi whr
do
    echo "snpid chr bp a1 a2 beta pval N" > ${DATA_DIR}/GeneticCorrelation/${i}.for-ldscr
    # zcat ${DATA_DIR}/data/${i}.giant-ukbb.meta-analysis.combined.23May2018.txt.gz | grep -v "CHR" | awk 'OFS="\t"{print $3,$1,$2,$4,$5,$7,$9,$10}' >> ${DATA_DIR}/GeneticCorrelation/${i}.for-ldscr
    cat ${DATA_DIR}/data/${i}.rsid.txt | grep -v "CHR" | awk 'OFS="\t"{print $3,$1,$2,$4,$5,$7,$9,$10}' >> ${DATA_DIR}/GeneticCorrelation/${i}.for-ldscr
    /lustre/groups/itg/shared/software/ldsc/munge_sumstats.py --sumstats ${DATA_DIR}/GeneticCorrelation/${i}.for-ldscr --out ${DATA_DIR}/GeneticCorrelation/${i}.ldsc --merge-alleles  ${DATA_DIR}/GeneticCorrelation/w_hm3.snplist
done

for i in whole_body_fat_mass body_fat_percentage
do
    echo "snpid chr bp a1 a2 beta pval N" > ${DATA_DIR}/GeneticCorrelation/${i}.for-ldscr
    zcat ${DATA_DIR}/data/${i}_neale_rsid.txt.gz | grep -v "CHR" | awk '{print $15,$1,$2,$4,$14,$10,$13,$7}' FS=',' OFS='\t' >> ${DATA_DIR}/GeneticCorrelation/${i}.for-ldscr
    /lustre/groups/itg/shared/software/ldsc/munge_sumstats.py --sumstats ${DATA_DIR}/GeneticCorrelation/${i}.for-ldscr --out ${DATA_DIR}/GeneticCorrelation/${i}.ldsc --merge-alleles  ${DATA_DIR}/GeneticCorrelation/w_hm3.snplist
done

#####################################################################################
#------------------------------ Run LDSC software ----------------------------------#
#####################################################################################
# Multiancestry with european LD structure
/lustre/groups/itg/shared/software/ldsc/ldsc.py --rg T2D.multiancestry.ldsc.sumstats.gz,SCZ.multiancestry.ldsc.sumstats.gz --ref-ld-chr ${DATA_DIR}/GeneticCorrelation/eur_w_ld_chr/ --w-ld-chr ${DATA_DIR}/GeneticCorrelation/eur_w_ld_chr/ --out multiancestry.correlation --pop-prev 0.1,0.01 --samp-prev 0.135,0.425

# EUR
/lustre/groups/itg/shared/software/ldsc/ldsc.py --rg T2D.eur.ldsc.sumstats.gz,SCZ.eur.ldsc.sumstats.gz --ref-ld-chr ${DATA_DIR}/GeneticCorrelation/eur_w_ld_chr/ --w-ld-chr ${DATA_DIR}/GeneticCorrelation/eur_w_ld_chr/ --out eur.correlation --pop-prev 0.1,0.01 --samp-prev 0.086,0.409

# EAS
/lustre/groups/itg/shared/software/ldsc/ldsc.py --rg T2D.asian.ldsc.sumstats.gz,SCZ.asian.ldsc.sumstats.gz --ref-ld-chr ${DATA_DIR}/GeneticCorrelation/eas_ldscores/ --w-ld-chr ${DATA_DIR}/GeneticCorrelation/eas_ldscores/ --out asian.correlation --pop-prev 0.1,0.01 --samp-prev 0.219,0.455

traits=(T2D.eur SCZ.eur T2D.multiancestry SCZ.multiancestry) 
pop_prev=(0.1 0.01 0.1 0.01)
sample_prev=(0.086 0.409 0.135 0.425)

for i in {0..3}
do 
    /lustre/groups/itg/shared/software/ldsc/ldsc.py \
    --rg ${traits[$i]}.ldsc.sumstats.gz,T2D.eur.ldsc.sumstats.gz,SCZ.eur.ldsc.sumstats.gz,T2D.multiancestry.ldsc.sumstats.gz,SCZ.multiancestry.ldsc.sumstats.gz,2hGlu.eur.ldsc.sumstats.gz,FG.eur.ldsc.sumstats.gz,FI.eur.ldsc.sumstats.gz,HbA1c.eur.ldsc.sumstats.gz,ADHD.ldsc.sumstats.gz,AnorexiaNervosa.ldsc.sumstats.gz,Autism.ldsc.sumstats.gz,BipolarDisorder.ldsc.sumstats.gz,Depression.ldsc.sumstats.gz,PanicDisorder.ldsc.sumstats.gz,PTSD_eur.ldsc.sumstats.gz \
    --ref-ld-chr ${DATA_DIR}/GeneticCorrelation/eur_w_ld_chr/ \
    --w-ld-chr ${DATA_DIR}/GeneticCorrelation/eur_w_ld_chr/ \
    --out ${traits[$i]}.correlation \
    --pop-prev ${pop_prev[$i]},0.1,0.01,0.1,0.01,nan,nan,nan,nan,0.03,0.02,0.01,0.024,0.06,0.3,0.04  \
    --samp-prev ${sample_prev[$i]},0.086,0.409,0.135,0.425,nan,nan,nan,nan,0.364,0.234,0.397,0.101,0.341,0.22,0.157
done

# $i.ldsc.sumstats.gz
# T2D.eur.ldsc.sumstats.gz
# SCZ.eur.ldsc.sumstats.gz
# T2D.multiancestry.ldsc.sumstats.gz
# SCZ.multiancestry.ldsc.sumstats.gz
# 2hGlu.eur.ldsc.sumstats.gz
# FG.eur.ldsc.sumstats.gz
# FI.eur.ldsc.sumstats.gz
# HbA1c.eur.ldsc.sumstats.gz
# ADHD.ldsc.sumstats.gz
# AnorexiaNervosa.ldsc.sumstats.gz
# Autism.ldsc.sumstats.gz
# BipolarDisorder.ldsc.sumstats.gz
# Depression.ldsc.sumstats.gz
# PanicDisorder.ldsc.sumstats.gz
# PTSD_eur.ldsc.sumstats.gz

# ------------------- adiposity ------------------- 
traits=(bmi whr whole_body_fat_mass body_fat_percentage) 
pop_prev=(nan nan nan nan)
sample_prev=(nan nan nan nan)

for i in {0..3}
do 
    /lustre/groups/itg/shared/software/ldsc/ldsc.py \
    --rg ${traits[$i]}.ldsc.sumstats.gz,T2D.eur.ldsc.sumstats.gz,SCZ.eur.ldsc.sumstats.gz,T2D.multiancestry.ldsc.sumstats.gz,SCZ.multiancestry.ldsc.sumstats.gz,2hGlu.eur.ldsc.sumstats.gz,FG.eur.ldsc.sumstats.gz,FI.eur.ldsc.sumstats.gz,HbA1c.eur.ldsc.sumstats.gz,ADHD.ldsc.sumstats.gz,AnorexiaNervosa.ldsc.sumstats.gz,Autism.ldsc.sumstats.gz,BipolarDisorder.ldsc.sumstats.gz,Depression.ldsc.sumstats.gz,PanicDisorder.ldsc.sumstats.gz,PTSD_eur.ldsc.sumstats.gz,bmi.ldsc.sumstats.gz,whr.ldsc.sumstats.gz,whole_body_fat_mass.ldsc.sumstats.gz,body_fat_percentage.ldsc.sumstats.gz \
    --ref-ld-chr ${DATA_DIR}/GeneticCorrelation/eur_w_ld_chr/ \
    --w-ld-chr ${DATA_DIR}/GeneticCorrelation/eur_w_ld_chr/ \
    --out ${traits[$i]}.correlation \
    --pop-prev ${pop_prev[$i]},0.1,0.01,0.1,0.01,nan,nan,nan,nan,0.03,0.02,0.01,0.024,0.06,0.3,0.04,nan,nan,nan,nan  \
    --samp-prev ${sample_prev[$i]},0.086,0.409,0.135,0.425,nan,nan,nan,nan,0.364,0.234,0.397,0.101,0.341,0.22,0.157,nan,nan,nan,nan
done

#####################################################################################
#---------------------------- Export results nicely --------------------------------#
#####################################################################################
for i in T2D.eur SCZ.eur T2D.multiancestry SCZ.multiancestry
do
    sed -e '1,/Summary of Genetic Correlation Results/d' $i.correlation.log | sed -e '/Analysis finished/,$d' > $i.correlation.result
done


for i in bmi whr whole_body_fat_mass body_fat_percentage
do
    sed -e '1,/Summary of Genetic Correlation Results/d' $i.correlation.log | sed -e '/Analysis finished/,$d' > $i.correlation.result
done




