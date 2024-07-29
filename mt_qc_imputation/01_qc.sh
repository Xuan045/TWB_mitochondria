#!/usr/bin/bash

#twb1
# batch='1'
# datadir=/staging/reserve/jacobhsu/TWB/TWB/microarray/TWB1/27719/Genotyped/
# prefix=TWB1.hg19

# twb2
batch='2'
datadir=/staging/reserve/jacobhsu/TWB/TWB/microarray/TWB2_genotyped_120163
prefix=TWB2.hg38.v4

mtdir=/staging/biology/u4432941/TWB1492_mt/microarray/mt_qc
mtqcdir=/staging/biology/u4432941/TWB1492_mt/microarray/mt_qc/qc_twb${batch}_final
mkdir -p $mtqcdir

PLINK="/opt/ohpc/Taiwania3/pkg/biology/PLINK/PLINK_v1.90/plink"
PLINK2="/opt/ohpc/Taiwania3/pkg/biology/PLINK2/PLINK_v2.00a2.3_AVX2/plink2"

TIME=`date +%Y%m%d%H%M`
logfile=${mtqcdir}/${TIME}_twb${batch}_callrate.log
# Redirect standard output and error to the log file
exec > "$logfile" 2>&1


# 1. Output a list of SNPs with call rate >0.95
$PLINK \
    --bfile ${mtqcdir}/twb${batch}_mt_recode \
    --chr 26 \
    --geno 0.05 \
    --write-snplist \
    --out ${mtqcdir}/twb${batch}_geno05

# 2. Output a list of samples with call rate >0.98
$PLINK2 \
    --bfile ${mtqcdir}/twb${batch}_mt_recode \
    --extract ${mtqcdir}/twb${batch}_geno05.snplist \
    --missing \
    --out ${mtqcdir}/twb${batch}_geno05_indv
awk 'NR>1{if($6<0.02){print $1,$2}}' ${mtqcdir}/twb${batch}_geno05_indv.smiss > ${mtqcdir}/twb${batch}_geno05_mind02.indvlist

# Plot the distribution of individual missing rate
module load compiler/gcc/9.4.0
module load biology/R/4.1.0
Rscript "00_idv_missing.r" ${mtqcdir} twb${batch}_geno05_indv.smiss


# 3. Filter to samples in the previous step, output a list of SNPs with call rate >0.98
$PLINK2 \
    --bfile ${mtqcdir}/twb${batch}_mt_recode \
    --extract ${mtqcdir}/twb${batch}_geno05.snplist \
    --keep ${mtqcdir}/twb${batch}_geno05_mind02.indvlist \
    --geno 0.02 \
    --write-snplist \
    --out ${mtqcdir}/twb${batch}_geno05_mind02_geno02


# 4. Keep samples and SNPs passing previous filters, calculate SNP-level missing rate
$PLINK2 \
    --bfile ${mtqcdir}/twb${batch}_mt_recode \
    --extract ${mtqcdir}/twb${batch}_geno05_mind02_geno02.snplist \
    --keep ${mtqcdir}/twb${batch}_geno05_mind02.indvlist \
    --missing \
    --out ${mtqcdir}/twb${batch}_geno05_mind02_geno02_miss

$PLINK2 \
    --bfile ${mtqcdir}/twb${batch}_mt_recode \
    --extract ${mtqcdir}/twb${batch}_geno05_mind02_geno02.snplist \
    --keep ${mtqcdir}/twb${batch}_geno05_mind02.indvlist \
    --make-bed \
    --out ${mtqcdir}/twb${batch}_geno05_mind02_geno02

