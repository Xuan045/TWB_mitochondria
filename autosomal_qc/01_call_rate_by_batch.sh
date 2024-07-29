#!/usr/bin/bash

batch='2'
datadir=/staging/reserve/jacobhsu/TWB/TWB/microarray/TWB2_genotyped_120163
prefix=TWB2.hg38.v4
preqcdir=/staging/biology/u4432941/TWB1492_mt/microarray/autosomal_qc/preqc_twb${batch}

PLINK="/opt/ohpc/Taiwania3/pkg/biology/PLINK/PLINK_v1.90/plink"

# Settings
mkdir -p $preqcdir

TIME=`date +%Y%m%d%H%M`
logfile=${preqcdir}/${TIME}_call_rate.log
# Redirect standard output and error to the log file
exec > "$logfile" 2>&1

# 1. Output a list of SNPs with call rate >0.95
$PLINK \
    --bfile ${datadir}/${prefix} \
    --autosome \
    --geno 0.05 \
    --write-snplist \
    --out ${preqcdir}/twb${batch}_geno05


# 2. Filter to SNPs in the previous step, output a list of samples with call rate >0.98
$PLINK \
    --bfile ${datadir}/${prefix} \
    --extract ${preqcdir}/twb${batch}_geno05.snplist \
    --missing \
    --out ${preqcdir}/twb${batch}_geno05_indv
awk 'NR>1{if($6<0.02){print $1,$2}}' ${preqcdir}/twb${batch}_geno05_indv.smiss > ${preqcdir}/twb${batch}_geno05_mind02.indvlist


# 3. Filter to samples in the previous step, output a list of SNPs with call rate >0.98
$PLINK \
    --bfile ${datadir}/${prefix} \
    --extract ${preqcdir}/twb${batch}_geno05.snplist \
    --keep ${preqcdir}/twb${batch}_geno05_mind02.indvlist \
    --geno 0.02 \
    --write-snplist \
    --out ${preqcdir}/twb${batch}_geno05_mind02_geno02


# 4. Keep samples and SNPs passing previous filters, calculate SNP-level missing rate
$PLINK \
    --bfile ${datadir}/${prefix} \
    --extract ${preqcdir}/twb${batch}_geno05_mind02_geno02.snplist \
    --keep ${preqcdir}/twb${batch}_geno05_mind02.indvlist \
    --missing \
    --out ${preqcdir}/twb${batch}_geno05_mind02_geno02_miss


# Summarize number of samples and SNPs removed at each step
nindv=$(cat ${datadir}/${prefix}.fam | wc -l)
nsnps=$(cat ${datadir}/${prefix}.bim | wc -l)
echo "Nindv "$nindv > qc_summary.txt
echo "Nsnps "$nsnps >> qc_summary.txt
echo -en "\n" >> qc_summary.txt

echo "Nsnps with call rate >0.95 "$(cat ${preqcdir}/twb${batch}_geno05.snplist | wc -l) >> qc_summary.txt
echo "Nindv with call rate >0.98 "$(cat ${preqcdir}/twb${batch}_geno05_mind02.indvlist | wc -l) >> qc_summary.txt
echo "Nsnps with call rate >0.98 "$(cat ${preqcdir}/twb${batch}_geno05_mind02_geno02.snplist | wc -l) >> qc_summary.txt
