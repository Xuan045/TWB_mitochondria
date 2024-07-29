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
mt_fasta="/staging/biology/u4432941/TWB1492_mt/chrM_ref/hg38_v0_chrM_Homo_sapiens_assembly38.chrM.fasta"

TIME=`date +%Y%m%d%H%M`
logfile=${mtqcdir}/${TIME}_twb${batch}_extractMT.log
# Redirect standard output and error to the log file
exec > "$logfile" 2>&1

# 1. Extract SNVs on MT
$PLINK2 \
    --bfile ${datadir}/${prefix} \
    --chr 26 \
    --make-bed \
    --missing \
    --snps-only \
    --out ${mtqcdir}/twb${batch}_mt

# 2. Reorder the ref, alt allele
$PLINK2 \
    --bfile ${mtqcdir}/twb${batch}_mt \
    --normalize \
    --ref-from-fa -fa ${mt_fasta} \
    --make-bed \
    --out ${mtqcdir}/twb${batch}_mt_recode

# 2-2. Calculate the AF
$PLINK2 \
    --bfile ${mtqcdir}/twb${batch}_mt_recode \
    --chr 26 \
    --geno-counts \
    --freq \
    --out ${mtqcdir}/twb${batch}_mt_recode

# 3. Export VCF
$PLINK2 \
    --bfile ${mtqcdir}/twb${batch}_mt_recode \
    --export vcf \
    --out ${mtqcdir}/twb${batch}_mt_recode
bcftools +fill-tags ${mtqcdir}/twb${batch}_mt_recode.vcf -- -t AN,AC,AF > ${mtqcdir}/twb${batch}_mt_recode_tmp.vcf
mv ${mtqcdir}/twb${batch}_mt_recode_tmp.vcf ${mtqcdir}/twb${batch}_mt_recode.vcf

# QC summary for samples and SNPs
nindv=$(cat ${datadir}/${prefix}.fam | wc -l)
nsnps=$(cat ${datadir}/${prefix}.bim | wc -l)
echo "Nindv "$nindv > ${mtdir}/twb${batch}_qc_summary.txt
echo "Nsnps "$nsnps >> ${mtdir}/twb${batch}_qc_summary.txt
echo -en "\n" >> ${mtdir}/twb${batch}_qc_summary.txt

echo "Nsnps on chrM "$(cat ${mtqcdir}/twb${batch}_mt_recode.bim | wc -l) >> ${mtdir}/twb${batch}_qc_summary.txt
echo -en "\n" >> ${mtdir}/twb${batch}_qc_summary.txt
echo "snp missing rate:" >> ${mtdir}/twb${batch}_qc_summary.txt
awk 'NR > 1 { sum += $5; if (NR == 2 || $5 < min) min = $5; if (NR == 2 || $5 > max) max = $5; } END { print "Min:", min; print "Max:", max; print "Average:", sum / (NR - 1); }' ${mtqcdir}/twb${batch}_mt.vmiss >> ${mtdir}/twb${batch}_qc_summary.txt

echo -en "\nindiviual missing rate:\n" >> ${mtdir}/twb${batch}_qc_summary.txt
awk 'NR > 1 { sum += $6; if (NR == 2 || $6 < min) min = $6; if (NR == 2 || $6 > max) max = $6; } END { print "Min:", min; print "Max:", max; print "Average:", sum / (NR - 1); }' ${mtqcdir}/twb${batch}_mt.smiss >> ${mtdir}/twb${batch}_qc_summary.txt

# Plot the distribution of individual missing rate
module load compiler/gcc/9.4.0
module load biology/R/4.1.0
Rscript "00_idv_missing.r" ${mtqcdir} twb${batch}_mt.smiss

