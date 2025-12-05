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
logfile=${mtqcdir}/${TIME}_twb${batch}_rm_mono.log
# Redirect standard output and error to the log file
exec > "$logfile" 2>&1


# Remove monomorhpic SNPs
$PLINK2 \
    --bfile ${mtqcdir}/twb${batch}_geno05_mind02_geno02 \
    --not-chr 0 \
    --mac 1 \
    --make-bed \
    --out ${mtqcdir}/twb${batch}_geno05_mind02_geno02_mono

# Convert to VCF
$PLINK2 \
    --bfile ${mtqcdir}/twb${batch}_geno05_mind02_geno02_mono \
    --ref-from-fa --fa ${mt_fasta} \
    --recode vcf \
    --out ${mtqcdir}/twb${batch}_geno05_mind02_geno02_mono
