#!/bin/bash

batch='2'
datadir=/staging/reserve/jacobhsu/TWB/TWB/microarray/TWB2_genotyped_120163
prefix=TWB2.hg38.v4
wkdir=/staging/biology/u4432941/TWB1492_mt/microarray
preqcdir=$wkdir/preqc_twb${batch}
pcadir=$wkdir/pca_twb${batch}

PLINK="/opt/ohpc/Taiwania3/pkg/biology/PLINK/PLINK_v1.90/plink"
PLINK2="/opt/ohpc/Taiwania3/pkg/biology/PLINK2/PLINK_v2.00a2.3_AVX2/plink2"


# Perform LD pruning: chrX
$PLINK2 \
    --bfile ${datadir}/${prefix} \
    --chr 23 \
    --keep ${pcadir}/twb${batch}_1kg_predpop0.8_EAS_indvlist \
    --geno 0.02 \
    --maf 0.05 \
    --snps-only just-acgt \
    --indep-pairwise 200 100 0.1 \
    --out ${preqcdir}/twb${batch}_EAS_chrx


# Check sex
$PLINK \
    --bfile ${datadir}/${prefix} \
    --chr 23 \
    --keep ${pcadir}/twb${batch}_1kg_predpop0.8_EAS_indvlist \
    --extract ${preqcdir}/twb${batch}_EAS_chrx.prune.in \
    --check-sex \
    --out ${preqcdir}/twb${batch}_EAS_sex

# Plot
module load compiler/gcc/9.4.0
module load biology/R/4.1.0
Rscript /staging/biology/u4432941/TWB1492_mt/microarray/autosomal_qc/07_pop_check_sex.R $preqcdir

