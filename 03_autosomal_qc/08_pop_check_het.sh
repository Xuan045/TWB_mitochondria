#!/bin/bash

batch='2'
datadir=/staging/reserve/jacobhsu/TWB/TWB/microarray/TWB2_genotyped_120163
prefix=TWB2.hg38.v4
wkdir=/staging/biology/u4432941/TWB1492_mt/microarray
preqcdir=$wkdir/preqc_twb${batch}
pcadir=$wkdir/pca_twb${batch}

PLINK="/opt/ohpc/Taiwania3/pkg/biology/PLINK/PLINK_v1.90/plink"
PLINK2="/opt/ohpc/Taiwania3/pkg/biology/PLINK2/PLINK_v2.00a2.3_AVX2/plink2"


# Estimate autosomal heterozygosity rate/inbreeding coef
# also use pruned SNPs! see: http://zzz.bwh.harvard.edu/plink/ibdibs.shtml
$PLINK \
    --bfile ${preqcdir}/twb${batch}_btqc_EAS \
    --extract ${pcadir}/twb${batch}_ref_ldpr.prune.in \
    --het \
    --out ${preqcdir}/twb${batch}_btqc_EAS_inbr


# Calculate hetorozygosity rate (i.e., the proportion of heterozygous genotypes for a given individual.)
awk 'NR>1{print ($5-$3)/$5}' ${preqcdir}/twb${batch}_btqc_EAS_inbr.het | sed '1i HetRate' | \
    paste ${preqcdir}/twb${batch}_btqc_EAS_inbr.het -> ${preqcdir}/twb${batch}_btqc_EAS_inbr.hetrate

# Plot
module load compiler/gcc/9.4.0
module load biology/R/4.1.0
Rscript /staging/biology/u4432941/TWB1492_mt/microarray/autosomal_qc/09_pop_check_het.R $preqcdir
