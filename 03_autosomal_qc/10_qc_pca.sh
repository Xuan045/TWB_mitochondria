#!/usr/bin/bash

batch='2'
datadir=/staging/reserve/jacobhsu/TWB/TWB/microarray/TWB2_genotyped_120163
prefix=TWB2.hg38.v4
wkdir=/staging/biology/u4432941/TWB1492_mt/microarray
preqcdir=$wkdir/preqc_twb${batch}
pcadir=$wkdir/pca_twb${batch}
fhighld_region=$wkdir/misc/long_range_LD_intervals.txt

PLINK="/opt/ohpc/Taiwania3/pkg/biology/PLINK/PLINK_v1.90/plink"
PLINK2="/opt/ohpc/Taiwania3/pkg/biology/PLINK2/PLINK_v2.00a2.3_AVX2/plink2"


# Settings
mkdir -p $pcadir

TIME=`date +%Y%m%d%H%M`
logfile=${pcadir}/${TIME}_pca.log
# Redirect standard output and error to the log file
exec > "$logfile" 2>&1

# PCA on study sample
$PLINK2 \
    --bfile ${preqcdir}/twb${batch}_btqc_EAS \
    --remove ${preqcdir}/twb${batch}_EAS_sexcheck_het_king_remove.indvlist \
    --indep-pairwise 200 100 0.1 \
    --out ${pcadir}/twb${batch}_EAS_qc_ldpr

$PLINK2 \
    --bfile ${preqcdir}/twb${batch}_btqc_EAS \
    --remove ${preqcdir}/twb${batch}_EAS_sexcheck_het_king_remove.indvlist \
    --extract ${pcadir}/twb${batch}_EAS_qc_ldpr.prune.in \
    --pca 20 approx \
    --out ${pcadir}/twb${batch}_EAS_pca_qc


# Plot PCA
module load compiler/gcc/9.4.0
module load biology/R/4.1.0
Rscript /staging/biology/u4432941/TWB1492_mt/microarray/autosomal_qc/10_qc_pca.R $pcadir
