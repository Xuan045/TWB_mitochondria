#!/usr/bin/bash

batch='2'
datadir=/staging/reserve/jacobhsu/TWB/TWB/microarray/TWB2_genotyped_120163
prefix=TWB2.hg38.v4
wkdir=/staging/biology/u4432941/TWB1492_mt/microarray
preqcdir=$wkdir/preqc_twb${batch}
pcadir=$wkdir/pca_twb${batch}

PLINK="/opt/ohpc/Taiwania3/pkg/biology/PLINK/PLINK_v1.90/plink"
PLINK2="/opt/ohpc/Taiwania3/pkg/biology/PLINK2/PLINK_v2.00a2.3_AVX2/plink2"


$PLINK2 \
    --bfile ${preqcdir}/twb${batch}_btqc \
    --keep ${pcadir}/twb${batch}_1kg_predpop0.8_EAS_indvlist \
    --make-bed \
    --out ${preqcdir}/twb${batch}_btqc_EAS

$PLINK2 \
    --bfile ${preqcdir}/twb${batch}_btqc_EAS \
    --keep ${pcadir}/twb${batch}_1kg_predpop0.8_EAS_indvlist \
    --freq \
    --out ${preqcdir}/twb${batch}_btqc_EAS_frq


