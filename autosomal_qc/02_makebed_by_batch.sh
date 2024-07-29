#!/usr/bin/bash

batch='2'
datadir=/staging/reserve/jacobhsu/TWB/TWB/microarray/TWB2_genotyped_120163
prefix=TWB2.hg38.v4
preqcdir=/staging/biology/u4432941/TWB1492_mt/microarray/autosomal_qc/preqc_twb${batch}

PLINK="/opt/ohpc/Taiwania3/pkg/biology/PLINK/PLINK_v1.90/plink"

# Settings
mkdir -p $preqcdir

TIME=`date +%Y%m%d%H%M`
logfile=${preqcdir}/${TIME}_make_bed.log
# Redirect standard output and error to the log file
exec > "$logfile" 2>&1


$PLINK \
    --bfile ${datadir}/${prefix} \
    --extract ${preqcdir}/twb${batch}_geno05_mind02_geno02.snplist \
    --keep ${preqcdir}/twb${batch}_geno05_mind02.indvlist \
    --mind 0.02 \
    --geno 0.02 \
    --make-bed \
    --out ${preqcdir}/twb${batch}_geno05_mind02_geno02
    