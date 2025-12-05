#!/usr/bin/bash

VCF="/staging/biology/u4432941/TWB1492_mt/microarray/mt/imputation/wgs_split/recode_vcf/TWB1465_wgs_recode.vcf.gz"
OUTDIR="/staging/biology/u4432941/TWB1492_mt/ld/hom_only_af01"
PARA="wgs_hom_only.af01"

PLINK="/opt/ohpc/Taiwania3/pkg/biology/PLINK/PLINK_v1.90/plink"
PLINK2="/opt/ohpc/Taiwania3/pkg/biology/PLINK2/PLINK_v2.00a2.3_AVX2/plink2"

mkdir -p $OUTDIR
cd $OUTDIR

# Convert VCF to plink bfile
$PLINK2 \
    --vcf $VCF \
    --make-bed \
    --out $OUTDIR/wgs_hom_only

# Filter by AF
$PLINK2 \
    --bfile $OUTDIR/wgs_hom_only \
    --maf 0.01 \
    --rm-dup list \
    --make-bed \
    --out $OUTDIR/$PARA

# Calculate pairwise LD
$PLINK \
    --bfile $OUTDIR/$PARA \
    --r2 square \
    --out $OUTDIR/$PARA
