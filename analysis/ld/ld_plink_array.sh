#!/usr/bin/bash

BFILE="/staging/biology/u4432941/TWB1492_mt/microarray/mt_qc/qc_twb2/twb2_mind02_geno02_mono_king01"
OUTDIR="/staging/biology/u4432941/TWB1492_mt/ld/array_af01"
PARA="array.af01"

PLINK="/opt/ohpc/Taiwania3/pkg/biology/PLINK/PLINK_v1.90/plink"
PLINK2="/opt/ohpc/Taiwania3/pkg/biology/PLINK2/PLINK_v2.00a2.3_AVX2/plink2"

mkdir -p $OUTDIR
cd $OUTDIR

# Filter by AF
$PLINK2 \
    --bfile $BFILE \
    --maf 0.01 \
    --rm-dup list \
    --make-bed \
    --out $OUTDIR/$PARA

# Calculate pairwise LD
$PLINK \
    --bfile $OUTDIR/$PARA \
    --keep-allele-order \
    --r2 square \
    --out $OUTDIR/$PARA
