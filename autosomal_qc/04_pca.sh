#!/usr/bin/bash

batch='2'
datadir=/staging/reserve/jacobhsu/TWB/TWB/microarray/TWB2_genotyped_120163
prefix=TWB2.hg38.v4
wkdir=/staging/biology/u4432941/TWB1492_mt/microarray
preqcdir=$wkdir/preqc_twb${batch}
pcadir=$wkdir/pca_twb${batch}
ref=/staging/biology/u4432941/1000GP/1000GP.GRCh38.snv.maf005
REF=/staging/reserve/paylong_ntu/AI_SHARE/reference/From_1000G/GRCh38_full_analysis_set_plus_decoy_hla.fa
fhighld_region=$wkdir/misc/long_range_LD_intervals.txt

PLINK="/opt/ohpc/Taiwania3/pkg/biology/PLINK/PLINK_v1.90/plink"
PLINK2="/opt/ohpc/Taiwania3/pkg/biology/PLINK2/PLINK_v2.00a2.3_AVX2/plink2"


# Settings
mkdir -p $pcadir

TIME=`date +%Y%m%d%H%M`
logfile=${pcadir}/${TIME}_pca.log
# Redirect standard output and error to the log file
exec > "$logfile" 2>&1

# Find SNPs in common between study sample and ref sample
awk '{seen[$2]++; if (seen[$2] > 1) print $2}' ${preqcdir}/twb${batch}_btqc.bim $ref.bim > ${pcadir}/twb${batch}_ref_comm.snplist

# Retain only overlapping variants
$PLINK2 \
    --bfile ${preqcdir}/twb${batch}_btqc_tmp \
    --extract ${pcadir}/twb${batch}_ref_comm.snplist \
    --make-bed \
    --out ${pcadir}/twb${batch}_comm

$PLINK2 \
    --bfile ${ref} \
    --extract ${pcadir}/twb${batch}_ref_comm.snplist \
    --make-bed \
    --out ${pcadir}/ref_comm

# Merge two panels
$PLINK \
    --bfile ${pcadir}/twb${batch}_comm \
    --keep-allele-order \
    --bmerge ${pcadir}/ref_comm \
    --make-bed \
    --out ${pcadir}/twb${batch}_ref

# Find strand ambiguous SNPs
python $wkdir/misc/find_atgc_snps.py ${pcadir}/twb${batch}_ref.bim > ${pcadir}/twb${batch}_ref_atgc.snplist

# Write a list of non-strand ambiguous SNPs to keep
awk 'NR==FNR{a[$1];next} !($2 in a) {print $2}' ${pcadir}/twb${batch}_ref_atgc.snplist ${pcadir}/twb${batch}_ref.bim > ${pcadir}/twb${batch}_ref_nonatgc.snplist

# Perform LD pruning
$PLINK \
    --bfile ${pcadir}/twb${batch}_ref \
    --autosome \
    --geno 0.02 \
    --maf 0.05 \
    --snps-only just-acgt \
    --extract ${pcadir}/twb${batch}_ref_nonatgc.snplist \
    --exclude range $fhighld_region \
    --indep-pairwise 200 100 0.1 \
    --out ${pcadir}/twb${batch}_ref_ldpr


# Run PCA
$PLINK2 \
    --bfile ${pcadir}/twb${batch}_ref \
    --extract ${pcadir}/twb${batch}_ref_ldpr.prune.in \
    --pca 20 approx \
    --out ${pcadir}/twb${batch}_ref_pca

# PCA on study sample
$PLINK2 \
    --bfile ${preqcdir}/twb${batch}_btqc \
    --extract ${pcadir}/twb${batch}_ref_ldpr.prune.in \
    --pca 20 approx \
    --out ${pcadir}/twb${batch}_pca

# Plot PCA
module load compiler/gcc/9.4.0
module load biology/R/4.1.0
ped_pop=/staging/biology/u4432941/TWB1492_mt/microarray/misc/ped_pop.tsv
Rscript /staging/biology/u4432941/TWB1492_mt/microarray/autosomal_qc/05_classify_ancestry.R $pcadir $ped_pop
