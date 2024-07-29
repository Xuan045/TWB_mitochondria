#!/usr/bin/bash

batch='2'
datadir=/staging/reserve/jacobhsu/TWB/TWB/microarray/TWB2_genotyped_120163
prefix=TWB2.hg38.v4
preqcdir=/staging/biology/u4432941/TWB1492_mt/microarray/preqc_twb${batch}
ref=/staging/biology/u4432941/1000GP/1000GP.GRCh38.snv.maf05.ldpr.bim
REF="/staging/reserve/paylong_ntu/AI_SHARE/reference/From_1000G/GRCh38_full_analysis_set_plus_decoy_hla.fa"

PLINK="/opt/ohpc/Taiwania3/pkg/biology/PLINK/PLINK_v1.90/plink"
PLINK2="/opt/ohpc/Taiwania3/pkg/biology/PLINK2/PLINK_v2.00a2.3_AVX2/plink2"


# Settings
mkdir -p $preqcdir

TIME=`date +%Y%m%d%H%M`
logfile=${preqcdir}/${TIME}_call_rate.log
# Redirect standard output and error to the log file
exec > "$logfile" 2>&1

# Find out duplicated SNPs (same chr, pos, ref, alt)
$PLINK \
    --bfile ${preqcdir}/twb${batch}_geno05_mind02_geno02 \
    --list-duplicate-vars ids-only suppress-first \
    --out ${preqcdir}/twb${batch}_dupcheck
$PLINK \
    --bfile ${preqcdir}/twb${batch}_geno05_mind02_geno02 \
    --exclude ${preqcdir}/twb${batch}_dupcheck.dupvar \
    --make-bed \
    --out ${preqcdir}/twb${batch}_geno05_mind02_geno02_dedup

# set ref/alt allels
# normalize the alleles
# set all the variant IDs to chr:pos:ref:alt and remove indels with length > 10
$PLINK2 \
    --bfile ${preqcdir}/twb${batch}_geno05_mind02_geno02_dedup \
    --ref-from-fa --fa $REF \
    --snps-only just-acgt \
    --output-chr chrM \
    --set-all-var-ids @:#:\$r:\$a \
    --make-bed \
    --out ${preqcdir}/twb${batch}_geno05_mind02_geno02_snps


$PLINK \
    --bfile ${preqcdir}/twb${batch}_geno05_mind02_geno02_snps \
    --not-chr 0 \
    --mac 1 \
    --make-bed \
    --out ${preqcdir}/twb${batch}_btqc
