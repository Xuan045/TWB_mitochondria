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

PLINK2="/opt/ohpc/Taiwania3/pkg/biology/PLINK2/PLINK_v2.00a2.3_AVX2/plink2"

TIME=`date +%Y%m%d%H%M`
logfile=${mtqcdir}/${TIME}_twb${batch}_king.log
# Redirect standard output and error to the log file
exec > "$logfile" 2>&1

# merge two lists
eas_indv="/staging/biology/u4432941/TWB1492_mt/microarray/pca_twb2/twb2_1kg_predpop0.8_non_EAS_indvlist"
related_indv="/staging/biology/u4432941/TWB1492_mt/microarray/preqc_twb2/twb2_EAS_king.indvlist"
cat "$related_indv" "$eas_indv" | sort | uniq > ${mtqcdir}/twb${batch}_rm_EAS_king.indv

PARA="mind02_geno02_mono_EAS_king"
$PLINK2 \
    --bfile ${mtqcdir}/twb${batch}_geno05_mind02_geno02_mono \
    --remove ${mtqcdir}/twb${batch}_rm_EAS_kingMT.indv \
    --make-bed \
    --out ${mtqcdir}/twb${batch}_${PARA}

$PLINK2 \
    --bfile ${mtqcdir}/twb${batch}_${PARA} \
    --recode vcf id-paste=iid \
    --out ${mtqcdir}/twb${batch}_${PARA}

# PCA on study sample
$PLINK2 \
    --bfile ${mtqcdir}/twb${batch}_${PARA} \
	--allow-extra-chr \
	--make-bed \
	--output-chr 26 \
	--chr 26 \
	--out ${mtqcdir}/twb${batch}_${PARA}_pca

$PLINK2 \
    --bfile ${mtqcdir}/twb${batch}_${PARA}_pca \
    --chr-set 30 \
    --pca 20 \
    --out ${mtqcdir}/twb${batch}_${PARA}
