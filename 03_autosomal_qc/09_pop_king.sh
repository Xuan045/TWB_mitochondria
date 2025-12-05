#!/bin/bash

batch='2'
datadir=/staging/reserve/jacobhsu/TWB/TWB/microarray/TWB2_genotyped_120163
prefix=TWB2.hg38.v4
wkdir=/staging/biology/u4432941/TWB1492_mt/microarray
preqcdir=$wkdir/preqc_twb${batch}
pcadir=$wkdir/pca_twb${batch}

PLINK="/opt/ohpc/Taiwania3/pkg/biology/PLINK/PLINK_v1.90/plink"
PLINK2="/opt/ohpc/Taiwania3/pkg/biology/PLINK2/PLINK_v2.00a2.3_AVX2/plink2"
KING="/opt/ohpc/Taiwania3/pkg/biology/KING/KING_v2.2.7/king"

TIME=`date +%Y%m%d%H%M`
logfile=${preqcdir}/${TIME}_king.log
# Redirect standard output and error to the log file
exec > "$logfile" 2>&1

# $KING \
#     -b ${preqcdir}/twb${batch}_btqc_EAS.bed \
#     --related --degree 2 \
#     --cpus 14 \
#     --prefix $preqcdir/twb${batch}_btqc_EAS_king_related

# make remove list for related individuals (mitochondria)
module load compiler/gcc/9.4.0
module load biology/R/4.1.0

indv_df="/staging/biology/u4432941/TWB1492_mt/microarray_association/pheno_data/twb2_survey.final.txt"
Rscript /staging/biology/u4432941/TWB1492_mt/microarray/autosomal_qc/09_pop_king.R $preqcdir $indv_df

cat ${preqcdir}/twb${batch}_EAS_sex_mismatch_F025_M075.indvlist ${preqcdir}/twb${batch}_EAS_het_outlier_6sd.indvlist ${preqcdir}/twb${batch}_EAS_king_mt.indvlist | \
    sort | uniq > ${preqcdir}/twb${batch}_EAS_sexcheck_het_king_remove_mt.indvlist

# make remove list for related individuals (autosomal)
cat ${preqcdir}/twb${batch}_EAS_sex_mismatch_F025_M075.indvlist ${preqcdir}/twb${batch}_EAS_het_outlier_6sd.indvlist ${preqcdir}/twb${batch}_EAS_king.indvlist | \
    sort | uniq > ${preqcdir}/twb${batch}_EAS_sexcheck_het_king_remove.indvlist

