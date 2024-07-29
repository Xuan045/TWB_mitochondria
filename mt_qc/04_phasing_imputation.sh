#!/usr/bin/bash

mtqcdir=/staging/biology/u4432941/TWB1492_mt/microarray/mt_qc/qc_twb2_final
WKDIR=/staging/biology/u4432941/TWB1492_mt/microarray/mt_qc/imputation/wgs1465.shapeit2.impute2.pre_impute_filter_king_eas_2
mkdir -p $WKDIR
cd $WKDIR

INPUT="${mtqcdir}/twb2_mind02_geno02_mono_EAS_kingMT_2"
PARA="twb2_study"
gmap="/staging/biology/u4432941/TWB1492_mt/microarray/gmap/chrMT.NC012920.gmap"
mtmap="/staging/biology/u4432941/TWB1492_mt/microarray/mt_qc/mtgeneticmap.map"
ref_panel="/staging/biology/u4432941/TWB1492_mt/microarray/mt_qc/validate_imputation/wgs_split/recode_vcf/TWB1465_wgs_recode.vcf.gz"
ref_para="wgs1465"

PLINK="/opt/ohpc/Taiwania3/pkg/biology/PLINK/PLINK_v1.90/plink"
PLINK2="/opt/ohpc/Taiwania3/pkg/biology/PLINK2/PLINK_v2.00a2.3_AVX2/plink2"
IMPUTE2="/work/opt/ohpc/Taiwania3/pkg/biology/IMPUTE2/impute_v2.3.2/impute2"
SHAPEIT2="/work/opt/ohpc/Taiwania3/pkg/biology/SHAPEIT/shapeit_v2.904.glibcv2.17/bin/shapeit"

TIME=`date +%Y%m%d%H%M`
logfile=${WKDIR}/${TIME}_twb_phasing_imputation.log
# Redirect standard output and error to the log file
exec > "$logfile" 2>&1
module load pkg/Anaconda3

# 1. Phasing of study genotypes with SHAPEIT2 (針對study GT (要imputation) 的做phasing)
# Filter out variants with AF deviated from WGS (AF_diff > 0.2)
# First get gcount file using PLINK
$PLINK2 \
    --bfile ${INPUT} \
    --geno-counts \
    --make-bed \
    --out ${WKDIR}/${PARA}

# Filter out variants with AF difference with WGS greater than 0.2 using python script
WGS_VCF="/staging/biology/u4432941/TWB1492_mt/microarray/mt_qc/validate_imputation/wgs_split/recode_vcf/TWB1465_wgs_recode.vcf"
python /staging/biology/u4432941/TWB1492_mt/microarray/misc/remove_deviate.py ${WGS_VCF} ${WKDIR}/${PARA} > ${WKDIR}/non_deviate.snplist

# # Find strand ambiguous SNPs
# python /staging/biology/u4432941/TWB1492_mt/microarray/misc/find_atgc_snps.py ${INPUT}.bim > ${WKDIR}/atcg.snplist
# # Write a list of non-strand ambiguous SNPs to keep
# awk 'NR==FNR{a[$1];next} !($2 in a) {print $2}' ${WKDIR}/atcg.snplist ${INPUT}.bim > ${WKDIR}/ref_nonatgc.snplist

# # Combine two lists
# sort ${WKDIR}/non_deviate.snplist > ${WKDIR}/sorted_non_deviate.snplist
# sort ${WKDIR}/ref_nonatgc.snplist > ${WKDIR}/sorted_ref_nonatgc.snplist
# comm -12 ${WKDIR}/sorted_non_deviate.snplist ${WKDIR}/sorted_ref_nonatgc.snplist > ${WKDIR}/keep.snplist

$PLINK2 \
    --bfile ${INPUT} \
    --output-chr chrM \
    --extract ${WKDIR}/non_deviate.snplist \
    --make-bed \
    --out ${WKDIR}/${PARA}


$SHAPEIT2 \
    --input-bed ${WKDIR}/${PARA}.bed ${WKDIR}/${PARA}.bim ${WKDIR}/${PARA}.fam \
    --rho 4.0E-12 \
    --thread 14 \
    --force \
    -O ${WKDIR}/twb2.shapeit2

# 2. Transform reference panel to OXFORD format
$PLINK2 \
    --vcf ${ref_panel} \
    --recode oxford \
    --out ${WKDIR}/${ref_para}.ref.oxf
sed -i '2d' ${WKDIR}/${ref_para}.ref.oxf.sample
cut -d' ' -f2,3,4,5 < ${WKDIR}/${ref_para}.ref.oxf.gen > ${WKDIR}/${ref_para}.ref.legend
echo "rsID position a0 a1" > ${WKDIR}/header1.txt
cat ${WKDIR}/header1.txt ${WKDIR}/${ref_para}.ref.legend > ${WKDIR}/fin.${ref_para}.ref.legend

# 3. Imputation
${IMPUTE2} \
    -merge_ref_panels \
    -m $mtmap \
    -h ${WKDIR}/${ref_para}.ref.oxf.gen \
    -l ${WKDIR}/fin.${ref_para}.ref.legend \
    -known_haps_g $WKDIR/twb2.shapeit2.haps \
    -int 1 16579 \
    -Ne 20000 \
    -o ${WKDIR}/twb2.EAS_kingMT.shapeit2.impute2.diploid
