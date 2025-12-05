#!/usr/bin/bash

vcf=/staging/biology/u4432941/TWB1492_mt/add_annotation_output/sample_vcf.vcf.bgz
SAMPLE_ID=TWB1492_PCA_gnomad_mt_pass
wkdir="/staging/biology/u4432941/TWB1492_mt/plink_pca"
QC_dir=$wkdir/gnomad_mt

REF="/staging/reserve/paylong_ntu/AI_SHARE/reference/GATK_bundle/2.8/hg38/Homo_sapiens_assembly38.fasta"
PLINK="/opt/ohpc/Taiwania3/pkg/biology/PLINK/PLINK_v1.90/plink"
PLINK2="/opt/ohpc/Taiwania3/pkg/biology/PLINK2/PLINK_v2.00a2.3_AVX2/plink2"
BCFTOOLS="/opt/ohpc/Taiwania3/pkg/biology/BCFtools/bcftools_v1.13/bin/bcftools"

mkdir -p ${QC_dir}
cd ${QC_dir}

# --vcf-filter -> keep only the variants with "PASS"
# --const-fid 0 -> prevent PLINK from converting SAMPLE_ID to FAM_ID
$PLINK --vcf $vcf \
	--vcf-filter \
	--make-bed \
	--const-fid 0 \
	--output-chr 26 \
	--out ${SAMPLE_ID}

$PLINK --bfile ${SAMPLE_ID} \
	--chr-set 26 \
	--pca 10 \
	--out ${SAMPLE_ID}

