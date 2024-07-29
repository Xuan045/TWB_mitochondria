#!/usr/bin/bash

vcf="/staging/biology/u4432941/TWB1492_mt/plink_pca/TWB1492_PCA_decom_norm_vqsr.vcf.gz"
SAMPLE_ID=TWB1492_PCA
wkdir="/staging/biology/u4432941/TWB1492_mt/plink_pca"

# for PLINK setting
mr_geno=0.02
mr_indv=0.02
# for output file name
geno=02
mind=02

REF="/staging/reserve/paylong_ntu/AI_SHARE/reference/GATK_bundle/2.8/hg38/Homo_sapiens_assembly38.fasta"
PLINK="/opt/ohpc/Taiwania3/pkg/biology/PLINK/PLINK_v1.90/plink"
PLINK2="/opt/ohpc/Taiwania3/pkg/biology/PLINK2/PLINK_v2.00a2.3_AVX2/plink2"
BCFTOOLS="/opt/ohpc/Taiwania3/pkg/biology/BCFtools/bcftools_v1.13/bin/bcftools"

QC_dir=${wkdir}/vqsr_geno${geno}_indv${mind}_snp_ldpr
mkdir -p ${QC_dir}
cd ${QC_dir}

TIME=`date +%Y%m%d%H%M`
logfile=${QC_dir}/${TIME}_${SAMPLE_ID}.log
# Redirect standard output and error to the log file
exec > "$logfile" 2>&1

# VCF to PLINK format
realated_list="/staging/biology/u4432941/TWB1492_mt/plink_pca/related_list"
$PLINK2 --vcf $vcf \
	--make-bed \
	--remove ${realated_list} \
	--allow-extra-chr \
	--out ${QC_dir}/${SAMPLE_ID}

# Filter out vaiants with missing rate > mr_geno or sample missing rate > mr_indv
# Remove indels
${PLINK2} --bfile ${QC_dir}/${SAMPLE_ID} \
	--autosome \
	--make-bed \
	--not-chr 0 \
	--allow-extra-chr \
	--geno ${mr_geno} \
	--mind ${mr_indv} \
	--snps-only \
	--out ${QC_dir}/${SAMPLE_ID}_vqsr_geno${geno}_mind${mind}_snp

# Filter out singleton with --mac 2 (minimum count)
${PLINK2} --bfile ${QC_dir}/${SAMPLE_ID}_vqsr_geno${geno}_mind${mind}_snp \
	--allow-extra-chr \
	--mac 2 \
	--make-bed \
	--out ${QC_dir}/${SAMPLE_ID}_vqsr_geno${geno}_mind${mind}_sing_snp

# Perform LD pruning
${PLINK2} --bfile ${QC_dir}/${SAMPLE_ID}_vqsr_geno${geno}_mind${mind}_sing_snp \
	--indep-pairwise 200 100 0.1 \
	--maf 0.01 \
	--out ${QC_dir}/${SAMPLE_ID}_vqsr_geno${geno}_mind${mind}_sing_snp_ldpr

# RUn PCA
${PLINK2} --bfile ${QC_dir}/${SAMPLE_ID}_vqsr_geno${geno}_mind${mind}_sing_snp \
	--extract ${QC_dir}/${SAMPLE_ID}_vqsr_geno${geno}_mind${mind}_sing_snp_ldpr.prune.in \
	--pca 20 \
	--out ${QC_dir}/${SAMPLE_ID}_vqsr_geno${geno}_mind${mind}_sing_snp_ldpr

# Output number of variants left to qc.txt
touch qc_vqsr_geno${geno}_mind${mind}_sing_snp.txt

n_total=$(cat ${QC_dir}/${SAMPLE_ID}.bim | wc -l)
n_var=$(cat ${QC_dir}/${SAMPLE_ID}_vqsr_geno${geno}_mind${mind}_snp.bim | wc -l)
n_sing==$(cat ${QC_dir}/${SAMPLE_ID}_vqsr_geno${geno}_mind${mind}_sing_snp.bim | wc -l)

echo "N_var in the VCF: "${n_total} > qc_vqsr_geno${geno}_mind${mind}_sing_snp.txt
echo "N_var after removing INDELs: "${n_var} >> qc_vqsr_geno${geno}_mind${mind}_sing_snp.txt
echo "N_var after removing singletons: "${n_sing} >> qc_vqsr_geno${geno}_mind${mind}_sing_snp.txt
