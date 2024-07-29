#!/usr/bin/bash

INPUT_VCF_PATH=/staging/biology/u4432941/TWB1492_mt/custom_annotation/TWB1492_woInfoGT.vcf
OUTPUT_VCF_PATH=/staging/biology/u4432941/TWB1492_mt/custom_annotation
SAMPLE_ID=$(basename ${INPUT_VCF_PATH})

# Files for custom annotation
ref_dir=/staging/biology/u4432941/TWB1492_mt/custom_annotation
gnomad_mt=${ref_dir}/gnomad_af.vcf.bgz
twb_mt=${ref_dir}/twb1465_af.vcf.bgz

VEP_PATH=/opt/ohpc/Taiwania3/pkg/biology/Ensembl-VEP/ensembl-vep/vep
VEP_CACHE_DIR=/opt/ohpc/Taiwania3/pkg/biology/DATABASE/VEP/Cache
VEP_FASTA=/staging/reserve/paylong_ntu/AI_SHARE/reference/VEP_ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
BCFTOOLS=/opt/ohpc/Taiwania3/pkg/biology/BCFtools/bcftools_v1.13/bin/bcftools

module load biology
module load Perl/5.28.1
export PATH=${PATH}:/opt/ohpc/Taiwania3/pkg/biology/HTSLIB/htslib_v1.13/bin:/opt/ohpc/Taiwania3/pkg/biology/SAMTOOLS/samtools_v1.15.1/bin
set -euo pipefail

cd $OUTPUT_VCF_PATH

# Log file settings
TIME=`date +%Y%m%d%H%M`
logfile=./${TIME}_run_${SAMPLE_ID}.log

# Redirect standard output and error to the log file
exec > "$logfile" 2>&1

echo "$(date '+%Y-%m-%d %H:%M:%S') Job started" >> ${logfile}

custom=${gnomad_mt},gnomAD_MT,vcf,exact,0,FILTER,AF_hom,AF_het,AF_hom_afr,AF_het_afr,AF_hom_ami,AF_het_ami,AF_hom_amr,AF_het_amr,AF_hom_asj,AF_het_asj,AF_hom_eas,AF_het_eas,AF_hom_sas,AF_het_sas,AF_hom_fin,AF_het_fin,AF_hom_nfe,AF_het_nfe,AF_hom_mid,AF_het_mid,AF_hom_oth,AF_het_oth
custom_twb=${twb_mt},TWB1465_MT,vcf,exact,0,FILTER,AF_hom,AF_het_vaf05,AF_het_vaf10

OUTPUT_PARA=${SAMPLE_ID}_allAF
$VEP_PATH --cache --offline \
	-i $INPUT_VCF_PATH \
	--format vcf \
  	--everything \
	--hgvsg \
	--allele_number \
	--force_overwrite \
	--no_stats \
	--minimal \
	--distance 0 \
	--dir_cache $VEP_CACHE_DIR \
	--assembly GRCh38 \
	--custom ${custom} \
	--custom ${custom_twb} \
	--merged \
	--fasta $VEP_FASTA \
	--vcf \
	-o ${OUTPUT_PARA}.vcf
    
date >> ${logfile}
echo "Finished vcf output"  >> ${logfile}

grep -v -E '(^##contig=<ID=chr[[:digit:]]_|^##contig=<ID=chr[[:digit:]][[:digit:]]_|^##contig=<ID=HLA|^##contig=<ID=chrUn|^##contig=<ID=chr[XY]_)' ${OUTPUT_PARA}.vcf > ${OUTPUT_PARA}_clean_header.vcf

echo -e "CHROM\tPOS\tREF\tALT\t$(${BCFTOOLS} +split-vep -l ${OUTPUT_PARA}.vcf | cut -f 2 | tr '\n' '\t' | sed 's/\t$//')" > ${OUTPUT_PARA}.tsv
${BCFTOOLS} +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%CSQ\n' -d -A tab ${OUTPUT_PARA}.vcf >> ${OUTPUT_PARA}.tsv

echo "$(date '+%Y-%m-%d %H:%M:%S') Job finished" >> ${logfile}
