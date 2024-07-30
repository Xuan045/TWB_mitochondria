#!/usr/bin/bash

VCF="/staging/biology/u4432941/TWB1492_mt/mitomap/mitomap_replaceMT.vcf"
BCFTOOLS="/work/opt/ohpc/Taiwania3/pkg/biology/BCFtools/bcftools_v1.13/bin/bcftools"
REF="/staging/biology/u4432941/TWB1492_mt/chrM_ref/hg38_v0_chrM_Homo_sapiens_assembly38.chrM.fasta"

# Split multiallelic sites and left-align
${BCFTOOLS} norm -m -any ${VCF} | ${BCFTOOLS} norm -f $REF - -o mitomap.split.leftaligned.vcf

