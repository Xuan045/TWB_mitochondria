#!/usr/bin/bash

vcf=/staging/biology/u4432941/TWB1492_mt/haplogroup_check/microarray/twb1492_mt_pass.vcf.gz
HAPLOGREP=/staging/biology/u4432941/apps/haplogrep

# Default tree: 17_FU1
$HAPLOGREP classify \
	--in $vcf \
	--format vcf \
	--extend-report \
	--hits 3 \
	--lineage 1 \
	--out TWB1492_hg_allPASS
