# Mitochondrial Variant Calling

This pipeline is adjusted based on the MT pipeline provided by gnomAD ([Mitochondria-SNPs-Indels-hg38 on Terra](https://app.terra.bio/#workspaces/help-gatk/Mitochondria-SNPs-Indels-hg38)), and the version used here is v.4.1.8.0. After running the pipeline, the output file can be further analyzed by the scripts provided by gnomAD. The files can be found through [broadinstitute/gnomad-mitochondria](https://github.com/broadinstitute/gnomad-mitochondria).

## Required Reference Files

The chrM reference genome and other related files can be downloaded from [gnomAD Google Cloud Storage](https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/chrM?prefix=&forceOnObjectsSortingFiltering=false), the following files are needed
- hg38_v0_chrM_Homo_sapiens_assembly38.chrM.fasta
- hg38_v0_chrM_Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta
- hg38_v0_chrM_ShiftBack.chain
- hg38_v0_chrM_non_control_region.chrM.interval_list
- hg38_v0_chrM_control_region_shifted.chrM.interval_list
- hg38_v0_chrM_blacklist_sites.hg38.chrM.bed

## Software dependencies
- [haplocheckCLI](https://github.com/leklab/haplocheckCLI)
    This tool can detect in-sample contamination in mtDNA and estimate MT haplogroup.

## Run the pipeline

1.	Run the Pipeline:
  - Execute the `single_sample_mt_pipeline.sh` script.
	- WGS hg38 aligned BAM or CRAM files are required as inputs.

2. Variant Annotation:
   - Use the `vep_annotation.sh` script to perform variant annotation using VEP. Adjustments include:
     - `--distance 0` to eliminate upstream and downstream annotations.
     - VCF files from gnomAD and TWB are prepared to conduct the custom annotation during the process (files in **custom_annotation/** folder).
	- gnomAD VCFs can be downloaded from the [official website](https://gnomad.broadinstitute.org/downloads#v3-mitochondrial-dna).

