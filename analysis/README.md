# Analysis Scripts

This folder contains various scripts for analysis, including allele frequency (AF) comparison, linkage disequilibrium (LD) analysis, and principal component analysis (PCA) for both genotyping array and whole genome sequencing (WGS) data.

## Directory Structure

- **af_compare/**: Directory containing scripts for allele frequency comparison.
  - `af_compare.R`: Script for comparing AFs with gnomAD.
  - `mitomap_af_compare.R`: Script for comparing AFs with MITOMAP.
- **ld/**: Directory containing scripts for LD analysis.
  - `ld_plink_array.sh`: Shell script for running LD analysis starting with PLINK bfile.
  - `ld_plink_wgs.sh`: Shell script for running LD analysis starting with a VCF file.
  - `plot_ld.R`: R script for generating a heatmap based on the LD R² values.
- **pca/**: Directory containing scripts for PCA analysis.
  - `autosomal_pca.sh`: Bash script for preprocessing autosomal VCF data and running PCA using PLINK.
  - `mt_pca.sh`: Bash script for preprocessing mitochondrial VCF data and running PCA using PLINK.
  - `pca_plot_array.R`: R script for generating PCA plots from genotyping array data.
  - `pca_plot_wgs.R`: R script for generating PCA plots from WGS data.

## AF Comparison Scripts

### `af_compare.R`

This script compares allele frequencies between the Taiwan Biobank (TWB) and the gnomAD database.

#### Input Files

- Annotated VCF file with gnomAD AFs:
  - `POS`
  - `REF`
  - `ALT`
  - `TWB1465_MT_FILTER`
  - `gnomAD_MT_FILTER`
  - `TWB1465_MT_AF_hom`
  - `TWB1465_MT_AF_het_vaf10`
  - Columns for gnomAD AFs (e.g., `gnomAD_MT_AF_hom_eas`, `gnomAD_MT_AF_het_eas`)

### `mitomap_af_compare.R`

This script compares allele frequencies with the MITOMAP database.

#### Input Files

- MITOMAP VCF file (downloaded and indels left-aligned using BCFtools)

## LD Analysis Scripts

### `ld_plink_array.sh`

This script runs LD analysis on genotyping array data using PLINK bfile.

#### Input Files

- PLINK bfile

### `ld_plink_wgs.sh`

This script runs LD analysis on WGS data starting with a VCF file.

#### Input Files

- VCF file

### `plot_ld.R`

This script generates a heatmap based on the LD R² values.

## PCA Analysis Scripts

### `autosomal_pca.sh`

This script preprocesses autosomal VCF data and performs PCA using PLINK.

#### Input Files

- Autosomal VCF file
- Reference genome (FASTA)
- List of related samples to remove

### `mt_pca.sh`

This script preprocesses mitochondrial VCF data and performs PCA using PLINK.

#### Input Files

- Mitochondrial VCF file
- Reference genome (FASTA)

### `pca_plot_array.R`

This script generates PCA plots from genotyping array data.

#### Input Files

- Eigenvalue file
- Eigenvector file
- Sample annotation file
- Ancestry information file

### `pca_plot_wgs.R`

This script generates PCA plots from WGS data.

#### Input Files

- Eigenvalue file
- Eigenvector file
- Sample annotation file
- Ancestry information file
