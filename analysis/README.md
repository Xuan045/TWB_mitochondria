# Analysis Scripts

This folder contains scripts for the analysis of both genotyping array and WGS data used in the thesis.

## Directory Structure

- **fig2_cov_cn_distribution.R**: Script for coverage and copy number distribution plots.
- **fig3a_pie_chart.R**: Script for generating pie charts.
- **fig3b_af_bin_proportion.R**: Script for generating AF bin proportion plots.
- **fig5_var_distribution_plot.R**: Script for generating variant distribution plots.
- **fig6_var_consequence_distribution.R**: Script for generating variant consequence distribution plots.
- **fig7_table2_cfrm.R**: Script for generating confirm pathogenic variants plot.

### `fig4_af_compare/`

- **af_compare.R**: Script for comparing AFs with gnomAD.
- **mitomap_af_compare.R**: Script for comparing AFs with MITOMAP.
- **mitomap/**: Directory containing MITOMAP VCF files and normalization script.
  - **mitomap.split.leftaligned.vcf**
  - **mitomap_polymorphisms.vcf**
  - **normalization.sh**

### `fig8_hg/`

- **run_haplogrep.sh**: Script for assigning haplogroups.
- **fig8a_hg_distribution.R**: Script for haplogroup distribution plots.
- **fig8c_phylotree.py**: Script for generating phylotrees.

### `fig8_fig16_pca/`

- **pca_autosomal.sh**: Bash script for preprocessing autosomal VCF data and running PCA using PLINK.
- **pca_mt.sh**: Bash script for preprocessing mitochondrial VCF data and running PCA using PLINK.
- **pca_plot_array.R**: R script for generating PCA plots from genotyping array data.
- **pca_plot_wgs.R**: R script for generating PCA plots from WGS data.

### `fig9_ld/`

- **ld_plink_array.sh**: Shell script for running LD analysis starting with PLINK bfile.
- **ld_plink_wgs.sh**: Shell script for running LD analysis starting with a VCF file.
- **plot_ld.R**: R script for generating a heatmap based on the LD R² values.

### `seq_info/`

- **NC_012920.1_sequence.txt**: Reference sequence file.
- **get_seq_info.py**: Python script for extracting sequence information and generating combined_positions.tsv.
- **combined_positions.tsv**: File containing combined positions information.

## Detailed Description of Scripts

### Fig. 2. Distribution of mtDNA Coverage and Copy Number in WGS Analyses

**fig2_cov_cn_distribution.R**

**Input Files**:
- `sample_annotations.txt`, `combined_sites_only.vcf`, and `output.tsv` generated using the GATK pipeline.
- `combined_positions.tsv` generated using `seq_info/get_seq_info.py`.

### Fig. 3. Distribution and Classification of mtDNA Variants in WGS Analyses

**fig3a_pie_chart.R**

**Input Files**:
- `stats_pass.txt` and `combined_sites_only.txt` generated using the GATK pipeline.

**fig3b_af_bin_proportion.R**

**Input Files**:
- `combined_sites_only.txt` generated using the GATK pipeline.

### Fig. 4. AF Comparison Scripts

**af_compare.R**

This script compares allele frequencies between the Taiwan Biobank (TWB) and the gnomAD database.

**Input Files**:
- Annotated VCF file with gnomAD AFs:
  - `POS`
  - `REF`
  - `ALT`
  - `TWB1465_MT_FILTER`
  - `gnomAD_MT_FILTER`
  - `TWB1465_MT_AF_hom`
  - `TWB1465_MT_AF_het_vaf10`
  - Columns for gnomAD AFs (e.g., `gnomAD_MT_AF_hom_eas`, `gnomAD_MT_AF_het_eas`).

**mitomap_af_compare.R**

This script compares allele frequencies with the MITOMAP database.

**Input Files**:
- MITOMAP VCF file (downloaded and indels left-aligned using BCFtools).

### Fig. 5. Analysis of mtDNA Variants and Age-Related Trends in Heteroplasmies

**fig5_var_distribution_plot.R**

**Input Files**:
- `sample_annotations.txt` generated using the GATK pipeline.

### Fig. 6. Variant Consequences of mtDNA Variants

**fig6_var_consequence_distribution.R**

**Input Files**:
- `combined_sites_only.txt` generated using the GATK pipeline.

### Fig. 7. Known Pathogenic Variants in TWB WGS Samples

**fig7_table2_cfrm.R**

**Input Files**:
- `reconstructed_vcf.txt` generated using the `reformat_vcf.R` script.

### Fig. 8. mtDNA Haplogroups Analysis

**run_haplogrep.sh**

This script assigns haplogroups for each individual using only the PASS variants.

**fig8a_hg_distribution.R**

**Input Files**:
- `sample_annotations.txt` generated using the GATK pipeline.
- `TWB1492_hg_allPASS_extend.txt` generated using the `run_haplogrep.sh`.

**fig8c_phylotree.py**

The script generates the phylotree with colored haplogroups presented in the TWB.

### Fig. 9. LD Analysis

**ld_plink_array.sh**

This script runs LD analysis on genotyping array data using PLINK bfile.

**Input Files**:
- PLINK bfile.

**ld_plink_wgs.sh**

This script runs LD analysis on WGS data starting with a VCF file.

**Input Files**:
- VCF file.

**plot_ld.R**

This script generates a heatmap based on the LD R² values.

### Fig. 8 and Fig. 16. PCA Analysis

**autosomal_pca.sh**

This script preprocesses autosomal VCF data and performs PCA using PLINK.

**Input Files**:
- Autosomal VCF file.
- Reference genome (FASTA).
- List of related samples to remove.

**mt_pca.sh**

This script preprocesses mitochondrial VCF data and performs PCA using PLINK.

**Input Files**:
- Mitochondrial VCF file.
- Reference genome (FASTA).

**pca_plot_array.R**

This script generates PCA plots from genotyping array data.

**Input Files**:
- Eigenvalue file.
- Eigenvector file.
- Sample annotation file.
- Ancestry information file.

**pca_plot_wgs.R**

This script generates PCA plots from WGS data.

**Input Files**:
- Eigenvalue file.
- Eigenvector file.
- Sample annotation file.
- Ancestry information file.
