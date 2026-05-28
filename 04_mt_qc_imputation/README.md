# Mitochondrial Quality Control & Imputation

This folder contains scripts and configurations for the mitochondrial Quality Control (QC) and genotype imputation pipeline.

## Pipeline Configuration

Before running any script, copy/edit and check the settings in [config.sh](file:///staging/biology/u4432941/Github/TWB_mitochondria/04_mt_qc_imputation/config.sh). It centralizes environment management, project directories, tool executables, and input paths:
- **`ENV_MANAGER`**: Choose `"micromamba"` or `"conda"` for environment activation.
- **`PROJECT_DIR`**: Working directory where outputs, logs, and imputed files will be stored.
- **`INPUT_PLINK`**: Path to the base PLINK files (excluding prefix).
- **`PREFIX`**: The dataset prefix (e.g. `twb1`, `twb2`).

## Steps in the Pipeline

All scripts support SLURM job submission headers and source `config.sh` dynamically to set up environment paths.

1. **`01_extract_mt.sh`**
   - Extracts mitochondrial DNA (MT) SNVs from the main genotype database (`--chr 26`).
   - Reorders and normalizes the alleles matching the FASTA reference (`--ref-from-fa`).
   - Calculates the allele frequencies and exports to VCF format.
   - Plots the sample missing rate distribution using `01_idv_missing.r`.

2. **`02_qc.sh`**
   - First filter: removes SNPs with missing rate > 5% (`--geno 0.05`).
   - Second filter: removes samples with missing rate > 5% (`--mind 0.05`).
   - Third filter: keeps remaining samples and removes SNPs with missing rate > 2% (`--geno 0.02`).
   - Plots individual missing rates at each step.

3. **`03_rm_mono.sh`**
   - Removes monomorphic SNPs (`--mac 1`) to focus on polymorphic variants.
   - Excludes non-mitochondrial chromosomes if any, and outputs the final quality-controlled PLINK bed/bim/fam and VCF files.

4. **`04_exclude_relate.sh`**
   - Excludes related individuals (e.g., from KING analysis) and non-East Asian ancestry samples based on exclusion files configured in `config.sh`.
   - Runs Principal Component Analysis (PCA) on the study sample to calculate the top 20 mitochondrial PCs.

5. **`05_phasing_imputation.sh`**
   - Performs genotype pre-phasing using **SHAPEIT2**.
   - Filters out variants whose allele frequencies deviate significantly from the Whole Genome Sequencing (WGS) panel (configured via `remove_deviate.py` in the `config.sh` settings).
   - Runs imputation using **IMPUTE2** with the pre-phased haplotypes and the WGS reference panel.

## Additional Helpers

- **`01_idv_missing.r`**: Generates diagnostic histogram plots of sample missing rate distributions.
- **`config.sh`**: Central configurations for environments, tools, reference panels, and path mappings.
