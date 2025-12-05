# Mitochondrial Variation in Taiwan Biobank

Code and workflows to reproduce the analyses in:

> **Chou, T.-H.**, Chien, P.-M., Chen, P.-X., Lee, N.-C., Wang, H.-Y., Chang, Y.-C., Chen, C.-Y., Chen, P.-L., & Hsu, J. S.-J. (2025).  
> *Mitochondrial variation in Taiwan Biobank reveals ancestry structure and trait associations.*  
> iScience, 28(11), 113746. https://doi.org/10.1016/j.isci.2025.113746 

Please cite the paper above if you use this code in your work.

---

## Overview

This repository contains scripts and workflows used to:

1. Call mitochondrial DNA (mtDNA) variants from whole-genome sequencing (WGS) data in the Taiwan Biobank (TWB).
2. Assign mitochondrial haplogroups and analyze mtDNA population structure.
3. Impute mtDNA variants into genotyping array samples using a TWB-specific mtDNA reference panel.
4. Perform mitochondrial genome-wide association studies (mtGWAS) across 86 traits.
5. Test haplogroup-based associations with renal function biomarkers.
6. Generate the main figures and tables in the manuscript.

**Important:**  
Raw TWB data (WGS, genotyping arrays, phenotypes) are **not** included in this repository due to data sharing restrictions.  
Researchers interested in accessing TWB data must apply through the Taiwan Biobank.
<!--
---

## Repository Structure

```
env/                    # Conda/R environment files
config/                 # Example path config, phenotype definitions
workflow/               # High-level sequential workflows
scripts/                # Individual analysis and plotting scripts
example_data/           # Toy VCF/phenotype for demonstration (no TWB data)
docs/                   # Figure mapping and extra notes
```

---

## Data Availability

### Taiwan Biobank
All WGS, array genotype, and phenotype data used in this study are available upon application from the **Taiwan Biobank**.  
Researchers must apply directly to TWB for access.

### External Data Sources
The workflow requires several public datasets (not included here), such as:  
- gnomAD v3.1 mitochondrial DNA calls  
- 1000 Genomes Project Phase 3 (mtDNA and nuclear variants)  
- HaploGrep2 reference files  

Configure paths in `config/paths_example.yaml`.

---

## Software Requirements

This analysis was performed on a Linux HPC environment.

### Core Tools
- **GATK** 4.2.3.0 – Mitochondrial SNP/Indel calling pipeline  
- **BWA-MEM** 0.7.17 – Alignment  
- **bcftools / htslib / vcftools** – VCF processing  
- **PLINK** 1.9 / 2.0 – Covariates & QC  
- **HaploGrep2** – Haplogroup assignment  
- **R** ≥ 4.0 – Statistical analysis and visualization  
  - Required packages: tidyverse, data.table, pegas, ade4, vcfr, poppr, etc.

-->