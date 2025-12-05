# Mitochondrial PheWAS Analysis

This folder contains scripts and resources for performing Mitochondrial Phenome-Wide Association Analysis (PheWAS) using generalized linear models and generating relevant plots.

## Folder Contents

- **run_assoc_r.sh**: A shell script for running association analyses with the provided R scripts.
- **linear_glm.R**: This script performs linear regression analyses for continuous phenotypes.
- **logistic_glm.R**: This script conducts logistic regression analyses for binary phenotypes.
- **phewas_plotting.R**: This script generates plots for visualizing PheWAS results.
- **phewas_analysis.R**: This comprehensive script orchestrates the PheWAS analysis, incorporating both linear and logistic models and generating summary statistics.

## Prerequisites

Before running the scripts, ensure that you have the following installed:

- R (version 4.0.0 or higher)
- The following R packages:
  - `dplyr`
  - `ggplot2`
  - `data.table`
  - `broom`
  - `tidyverse`

## Usage

1. **Running All Analyses**:
    - The `run_assoc_r.sh` script can be used to run all the association analyses in one go. Ensure the script has execute permissions:
      ```sh
      chmod +x run_assoc_r.sh
      ./run_assoc_r.sh
      ```
    - This script will perform all the necessary analysis steps, including linear and logistic regressions, and generate summary statistics.

2. **Linear Regression Analysis**:
    - The `linear_glm.R` script can perform linear regression analysis for continuous phenotypes.

3. **Logistic Regression Analysis**:
    - The `logistic_glm.R` script can conduct logistic regression analysis for binary phenotypes.

4. **PheWAS Analysis**:
    - The `phewas_analysis.R` script can run a full PheWAS analysis. This script integrates both linear and logistic regression analyses.

5. **Plotting Results**:
    - The `phewas_plotting.R` script can generate plots for visualizing PheWAS results.
