library(tidyverse)
library(ggplot2)
library(ggrepel)

# Parse command line arguments
# args[1] is expected to be the path to the genotype file
# args[2] is expected to be the path to the phenotype file
# args[3] is expected to be the path to the PCA file
# args[4] is expected to be the path to the output directory
args <- commandArgs(trailingOnly = TRUE)

# Assign each argument to a variable
genotype_file <- args[1]
phenotype_file <- args[2]
pca_file <- args[3]
outdir <- args[4]
file_prefix <- args[5]
min_cases <- args[6]

source("/staging/biology/u4432941/TWB1492_mt/microarray_association/phewas_analysis.R")
source("/staging/biology/u4432941/TWB1492_mt/microarray_association/phewas_plotting.R")

# Set the random seed so it is replicable
set.seed(2005)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# Load genotype file
message("Get genotype data")
genotype <- read.delim(genotype_file, header = TRUE, stringsAsFactors = FALSE, na.strings = ".")
colnames(genotype) <- gsub("^X", "m.", colnames(genotype))
names(genotype)[1] = "ID"

# Add PCA eigenvec
# Store proper column names
message("Get PCA data")
pca_df <- read.table(pca_file, header = FALSE)[,-1]
col_names <- c("ID", paste0("PC", 1:(ncol(pca_df) - 1)))
colnames(pca_df) <- col_names
pca_df$ID <- sub("SM_", "", pca_df$ID)

####################
# Sort phenotypes
####################
message("Get phenotype data")
select_df <- read.delim(phenotype_file, header = TRUE)
select_df$ID <- gsub("-", "_", select_df$ID)
dz_cols <- colnames(select_df)[-c(1:3)]

## Join with genotype
select_df <-  inner_join(select_df, genotype, by="ID")
## Combine with pca_df
select_df <- inner_join(select_df, pca_df, by = "ID")

#############################
# Conduct association study
#############################
outfull <- paste0(outdir, file_prefix)

message("Perform association tests...")
geno_cols <- names(genotype)[-1]
covariate_cols <- c("AGE", "SEX", "TWB2_Batch", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
final_df <- perform_phewas_analysis(select_df, dz_cols, geno_cols, covariates = covariate_cols, min_cases = min_cases, age_sex_cov = TRUE)

# df to output
message(paste0("Writing summary statistics to ", outfull, "_logistic_summary.csv"))
out_df <- final_df %>% arrange(p, desc(OR))
write_csv(out_df, file = paste0(outfull, "_logistic_summary.csv"))

# Make pheWAS manhatton plot
plot_phewas_logistic(final_df, dz_cols, geno_cols, paste0(outfull, "_logistic"))
