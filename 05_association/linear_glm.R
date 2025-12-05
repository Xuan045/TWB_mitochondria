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

# Create the output directory if it does not exist
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

source("/staging/biology/u4432941/TWB1492_mt/microarray_association/phewas_analysis.R")
source("/staging/biology/u4432941/TWB1492_mt/microarray_association/phewas_plotting.R")

#Set the random seed so it is replicable
set.seed(2005)

# Load genotype file
message("Get data")
genotype <- read.delim(genotype_file, header = TRUE, stringsAsFactors = FALSE, na.strings = ".")
colnames(genotype) <- gsub("^X", "m.", colnames(genotype))
names(genotype)[1] = "ID"

# Load phenotype file
full_df <- read.table(phenotype_file, header = TRUE)
colnames(full_df)[colnames(full_df) == "Sample_Name"] <- "ID"

# Add PCA eigenvec
# Store proper column names
pca_df <- read.table(pca_file, header = FALSE)[,-1]
col_names <- c("ID", paste0("PC", 1:(ncol(pca_df) - 1)))
colnames(pca_df) <- col_names
pca_df$ID <- sub("SM_", "", pca_df$ID)

# Sort phenotypes
## Select target columns
select_col <- c("BODY_HEIGHT", "BODY_WEIGHT", "BMI", "BODY_FAT_RATE", "BODY_WAISTLINE", "BODY_BUTTOCKS", "WHR", "T_SCORE",
                "Z_SCORE", "FVC", "FEV10", "FEV10_FVC", "SYSTOLIC_PRESSURE", "DIASTOLIC_PRESSURE", "HEARTBEAT_SPEED", 
                "RBC", "WBC", "HB", "HCT", "PLATELET", "HBA1C", "AC_GLUCOSE", "T_CHO", "TG", "HDL_C", "LDL_C", 
                "T_BILIRUBIN", "ALBUMIN", "SGOT", "SGPT", "GAMMA_GT", "AFP", "BUN", "CREATININE", "URIC_ACID",
                "MICROALB", "CREATININE_URINE", "egfr_cr")
select_df <- full_df#  %>% select(ID, AGE, SEX, select_col)
genotype <- genotype %>% filter(ID %in% select_df$ID)

## For "SEX" column, convert 1 to M; 2 to F
## Join with genotype
select_df <- select_df %>% 
  mutate(SEX = recode(SEX, "1" = "M", "2" = "F")) 
select_df <-  inner_join(select_df, genotype, by="ID") %>% 
  mutate_at(vars(select_col), as.numeric)
select_df <- inner_join(select_df, pca_df, by = "ID")

geno_cols <- names(genotype)[-1]

###################
# Perform PheWAS
###################
outfull <- paste0(outdir, file_prefix)

message("Perform association tests...")
final_df <- perform_phewas_analysis(select_df, select_col, geno_cols, c("AGE", "SEX", "TWB2_Batch", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), min_cases = 1000, age_sex_cov = TRUE)

# Make pheWAS manhatton plot
message("Make PheWAS plots...")
plot_phewas_linear(final_df, select_col, geno_cols, paste0(outfull, "_linear"))

# df to output
message(paste0("Writing summary statistics to ", outfull, "_linear_summary.csv"))
out_df <- final_df %>% arrange(p)
write_csv(out_df, file = paste0(outfull, "_linear_summary.csv"))
