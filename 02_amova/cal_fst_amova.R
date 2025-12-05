library(ape)
library(adegenet)
library(poppr)
library(pegas)
library(dplyr)
library(vcfR)
library(pheatmap)


setwd("/Users/tingxuan/Documents//TWB_1492/amova/")
set.seed(2025)

# =================================
# AMOVA analysis
# =================================
# Load and filter
pop_df <- read.csv("/Users/tingxuan/Documents/Github/twb_data_cleaning/twb_ancestry/maternal_origin.csv") %>%
  filter(!is.na(Sample_Name) & Sample_Name != "" &
           Platform != "Ion Proton" & NATIVE_MOM != "Undefined China") %>%
  select(Sample_Name, NATIVE_MOM) %>%
  mutate(NATIVE_MOM = recode(NATIVE_MOM,
                             "Southern China" = "Southern Chinese",
                             "Northern China" = "Northern Chinese"))

# Object to store pop assignments
pop_assignments <- pop_df$NATIVE_MOM
names(pop_assignments) <- pop_df$Sample_Name

pop_table_ids <- pop_df$Sample_Name

# ############ FASTA solution ############
# seq <- read.FASTA("mtDNA_aligned.fasta")
# 
# # Filter samples present in only one data type
# fasta_ids <- labels.DNAbin(seq)
# 
# ## Identify IDs specific in FASTA or pop_df
# common_ids <- intersect(fasta_ids, pop_table_ids)
# ids_in_fasta_only <- setdiff(fasta_ids, pop_table_ids)
# ids_in_pop_only <- setdiff(pop_table_ids, fasta_ids)
# 
# ## Keep only the IDs that have matching information
# seq_filtered <- seq[common_ids]
# pop_df_filtered <- pop_assignments[common_ids]
# 
# 
# # Convert DNAbin to genind object
# # For mitochondrial DNA, ploidy is 1 (haploid).
# genind_data <- DNAbin2genind(seq_filtered, pop = pop_df_filtered)
# ############ FASTA solution ############

############ VCF solution ############
vcf <- read.vcf("TWB1465_wgs_recode.snps.vcf.gz") 
vcfR <- read.vcfR("TWB1465_wgs_recode.snps.vcf.gz") #vcfR package

# Annotate POS
colnames(vcf) <- getPOS(vcfR)

# Keep only the IDs that have matching information
common_ids <- intersect(rownames(vcf), pop_table_ids)
ids_in_fasta_only <- setdiff(rownames(vcf), pop_table_ids)
ids_in_pop_only <- setdiff(pop_table_ids, rownames(vcf))

# Filter to the vcf and pop df
vcf_filtered <- vcf[rownames(vcf) %in% common_ids,]
pop_df_filtered <- pop_assignments[common_ids]

# Convert to genind object and assign populations
genind_data <- loci2genind(vcf_filtered, ploidy = 1) #pegas package
pop(genind_data) <- pop_df_filtered

# Convert to genclone object (recommended for poppr analyses)
gc_data <- as.genclone(genind_data) #poppr package

print(gc_data)

# Calculate AMOVA within and between populations
amova_result <- poppr.amova(gc_data, ~pop)

# Print the AMOVA results
print(amova_result)

# To test the significance of the results, perform a permutation test.
# nrepet = 1000 is a common number of permutations for statistical significance.
amova_test <- randtest(amova_result, nrepet = 1000)

# Print the permutation test results
print(amova_test)

# Plot the permutation test results to visualize the observed statistic
# against the distribution of randomized statistics.
plot(amova_test)

# Create an AMOVA output
report <- c(
  "AMOVA Analysis Report",
  "====================",
  paste("Date:", Sys.Date()),
  paste("Seed:", 2025),
  "",
  "Data Summary:",
  paste("Total samples analyzed:", nInd(gc_data)),
  paste("Number of populations:", nPop(gc_data)),
  paste("Population names:", paste(popNames(gc_data), collapse = ", ")),
  "",
  "AMOVA Results:",
  capture.output(print(amova_result)),
  "",
  "Permutation Test Results:",
  capture.output(print(amova_test)),
  "",
  "Statistical Significance:",
  paste("P-value:", amova_test$pvalue),
  ifelse(amova_test$pvalue < 0.05, "Result is statistically significant (p < 0.05)", 
         "Result is not statistically significant (p >= 0.05)")
)

writeLines(report, "amova_report.txt")


# =====================
# Pariwise calculation
# =====================
vcfR <- read.vcfR("TWB1465_wgs_recode.snps.haploid.vcf.gz") #vcfR package
pairwise_gst <- pairwise_genetic_diff(vcfR, pops = as.factor(pop_df_filtered))
gst_values <- colMeans(pairwise_gst[, c(4:ncol(pairwise_gst))], na.rm = TRUE)

# Define pop order
pop_order <- c("Holo", "Hakka", "Southern Chinese", "Northern Chinese")

# Create gst matrix for plotting
gst_matrix <- matrix(0, nrow = 4, ncol = 4, dimnames = list(pop_order, pop_order))
gst_matrix["Holo", "Southern Chinese"] <- gst_values["Gst_Holo_Southern Chinese"]
gst_matrix["Holo", "Hakka"] <- gst_values["Gst_Hakka_Holo"]
gst_matrix["Holo", "Northern Chinese"] <- gst_values["Gst_Holo_Northern Chinese"]
gst_matrix["Hakka", "Southern Chinese"] <- gst_values["Gst_Hakka_Southern Chinese"]
gst_matrix["Southern Chinese", "Northern Chinese"] <- gst_values["Gst_Northern Chinese_Southern Chinese"]
gst_matrix["Hakka", "Northern Chinese"] <- gst_values["Gst_Hakka_Northern Chinese"]

# Mirror to lower triangle
gst_matrix[lower.tri(gst_matrix)] <- t(gst_matrix)[lower.tri(gst_matrix)]

# Plot heatmap
pdf("mt_pairwise_gst.pdf", width = 6, height = 4)
pheatmap(gst_matrix,
         display_numbers = TRUE,
         color = colorRampPalette(c("lightblue", "white", "salmon"))(100),
         number_format = "%.4f",
         number_color = "black",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = expression(bold("Pairwise Population G"[ST])),
         fontsize_row = 10,
         fontsize_col = 10,
         fontsize_number = 8,
         angle_col = 45,
         border_color = "grey60")
dev.off()

