#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
preqcdir <- args[1]
indv_info <- args[2] # tsv with ind_id, age and sex info

# preqcdir <- "/Users/xuanchou/Documents/TWB_1492/microarray/autosomal_qc/king/"
# indv_info <- "/Users/xuanchou/Documents/TWB_1492/microarray/assoc/twb2_survey.final.txt"

library(tidyverse)
setwd(preqcdir)

# Function to remove individuals from the list based on frequency
remove_individuals <- function(data) {
  left_keep_indv <- data[, c("ID1", "ID2")]
  
  final_remove <- c() # Initiate final remove list
  
  while (TRUE) {
    freq_df <- table(c(left_keep_indv$ID1, left_keep_indv$ID2))
    
    if (all(freq_df == 1)) {
      break  # If all the frequency is 1, then break the loop
    }
    
    max_freq_element <- names(freq_df)[which.max(freq_df)]
    final_remove <- c(final_remove, max_freq_element)
    left_keep_indv <- left_keep_indv[!(left_keep_indv$ID1 == max_freq_element | left_keep_indv$ID2 == max_freq_element), ]
  }
  
  final_remove <- c(final_remove, left_keep_indv$ID2)
  
  if (any(duplicated(final_remove))) {
    warning("There are duplicated elements in final_remove.")
  }
  return(final_remove)
}

# Output list for mitochondria =====
# Read individual information
indv_df <- read.table(indv_info, sep = "\t", header = TRUE) %>% 
  select(ID, SEX, AGE)

# Read king table
king_df <- read.table("twb2_btqc_EAS_king_related.kin0", sep = "\t", header = TRUE)

# Merge individual and king dataframes
king_merge <- merge(king_df, indv_df, by.x = "ID1", by.y = "ID", all.x = TRUE) %>%
  merge(., indv_df, by.x = "ID2", by.y = "ID", all.x = TRUE, suffixes = c("_1", "_2"))

# Filter out father-offspring pairs
king_fo <- king_merge %>%
  filter(InfType == "PO") %>%
  filter((AGE_1 > AGE_2 & SEX_1 == "M") | (AGE_1 < AGE_2 & SEX_2 == "M"))
keep_indv <- anti_join(king_merge, king_fo)

# Remove "3rd" from the list
keep_indv <- keep_indv %>%
  filter(InfType != "3rd")

# Write individual list to be removed
final_remove <- remove_individuals(keep_indv)
output_string <- paste(final_remove, final_remove, sep = "\t")
write.table(output_string, "twb2_EAS_king_mt.indvlist", quote = FALSE, col.names = FALSE, row.names = FALSE)

# Output list for autosomes =====
keep_indv <- king_df %>%
  filter(InfType != "3rd")

# Write individual list to be removed separately for autosomes
final_remove <- remove_individuals(keep_indv)
output_string <- paste(final_remove, final_remove, sep = "\t")
write.table(output_string, "twb2_EAS_king.indvlist", quote = FALSE, col.names = FALSE, row.names = FALSE)
