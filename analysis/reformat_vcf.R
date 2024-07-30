# Split VEP annotations and store to a new table
# Extract GT info in every sample, so that they can be analyzed based on variant

library(tidyverse)
library(jsonlite)

setwd("/Users/xuanchou/Documents/TWB_1492/Analysis/")
# quote = "" to keep "" around character fields
variant_table <- read.delim("/Users/xuanchou/Documents/TWB_1492/add_annotation/output_age/combined_sites_only.txt", header = TRUE, quote = "")

# Import sample VCF
vcf_lines <- readLines("/Users/xuanchou/Documents/TWB_1492/add_annotation/output_age/sample_vcf.vcf")
vcf_lines <- vcf_lines[grep("^#CH", vcf_lines)] %>% strsplit("\t") %>% unlist()
sample_vcf <- read.table("/Users/xuanchou/Documents/TWB_1492/add_annotation/output_age/sample_vcf.vcf", header = FALSE, 
                         sep = "\t", col.names = vcf_lines)
sample_list <- vcf_lines[grep("^NGS", vcf_lines)]

# Create a new dataframe to store VEP information
variant_vep <- fromJSON(variant_table[1, 'vep'], simplifyDataFrame = TRUE)$transcript_consequences %>% as.data.frame()

# Extract VEP info to new columns
for(i in 2:nrow(variant_table)){
  # Parse the JSON string
  json_str <- fromJSON(variant_table[i, "vep"], simplifyDataFrame = TRUE)$transcript_consequences
 
  # JSON to dataframe
  vep_table <- as.data.frame(json_str)
  
  # Cobine with the new dataframe
  variant_vep <- rbind(variant_vep, vep_table)
} 

# Combine with original variant table to keep important information
variant_vep <- cbind(variant_table, variant_vep) %>% select(-c(vep, allele))

write.table(variant_vep, file = "/Users/xuanchou/Documents/TWB_1492/Analysis/twb1465_variant_vep.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#######################################################################
# Reconstruct the VCF table. Separate every GT info into a new column.
#######################################################################
sample_df <- sample_vcf[, sample_list]

extract_and_add_values <- function(df, extract_index, new_col_name) {
  # :param extract_index: which value should be extracted
  #: param new_col_name: the column name that extracted values will be added
  
  df[[new_col_name]] <- apply(df, 1, function(x){
    # Split fields by ":"
    field_values <- strsplit(x, ":")
    # Extract the desired value
    extracted_value <- sapply(field_values, function(y) y[extract_index])
    # Concatenate the extracted values into a string separate by coma
    paste(extracted_value, collapse = ",")
  })
  return(df)
}

# Extract the values based on the information in the FORMAT field and add them to a new column.
format_info <- sample_vcf[1, "FORMAT"] %>% strsplit(":") %>% unlist()
format_info <- format_info[1:length(format_info)-1]

for (i in seq_along(format_info)) {
  # Get the format info that will be used as new column name
  info <- format_info[i]
  # Add the new column using the function "extract_and_add_values"
  sample_df <- extract_and_add_values(sample_df, i, info)
}

# Combine to the original sample_vcf
sample_vcf <- cbind(sample_vcf, sample_df[, (ncol(sample_df) - length(format_info) + 1):ncol(sample_df)])

# Confirm that if the number of combined values equal to the sample count
num_fields <- sample_vcf[440, "GT"] %>% strsplit(",") %>% unlist() %>% length()
num_sample <- length(sample_list)
if(num_fields != num_sample) {
  stop("Error: Number of combined values is not equal to the sample count")
}

# Combine with VEP table
sample_vcf <- sample_vcf %>% 
  mutate(variant_collapsed = paste(REF, POS, ALT, sep = "")) %>% 
  select(-c(X.CHROM, ID, QUAL, FILTER, INFO, REF, ALT))
final_combined <- merge(variant_vep, sample_vcf, by = "variant_collapsed")
final_combined <- final_combined[order(final_combined$POS),]

# Clean the table
colnames(final_combined) <- gsub("X.", "", colnames(final_combined))
# Move the columns that column name start with "NGS" to the last lines
new_order <- c(grep("^NGS", names(final_combined), value = TRUE, invert = TRUE), grep("^NGS", names(final_combined), value = TRUE))
final_combined <- select(final_combined, all_of(new_order))

# Output sample_order.txt and reconstructed.txt
cat(paste("sample order: ", paste(sample_list, collapse = ",")), file = "/Users/xuanchou/Documents/TWB_1492/Analysis/reconstructed_vcf_sample_order.txt")
write.table(final_combined, file = "/Users/xuanchou/Documents/TWB_1492/Analysis/reconstructed_vcf.txt", sep = "\t", row.names = FALSE, quote = FALSE)
