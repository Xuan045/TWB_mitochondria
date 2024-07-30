library(tidyverse)

# Set the random seed for reproducibility
set.seed(2005)

# Set the working directory
setwd("/Users/xuanchou/Documents/TWB_1492/microarray/assoc/kidney/conditional_analysis/")

# Load the data
kidney_df <- read_delim("/Users/xuanchou/Documents/TWB_1492/microarray/assoc/kidney/twb2_kidney_var.csv", delim = ",") 
hg_df <- read_delim("/Users/xuanchou/Documents/TWB_1492/microarray/haplogroup/twb2_hg_sorted.csv", delim = ",") %>%
  select(-c(Haplogroup, Macrohaplogroup, Haplogroup_digit2)) %>%
  rename(Haplogroup = Haplogroup_digit3)

# Merge and prepare the data
merged_df <- merge(kidney_df, hg_df, by.x = "ID", by.y = "SampleID") %>%
  select(-all_of(setdiff(grep("^m\\.", colnames(kidney_df), value = TRUE), c("m.16217C", "m.827G", "m.15535T", "m.6216C", "m.6023A", "m.6413C")))) %>%
  mutate(across(grep("^m\\.", colnames(.), value = TRUE), ~as.integer(substr(., 1, 1)) + as.integer(substr(., 3, 3))))

# Define phenotypes and covariates, convert factors
merged_df$SEX <- as.factor(merged_df$SEX)
merged_df$TWB2_Batch <- as.factor(merged_df$TWB2_Batch)
merged_df$Haplogroup <- as.factor(merged_df$Haplogroup)
merged_df$is_B4b <- ifelse(merged_df$Haplogroup == "B4b", 1, 0)

# Define a function to fit models for both eGFRcr and CREATININE
fit_models <- function(phenotype) {
  variant_results <- lapply(names(merged_df)[grepl("^m\\.", names(merged_df))], function(variant) {
    print(paste("Processing:", variant, "for", phenotype))
    formula <- as.formula(paste(phenotype, "~", variant, "+ is_B4b + TWB2_Batch + AGE + SEX"))
    model <- glm(formula, data = merged_df, family = gaussian())
    summary_model <- summary(model)
    
    # Return the coefficients for the variant
    list(
      Variant = variant,
      Phenotype = phenotype,
      Coeff = summary_model$coefficients[variant, "Estimate"],
      SE = summary_model$coefficients[variant, "Std. Error"],
      P_value = summary_model$coefficients[variant, "Pr(>|t|)"]
    )
  })
  return(variant_results)
}

# Fit models and collect results
results_eGFRcr <- fit_models("eGFRcr")
results_creatinine <- fit_models("CREATININE")

# Combine and adjust results
results_df <- bind_rows(lapply(c(results_eGFRcr, results_creatinine), function(x) {
  data.frame(
    Variant = x$Variant,
    Phenotype = x$Phenotype,
    Coeff = x$Coeff,
    SE = x$SE,
    P_value = x$P_value
  )
}))
results_df$P_adjusted <- p.adjust(results_df$P_value, method = "fdr")

# Write the results to a file
write_csv(results_df, "combined_results.csv")

# Print a summary to check the table structure
print(head(results_df))

formatted_sorted_results_df <- results_df %>%
  mutate(
    Coeff = round(Coeff, 4),
    SE = round(SE, 4),
    P_value = round(P_value, 4)
  ) %>%
  arrange(P_value, Coeff, SE)

# Write the formatted and sorted results to a new CSV file
write_csv(formatted_sorted_results_df, "formatted_sorted_results.csv")

# Print a summary to check the formatted and sorted table structure
print(head(formatted_sorted_results_df))
