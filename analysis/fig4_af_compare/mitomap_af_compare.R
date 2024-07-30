library(tidyverse)

setwd("/Users/xuanchou/Documents/TWB_1492/Analysis/af_compare/mitomap/")

mitomap_df <- readr::read_delim("mitomap.split.leftaligned.vcf", comment = "##")
twb_df <- readr::read_delim("/Users/xuanchou/Documents/TWB_1492/add_annotation/output_age/combined_sites_only.txt")
twb_df <- twb_df %>% filter(filters == "[]")

# Preprocess MITOMAP df
mitomap_df <- mitomap_df %>% 
  mutate(variant_collapsed = paste0(REF, POS, ALT)) %>% 
  rename(Mitomap_info = INFO)

merged_df <- merge(twb_df, mitomap_df[, c("variant_collapsed", "Mitomap_info")], by = "variant_collapsed", all = TRUE)

# Separate the data
merged_df <- merged_df %>% 
  # Add a row identifier
  mutate(row_id = row_number()) %>% 
  # Separate rows into multiple rows at the semicolon (;)
  separate_rows(Mitomap_info, sep = ";") %>%
  # Separate the keys and values into two new columns based on the equal sign (=)
  separate(Mitomap_info, into = c("key", "value"), sep = "=") %>%
  # Spread the key-value pairs into separate columns
  spread(key = key, value = value, convert = TRUE) %>%
  select(-row_id)

# Recalculate AF, AN is 61,168
merged_df <- merged_df %>% 
  mutate(AF = AC / 61168)

# Remove rows if both AFs in TWB and MITOMAP are 0
merged_df <- merged_df %>% filter(!(AF == 0 & AF_hom == 0))

# Delete GCACA513G cause it had duplicated info
merged_df <- merged_df %>% 
  filter(!((variant_collapsed == "GCACA513G") & (AF == 0)))

# Variants shared between two datasets
shared_df <- merged_df %>% 
  filter(as.numeric(AF_hom) > 0 & as.numeric(AF) > 0)
p_cor <- cor(shared_df$AF_hom, shared_df$AF, method = "pearson", use = "pairwise.complete.obs")
s_cor <- cor(shared_df$AF_hom, shared_df$AF, method = "spearman", use = "pairwise.complete.obs")

# Annotation of merged_df for plotting
merged_df <- merged_df %>% 
  mutate(shared_var = ifelse(variant_collapsed %in% shared_df$variant_collapsed, TRUE, FALSE),
         only_mito = ifelse(AF > 0 & (AF_hom == 0 | is.na(AF_hom)), TRUE, FALSE),
         only_twb = ifelse(AF_hom > 0 & (AF == 0 | is.na(AF)), TRUE, FALSE),
         true_count = shared_var + only_mito + only_twb)
unexpected_rows <- merged_df %>% filter(true_count != 1)
num_unexpected_rows <- nrow(unexpected_rows)
print(unexpected_rows)

# Make point plot for shared and only_twb
plt_df <- merged_df %>% filter(shared_var == TRUE | only_twb == TRUE)
# Replace NA with 0 in AF column when only_twb is TRUE
plt_df <- plt_df %>% 
  mutate(AF = if_else(only_twb & is.na(AF), 0, AF))

# Count the number of each category
category_count <- plt_df %>%
  group_by(only_twb) %>%
  summarise(count = n()) %>%
  mutate(label = if_else(only_twb, paste("Only in TWB (", count, ")", sep = ""), 
                         paste("Shared variants (", count, ")", sep = "")))

# Merge the count data back into the original data frame
plt_df <- plt_df %>%
  left_join(category_count, by = "only_twb")

# Plotting
p <- ggplot(data = plt_df, aes(x = AF_hom, y = AF)) +
  geom_point(aes(color = label)) +
  geom_abline(slope = 1, linetype = "dashed") +
  theme_classic() +
  xlim(0, 1) +
  ylim(0, 1) +
  scale_color_manual(values = c("#9B59B6", "#99A3A4")) +
  theme(legend.position = c(0.05, 1), legend.justification = c(0, 1),
        legend.title = element_blank()) +
  labs(x = "TWB homoplasmic AF", y = "MITOMAP homoplasmic AF", color = "Category")
# Add correlation results to the plot
text_to_add <- paste(
  "Correlation of shared variants\n",
  "Pearson's r: ", round(p_cor, 3), "\n",
  "Spearman's r: ", round(s_cor, 3), "\n",
  sep = ""
)
p <- p + annotate("text", x = 0.7, y = 0.05, hjust = 0, vjust = 0, label = text_to_add, size = 3)
p

ggsave(p, filename = "af_compare.png", dpi = 300, width = 8, height = 4.5)

###################################
# Make point plot for all variants
###################################
# Replace NA values with 0 in AF and AF_hom
plt_df <- merged_df %>% 
  mutate(AF = ifelse(is.na(AF), 0, AF),
         AF_hom = ifelse(is.na(AF_hom), 0, AF_hom),
         category = case_when(
           only_twb ~ "Only in TWB",
           only_mito ~ "Only in MITOMAP",
           shared_var ~ "Shared variants"
         ))

# Count the number of each category
category_count <- plt_df %>%
  group_by(category) %>%
  summarise(count = n()) %>%
  mutate(label = paste(category, " (", count, ")", sep = ""))

# Merge the count data back into the original data frame
plt_df <- plt_df %>%
  left_join(category_count, by = "category")

# Plotting
p <- ggplot(data = plt_df, aes(x = AF_hom, y = AF)) +
  geom_point(aes(color = label)) +
  geom_abline(slope = 1, linetype = "dashed") +
  theme_classic() +
  xlim(0, 1) +
  ylim(0, 1) +
  scale_color_manual(values = c("#5DADE2", "#9B59B6", "#99A3A4")) +
  theme(legend.position = c(0.05, 1), legend.justification = c(0, 1),
        legend.title = element_blank()) +
  labs(x = "TWB homoplasmic AF", y = "MITOMAP homoplasmic AF", color = "Category")
# Add correlation results to the plot
text_to_add <- paste(
  "Correlation of shared variants\n",
  "Pearson's r: ", round(p_cor, 3), "\n",
  "Spearman's r: ", round(s_cor, 3), "\n",
  sep = ""
)
p <- p + annotate("text", x = 0.7, y = 0.05, hjust = 0, vjust = 0, label = text_to_add, size = 3)
p

ggsave(p, filename = "af_compare_all.png", dpi = 300, width = 8, height = 4.5)
