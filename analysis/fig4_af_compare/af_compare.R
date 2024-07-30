library(ggpubr)
library(tidyverse)

setwd("/Users/xuanchou/Documents/TWB_1492/Analysis/af_compare/")
population_col <- "eas"
population_hom_col <- paste0("gnomAD_MT_AF_hom_", population_col)
population_het_col <- paste0("gnomAD_MT_AF_het_", population_col)
annotate_df <- read.delim("/Users/xuanchou/Documents/TWB_1492/Analysis/af_compare/TWB1492_woInfoGT.vcf_allAF.tsv")

# Remove duplicated orws
duplicated_rows <- duplicated(annotate_df[, c("POS", "REF", "ALT")])
annotate_df <- annotate_df[!duplicated_rows, ]

# Select pass variants
annotate_df <- annotate_df %>% 
  filter(TWB1465_MT_FILTER == "PASS") %>% 
  filter(gnomAD_MT_FILTER == "PASS")

# Plot the AF compare with all
long_table_hom <- annotate_df %>%
  select(CHROM:Allele, population_hom_col, TWB1465_MT_AF_hom) %>% 
  mutate(type = "Homoplasmy") %>% 
  rename(gnomAD_eas = population_hom_col, TWB = TWB1465_MT_AF_hom)

long_table_het <- annotate_df %>%
  select(CHROM:Allele, population_het_col, TWB1465_MT_AF_het_vaf10) %>% 
  mutate(type = "Heteroplasmy") %>% 
  rename(gnomAD_eas = population_het_col, TWB = TWB1465_MT_AF_het_vaf10)

long_table <- rbind(long_table_hom, long_table_het)

# Remove rows if gnomAD_eas is equal to 0
long_table <- long_table %>% 
  filter(as.numeric(gnomAD_eas) != 0) %>% 
  filter(as.numeric(TWB) != 0)

# Make the AF scatter plot
unique_combinations <- unique(long_table[, c("POS", "REF", "ALT")])
unique_count <- nrow(unique_combinations)
n_hom <- sum(long_table$type == "Homoplasmy")
n_het <- sum(long_table$type == "Heteroplasmy")

long_table <- long_table %>% 
  mutate(label = ifelse(type == "Homoplasmy", paste0("Homplasmy (n=", n_hom, ")"), paste0("Heteroplasmy (n=", n_het, ")")))

af_plt <- ggplot(long_table, aes(x = TWB, y = as.numeric(gnomAD_eas))) +
  geom_point(aes(color = label)) +
  theme_classic() +
  geom_abline(slope = 1, linetype = "dashed") +
  labs(x = "TWB1465", y = paste0("gnomAD_", population_col),
       title = paste0("TWB1465 vs. gnomAD_", population_col),
       subtitle = paste0(unique_count, " shared variants, ", "r = ", round(cor(long_table$TWB, as.numeric(long_table$gnomAD_eas)), 3))) +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.2),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 15),
        axis.title = element_text(size = 15, face = "bold"))

# Calculate R-square value for each 'type'
correlations <- long_table %>%
  group_by(label) %>%
  summarise(R = cor(TWB, as.numeric(gnomAD_eas)))

# Make the AF scatter plot and add correlation text
af_facet <- ggplot(long_table, aes(x = TWB, y = as.numeric(gnomAD_eas))) +
  geom_point(aes(color = type)) +
  geom_abline(slope = 1, linetype = "dashed") +
  theme_classic() +
  theme(legend.position = "none",
        legend.key = element_blank(),
        axis.title = element_blank()) +
  facet_wrap(~label) +
  geom_text(data = correlations, aes(label = paste("r = ", round(R, 3)), x = Inf, y = -Inf, hjust = 1, vjust = -0.5))

af_combined <- ggarrange(af_plt, af_facet, nrow = 2)
ggsave(af_combined, filename = paste0("af_combined_", population_col, ".pdf"), dpi = 300, width = 8, height = 8)


######################################
# Categorize into three groups in EAS
######################################
annotate_df <- read.delim("/Users/xuanchou/Documents/TWB_1492/custom_annotation/TWB1492_gnomAD_twb_mt.tsv")

# Remove duplicated orws
duplicated_rows <- duplicated(annotate_df[, c("POS", "REF", "ALT")])
annotate_df <- annotate_df[!duplicated_rows, ]

# Select pass variants
annotate_df <- annotate_df %>% 
  filter(TWB1465_MT_FILTER == "PASS") %>% 
  filter(gnomAD_MT_FILTER %in% c("PASS", "."))

# Group into three categories and compare each af
# 1. Separate hom and het into two tables and then combine together
long_table_hom <- annotate_df %>%
  select(CHROM:Allele, population_hom_col, gnomAD_MT_FILTER, TWB1465_MT_AF_hom) %>% 
  mutate(type = "Homoplasmy") %>% 
  rename(gnomAD = population_hom_col, TWB = TWB1465_MT_AF_hom) %>% 
  filter(TWB != 0)

long_table_het <- annotate_df %>%
  select(CHROM:Allele, population_het_col, gnomAD_MT_FILTER, TWB1465_MT_AF_het_vaf10) %>% 
  mutate(type = "Heteroplasmy") %>% 
  rename(gnomAD = population_het_col, TWB = TWB1465_MT_AF_het_vaf10) %>% 
  filter(TWB != 0)

long_table <- rbind(long_table_hom, long_table_het)

# 2. Categorize the variants into differenct groups
# a. Shared variants (grey)
# b. Only in TWB (hom) (blue)
# c. Only in TWB (het). (purple)
long_table$gnomAD <- replace(long_table$gnomAD, long_table$gnomAD == ".", 0) %>% as.numeric()
long_table$gnomAD_MT_FILTER <- replace(long_table$gnomAD_MT_FILTER, long_table$gnomAD_MT_FILTER == ".", NA)

long_table <- long_table %>% 
  mutate(group = case_when(
    TWB > 0 & gnomAD > 0 & gnomAD_MT_FILTER == "PASS" ~ "Shared variants",
    type == "Homoplasmy" ~ "Only in TWB (hom)",
    TRUE ~ "Only in TWB (het)"
  ))
# Add count of each group for plotting
# For shared variants are both hom and het, those should be counted one time only
category_count <- long_table %>%
  filter(group == "Shared variants") %>% 
  distinct(POS, REF, ALT, .keep_all = TRUE) %>% # Remove duplicated combinations
  bind_rows(long_table %>% filter(group != "Shared variants")) %>% # Add other groups back to the df
  group_by(group) %>%
  summarise(count = n()) %>%
  mutate(group_label = paste(group, " (", count, ")", sep = ""))
plot_table <- long_table %>% 
  left_join(category_count, by = "group")

# 3. Plot the point plot to compare the AFs between two datasets
af_plt <- ggplot() +
  geom_point(data = filter(plot_table, group == "Shared variants"), aes(x = TWB, y = gnomAD, color = group_label)) +
  geom_point(data = filter(plot_table, group != "Shared variants"), aes(x = TWB, y = gnomAD, color = group_label)) +
  geom_abline(slope = 1, linetype = "dashed") +
  theme_classic() +
  scale_color_manual(values = c("#5DADE2", "#9B59B6", "#99A3A4")) +
  scale_alpha_manual(values = c(1, 1, 0.5)) +
  labs(y = "gnomAD") +
  theme(legend.position = c(0.05, 1), legend.justification = c(0, 1),
        legend.title = element_blank())
af_plt

# Add correlation results to the plot
shared_df <- long_table %>% filter(group == "Shared variants")
text_to_add <- paste(
  "Correlation of shared variants\n",
  "Pearson's r: ", round(cor(shared_df$gnomAD, shared_df$TWB, method = "pearson", use = "pairwise.complete.obs"), 3), "\n",
  sep = ""
)
af_plt <- af_plt + annotate("text", x = 0.7, y = 0.05, hjust = 0, vjust = 0, label = text_to_add, size = 3.5)
af_plt

category_facet <- ggplot(plot_table, aes(x = TWB, y = gnomAD, color = group_label)) +
  geom_point() +
  scale_color_manual(values = c("#5DADE2", "#9B59B6", "#99A3A4")) +
  facet_wrap(~group) +
  theme_classic()
af_combine <- ggarrange(af_plt, category_facet, ncol = 1, common.legend = TRUE)
af_combine

category_combine_var_type <- ggarrange(af_plt, af_facet, ncol = 1)
category_combine_var_type
ggsave(category_combine_var_type, filename = "category_combine_var_type.png", dpi = 300, width = 6.5, height = 4.5)
