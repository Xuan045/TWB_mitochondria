library(tidyverse)
library(ggpubr)

# Set working directory and read data
setwd("/Users/xuanchou/Documents/TWB_1492/Analysis/variant_stat")
sample_annote_vaf10 <- read.delim("/Users/xuanchou/Documents/TWB_1492/add_annotation/output_age/sample_annotations.txt")

plot_a_theme <- theme(axis.title = element_text(colour = "black", size = 14, face = "bold"),
                      axis.text = element_text(colour = "black", size = 11),
                      axis.ticks.x = element_blank(),
                      axis.line.x = element_blank(),
                      axis.line.y = element_line(size = 1))
other_theme <- theme(axis.title = element_text(colour = "black", size = 14, face = "bold"),
                     axis.text = element_text(colour = "black", size = 11),
                     axis.ticks.x = element_blank(),panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) 

# Create a function to generate plots
generate_plots <- function(df, para) {
  # Add the number of heteroplasmic SNPs and indels to get the total number of heteroplasmies
  df <- mutate(df, n_total_het = n_snp_het + n_indel_het)
  
  # Count the number of occurrences for each total number of heteroplasmies
  twb_hets <- df %>% group_by(n_total_het) %>% summarise(count = n()) %>% ungroup()
  twb_hets <- mutate(twb_hets, percent = count / sum(count) * 100)

  # Add the number of homoplasmic SNPs and indels to get the total number of homoplasmies
  df <- mutate(df, n_total_hom = n_snp_hom + n_indel_hom)

  # Count the number of occurrences for each total number of homoplasmies
  twb_homs <- df %>% group_by(n_total_hom) %>% summarise(count = n()) %>% ungroup()
  twb_homs <- mutate(twb_homs, percent = count / sum(count) * 100)
  
  # Combine heteroplasmic and homoplasmic data
  twb_homs <- mutate(twb_homs, type = "Homoplasmic") %>% rename(n_total = n_total_hom)
  twb_hets <- mutate(twb_hets, type = "Heteroplasmic") %>% rename(n_total = n_total_het) 
  twb_table <- rbind(twb_homs, twb_hets)
  
  hom_het_hist <- ggplot(twb_table, aes(x = n_total, y = percent, fill = type)) +
    facet_grid(type ~ .) +
    geom_histogram(stat = "identity") +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.85, 0.85),
          legend.key.size = unit(16, "pt"), legend.text = element_text(size = 15)) +
    plot_a_theme +
    labs(x = "# variants / sample", y = "Percent samples", title = "") 
  
  return(list(hom_het_hist = hom_het_hist, sort_table = twb_table))
}

# Generate plots for vaf10 data
plots_vaf10 <- generate_plots(sample_annote_vaf10, "vaf10_")

# Data preparation function
prepare_data <- function(df) {
  df %>%
    mutate(age_bin = cut(age, breaks = seq(min(age) - 1, max(age) + 1, by = 3), labels = FALSE)) %>% 
    group_by(age_bin) %>% 
    summarise(mean_het_snp = mean(n_snp_het, na.rm = TRUE),
              sem_het_snp = sd(n_snp_het, na.rm = TRUE) / sqrt(n()),
              mean_het_indel = mean(n_indel_het, na.rm = TRUE),
              sem_het_indel = sd(n_indel_het, na.rm = TRUE) / sqrt(n()),
              age_label = mean(unique(age), na.rm = TRUE))
}

# Plot function
plot_data <- function(data, y_value, y_error, title, y_label) {
  ggplot(data, aes(x = age_label, y = y_value)) +
    geom_point() +
    geom_errorbar(aes(ymin = y_value - y_error, ymax = y_value + y_error), width = 0.2) +
    theme_classic() +
    other_theme +
    scale_y_continuous(limits = c(0, max(y_value) + 0.5)) +
    labs(x = "Age bin", y = y_label, title = title)
}

# Prepare data for age bins
age_bin_vaf10 <- prepare_data(sample_annote_vaf10)

# Generate plots for SNV and Indel
snv_plot <- plot_data(age_bin_vaf10, age_bin_vaf10$mean_het_snp, age_bin_vaf10$sem_het_snp, "SNV", "Mean heteroplasmy count (VAF ≥ 10%)")
indel_plot <- plot_data(age_bin_vaf10, age_bin_vaf10$mean_het_indel, age_bin_vaf10$sem_het_indel, "Indel", "Mean heteroplasmy count (VAF ≥ 10%)")

# Combine the plots and add labels
combined_plots <- ggarrange(plots_vaf10$hom_het_hist, 
                            ggarrange(snv_plot, indel_plot, ncol = 2, labels = c("b", "c"), font.label = list(size = 20, family = "Times New Roman")),
                            labels = c("a"), ncol = 1, nrow = 2, font.label = list(size = 20, family = "Times New Roman"))

print(combined_plots)

# Create the thesis folder if it doesn't exist
dir.create("thesis", showWarnings = FALSE)

# Save the combined plot to the thesis folder
ggsave("thesis/combined_plots.png", combined_plots, width = 12, height = 9, dpi = 400)

# Variant count df
var_count_df <- plots_vaf10$sort_table %>% 
  mutate(percent = paste0(round(percent, 2), "%")) %>% 
  rename("Variant count" = n_total,
         n_indv = count,
         Percent = percent,
         Type = type)
write_csv(var_count_df, file = "thesis/var_count.csv")

# Output text file documenting mean of variant counts
mean_het <- mean(join_hg$n_total_het)
median_het <- median(join_hg$n_indel_het)

mean_hom <- mean(join_hg$n_total_hom)
median_hom <- median(join_hg$n_total_hom)

# Prepare the output text
output_text <- paste(
  "Documentation of Mean Variant Counts",
  "-------------------------------------",
  paste("Mean number of heteroplasmic variants:", mean_het),
  paste("Median number of heteroplasmic variants:", median_het),
  paste("Mean number of homoplasmic variants:", mean_hom),
  paste("Median number of homoplasmic variants:", median_hom),
  sep = "\n"
)

# Define the output file path
output_file <- "thesis/mean_variant_counts.txt"

# Write the output text to the file
writeLines(output_text, output_file)

