library(tidyverse)
library(grid)
library(gridExtra)

setwd("/Users/xuanchou/Documents/TWB_1492/microarray/assoc/kidney_related/hg_assoc")
summary_df <- read_delim("kidney_hg_assoc.tsv", delim = "\t")
dir.create("plot/", showWarnings = FALSE, recursive = TRUE)

# Define haplogroup orders
African <- c("L0", "L1", "L5", "L2", "L6", "L4", "L3")
Asian <- c("M", "C", "Z", "E", "G", "Q", "D", "N", "Y", "A", "O", "S", "F", "B", "P")
European <- c("I", "W", "X", "R", "HV", "V", "H", "J", "T", "U", "K")

# Function to sort haplogroups
sort_haplogroups <- function(haplogroups) {
  order <- c(Asian, European, African)
  haplogroups <- unique(haplogroups)
  haplogroups <- haplogroups[order(match(substr(haplogroups, 1, 1), order))]
  return(haplogroups)
}

# Function to plot the results from the summary table
plot_summary_table <- function(summary_df, model_keys, pheno, angle = 45, significance_level = 0.05) {
  plots <- list()
  significant_color <- "#AE2012"
  non_significant_color <- "black"
  
  for (model_key in model_keys) {
    plot_data <- summary_df %>% 
      filter(model_type == model_key & phenotype == pheno)
    
    # Bonferroni's correction threshold
    bonferroni_threshold <- significance_level / nrow(plot_data)
    
    # Filter to haplogroup with count > 500
    plot_data <- plot_data[plot_data$count > 500, ]
    
    # Sort haplogroups
    sorted_haplogroups <- sort_haplogroups(plot_data$hap)
    
    plot_data$hap <- factor(plot_data$hap, levels = sorted_haplogroups)
    
    # Determine x-axis text colors based on significance
    x_text_colors <- sapply(levels(plot_data$hap), function(h) {
      ifelse(any(plot_data$hap == h & plot_data$`Pr(>|t|)` < bonferroni_threshold), significant_color, non_significant_color)
    })
    
    # Ensure the first label is black
    if (length(x_text_colors) > 0) {
      x_text_colors[1] <- non_significant_color
    }
    
    # Generate the plot using ggplot2 with rotated x-axis labels
    p <- ggplot(plot_data, aes(x = hap, y = Estimate)) +
      geom_bar(stat = "identity", aes(fill = `Pr(>|t|)` < bonferroni_threshold), color = "black") +
      geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), width = 0.2) +
      theme_bw() +
      scale_fill_manual(values = c("TRUE" = significant_color, "FALSE" = "grey"), guide = "none") +
      labs(title = paste("Association between", pheno, "and", model_key),
           x = "Haplogroup",
           y = "Effect Size") +
      theme(axis.text.x = element_text(angle = angle, hjust = 1, color = x_text_colors)) +
      geom_text(aes(y = max(`97.5 %`, na.rm = TRUE) + 0.1 * max(abs(`97.5 %`), na.rm = TRUE),
                    label = ifelse(`Pr(>|t|)` < bonferroni_threshold,
                                   paste0("p = ", formatC(`Pr(>|t|)`, format = "e", digits = 2)),
                                   ""),
                    color = `Pr(>|t|)` < bonferroni_threshold),
                vjust = 0.5, size = 5) +
      scale_color_manual(values = c("TRUE" = significant_color, "FALSE" = "black"), guide = "none")
    
    ggsave(plot = p, filename = paste0("plot/", pheno, "_", model_key, ".pdf"), width = 10, height = 5.5)
    plots[[model_key]] <- p
  }
  
  return(plots)
}

cre_model_keys <- c("Macrohaplogroup", "Haplogroup_digit2", "Haplogroup_digit3")
egfr_model_keys <- c("Macrohaplogroup", "Haplogroup_digit2", "Haplogroup_digit3")

# Run the plotting function
cre_plots <- plot_summary_table(summary_df, cre_model_keys, "CREATININE", angle = 45)
egfr_plots <- plot_summary_table(summary_df, egfr_model_keys, "eGFRcr", angle = 45)

# Save as PNG
ggsave(filename = "plot/cre_hg_sorted.png", plot = grid.arrange(grobs = cre_plots, ncol = 1), 
       dpi = 300, width = 14, height = 8)
ggsave(filename = "plot/egfr_hg_sorted.png", plot = grid.arrange(grobs = egfr_plots, ncol = 1), 
       dpi = 300, width = 14, height = 8)

# Save as PDF
pdf("plot/hg_assoc_sorted.pdf", width = 9, height = 7)
for (p in cre_plots) {
  grid.draw(p)
}
for (p in egfr_plots) {
  grid.draw(p)
}
dev.off()


