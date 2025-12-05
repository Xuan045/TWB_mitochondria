library(ggplot2)
library(dplyr)
library(ggrepel)

# Function to create a cumulative data frame for plotting
prepare_cumulative_df <- function(df, phenotype_order) {
  df %>% 
    group_by(phenotype) %>% 
    summarise(num_pheno = n()) %>%
    mutate(phenotype = factor(phenotype, levels = phenotype_order)) %>% 
    arrange(phenotype) %>% 
    mutate(pheno_add = lag(cumsum(num_pheno), default = 0)) %>% 
    select(phenotype, pheno_add) %>% 
    mutate(next_pheno_add = lead(pheno_add, default = max(pheno_add) + length(unique(df$snp))),
           x_label = (pheno_add + next_pheno_add) / 2)
}

# Plot features
plot_feature <- theme(legend.position = "none",
                      panel.grid.minor.x = element_blank(),
                      axis.ticks=element_blank(),
                      axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
                      axis.text.y = element_text(size = 15),
                      axis.title = element_text(size = 18, face = "bold"))

# PheWAS plotting function
plot_phewas <- function(plot_df, df_cum, select_cols, gene_cols, outdir) {
  sig <- 0.05 / length(geno_cols)
  strict_sig <- sig / length(select_cols)
  
  annotation_data <- data.frame(
    x = c(max(plot_df$cum_add), max(plot_df$cum_add)),
    y = c(-log10(sig) - 0.1, -log10(strict_sig) - 0.1),
    label = c(sprintf("Significance level: %.1e", sig), sprintf("Significance level: %.1e", strict_sig))
  )
  
  # Y-axis limit
  if (strict_sig < min(plot_df$p)) {
    y_max <- strict_sig
  } else {
    y_max <- min(plot_df$p)
  }
  
  p <- ggplot(plot_df, aes(x = cum_add, y = -log10(p))) + 
    geom_point(aes(color = phenotype, size = ifelse(is.na(OR), 2.5, OR))) + 
    theme_classic() + 
    scale_x_continuous(labels = select_cols, breaks = df_cum$x_label) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, -log10(y_max)) + 0.5) +
    plot_feature +
    labs(color = "Phenotype", x = "Phenotypes", y = "-log(p-value)") +
    geom_text_repel(data=. %>% mutate(label = ifelse(p < strict_sig , snp, "")),
                    aes(label=label), size = 5, box.padding = unit(0.7, "lines")) +
    geom_text_repel(data = annotation_data, 
                    aes(x = x, y = y, label = label), 
                    size = 6, color = "black") +
    geom_hline(yintercept = -log10(sig), color = "red", size = 1, alpha = 0.5) +
    geom_hline(yintercept = -log10(strict_sig), color = "red", size = 2, alpha = 0.7)
  
  if (all(!is.na(plot_df$OR))) {
    p <- p + labs(size = "OR")
  } else {
    p <- p + theme(legend.position = "none")
  }

  # Save plot
  png(outdir, width = 20, height = 12, units = "in", res = 300, type = "cairo")
  print(p)
  # ggsave(filename = paste0(outdir), plot = p, dpi = 300, width = 20, height = 12)
  return(p)
}

# Function to plot PheWAS for linear regression results
plot_phewas_linear <- function(final_df, select_cols, geno_cols, outdir) {
  df_cum <- prepare_cumulative_df(final_df, select_cols)
  
  # The points are sorted by POS
  final_df$POS <- as.numeric(gsub("[^0-9]", "", final_df$snp))
  plot_df <- final_df %>%
    group_by(phenotype) %>% 
    inner_join(df_cum) %>%
    mutate(cum_add = row_number() + pheno_add) %>%
    arrange(POS)
  plot_df$phenotype <- factor(plot_df$phenotype, levels = select_cols)
  ori_plot <- plot_phewas(plot_df, df_cum, select_cols, gene_cols, paste0(outdir, ".png"))
  
  # Sorted by p value
  plot_df_sorted <- final_df %>% 
    group_by(phenotype) %>% 
    arrange(desc(p)) %>% 
    inner_join(df_cum) %>% 
    mutate(cum_add = row_number() + pheno_add)
  plot_df_sorted$phenotype <- factor(plot_df_sorted$phenotype, levels = select_cols)
  sorted_plot <- plot_phewas(plot_df_sorted, df_cum, select_cols, gene_cols, paste0(outdir, "_sorted.png"))
  sorted_loose_plot <- sorted_plot + geom_text_repel(data=. %>% mutate(label = ifelse((p < (0.05 / length(geno_cols))) & (p > (0.05 / length(geno_cols) / length(select_cols))) , snp, "")),
                                                     aes(label=label), size = 5, box.padding = unit(0.7, "lines"))

  # Save plot
  png(paste0(outdir, "_sorted_loose.png"), width = 20, height = 12, units = "in", res = 300, type = "cairo")
  print(sorted_loose_plot)
  # ggsave(filename = paste0(outdir, "_sorted_loose.png"), plot = sorted_loose_plot, dpi = 300, width = 20, height = 12)
}


# Function to plot PheWAS for logistic regression results
plot_phewas_logistic <- function(final_df, select_cols, geno_cols, outdir) {
  # Remove rows which OR is NA
  plot_df <- final_df %>% filter(!is.na(OR) & Beta < 2)
  select_cols <- select_cols[select_cols %in% unique(plot_df$phenotype)]
  
  df_cum <- prepare_cumulative_df(plot_df, select_cols)
  
  # The points are sorted by POS
  plot_df$POS <- as.numeric(gsub("[^0-9]", "", plot_df$snp))
  plot_df <- plot_df %>%
    group_by(phenotype) %>% 
    inner_join(df_cum) %>%
    mutate(cum_add = row_number() + pheno_add) %>%
    arrange(POS)
  plot_df$phenotype <- factor(plot_df$phenotype, levels = select_cols)
  ori_plot <- plot_phewas(plot_df, df_cum, select_cols, gene_cols, paste0(outdir, ".png"))

  # Sort by p value
  plot_df_sorted <- plot_df %>%
    group_by(phenotype) %>%
    arrange(desc(p)) %>%
    inner_join(df_cum) %>%
    mutate(cum_add = row_number() + pheno_add)
  plot_df_sorted$phenotype <- factor(plot_df_sorted$phenotype, levels = select_cols)
  sorted_plot <- plot_phewas(plot_df_sorted, df_cum, select_cols, gene_cols, paste0(outdir, "_sorted.png"))

  sorted_loose_plot <- sorted_plot + geom_text_repel(data=. %>% mutate(label = ifelse((p < (0.05 / length(geno_cols))) & (p > (0.05 / length(geno_cols) / length(select_cols))) , snp, "")),
                                                     aes(label=label), size = 5, box.padding = unit(0.7, "lines"))

  # Save plot
  png(paste0(outdir, "_sorted_loose.png"), width = 20, height = 12, units = "in", res = 300, type = "cairo")
  print(sorted_loose_plot)
  # ggsave(filename = paste0(outdir, "_sorted_loose.png"), plot = sorted_loose_plot, dpi = 300, width = 20, height = 12)
}
