library(tidyverse)
library(patchwork)
library(gggenes)

setwd("/Users/xuanchou/Documents/TWB_1492/microarray/assoc/")
outdir <- "kidney_related/"
filedir <- "/Users/xuanchou/Documents/TWB_1492/microarray/assoc/twb2_maf.01_info0.7_refractor/"

linear_df <- read.csv(paste0(filedir, "twb2_maf.01_info0.7.rmOut.rint_linear_summary.csv"))
linear_sig <- 0.05 / length(unique(linear_df$snp))
logistic_df <- read.csv(paste0(filedir, "twb2_maf.01_info0.7.case1000_logistic_summary.csv"))
logistic_sig <- 0.05 / length(unique(logistic_df$snp))

linear_kidney <- c("BUN", "CREATININE", "egfr_cr", 
                   "URIC_ACID", "MICROALB", "CREATININE_URINE")
logistic_kidney <- c("KIDNEY_STONE", "DIABETES", "DIABETES_Type2", "DIABETES_Gestational", "RENAL_FAILURE")

linear_kidney_df <- linear_df %>% 
  filter(phenotype %in% linear_kidney) %>% 
  mutate(phenotype = ifelse(phenotype == "egfr_cr", "EGFRcr", phenotype))
linear_kidney <- c("BUN", "CREATININE", "EGFRcr", "URIC_ACID", "MICROALB", "CREATININE_URINE")
logistic_kidney_df <- logistic_df %>% 
  filter(phenotype %in% logistic_kidney)
kidney_df <- rbind(linear_kidney_df, logistic_kidney_df)

# add corresponding significant values to a new column
kidney_df <- kidney_df %>% 
  mutate(sig = ifelse(phenotype %in% linear_kidney, linear_sig, logistic_sig),
         if_sig = ifelse(p < sig, TRUE, FALSE))

# plot all the kidney-related phenotypes in a single plot
pheno_order <- c(linear_kidney, logistic_kidney)
kidney_df <- kidney_df %>%
  separate(snp, into = c("CHROM", "POS", "EA", "Other Allele"), sep = "\\.|\\.") %>% 
  mutate(CHROM = paste0("chr", toupper(CHROM)))
kidney_df$phenotype <- factor(kidney_df$phenotype, levels = pheno_order)

# -log10(p-val)
kidney_df$log_p <- -log10(kidney_df$p)
y_lim <- ceiling(max(kidney_df$log_p))

create_kidney_plot <- function(data) {
  plt <- ggplot(data, aes(x = as.numeric(POS), y = log_p, color = phenotype)) +
    geom_point() +
    labs(x = "POS", y = "-log10(p)", color = "Phenotype") +
    theme_minimal() +
    facet_wrap(~phenotype, scales = "free_y", ncol = 1) +
    geom_hline(yintercept = -log10(linear_sig), color = "red", linetype = "dashed") +
    ylim(0, y_lim) +
    guides(color = "none") +
    coord_cartesian(xlim = c(1, 16569), expand = FALSE) +
    theme(axis.title.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank())
  
  return(plt)
}

kidney_plt <- create_kidney_plot(kidney_df)
kidney_plt

# combine with gene info =====
# load gene information
gene_info <- read.delim("/Users/xuanchou/Documents/TWB_1492/Analysis/seq_info/combined_positions.tsv", sep = "\t")

# add strand orientation
gene_info <- gene_info %>% 
  mutate(Strand = ifelse(Gene.Name == "ND6", FALSE, TRUE))

# plot a rectangle representing the gene position and name
arrow_label <- c("rRNA", "Protein coding")
block_label <- c("Control region", "tRNA", "Intergenic region")

gene_label <- ggplot(data = gene_info,
                     aes(xmin = Start.Position, xmax = End.Position, y = "", label = Gene.Name)) +
  geom_gene_arrow(data = filter(gene_info, Function %in% arrow_label), # those genes need arrow representing the direction
                  aes(forward = Strand, fill = Function),
                  arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  geom_gene_arrow(data = filter(gene_info, Function %in% block_label),
                  aes(forward = Strand, fill = Function),
                  arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm")) +
  geom_gene_label(data = filter(gene_info, Complement != "")) +
  scale_fill_brewer(palette = "Set2") +
  theme_genes() +
  coord_cartesian(xlim = c(1, 16569), expand = FALSE) +
  theme(axis.title.y = element_blank(),
        legend.title = element_blank()) +
  scale_x_discrete(limits = c(1, 3000, 6000, 9000, 12000, 15000, 16569))
gene_label

combine_gene <- 
  kidney_plt / gene_label + plot_layout(nrow = 2, heights = c(0.7, 0.07),
                                               guides = "collect")
combine_gene
ggsave(filename = paste0(outdir, "all_traits.pdf"), plot = combine_gene, width = 8.1, height = 10)

# remove phenotypes without significant signals =====
sig_kidney_pheno <- c("CREATININE", "EGFRcr")

kidney_sig_df <- kidney_df %>% 
  filter(phenotype %in% sig_kidney_pheno)
kidney_sig_df$POS <- as.numeric(kidney_sig_df$POS)

# Map Gene.Name and Function based on position
kidney_sig_df <- within(kidney_sig_df, {
  Gene.Name <- sapply(POS, function(x) {
    gene_subset <- gene_info[gene_info$Start.Position <= x & gene_info$End.Position >= x, ]
    if (nrow(gene_subset) == 0) return(NA)
    gene_subset$Gene.Name[1]
  })
  Function <- sapply(POS, function(x) {
    func_subset <- gene_info[gene_info$Start.Position <= x & gene_info$End.Position >= x, ]
    if (nrow(func_subset) == 0) return(NA)
    func_subset$Function[1]
  })
})

# Define consistent color mapping
function_order <- c("Control region", "Intergenic region", "Protein coding", "rRNA", "tRNA")
color_palette <- scale_color_brewer(palette = "Set2", name = "Function", labels = function_order)

# Make manhattan plot
# Determine the limits for the y-axis, adding a buffer
y_min <- min(kidney_sig_df$log_p, na.rm = TRUE) - 0.1
y_max <- max(kidney_sig_df$log_p, na.rm = TRUE) + 0.1

# Create the plot
kidney_sig_plt <- ggplot(kidney_sig_df, aes(x = as.numeric(POS), y = log_p, color = Function)) +
  geom_point(alpha = 0.7) +
  labs(x = "POS", y = "-log10(p)", color = "Phenotype") +
  theme_bw() +
  geom_hline(yintercept = -log10(linear_sig), color = "red", linetype = "dashed") +
  guides(color = "none") +
  facet_wrap(~phenotype, scales = "free_y", ncol = 1) +
  coord_cartesian(xlim = c(1, 16569), ylim = c(y_min, y_max), expand = FALSE) +
  color_palette +
  theme(axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank())
kidney_sig_plt
kidney_sig_plt / gene_label + plot_layout(nrow = 2, heights = c(0.7, 0.07),
                                      guides = "collect")

ggsave(filename = paste0(outdir, "sig_traits_poster.pdf"), width = 8, height = 6)

merged_df <- kidney_sig_df %>%
  filter(if_sig == TRUE) %>% 
  mutate(POS = as.numeric(POS)) %>% 
  arrange(POS) %>% 
  select(-c("OR", "OR_CI_lower", "OR_CI_upper", "type", "n_cases", "n_controls", "formula", "note", "sig", "if_sig", "log_p"))
  

# Add gene information
merge_df_with_gene <- merged_df %>%
  mutate(
    gene_name = sapply(POS, function(pos) {
      gene <- gene_info$Gene.Name[which(pos >= gene_info$Start.Position & pos <= gene_info$End.Position)]
      if (length(gene) == 0) {
        gene <- NA
      }
      return(gene)
    })
  )

# Add AF information
var_af <- read.table("/Users/xuanchou/Documents/TWB_1492/microarray/replicate/twb2_assoc_var_summary.tsv", header = T)
final_merge <- merge(merge_df_with_gene, var_af, 
                     by.x = c("POS", "EA"),
                     by.y = c("Position", "EA"),
                     all.x = TRUE) %>% 
  arrange(POS) %>% 
  select(-c("Other Allele"))

# Rearrange the column order
first_columns <- c("ID", "CHROM", "POS", "EA", "ALT", "REF")
remaining_columns <- setdiff(colnames(final_merge), first_columns) # Get the names of the remaining columns

new_column_order <- c(first_columns, remaining_columns)
final_merge <- final_merge[, new_column_order]

write.table(final_merge, file = paste0(outdir, "sig_w_gene.tsv"), sep = "\t", 
            row.names = FALSE, quote = FALSE)  
