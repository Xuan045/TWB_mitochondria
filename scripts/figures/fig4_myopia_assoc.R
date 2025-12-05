library(tidyverse)
library(gggenes)
library(patchwork)

setwd("/Users/xuanchou/Documents/TWB_1492/microarray/assoc/myopia/")

filedir <- "/Users/xuanchou/Documents/TWB_1492/microarray/assoc/twb2_maf.01_info0.7_refractor/"
logistic_df <- read.csv(paste0(filedir, "twb2_maf.01_info0.7.case1000_logistic_summary.csv"))
logistic_sig <- 0.05 / length(unique(logistic_df$snp))

# Gene label beneath the plot
gene_info <- read.delim("/Users/xuanchou/Documents/TWB_1492/Analysis/seq_info/combined_positions.tsv", sep = "\t")

# add strand orientation
gene_info <- gene_info %>% 
  mutate(Strand = ifelse(Gene.Name == "ND6", FALSE, TRUE))

# plot a rectangle representing the gene position and name
arrow_label <- c("rRNA", "Protein coding")
block_label <- c("Control region", "tRNA", "Intergenic region")

# Order the function
function_order <- c("Control region", "Intergenic region",
                    "Protein coding", "rRNA", "tRNA")
gene_info$Function <- factor(gene_info$Function, levels = function_order)

# Define consistent color mapping
color_palette <- scale_color_brewer(palette = "Set2", name = "Function", labels = function_order)

# Make gene label plot
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

# Filter to phenotype == "MYOPIA_1000"
# split the snp column
myopia_df <- logistic_df %>%
  filter(phenotype == "MYOPIA_1000") %>% 
  mutate(
    chr = "MT",  # Assuming mitochondrial data "m"
    pos = as.numeric(gsub("m\\.([0-9]+).*", "\\1", snp)),
    effect = gsub("m\\.\\d+\\.([A-Z]+)\\..*", "\\1", snp),
    other = gsub("m\\.\\d+\\.[A-Z]+\\.([A-Z]+)", "\\1", snp)
  )

# Map Gene.Name and Function based on position
myopia_df <- within(myopia_df, {
  Gene.Name <- sapply(pos, function(x) {
    gene_subset <- gene_info[gene_info$Start.Position <= x & gene_info$End.Position >= x, ]
    if (nrow(gene_subset) == 0) return(NA)
    gene_subset$Gene.Name[1]
  })
  Function <- sapply(pos, function(x) {
    func_subset <- gene_info[gene_info$Start.Position <= x & gene_info$End.Position >= x, ]
    if (nrow(func_subset) == 0) return(NA)
    func_subset$Function[1]
  })
})

# Make manhattan plot
# First, ensure there are no zero or NA values in the p-values
myopia_df <- myopia_df %>% 
  filter(p > 0) %>%  # Remove rows where p-value is zero or negative
  mutate(log_p = -log10(p))  # Calculate the log10 transformation

# Determine the limits for the y-axis, adding a buffer
y_min <- min(myopia_df$log_p, na.rm = TRUE) - 0.1
y_max <- max(myopia_df$log_p, na.rm = TRUE) + 0.1

# Create the plot
myopia_plt <- ggplot(myopia_df, aes(x = as.numeric(pos), y = log_p, color = Function)) +
  geom_point(alpha = 0.7) +
  labs(x = "POS", y = "-log10(p)", color = "Phenotype") +
  theme_bw() +
  geom_hline(yintercept = -log10(logistic_sig), color = "red", linetype = "dashed") +
  guides(color = "none") +
  coord_cartesian(xlim = c(1, 16569), ylim = c(y_min, y_max), expand = FALSE) +
  color_palette +
  theme(axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank())

combine_gene <- 
  myopia_plt / gene_label + plot_layout(nrow = 2, heights = c(0.7, 0.07),
                                        guides = "collect")
combine_gene
ggsave(filename = "myopia_plt.png", width = 8.1, height = 4, dpi = 300)
ggsave(filename = "myopia_plt.pdf", width = 8.1, height = 4)
ggsave(filename = "myopia_plt_poster.pdf", width = 8.1, height = 3.5)

# Sort the table
var_af <- read.table("/Users/xuanchou/Documents/TWB_1492/microarray/replicate/twb2_assoc_var_summary.tsv", header = T)
sort_df <- merge(myopia_df, var_af, by.x = c("pos"), by.y = c("Position"), all.x = TRUE)

column_order <- c("ID", "pos", "EA", "ALT", "REF", "phenotype", "Beta", "OR", "SE", 
                  "p", "Gene.Name", "status", "info", "EAF", "AF_array", "AF_wgs")
final_df <- sort_df[sort_df$p < logistic_sig, column_order] %>%
  rename(Position = pos,
         INFO = info,
         Status = status) %>%
  mutate(p = format(p, scientific = TRUE, digit = 4)) %>%
  mutate(across(.cols = c("Beta", "OR", "SE"), .fns = ~round(., 3))) %>%
  mutate(`Beta (SE)` = paste0(Beta, " (", SE, ")")) %>%
  select(-Beta, -SE) %>%
  select(ID, Position, EA, ALT, REF, phenotype, `Beta (SE)`, everything())
write_delim(final_df, file = "sig_w_gene.tsv", delim = "\t")
