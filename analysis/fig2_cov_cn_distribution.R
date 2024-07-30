library(tidyverse)
library(gggenes)
library(patchwork)

setwd("/Users/xuanchou/Documents/TWB_1492/Analysis/cov_cn/")
sample_df <- read.delim("/Users/xuanchou/Documents/TWB_1492/add_annotation/output_age/sample_annotations.txt")
pos_cov <- read_tsv("/Users/xuanchou/Documents/TWB_1492/annotate_coverage/output.tsv")
pos_cov$locus <- as.numeric(gsub("chrM:", "", pos_cov$locus))
sites_df <- read.table("/Users/xuanchou/Documents/TWB_1492/add_annotation/output_age/combined_sites_only.vcf", comment.char = "#", sep = "\t", 
                       col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
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
  scale_fill_brewer(palette = "Set3") +
  theme_genes() +
  coord_cartesian(xlim = c(1, 16569), expand = FALSE) +
  theme(axis.title.y = element_blank(),
        legend.title = element_blank()) +
  scale_x_discrete(limits = c(1, 3000, 6000, 9000, 12000, 15000, 16569))
gene_label

# Collcet pos info by sites_df
pos_filter <- sites_df %>% 
  group_by(POS) %>%
  summarize(FILTERS = toString(unique(FILTER))) %>%
  mutate(FILTERS = ifelse(grepl("artifact_prone_site", FILTERS), "artifact_prone_site",
                          ifelse(grepl("indel_stack", FILTERS), "indel_stack", FILTERS)))
pos_cov <- merge(pos_cov, pos_filter, by.x = "locus", by.y = "POS", all.x = TRUE)
pos_cov <- pos_cov %>% 
  mutate(FILTERS = ifelse(as.numeric(locus) == 3107, "artifact_prone_site", FILTERS))

# Coverage at ever base
pos_cov_plt <- ggplot(pos_cov, aes(x = locus, y = mean)) +
  geom_line() +
  labs(x = "Position",
       y = "Mean coverage",
       subtitle = paste("Mean covergae:", round(mean(pos_cov$mean), 2))) +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  theme(axis.title = element_text(size = 15, face = "bold"))
pos_cov_plt
#ggsave(filename = "pos_cov.png", dpi = 300)

# add filter label
pos_cov_label_plt <- pos_cov_plt +
  geom_vline(data = filter(pos_cov, FILTERS == "artifact_prone_site"),
             aes(xintercept = locus, color = "Artifact Prone Site"), 
             linetype = "dashed", show.legend = TRUE) +
  geom_vline(data = filter(pos_cov, FILTERS == "indel_stack"),
             aes(xintercept = locus, color = "Indel Stack"), 
             linetype = "dashed", show.legend = TRUE) +
  scale_color_manual(values = c("#5499C7", "#F0B27A"), 
                     breaks = c("Artifact Prone Site", "Indel Stack"),
                     labels = c("Artifact-prone site", "Indel stack")) +
  labs(color = "Filter") +
  theme(axis.title.y = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank())

pos_cov_label_plt_combine <- 
  pos_cov_label_plt / gene_label + plot_layout(nrow = 2, heights = c(0.7, 0.07),
                                               guides = "collect")

pos_cov_label_plt_combine
ggsave("cov_gene_annot.png", dpi = 300, width = 9, height = 4.5)

# Add annotations for the marked loci
# Create a table of annotated loci with FILTERS information
output_table <- pos_cov %>% 
  filter(FILTERS %in% c("artifact_prone_site", "indel_stack"))

# Get gene name from the gene_info df
get_gene_name <- function(position, gene_info) {
  gene_row <- gene_info %>%
    filter(Start.Position <= position & End.Position >= position)
  
  if (nrow(gene_row) > 0) {
    return(gene_row$Gene.Name)
  } else {
    return(NA)
  }
}

output_table <- output_table %>% 
  mutate(Gene = ifelse(locus < 576 | locus > 16024, "DLOOP",
         sapply(locus, get_gene_name, gene_info = gene_info))) %>% 
  select(-c(over_100, over_1000))
colnames(output_table) <- c("Position", "Mean coverage", "Median coverage", "Filter", "Gene")
write_delim(output_table, "artifact_indel.csv", delim = ",")

### Add table to the plot ###
annotations_table <- pos_cov %>%
  filter(FILTERS %in% c("artifact_prone_site", "indel_stack")) %>%
  group_by(FILTERS) %>%
  summarise("Annotated Loci" = toString(unique(locus))) %>%
  mutate(Filter = case_when(
    FILTERS == "artifact_prone_site" ~ "Artifact Prone Site",
    FILTERS == "indel_stack" ~ "Indel Stack",
    TRUE ~ FILTERS
  )) %>%
  select(Filter, "Annotated Loci")

tbl_plt <- tableGrob(annotations_table, rows = NULL, theme = ttheme_minimal())
pos_cov_label_locus_plt <- pos_cov_label_plt_combine / tbl_plt + 
  plot_layout(heights = c(2, 0.1, 1), nrow = 3, guides = "keep")

pos_cov_label_locus_plt
ggsave(pos_cov_label_locus_plt, file = "pos_cov_locus_label.png", width = 9, height = 5, dpi = 200)

pos_cov_tbl <- pos_cov_label_plt / tbl_plt +
  plot_layout(heights = c(2, 1), nrow = 2)
ggsave(pos_cov_tbl, file = "pos_cov_tbl.png", width = 9, height = 5, dpi = 200)


# mean coverage and CN per individual and the frequency =====
# Distribution of WGS coverage
distribution_plot <- ggplot(sample_df, aes(x = wgs_mean_coverage)) +
  geom_histogram(binwidth = 1) +
  labs(x = "nDNA mean coverage",
       y = "Frequency") +
  theme_classic() +
  theme(axis.title = element_text(size = 15, face = "bold"))
distribution_plot
mean(sample_df$wgs_mean_coverage)
ggsave(distribution_plot, filename = "wgs_cov.png", width = 10.5, height = 6.3, dpi = 300)

# Distribution of MT coverage
mt_distribution_plot <- ggplot(sample_df, aes(x = mt_mean_coverage)) +
  geom_histogram(binwidth = 100) +
  labs(x = "mtDNA mean coverage",
       y = "Frequency") +
  theme_classic() +
  theme(axis.title = element_text(size = 15, face = "bold"))
mt_distribution_plot
mean(sample_df$mt_mean_coverage)
ggsave(mt_distribution_plot, filename = "mt_cov.png", width = 10.5, height = 6.3, dpi = 300)

# Distribution of MtCN
mtcn_distribution_plot <- ggplot(sample_df, aes(x = mito_cn)) +
  geom_histogram(binwidth = 10) +
  labs(x = "mtDNA copy number estimate",
       y = "Frequency",
       subtitle = paste0("Mean CN: ", round(mean(sample_df$mito_cn), digits = 2))) +
  theme_classic() +
  theme(axis.title = element_text(size = 15, face = "bold"))
mtcn_distribution_plot
mean(sample_df$mito_cn)
median(sample_df$mito_cn)
ggsave(mtcn_distribution_plot, filename = "mt_cn.png", width = 10.5, height = 6.3, dpi = 300)

