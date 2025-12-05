library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(GGally)
library(colorspace)

# para <- "gnomad_mt_pass"
para <- "vqsr_geno02_mind05_sing_snp"
file_dir <- "/Users/xuanchou/Documents/TWB_1492/plink_pca/plink_file/"
out_dir <- "/Users/xuanchou/Documents/TWB_1492/plink_pca/"

if (grepl("mt", para)) {
  pc_label <- "mtPC"
  eigenval_df <- "twb2_imp.geno05.mind05.maf01.eigenval"
  eigenvec_df <- "twb2_imp.geno05.mind05.maf01.eigenvec"
} else {
  pc_label <- "nucPC"
  eigenval_df <- "twb2_EAS_pca_qc.eigenval"
  eigenvec_df <- "twb2_EAS_pca_qc.eigenvec"
}

# Check if the folder named para exists
if (!dir.exists(paste0(out_dir, para))) {
  # If not, create a new folder
  dir.create(paste0(out_dir, para))
}
setwd(paste0(out_dir, para))

# Load files
eigneval <- read.delim(paste0(file_dir, "TWB1492_PCA_", para, ".eigenval"), header = FALSE, nrows = 10)
eigenvec <- read.table(paste0(file_dir, "TWB1492_PCA_", para, ".eigenvec"), header = FALSE) %>% select(1:12)
eigenvec <- eigenvec %>% select(-1)
colnames(eigenvec) <- c("ID", paste0("PC", 1:10))
## sample annotation
hap_df <- read.delim("/Users/xuanchou/Documents/TWB_1492/add_annotation/output_age/sample_annotations.txt")
hap_df <- hap_df %>% select(s, hap) %>% rename(ID = s)
## ancestry info
pop_df <- read.csv("/Users/xuanchou/Documents/Github/twb_survey_cleaning/twb_ancestry/processed_origin.csv")
pop_df <- pop_df %>% 
  rename(ID = Sample_Name) %>%
  mutate(across(everything(), ~ str_replace_all(., "China", "Chinese")))

# Make PVE plot
total_variance <- sum(eigneval$V1)
eigneval <- eigneval %>% mutate(Proportion = V1 / total_variance * 100)

ggplot(eigneval, aes(x = factor(rev(Proportion)), y = Proportion)) +
  geom_col(fill = "#5499C7") +
  geom_text(aes(label = sprintf("%.1f%%", Proportion), vjust = -0.5), size = 4) +
  labs(title = "Proportion of Variance Explained",
       x = "Principal Components",
       y = "Proportion of Variance Explained (%)") +
  scale_x_discrete(labels = paste0("PC", 1:10)) +
  theme_classic()
ggsave(filename = paste0(para, "_pve.jpg"), dpi=400)

# Make PCA
## Merge eigenvec with hap_df
eigenvec$ID <- gsub("SM_", "", eigenvec$ID)
merge_df <- inner_join(hap_df, eigenvec, by = "ID") %>%
  inner_join(pop_df, by = "ID") %>% 
  mutate(hg_branch = case_when(
    hap %in% c("M", "C", "Z", "E", "G", "D") ~ "Branch M",
    hap %in% c("N", "Y", "A") ~ "Branch N",
    hap %in% c("F", "B", "R") ~ "Branch R",
    TRUE ~ hap
  ))

# Remove outliers (> 6SD) for PC1 to PC6
for (i in 1:6) {
  pc_col <- paste0("PC", i)
  mean_pc <- mean(merge_df[[pc_col]], na.rm = TRUE)
  sd_pc <- sd(merge_df[[pc_col]], na.rm = TRUE)
  merge_df <- merge_df %>%
    filter((!!sym(pc_col)) < mean_pc + 8 * sd_pc & (!!sym(pc_col)) > mean_pc - 8 * sd_pc)
}

# Add PVE information to column names
pve_labels <- eigneval$Proportion
new_colnames <- c(sprintf("%s1 (%.1f%%)", pc_label, pve_labels[1]), 
                  sprintf("%s2 (%.1f%%)", pc_label, pve_labels[2]),
                  sprintf("%s3 (%.1f%%)", pc_label, pve_labels[3]), 
                  sprintf("%s4 (%.1f%%)", pc_label, pve_labels[4]),
                  sprintf("%s5 (%.1f%%)", pc_label, pve_labels[5]), 
                  sprintf("%s6 (%.1f%%)", pc_label, pve_labels[6]),
                  sprintf("%s7 (%.1f%%)", pc_label, pve_labels[7]), 
                  sprintf("%s8 (%.1f%%)", pc_label, pve_labels[8]),
                  sprintf("%s9 (%.1f%%)", pc_label, pve_labels[9]), 
                  sprintf("%s10 (%.1f%%)", pc_label, pve_labels[10]))
pc1_col <- which(colnames(merge_df) == "PC1")
pc10_col <- which(colnames(merge_df) == "PC10")
colnames(merge_df)[pc1_col:pc10_col] <- new_colnames

# Ancestry
# Define base colors for each prefix
base_colors <- c(
  "Holo" = "#D6C078",
  "Hakka" = "#99B581",
  "Southern Chinese" = "#4C87B8",
  "Northern Chinese" = "#9C5CBF",
  "Other Chinese" = "#C46A9C",
  "Other" = "#797D7F"
)

# Generate gradient colors for each category
generate_gradient <- function(base_color) {
  colorRampPalette(c(base_color, lighten(base_color, 0.3)))(5)
}

ancestry_colors <- c(
  "Holo" = base_colors[["Holo"]],
  "Holo/Hakka" = generate_gradient(base_colors["Holo"])[2],
  "Holo/Southern Chinese" = generate_gradient(base_colors["Holo"])[3],
  "Holo/Northern Chinese" = generate_gradient(base_colors["Holo"])[4],
  "Holo/Other Chinese" = generate_gradient(base_colors["Holo"])[5],
  "Holo/Other" = generate_gradient(base_colors["Holo"])[5],
  
  "Hakka" = base_colors[["Hakka"]],
  "Hakka/Southern Chinese" = generate_gradient(base_colors["Hakka"])[2],
  "Hakka/Northern Chinese" = generate_gradient(base_colors["Hakka"])[3],
  "Hakka/Other Chinese" = generate_gradient(base_colors["Hakka"])[4],
  "Hakka/Other" = generate_gradient(base_colors["Hakka"])[5],
  
  "Southern Chinese" = base_colors[["Southern Chinese"]],
  "Southern Chinese/Northern Chinese" = generate_gradient(base_colors["Southern Chinese"])[2],
  "Southern Chinese/Other Chinese" = generate_gradient(base_colors["Southern Chinese"])[3],
  "Southern Chinese/Other" = generate_gradient(base_colors["Southern Chinese"])[4],
  
  "Northern Chinese" = base_colors[["Northern Chinese"]],
  "Northern Chinese/Other Chinese" = generate_gradient(base_colors["Northern Chinese"])[2],
  "Northern Chinese/Other" = generate_gradient(base_colors["Northern Chinese"])[3],
  
  "Other Chinese" = base_colors[["Other Chinese"]],
  
  "Other" = base_colors[["Other"]]
)

# Define the custom order
ancestry_custom_order <- c(
  "Holo", "Holo/Hakka", "Holo/Southern Chinese", "Holo/Northern Chinese", "Holo/Other Chinese", "Holo/Other",
  "Hakka", "Hakka/Southern Chinese", "Hakka/Northern Chinese", "Hakka/Other Chinese", "Hakka/Other",
  "Southern Chinese", "Southern Chinese/Northern Chinese", "Southern Chinese/Other Chinese", "Southern Chinese/Other",
  "Northern Chinese", "Northern Chinese/Other Chinese", "Northern Chinese/Other",
  "Other Chinese", "Other"
)

# Function to create and save plots
create_and_save_plot <- function(df, column, filename, legend_label, custom_order, ancestry_colors, para) {
  df[[column]] <- factor(df[[column]], levels = custom_order)
  plt <- ggplot(df, aes(x = !!sym(new_colnames[1]), y = !!sym(new_colnames[2]))) +
    geom_point(aes(color = !!sym(column)), size = 3, alpha = 0.7) +
    theme_bw() +
    scale_color_manual(values = ancestry_colors) +
    theme(panel.grid = element_blank()) +
    labs(color = legend_label)
  
  ggsave(plt, filename = paste0(para, "_", filename, ".jpg"), dpi=400, width = 10.8, height = 6.85)
  saveRDS(plt, file = paste0(para, "_", filename, ".rds"))
}

# Apply the function to different columns
# Platform
platform_colors <- c("#1A759F", "#52B788", "#BE95C4")
create_and_save_plot(merge_df, "Platform", "platform", "Platform", c("HiSeq 2500", "HiSeq 4000", "NovaSeq 6000"), platform_colors, para)

# Haplogroup
hg_order <- c("M", "C", "Z", "E", "G", "D", 
              "N", "Y", "A",
              "F", "B", "R",
              "H", "J", "K")
hg_colors <- c(c("#184E77", "#1E6091", "#1A759F", "#3479A3", "#168AAD", "#6BA1B3"),
               "#B7E4C7", "#95D5B2", "#52B788",
               "#BE95C4", "#B392AC", "#D1B3C4",
               "#797D7F", "#717D7E", "#7E846B")
create_and_save_plot(merge_df, "hap", "hap", "Haplogroup", hg_order, hg_colors, para)

# Haplogorup branch
hg_branch_order <- c("Branch M", "Branch N", "Branch R", "H", "J", "K")
hg_branch_colors <- c("#1A759F", "#52B788", "#BE95C4", "#797D7F", "#717D7E", "#7E846B")
create_and_save_plot(merge_df, "hg_branch", "hg_branch", "Haplogroup branch", hg_branch_order, hg_branch_colors, para)

create_and_save_plot(merge_df, "NATIVE_MOM", "native_mom", "Maternal ancestry", ancestry_custom_order, ancestry_colors, para)
create_and_save_plot(merge_df, "NATIVE_FA", "native_fa", "Paternal ancestry", ancestry_custom_order, ancestry_colors, para)
create_and_save_plot(merge_df, "NATIVE_COMBINE", "native_combined", "Parental ancestry", ancestry_custom_order, ancestry_colors, para)

# Pairplot from PC1 to PC10

pair_plt <- ggpairs(
  merge_df, 
  columns = pc1_col:pc10_col, 
  aes(color = NATIVE_MOM),
  upper = list(continuous = "blank")  # This line will blank out the upper right correlations
) + 
  scale_color_manual(values = ancestry_colors)
ggsave(filename = paste0(para, "_native_mom_pair.jpg"), plot = pair_plt, width = 10, height = 10)
