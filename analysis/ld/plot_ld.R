library(tidyverse)
library(reshape2)
library(RColorBrewer)

setwd("/Users/xuanchou/Documents/TWB_1492/LD")

file_dir <- "array_impute/"
para <- "info0.7_maf0.05"
ld_data <- read.table(paste0(file_dir, para, ".ld"))

ld_bim <- read.table(paste0(file_dir, para, ".bim"))
ld_bim <- ld_bim[order(ld_bim$V4), ] # Sort ld_bim by position

# Identify duplicated values in the 'V2' column
duplicated_rows <- duplicated(ld_bim$V2)

# Replace duplicated values in 'V2' with a new string
ld_bim$V2[duplicated_rows] <- paste(ld_bim$V2[duplicated_rows], ld_bim$V4[duplicated_rows], ld_bim$V5[duplicated_rows], ld_bim$V6[duplicated_rows], sep="_")

ld_id <- ld_bim$V2 # Extract sorted SNP IDs
colnames(ld_data) <- ld_id

# Create the ld_matrix
ld_matrix <- ld_data %>% 
  mutate(snp = ld_id) %>% 
  melt()
colnames(ld_matrix) <- c("snp1", "snp2", "r2")

# Reorder snp1 and snp2 based on ld_id
ld_matrix <- ld_matrix %>%
  mutate(snp1 = factor(snp1, levels = ld_id),
         snp2 = factor(snp2, levels = ld_id))

# Merge with ld_bim to get additional information
ld_pos <- unique(ld_bim$V4)
ld_matrix <- ld_matrix %>%
  left_join(ld_bim %>% select(V2, V4), by = c("snp1" = "V2")) %>%
  rename(snp1_info = V4) %>%
  left_join(ld_bim %>% select(V2, V4), by = c("snp2" = "V2")) %>%
  rename(snp2_info = V4) %>% 
  mutate(snp1_info = factor(snp1_info, levels = ld_pos),
         snp2_info = factor(snp2_info, levels = ld_pos))

# Make heatmap
pal_c <- colorRampPalette(brewer.pal(8, "Blues"))(25)
ggplot(ld_matrix, aes(x = snp1_info, y = snp2_info, fill = r2)) +
  geom_tile(color = "#E5E7E9") +
  scale_fill_gradientn(colors = pal_c) +
  labs(title = "Pairwise LD",
       y = "Position",
       fill = "r²") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.ticks.y = element_line()) +
  scale_y_discrete(breaks = ld_pos[seq(1, length(ld_pos), by = 3)])
ggsave(filename = paste0(file_dir, para, "_ld.png"), width = 8, height = 6, dpi = 200)

# Write LD matrix
write.csv(ld_matrix, file = paste0(file_dir, para, "_ld_matrix.csv"))

