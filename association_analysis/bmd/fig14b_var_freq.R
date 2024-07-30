library(tidyverse)
library(ggpubr)

setwd("/Users/xuanchou/Documents/TWB_1492/microarray/assoc/bone/")

pheno_df <- read_delim("/Users/xuanchou/Documents/TWB_1492/microarray/assoc/lab_test_cleaning/lab_test_rmOut.txt")
var_df <- read_delim("/Users/xuanchou/Documents/TWB_1492/microarray/replicate/wgs1465.twb2.imputed.filtered.info0.7.prob0")

var_df <- var_df %>% filter(RS_ID %in% c("rs41419246"))

# Transpose the df
var_t_df <- var_df %>% 
  select(RS_ID, starts_with("TV")) %>%
  pivot_longer(cols = -RS_ID, names_to = "ID", values_to = "value") %>%
  pivot_wider(names_from = RS_ID, values_from = value)

# Merge with phenotype
merge_df <- pheno_df %>% 
  select(ID, TWB2_Batch, AGE, SEX, T_SCORE, Z_SCORE) %>% 
  merge(var_t_df, by = "ID") %>% 
  filter(!(is.na(T_SCORE) | is.na(Z_SCORE)))

# Difference of two values in two genotype groups
val_diff <- merge_df %>%
  group_by(rs41419246, SEX) %>%
  summarize(
    mean_T_SCORE = mean(T_SCORE, na.rm = TRUE),
    mean_Z_SCORE = mean(Z_SCORE, na.rm = TRUE),
    mean_age = mean(AGE),
    .groups = "drop"
  ) %>% 
  ungroup() 

val_diff_plot <- val_diff %>% 
  pivot_longer(., c("mean_T_SCORE", "mean_Z_SCORE"),
               values_to = "BMD score", names_to = "BMD")
  
# Make bar plot
bmd_score_plt <- ggplot(val_diff_plot, aes(x = rs41419246, y = `BMD score`, fill = SEX)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "rs41419246",
       y = "Mean BMD Score") +
  scale_fill_discrete(
    name = "Sex",
    labels = c("Female", "Male")
  ) +
  facet_wrap(~BMD, labeller = labeller(BMD = c(
    "mean_T_SCORE" = "T-score",
    "mean_Z_SCORE" = "Z-score"
  ))) +
  theme_bw() +
  theme(legend.position = "bottom")
bmd_score_plt
ggsave("bmd_bar.png", plot = bmd_score_plt, dpi = 300, width = 8, height = 5)

# Output dataframe
write_delim(val_diff, "bmd_score.tsv", delim = "\t")

