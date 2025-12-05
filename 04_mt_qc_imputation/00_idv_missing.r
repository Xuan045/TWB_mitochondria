#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
filedir <- args[1]
filename <- args[2]

library(tidyverse)

setwd(filedir)
missing_df <- read.table(filename, header = TRUE, comment.char = "")

xmax <- max(0.06, max(missing_df$F_MISS))
p <- ggplot(missing_df, aes(x = F_MISS)) +
  geom_histogram(alpha = 0.6, binwidth = 0.001, color = "black") +
  geom_vline(xintercept = 0.05) +
  geom_vline(xintercept = 0.02) +
  xlim(0, xmax) +
  theme_bw() +
  labs(x = "Per-individual missing rate", y = "Count")

# add text indicating count of missing rate > 0.02 and 0.05
max_count <- max(ggplot_build(p)$data[[1]]$count)
text_y <- max_count / 3 * 2
p <- p +
  annotate("text", x = 0.055, y = text_y, label = paste("F_MISS > 0.05:", sum(missing_df$F_MISS > 0.05)), color = "#F07167", fontface = "bold") +
  annotate("text", x = 0.025, y = text_y, label = paste("F_MISS > 0.02:", sum(missing_df$F_MISS > 0.02)), color = "#0081A7", fontface = "bold")

ggsave(paste0(filename, ".pdf"), p, width=7, height=4)