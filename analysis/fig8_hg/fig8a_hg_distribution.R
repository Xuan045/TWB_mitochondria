library(scales)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(reshape2)

setwd("/Users/xuanchou/Documents/TWB_1492/Analysis/haplogroup/reanalyze_with_pass/")
sample_annot <- read.delim("/Users/xuanchou/Documents/TWB_1492/add_annotation/output_age/sample_annotations.txt")
# Hg classification using high-qual variants
pass_hg <- read.delim("/Users/xuanchou/Documents/TWB_1492/Analysis/haplogroup/hg_check/TWB1492_hg_allPASS_extend.txt") %>% 
  filter(Rank == 1) %>% 
  rename(major_haplogroup = Haplogroup) %>% 
  mutate(hap = substr(major_haplogroup, start = 1, stop = 1)) %>% 
  select(-Rank)
# Remove haplogroup info in the original df, and then retrieve more reliable version from the pass_hg
sample_annot <- sample_annot %>% 
  select(-c(hap, major_haplogroup)) %>% 
  merge(pass_hg, by.x = "s", by.y = "SampleID")

sample_full_df <- read.csv("/Users/xuanchou/Documents/TWB_1492/Analysis/TWB1492_full_table.csv")
hg_file <- sample_annot %>% group_by(hap) %>% summarise(n = n())
n_twb <- sum(hg_file$n)

# Define which colors to use for the high-level haplogroups
african_color <- "#955F90"
asian_color <- "#F4A259"
european_color <- "#F4E285"

# Group haplogroups into high level categories according to mitomap paper
African <- c("L0", "L1", "L5", "L2", "L6", "L4", "L3")
Asian <- c("M", "C", "Z", "E", "G", "Q", "D", "N", "Y", "A", "O", "S", "F", "B", "P")
European <- c("I", "W", "X", "R", "HV", "V", "H", "J", "T", "U", "K")

##### Haplogroups in TWB #####
# Sort and check if all Asian haplogroups in TW population, if not, add a new row with n=0
hg_sorted = tibble(hg_file[0,]) %>% add_column("place")
for(hg in Asian) {
  if(hg %in% hg_file$hap){
    hg_sorted[nrow(hg_sorted) + 1, ] = list(hg, as.integer(hg_file[hg_file$hap==hg,2]), "Asian")
  }else{hg_sorted[nrow(hg_sorted) + 1, ] = list(hg, 0, "Asian")}
}

# Sort and check if European hg in TW population
for (hg in European) {
  if(hg %in% hg_file$hap){
    hg_sorted[nrow(hg_sorted) + 1, ] = list(hg, as.integer(hg_file[hg_file$hap==hg,2]), "European")
  }
}
hg_sorted = hg_sorted %>% rename(place = `"place"`)

# Percentage of haplogroups
hg_sorted <- hg_sorted %>% mutate(percentage = percent(n/sum(n), accuracy = 0.01))
total <- sum(hg_sorted$n)

# Define colors
base_colors <- c("#184E77", "#1E6091", "#1A759F", "#3479A3", "#168AAD", "#6BA1B3",
                 "#B7E4C7", "#95D5B2", "#52B788",
                 "#BE95C4", "#B392AC", "#D1B3C4",
                 "#797D7F", "#717D7E", "#7E846B")
names(base_colors) <- c("M", "C", "Z", "E", "G", "D", 
                        "N", "Y", "A",
                        "F", "B", "R",
                        "H", "J", "K")

# Set a default color for unspecified levels
default_color <- "#D3D3D3"

# Generate a full color vector including a default color for unspecified levels
full_color_vector <- setNames(rep(default_color, length(levels(hg_sorted$hap))), levels(hg_sorted$hap))
full_color_vector[names(base_colors)] <- base_colors

# Plot
twb_hg_plot <- ggplot(hg_sorted, aes(x = hap, y = n, fill = hap)) +
  geom_bar(stat = "identity") +  # Use geom_bar for histograms with pre-summarized data
  scale_fill_manual(values = full_color_vector) +
  labs(x = "Haplogroups", y = "Counts",
       subtitle = paste("Total: ", sum(hg_sorted$n))) +
  geom_text(aes(label = n), vjust = -3, fontface = "bold", color = "black") +
  geom_text(aes(label = percentage), vjust = -1, size = 3, color = "black") +
  theme_bw() +
  ylim(0, max(hg_sorted$n) + 30) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
twb_hg_plot
ggsave("hg_color_order_submit.png", plot = twb_hg_plot, height = 5, width = 8.5, dpi = 400)
