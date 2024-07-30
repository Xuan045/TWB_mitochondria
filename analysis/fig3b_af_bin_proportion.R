library(ggplot2)
library(tidyverse)

setwd("/Users/xuanchou/Documents/TWB_1492/Analysis/vartype_af/")
# Read in the sites-only release data to obtain variant information
variant_data <- read.delim("/Users/xuanchou/Documents/TWB_1492/add_annotation/output_age/combined_sites_only.txt")
min_hom_threshold <- 0.95

# Filter to PASS-only variants
variant_data <- variant_data %>% filter(filters == "[]")

# Convert AFs to percentage
variant_data <- mutate(variant_data, ac = AC_het + AC_hom)
variant_data <- mutate(variant_data, af = ac/AN)
variant_data <- mutate(variant_data, af_percentage = af * 100)

# Group AFs into specific percentage categories to plot singletons, doubletons, and specific AFs
variant_data <- mutate(variant_data, 
                       variant_af_group = case_when(
                         ac == 1 ~ "singleton",
                         ac == 2 ~ "doubleton",
                         af_percentage < 1 ~ "<1%",
                         af_percentage >= 1 & af_percentage < 10 ~ "1-10%",
                         af_percentage >= 10 ~ paste0(floor(af_percentage/10)*10, " ~ ", floor(af_percentage/10)*10 + 10, "%"),
                         TRUE ~ "other"
                       ))

variant_data <- mutate(variant_data, count_variable = 1)
# Calculate proportions
proportions <- variant_data %>%
  group_by(variant_af_group) %>%
  summarize(prop = sum(count_variable) / nrow(variant_data))

# Sort the x-axis label
custom_order <- c("singleton", "doubleton", "<1%", "1-10%", "10 ~ 20%", "20 ~ 30%", "30 ~ 40%", "40 ~ 50%", "50 ~ 60%", "60 ~ 70%", "90 ~ 100%")
proportions$variant_af_group <- factor(proportions$variant_af_group, levels = custom_order)

# Make the bar plot
ggplot(proportions, aes(x = variant_af_group, y = prop)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(y = "Proportion",
       x = "AF bins") +
  theme(axis.text.x = element_text(size = 12, colour = "black", angle = 35, hjust = 1),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(colour = "black", size = 15, face = "bold"))

ggsave(filename = "af_bin_proportion.png", width = 8, height = 5, dpi = 200)
