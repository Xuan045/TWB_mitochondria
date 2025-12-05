library(tidyverse)
library(ggplot2)
library(patchwork)

setwd("/Users/xuanchou/Documents/TWB_1492/Analysis/pie_chart/")
data_text <- read_lines("/Users/xuanchou/Documents/TWB_1492/add_annotation/output_age/stats_pass.txt",
                        skip = 1)

data_list <- strsplit(data_text, "\n")
data_list <- unlist(data_list)

# Extract metrics and values from the text
metrics <- gsub("\"", "", sapply(strsplit(data_list, ": "), "[", 1))
values <- as.numeric(gsub("[[:alpha:] ]", "", sapply(strsplit(data_list, ": "), "[", 2)))

# Create the stat_pass data frame
stat_pass <- data.frame(metric = metrics, value = values)

# Calculate percentages for pie charts
total_bases_with_variation <- stat_pass$value[6]
total_bases <- 16569
percentage_bases_with_variants <- (total_bases_with_variation / total_bases) * 100

total_unique_variants <- stat_pass$value[5]
total_SNV_transition <- stat_pass$value[20]
total_SNV_transversion <- stat_pass$value[21]
total_indels <- stat_pass$value[19]
percentage_SNV_transition <- (total_SNV_transition / total_unique_variants) * 100
percentage_SNV_transversion <- (total_SNV_transversion / total_unique_variants) * 100
percentage_indel <- (total_indels / total_unique_variants) * 100

total_homoplasmic_only <- stat_pass$value[9]
total_hom_het_sites <- stat_pass$value[11]
total_heteroplasmic_only <- stat_pass$value[10]
percentage_homoplasmic_only <- (total_homoplasmic_only / total_unique_variants) * 100
percentage_hom_het_sites <- (total_hom_het_sites / total_unique_variants) * 100
percentage_heteroplasmic_only <- (total_heteroplasmic_only / total_unique_variants) * 100

total_variants <- stat_pass$value[8]
total_homoplasmic_variants <- stat_pass$value[12]
total_heteroplasmic_variants <- stat_pass$value[14]
percentage_homoplasmic_variants <- (total_homoplasmic_variants / total_variants) * 100
percentage_heteroplasmic_variants <- (total_heteroplasmic_variants / total_variants) * 100

# Create data frames for pie charts
df_bases <- data.frame(values = c(percentage_bases_with_variants, 100 - percentage_bases_with_variants), labels = c("Bases with Variants", "Bases without Variants"))
df_SNV <- data.frame(values = c(percentage_SNV_transition, percentage_SNV_transversion, percentage_indel), labels = c("SNV Transition", "SNV Transversion", "Indel"))
df_hom_het <- data.frame(values = c(percentage_homoplasmic_only, percentage_hom_het_sites, percentage_heteroplasmic_only), labels = c("Homoplasmic Only", "Homoplasmic/Heteroplasmic", "Heteroplasmic Only"))
df_hom_het_variants <- data.frame(values = c(percentage_homoplasmic_variants, percentage_heteroplasmic_variants), labels = c("Homoplasmic Variants", "Heteroplasmic Variants"))

# Function to create a pie chart with percentage labels
create_pie_chart <- function(data_frame, title, colors) {
  pie_chart <- ggplot(data_frame, aes(x = 1, y = values, fill = labels)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    labs(fill = title) +
    theme_void() +
    geom_text(aes(label = paste0(round(values), "%")), position = position_stack(vjust = 0.5),
              colour = "white", fontface = "bold") +
    scale_fill_manual(values = colors)
  return(pie_chart)
}

# Create pie charts
pie_chart_bases <- create_pie_chart(df_bases, "16569 mtDNA bases", c("#627264", "#A1CDA8"))
pie_chart_SNV <- create_pie_chart(df_SNV, paste(total_unique_variants, "variants"), c("#8ECAE6", "#219EBC", "#023047"))
pie_chart_hom_het <- create_pie_chart(df_hom_het, paste(total_unique_variants, "variants"), c("#006D77", "#42999B", "#83C5BE"))
pie_chart_hom_het_variants <- create_pie_chart(df_hom_het_variants, paste(total_variants, "variant calls"), c("#0081A7", "#00AFB9"))

# Combine pie charts using patchwork
combined_pie_sumbit <- pie_chart_SNV + pie_chart_hom_het
ggsave(combined_pie_sumbit, filename ="combined_pie_sumbit.png", dpi = 500, width = 9, height = 3.5)

combined_pie_charts <- pie_chart_bases + pie_chart_SNV + pie_chart_hom_het + pie_chart_hom_het_variants

# Display the combined pie charts
print(combined_pie_charts)

ggsave(combined_pie_charts, filename ="combined_pie.png", dpi = 500, width = 9, height = 5.5)

# Haplogroup-defining variants
var_df <- read.delim("/Users/xuanchou/Documents/TWB_1492/add_annotation/output_age/combined_sites_only.txt") %>% 
  filter(filters == "[]" & AC_hom != 0)
hap_marker_af <- var_df %>% 
  filter(hap_defining_variant == "false") %>% 
  arrange(desc(AF_hom)) %>% 
  select(locus, alleles, hap_defining_variant, AC_hom, AF_hom) %>% 
  filter(AF_hom > 0.01)
write_csv(hap_marker_af, file = "hap_marker_high_af.csv")

hap_marker <- sum(var_df$hap_defining_variant == "true") / nrow(var_df) *100
df_hap_marker <- data.frame(values = c(hap_marker, 100 - hap_marker), labels = c("Hap-defining variants", "Not hap-defining variants"))
pie_chart_marker <- create_pie_chart(df_hap_marker, "Homoplasmies", c("#627264", "#A1CDA8"))
ggsave(filename = "marker_pie.png", pie_chart_marker, width = 6.5, height = 3.5)
