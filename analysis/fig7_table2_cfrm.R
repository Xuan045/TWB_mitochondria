# Plot confirmed pathogenic variants , which are documented in MITOMAP, present in TWB
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(scales)

setwd("/Users/xuanchou/Documents/TWB_1492/Analysis/mito_cfrm/")
twb <- read.delim(file = "/Users/xuanchou/Documents/TWB_1492/Analysis/reconstructed_vcf.txt", 
                  sep = "\t", na.strings = "", comment.char = "#")
twb <- twb %>% filter(filters == "[]")

# Remove "chrM:" in hgvsg annotations
twb$hgvsg <- str_extract(twb$hgvsg, "(?<=chrM:)[^|]+")
# Standardize the annotations for hgvsg (some are "g." and some are "m.")
twb$hgvsg <- sub("g", "m", twb$hgvsg)

# Input the cfrm file downloaded from MITOMAP
mitomap_cfrm <- read.csv(file = "ConfirmedMutations_MITOMAP_20230728.csv", row.names = 1)
coding_cfrm <- read.csv(file = "cfrm_coding_20230728.csv")
coding_cfrm[c("hom", "het")] <- str_split_fixed(coding_cfrm$Plasmy.Reports.Homo.Hetero., "/", 2) 
coding_cfrm <- coding_cfrm %>% select(Allele, hom, het)

rna_cfrm <- read.csv(file = "cfrm_rna_20230728.csv") %>% 
  select(Locus, Disease, Homoplasmy, Heteroplasmy) %>% 
  rename(Associated.Diseases = Disease)
mitomap_cfrm <- merge(mitomap_cfrm, coding_cfrm, by = "Allele", all.x = TRUE) %>% 
  merge(., rna_cfrm, by = c("Locus", "Associated.Diseases") , all.x = TRUE)
mitomap_cfrm <- mitomap_cfrm %>% 
  mutate("Hom. reported" = coalesce(hom, Homoplasmy),
         "Het. reported" = coalesce(het, Heteroplasmy)) %>% 
  select(-hom, -Homoplasmy, -het, -Heteroplasmy)

names(mitomap_cfrm)[names(mitomap_cfrm) == "Allele"] <- "hgvsg"

# cfrm in TWB and count the AC (hom + het)
twb_cfrm <- inner_join(mitomap_cfrm, twb, by = "hgvsg")
twb_cfrm$count <- twb_cfrm$AC_het + twb_cfrm$AC_hom

# Calculate the heteroplasmy levels of these cfrms and save it to a new column named "HL_defined"
twb_cfrm <- twb_cfrm %>% 
  mutate(HL_list = strsplit(HL, ",")) %>%
  mutate(HL_defined = map(HL_list, ~ as.numeric(.x))) %>%
  mutate(HL_defined = map(HL_defined, ~ .x[!is.na(.x) & .x != 0])) %>%
  select(-HL_list)

# total carrier frequency - write to file
write.table(sub("^", "Total carrier frequency of pathogenic variants in TWB (>10% heteroplasmy): 1 in ", comma(ceiling(1  / ((sum(twb_cfrm$AC_hom) + sum(twb_cfrm$AC_het)) / 1465)), accuracy = 1)),
            file = "total_pathogenic_carrier_freq.txt", append = FALSE, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# dotplot of heteroplasmy level
# variant order for plotting (AC -> POS)
twb_cfrm$hgvsg <- factor(twb_cfrm$hgvsg, levels = c(unique(twb_cfrm[order(twb_cfrm$count, -twb_cfrm$POS), "hgvsg"]))) #order by count for plot
## text on graph
twb_cfrm <- twb_cfrm %>% mutate(Freq = paste0("1 in ", round(AN/count, digits = 0), "\n(", percent(count/AN, accuracy = 0.01), ")"))
twb_cfrm$Associated.Diseases <- gsub("Leigh Syndrome|Leigh Disease|Leigh syndrome|- LD", "LS", twb_cfrm$Associated.Diseases)
twb_cfrm$Associated.Diseases <- gsub("NARP-like disease", "NARP", twb_cfrm$Associated.Diseases)
twb_cfrm$Associated.Diseases <- gsub("Mitochondrial myopathy, lactic acidosis and sideroblastic anemia \\(MLASA)", "MLASA", twb_cfrm$Associated.Diseases)
twb_cfrm$Associated.Diseases <- gsub("\\/Deafness",", DEAF", twb_cfrm$Associated.Diseases)
twb_cfrm$Associated.Diseases <- ifelse(grepl("Multiorgan failure / |\\ / HCM|LDYT / |MICM\\+DEAF / |\\ / MILS|Encephalopathy / | Autism /| Other|\\ / Depressive mood disorder / leukoencephalopathy / HiCM|; autism spectrum intellectual disability; possibly antiatherosclerotic|\\ / Ataxia\\+Lipomas|\\ / carotid atherosclerosis risk|\\ / Progressive Dystonia|FBSN / |BSN / |\\ / other|\\ / FSGS / ASD / Cardiac\\+multi-organ dysfunction|\\ / dystonia|\\/ DMDF / MIDD |\\ / CPEO|\\ / IgG nephropathy", twb_cfrm$Associated.Diseases),
                                   paste(gsub("Multiorgan failure / |\\ / HCM|LDYT / |MICM\\+DEAF / |\\ / MILS|Encephalopathy / | Autism /| Other|\\ / Depressive mood disorder / leukoencephalopathy / HiCM|; autism spectrum intellectual disability; possibly antiatherosclerotic|\\ / Ataxia\\+Lipomas|\\ / carotid atherosclerosis risk|\\ / Progressive Dystonia|FBSN / |BSN / |\\ / other|\\ / FSGS / ASD / Cardiac\\+multi-organ dysfunction|\\ / dystonia|\\/ DMDF / MIDD |\\ / CPEO|\\ / IgG nephropathy", "", twb_cfrm$Associated.Diseases), ", other", sep = ""), 
                                   as.character(twb_cfrm$Associated.Diseases))
twb_cfrm$Associated.Diseases <- gsub("\\ /|;", ",", twb_cfrm$Associated.Diseases)
twb_cfrm$Associated.Diseases <- gsub("myopathy,", "MM,", twb_cfrm$Associated.Diseases)
twb_plot <- unnest(twb_cfrm, HL_defined) %>% 
  select(-starts_with("NGS"))
# Set duplicated as NA to prevent from overlapped text in graph
twb_plot$dis_text <- ifelse(duplicated(twb_plot$Associated.Diseases), NA, twb_plot$Associated.Diseases)

# Set seed so that the jitter plot would be the same every time
set.seed(3)
cfrm_plt <- ggplot(twb_plot, aes(x = HL_defined, y = hgvsg, label = Associated.Diseases)) +
  geom_jitter(cex = 2.5, position = position_jitter(h = 0.1)) +
  labs(x = "Heteroplasmy level") + 
  theme(axis.title.x = element_text(size = 13, face = "bold"), 
        axis.text.x  = element_text(size = 14), 
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size = 10),
        panel.grid.major = element_line(colour = "light grey", size = 0.2), 
        panel.grid.minor = element_line(colour = "light grey", size = 0.2),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1, 20, 1, 1), "lines"),
        plot.tag.position = c(1.27, 0.015)) + 
  xlim(0, 1) +
  scale_x_continuous(breaks = seq(0.1, 1.0, 0.1), labels = percent_format(accuracy = 1)) +
  geom_vline(xintercept = c(0.95), linetype = "solid", color = "dark grey") +
  geom_vline(xintercept = c(0.6), linetype = "dashed", color = "#F5B041") +
  geom_text(aes(label = Freq), x = 1.1, hjust = 0) +
  geom_text(x = 1.3, hjust = 0) +
  labs(tag = expression(bold("Freq. in TWB      Associated diseases"))) +
  coord_cartesian(clip = "off")
cfrm_plt

ggsave(cfrm_plt, filename = "hl_of_cfrm.png", dpi = 500, width = 10.8, height = 7.7)
 
