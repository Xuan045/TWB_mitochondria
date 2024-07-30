library(tidyverse)
library(ggpubr)

############################
# input files
############################
setwd("/Users/xuanchou/Documents/TWB_1492/Analysis/annotation/")
pass <- read.delim("/Users/xuanchou/Documents/TWB_1492/add_annotation/output_age/combined_sites_only.txt", stringsAsFactors=FALSE) 
pass[] <- lapply(pass, function(x) if(is.character(x)) gsub("\\[|\\]", "", x) else x)
pass <- pass %>% separate(alleles, into = c("REF", "ALT"), sep = ",", remove = FALSE)
pass$POS <- gsub("chrM:", "", pass$locus)

pass$count <- 1

######################################################################################
# 3E: stacked bar x variant type (CR, intergenic, tRNA, rRNA, coding) file=var_by_type.combo.pdf
# make normalized stacked bar chart by variant type (CR, intergenic, tRNA, rRNA, coding LOW,MODERATE,HIGH)#
######################################################################################

# step 1: read in annotations of all chrM positions, and get 2 columns: varType and subset
annot <- read.delim("stack_bar/chrM.pos2annot.txt", stringsAsFactors = FALSE)
annot[annot$Type=="coding","Type"]="protein-coding (all)"
annot$varType <- annot$Type
annot$subset <- "mtDNA\npositions\nn=16569"
annot$varType <- ifelse(annot$varType == "CR", "Control region", annot$varType)

# step 2: get all unique variants and exclude indel_stack
# assign these varType and subset too
m <- pass %>%
  filter(nchar(REF) == 1 & nchar(ALT) == 1 & filters == "")

# subset = homoplasmic-haplogroup, homoplasmic-non-haplogroup, heteroplasmic-only (add n counts)
m$subset <- "none"
homMarker <- unique(m[(m$max_hl >= 0.95) & (m$hap_defining_variant == "true"), "variant_collapsed"])
homNoMarker <- unique(m[(m$max_hl >= 0.95) & (m$hap_defining_variant=="false"),"variant_collapsed"])
hetOnly <- unique(m[!(m$variant_collapsed %in% c(homMarker,homNoMarker)),"variant_collapsed"])
#
commonLowHet <- unique(m[!(m$variant_collapsed %in% c(homMarker,homNoMarker)) & m$common_low_heteroplasmy == "true", "variant_collapsed"])
#
m[m$variant_collapsed %in% homMarker,"subset"] <- paste("homoplasmic\n(haplogroup)\nn=",length(homMarker),sep="")
m[m$variant_collapsed %in% homNoMarker,"subset"] <- paste("homoplasmic\n(non-haplogroup)\nn=",length(homNoMarker),sep="")
m[m$variant_collapsed %in% hetOnly,"subset"] <- paste("heteroplasmic\nonly\nn=",length(hetOnly),sep="")
m[m$variant_collapsed %in% commonLowHet,"subset"] <- paste0("common low\nheteroplasmy\nn=",length(commonLowHet))

# assign mu to be just unique variants
mu <- unique(m[,c("variant_collapsed","subset","count")])

# pull out the VEP annotations and assign to categories
vep <- read.delim("/Users/xuanchou/Documents/TWB_1492/Analysis/reformat_vcf/reformated.vcf", stringsAsFactors=FALSE)
vep$variant_collapsed <- paste0(vep$REF, vep$POS, vep$ALT)

mu$varType <- "none"
mu[mu$variant_collapsed %in% vep[vep$IMPACT %in% c("HIGH"),"variant_collapsed"],"varType"] <- "non-synonymous"
mu[mu$variant_collapsed %in% vep[vep$IMPACT %in% c("MODERATE"),"variant_collapsed"],"varType"] <- "non-synonymous"
mu[mu$variant_collapsed %in% vep[vep$IMPACT %in% c("LOW"),"variant_collapsed"],"varType"] <- "synonymous"
mu[mu$variant_collapsed %in% vep[vep$BIOTYPE %in% c("Mt_tRNA"),"variant_collapsed"],"varType"] <- "tRNA"
mu[mu$variant_collapsed %in% vep[vep$BIOTYPE %in% c("Mt_rRNA"),"variant_collapsed"],"varType"] <- "rRNA"
mu[mu$variant_collapsed %in% vep[(vep$BIOTYPE=="") & (vep$POS >= 577) & (vep$POS < 16024),"variant_collapsed"],"varType"] <- "intergenic"
mu[mu$variant_collapsed %in% vep[(vep$IMPACT=="MODIFIER") & (vep$BIOTYPE=="protein_coding"),"variant_collapsed"],"varType"] <- "intergenic"
mu[mu$variant_collapsed %in% vep[(vep$BIOTYPE=="") & !((vep$POS >= 577) & (vep$POS < 16024)),"variant_collapsed"],"varType"] <- "Control region"

# step 3: combine the mtDNA annotations plus the unique variant annotations
tmp1 <- annot[,c("subset","varType")]
tmp2 <- mu[,c("subset","varType")]
tmp3 <- rbind(tmp1,tmp2)

# create factors
tmp3$varType <- factor(tmp3$varType,levels=c("protein-coding (all)","non-synonymous","synonymous","rRNA","tRNA","intergenic","Control region"))
tmp3$subset <- factor(tmp3$subset,levels = rev(c("mtDNA\npositions\nn=16569", "homoplasmic\n(haplogroup)\nn=1715", "homoplasmic\n(non-haplogroup)\nn=465","heteroplasmic\nonly\nn=109", "common low\nheteroplasmy\nn=16")))

# now plot
mycolors=c("olivedrab","olivedrab3","olivedrab2","lightsteelblue","lightsteelblue1","gray85","black")

theme_custom <- function() {
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold")
  )
}

stack_bar <- ggplot(tmp3, aes(x = subset, fill = varType)) +
  geom_bar(position="fill") +
  scale_fill_manual(values=mycolors) +
  ggthemes::theme_economist_white(gray_bg=FALSE) +
  theme_custom() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(0, 1, 0, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +  # Adjust y-axis scale due to coord_flip()
  coord_flip() +
  guides(fill = guide_legend(nrow = 3, byrow = FALSE))

# ggsave(filename = "annotation_stack_add_common_low.png", dpi = 400, width = 12, height = 7)


twb <- read.delim(file = '/Users/xuanchou/Documents/TWB_1492/Analysis/reformat_vcf/reformated.vcf', header = TRUE, sep = "\t", na.strings = "")

twb$locate <- ifelse(twb$Consequence == "stop_gained", "Stop gain",
                     ifelse(twb$BIOTYPE == "Mt_tRNA" & !is.na(twb$BIOTYPE), "tRNA",
                            ifelse(twb$BIOTYPE == "Mt_rRNA" & !is.na(twb$BIOTYPE), "rRNA",
                                   ifelse(twb$Consequence == "missense_variant", "Missense",
                                          ifelse(twb$Consequence == "synonymous_variant", "Synonymous",
                                                 ifelse(twb$VARIANT_CLASS != "SNV" & twb$BIOTYPE == "protein_coding" & !is.na(twb$BIOTYPE), "Protein coding\nindel",
                                                        ifelse(twb$Consequence == "intergenic_variant" & (twb$POS < 577 | twb$POS > 16023), "Control region", 
                                                               ifelse(twb$Consequence == "intergenic_variant", "Intergenic region", "Other"))))))))

# Custom theme function for enhanced labels
theme_custom <- function() {
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold")
  )
}

twb$locate <- factor(twb$locate, levels = c("Control region", "Intergenic region", "Synonymous", "Missense", "Stop gain", "Protein coding\nindel", "rRNA", "tRNA", "Other"))

# Plot 'consequence' with custom theme
consequence_plot <- ggplot(twb, aes(x=locate)) +
  geom_bar(fill="skyblue") +
  theme_bw() +
  theme_custom() +
  labs(title="Distribution of Consequence", x="Consequence", y="Count") +
  geom_text(aes(label=after_stat(count)), stat='count', vjust=-0.5) +
  scale_y_continuous(limits = c(0, 1200), expand = c(0, 0))

# Save the combined plot
outdir <- "thesis/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Combine the plots and add labels
combined_plots <- ggarrange(consequence_plot, stack_bar, 
                            labels = c("a", "b"),
                            ncol = 1, nrow = 2, 
                            font.label = list(size = 20, family = "Times New Roman"))
combined_plots
ggsave(combined_plots, filename = "thesis/var_func.png", height = 11, width = 12)

