#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
pcadir <-args[1]
ped_pop <- args[2]  # tsv with ind_id, pop, and superpop info
pcplot0 <- TRUE     # TRUE or FALSE
npc <- 6
pred_prob <- 0.8


setwd(pcadir)
library(randomForest)
library(tidyverse)
library(ggsci)
library(RColorBrewer)


print("Read in data")
# Read in data
pca_w_ref <- read.table("twb2_ref_pca.eigenvec", header=TRUE, sep="\t", stringsAsFactors=F, comment.char = "")
names(pca_w_ref) <- gsub("X.", "", names(pca_w_ref))

# ped <- read.table("/Users/xuanchou/Documents/1KGP/20130606_g1k.ped", sep = "\t", header = TRUE)
# sup_pop <- read.table("/Users/xuanchou/Documents/1KGP/20131219.populations.tsv", sep = "\t", header = TRUE) %>% 
#   select(all_of(c("Population.Code", "Super.Population")))
# ped_pop <- merge(ped, sup_pop, by.x = "Population", by.y = "Population.Code", all.x = TRUE)
# write.table(ped_pop, file="/Users/xuanchou/Documents/1KGP/ped_pop.tsv", sep = "\t", row.names = FALSE)
ped_pop <- read.table(ped_pop, sep = "\t", header = TRUE) %>% 
  select(all_of(c("Population", "Individual.ID", "Super.Population")))

pca_w_ref <- merge(pca_w_ref, ped_pop, by.x = "IID", by.y = "Individual.ID", all.x = TRUE)
pca_w_ref <- rename(pca_w_ref, "pop" = "Population", "superpop" = "Super.Population")
pca_w_ref$superpop <- ifelse(is.na(pca_w_ref$superpop), "TWB2", pca_w_ref$superpop)

# pca_w_ref$pop <- pca_w_ref$FID
# 
# 
# pca_w_ref$superpop <- "TWB2"
# pca_w_ref$superpop <- ifelse(pca_w_ref$pop=="CHB" | pca_w_ref$pop=="JPT" | pca_w_ref$pop=="CHS" | pca_w_ref$pop=="CDX" | pca_w_ref$pop=="KHV", "EAS", pca_w_ref$superpop)
# pca_w_ref$superpop <- ifelse(pca_w_ref$pop=="CEU" | pca_w_ref$pop=="TSI" | pca_w_ref$pop=="FIN" | pca_w_ref$pop=="GBR" | pca_w_ref$pop=="IBS", "EUR", pca_w_ref$superpop)
# pca_w_ref$superpop <- ifelse(pca_w_ref$pop=="YRI" | pca_w_ref$pop=="LWK" | pca_w_ref$pop=="GWD" | pca_w_ref$pop=="MSL" | pca_w_ref$pop=="ESN" | pca_w_ref$pop=="ASW" | pca_w_ref$pop=="ACB", "AFR", pca_w_ref$superpop)
# pca_w_ref$superpop <- ifelse(pca_w_ref$pop=="MXL" | pca_w_ref$pop=="PUR" | pca_w_ref$pop=="CLM" | pca_w_ref$pop=="PEL", "AMR", pca_w_ref$superpop)
# pca_w_ref$superpop <- ifelse(pca_w_ref$pop=="GIH" | pca_w_ref$pop=="PJL" | pca_w_ref$pop=="BEB" | pca_w_ref$pop=="STU" | pca_w_ref$pop=="ITU", "SAS", pca_w_ref$superpop)

pca_w_ref$superpop <- factor(pca_w_ref$superpop, levels=c("TWB2", "AFR", "AMR", "EAS", "EUR", "SAS"))


print("Make PC plots: twb + ref samples, colored by super population")
# PC plots: twb + ref samples, colored by super population
if ( pcplot0 == "TRUE" ){
  for(ii in 1:5){
    p = ggplot(pca_w_ref, aes_string(x=paste0('PC',ii), y=paste0('PC',ii+1))) +
      geom_point(aes(color=superpop), alpha=0.57, size=0.97) +
      scale_color_manual(values = c("grey25", pal_d3("category20")(20)[c(1:5)])) + 
      theme_bw() +
      guides(color = guide_legend(override.aes = list(size=2))) +
      labs(x=paste0('PC',ii), y=paste0('PC',ii+1),
           color="POP")
    
    ggsave(paste0("TWB2_1KG_PC",ii,"_PC",ii+1,"_by_superpop.pdf"), p, width=7.5, height=6)
  }
}



print(paste0("Use RF to predict ancestry based on top ", npc, " PCs"))
# RF function to predict ancestry using PCs:
pop_forest <- function(training_data, data, ntree=100, seed=42, pcs=1:npc) {
  set.seed(seed)
  form <- formula(paste('as.factor(known_pop) ~', paste0('PC', pcs, collapse = ' + ')))
  forest <- randomForest(form, data = training_data, importance = T, ntree = ntree)
  print(forest)
  
  fit_data <- data.frame(predict(forest, data, type='prob'), sample = data$sample)
  fit_data %>%
    gather(predicted_pop, probability, -sample) %>%
    group_by(sample) %>%
    slice(which.max(probability))
}

# Prep data
# Explanation:
# `trdat` and `tedat` are training and testing data. training data has to have a 
# column `known_pop` and `PC1` to `PC10` or so. `tedat` is expected to have a 
# column `sample` which is just the sample ID, and also PC columns.

trdat <- pca_w_ref %>%
  filter(superpop != 'TWB2') %>%
  select(superpop, PC1:PC10) %>%
  rename(known_pop = superpop)
trdat$known_pop <- as.character(trdat$known_pop)

tedat <- pca_w_ref %>%
  filter(superpop == 'TWB2') %>%
  select(IID, PC1:PC10)
names(tedat)[1] <- "sample"
tedat$sample <- as.character(tedat$sample)


print("Prediction results:")
# Make prediction
pop_pred <- as.data.frame(pop_forest(training_data = trdat, data = tedat))

print("Overall ancestry assignment:")
table(pop_pred$predicted_pop)

print(paste0("Subset to prediction prob. > ", pred_prob, ":"))
pop_pred.sub <- pop_pred %>% filter(probability > pred_prob)
table(pop_pred.sub$predicted_pop)

# Remove non-EAS samples and samples with low pred prob
names(pop_pred)[1] <- "IID"
remove_indv <- pop_pred %>% filter(predicted_pop!="EAS" | probability<pred_prob)
pca_w_ref_sub <- pca_w_ref[ ! pca_w_ref$IID %in% remove_indv$IID, ]

print("Make PC plots: qc'ed twb + ref samples, colored by super population")
# PC plots: qc'ed twb + ref samples, colored by super population
if ( pcplot0 == "TRUE" ){
  for(ii in 1:5){
    p = ggplot(pca_w_ref_sub, aes_string(x=paste0('PC',ii), y=paste0('PC',ii+1))) +
      geom_point(aes(color=superpop), alpha=0.57, size=0.97) +
      scale_color_manual(values = c("grey25", pal_d3("category20")(20)[c(1:5)])) +
      theme_bw() +
      guides(color = guide_legend(override.aes = list(size=2))) +
      labs(x=paste0('PC',ii), y=paste0('PC',ii+1),
           color="POP")
    
    ggsave(paste0("TWB2_1KG_QC_PC",ii,"_PC",ii+1,"_by_superpop.pdf"), p, width=7.5, height=6)
  }
}


# Merge back with pca_w_ref to udpate PC plots (colored by pred. super population)
pca_w_ref.twb <- subset(pca_w_ref, superpop=="TWB2")
pca_w_ref.pred <- merge(pca_w_ref.twb, pop_pred, by="IID")
pca_w_ref.pred$predicted_pop <- factor(pca_w_ref.pred$predicted_pop, levels=c("AFR", "AMR", "EAS", "EUR", "SAS"))
pca_w_ref.pred.sub <- subset(pca_w_ref.pred, probability > pred_prob)

print("Export predicted ancestral groups for study samples:")
# Save predicted population labels and EUR family and individual IDs
write.table(pca_w_ref.pred, "twb2_1kg_predpop.tsv", quote=F, col.names=T, row.names=F, sep="\t")

eas_indv <- pca_w_ref.pred.sub[pca_w_ref.pred.sub$predicted_pop=="EAS", c("FID","IID")]
write.table(eas_indv, paste0("twb2_1kg_predpop",pred_prob,"_EAS_indvlist"), quote=F, col.names=F, row.names=F, sep="\t")

non_eas_indv <- pca_w_ref.pred.sub[pca_w_ref.pred.sub$predicted_pop!="EAS", c("FID","IID")]
write.table(non_eas_indv, paste0("twb2_1kg_predpop",pred_prob,"_non_EAS_indvlist"), quote=F, col.names=F, row.names=F, sep="\t")

