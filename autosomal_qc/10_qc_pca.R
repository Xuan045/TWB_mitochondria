#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
pcadir <- args[1]


setwd(pcadir)
library(ggplot2)
library(ggsci)


# List of outliers
# remove_indv <- read.table("twb02_EAS_pc_outlier_6sd.indvlist", header=F, sep="\t", stringsAsFactors=F)
# names(remove_indv)[1] <- "FID"
# names(remove_indv)[2] <- "IID"


# Read in data: pop-specific TWB samples PCA
pca <- read.table("twb2_EAS_pca_qc.eigenvec", header=T, sep="\t", stringsAsFactors=F, comment.char = "")
names(pca) <- gsub("X.", "", names(pca))
# pca <- pca[ ! pca$IID %in% remove_indv$IID, ]
pca$pop <- "TWB02"


# PC plots: pop-specific TWB samples
for(ii in 1:5){
    p = ggplot(pca, aes_string(x=paste0('PC',ii), y=paste0('PC',ii+1))) +
    geom_point(aes(color=pop), alpha=0.57, size=0.97) +
    scale_color_manual(values = pal_d3("category20")(20)[1]) + 
    theme_bw() +
    guides(color = guide_legend(override.aes = list(size=2))) +
    labs(x=paste0('PC',ii), y=paste0('PC',ii+1),
    color="Sample")

    ggsave(paste0("twb2_EAS_QC_PC",ii,"_PC",ii+1,".pdf"), p, width=7.5, height=6)
}


# # Read in data: pop-specific TWB samples + ref data PCA
# pca_w_ref <- read.table("twb02_EAS_ref_pca_qc.eigenvec", header=T, sep="\t", stringsAsFactors=F)
# # pca_w_ref <- pca_w_ref[ ! pca_w_ref$IID %in% remove_indv$IID, ]

# pca_w_ref$pop <- pca_w_ref$FID
# pca_w_ref$pop <- ifelse(pca_w_ref$FID!="CHB" & pca_w_ref$FID!="JPT" & pca_w_ref$FID!="CHS" & pca_w_ref$pop!="CDX" & pca_w_ref$pop!="KHV", "TWB02", pca_w_ref$pop)
# pca_w_ref$pop = factor(pca_w_ref$pop, levels=c("TWB02", "CHB", "JPT", "CHS", "CDX", "KHV"))
# pca_w_ref <- pca_w_ref[order(pca_w_ref$pop),]
# table(pca_w_ref$pop)


# # PC plots: pop-specific TWB + ref samples, colored by 1kg population
# for(ii in 1:5){
#     p = ggplot(pca_w_ref, aes_string(x=paste0('PC',ii), y=paste0('PC',ii+1))) +
#     geom_point(aes(color=pop), alpha=0.55) +
#     scale_color_manual(values = c("grey25",pal_d3("category20")(20)[c(1:nlevels(pca_w_ref$pop))])) + 
#     theme_bw() +
#     guides(color = guide_legend(override.aes = list(size=2))) +
#     labs(x=paste0('PC',ii), y=paste0('PC',ii+1),
#     color="POP")

#     ggsave(paste0("TWB02_EAS_1KG_QC_PC",ii,"_PC",ii+1,"_by_pop.pdf"), p, width=7.5, height=6)
# }


