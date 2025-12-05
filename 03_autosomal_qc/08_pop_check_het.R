#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
preqcdir <- args[1]

setwd(preqcdir)
library(ggplot2)


# Read in data
het <- read.table("twb2_btqc_EAS_inbr.hetrate", header=T)


# Plot F-stat distribution
p1 <- ggplot(het, aes(F)) +
    geom_histogram(alpha=0.6, binwidth=0.01, color="black", fill="#4DBBD5B2") +
    xlab("Inbreeding Coefficient (F-stat)") +
    ylab("Frequency") +
    geom_vline(xintercept=-0.2, linetype="dashed") +
    geom_vline(xintercept=0.2, linetype="dashed") +
    theme_bw()

ggsave("twb2_EAS_hetcheck_Fstat_thresh02.pdf", p1, width=7, height=4)


# Plot heterozygosity rate distribution
cutoff1 <- c(mean(het$HetRate)-3*sd(het$HetRate), mean(het$HetRate)+3*sd(het$HetRate))
cutoff2 <- c(mean(het$HetRate)-6*sd(het$HetRate), mean(het$HetRate)+6*sd(het$HetRate))

p2 <- ggplot(het, aes(HetRate)) +
    geom_histogram(alpha=0.6, binwidth=0.001, color="black", fill="#4DBBD5B2") +
    xlab("Heterozygosity Rate") +
    ylab("Frequency") +
    geom_vline(xintercept=cutoff1[1], linetype="dashed") +
    geom_vline(xintercept=cutoff1[2], linetype="dashed") +
    geom_vline(xintercept=cutoff2[1], linetype="dashed") +
    geom_vline(xintercept=cutoff2[2], linetype="dashed") +
    annotate('text', label="mean - 3SD", x=cutoff1[1], y=4000, angle=90, size=2.5, vjust=-0.5) +
    annotate('text', label="mean + 3SD", x=cutoff1[2], y=4000, angle=90, size=2.5, vjust=-0.5) +
    annotate('text', label="mean - 6SD", x=cutoff2[1], y=4000, angle=90, size=2.5, vjust=-0.5) +
    annotate('text', label="mean + 6SD", x=cutoff2[2], y=4000, angle=90, size=2.5, vjust=-0.5) +
    theme_bw()

ggsave("twb2_EAS_hetcheck_hetrate_dev_from_mean.pdf", p2, width=7, height=4)


# Write out a list of het outlier samples for removal: given a few different thresholds
het_fstat_outliers <- het[(het$F < -0.2 | het$F > 0.2), c("FID","IID")]
write.table(het_fstat_outliers, "twb2_EAS_het_outlier_f02.indvlist", quote=F, col.names=F, row.names=F)

# hetrate_cutoff1_outliers <- het[(het$HetRate < cutoff1[1] | het$HetRate > cutoff1[2]), c("FID","IID")]
hetrate_cutoff2_outliers <- het[(het$HetRate < cutoff2[1] | het$HetRate > cutoff2[2]), c("FID","IID")]
# write.table(hetrate_cutoff1_outliers, paste0(pop,"_pbk_het_outlier_3sd.indlist"), quote=F, col.names=F, row.names=F)
write.table(hetrate_cutoff2_outliers, "twb2_EAS_het_outlier_6sd.indvlist", quote=F, col.names=F, row.names=F)


