#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
preqcdir <- args[1]

setwd(preqcdir)
library(ggplot2)


# Read in data
imp_sex <- read.table("twb2_EAS_sex.sexcheck", header=T)

imp_sex$PEDSEX <- factor(imp_sex$PEDSEX)
levels(imp_sex$PEDSEX) <- c("Male", "Female")
imp_sex$SNPSEX <- factor(imp_sex$SNPSEX)
levels(imp_sex$SNPSEX) <- c("Unknown", "Male", "Female")


# Plot F-stat distribution: histogram
p <- ggplot(imp_sex, aes(x=F, fill=SNPSEX)) +
    geom_histogram(alpha=0.6, binwidth=0.02, color="black") +
    geom_vline(xintercept=0.8, linetype="dashed") + #default cutoff: >0.8: male
    geom_vline(xintercept=0.2, linetype="dashed") + #default cutoff: <0.2: female
    theme_bw() +
    labs(x='F-statistic',y='Frequency', fill='Imputed sex') +
    scale_fill_manual(values=c("grey30","#56B4E9","#FD8D3C"), labels=c("Unknown","Male","Female"))

ggsave("twb2_EAS_sexcheck_Fstat_F02_M08.pdf", p, width=7.5, height=4)


# Plot: imputed vs. reported gender
pdf("twb2_EAS_sexcheck_Fstat_F02_M08_pedsex_vs_snpsex.pdf", width=7.5, height=5)
ggplot(imp_sex, aes(y=F, x=factor(SNPSEX), color=factor(SNPSEX))) +
    geom_jitter(alpha=0.7) +
    labs(x='Imputed Gender', y='F-statistic', color='Imputed Gender') +
    scale_color_manual(values=c("#7F7F7FFF","#56B4E9","#FD8D3C")) +
    theme_bw()

ggplot(imp_sex, aes(y=F, x=factor(SNPSEX), color=factor(PEDSEX))) +
    geom_jitter(alpha=0.7) +
    labs(x='Imputed Gender', y='F-statistic', color='Reported Gender') +
    scale_color_manual(values=c("#56B4E9","#FD8D3C")) +
    theme_bw()

dev.off()


# Write out a list of sex discordant samples for removal
sex_mismatch_F02_M08_indv <- imp_sex[(imp_sex$STATUS=="PROBLEM"), c("FID","IID")]

imp_sex$PEDSEX <- as.character(imp_sex$PEDSEX)
imp_sex$SNPSEX <- as.character(imp_sex$SNPSEX)
imp_sex$SNPSEX <- ifelse(imp_sex$F > 0.75, "Male", imp_sex$SNPSEX)
imp_sex$SNPSEX <- ifelse(imp_sex$F < 0.25, "Female", imp_sex$SNPSEX)
sex_mismatch_F025_M075_indv <- imp_sex[imp_sex$SNPSEX != imp_sex$PEDSEX, c("FID","IID")]

write.table(sex_mismatch_F025_M075_indv, "twb2_EAS_sex_mismatch_F025_M075.indvlist", quote=F, col.names=F, row.names=F)


