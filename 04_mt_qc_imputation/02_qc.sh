#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J mt_qc
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92G
#SBATCH -o /dev/null
#SBATCH -e /dev/null

# ---- Source Configuration ----
if [ -f "./config.sh" ]; then
    source "./config.sh"
elif [ -n "$SLURM_SUBMIT_DIR" ] && [ -f "$SLURM_SUBMIT_DIR/config.sh" ]; then
    source "$SLURM_SUBMIT_DIR/config.sh"
else
    source "/staging/reserve/jacobhsu/TWB/TWBR11106-05/Phenotypes_TX/TWB1492_mt/microarray/mt_qc/config.sh"
fi

# ---- Environment ----
# Environment is initialized in config.sh
activate_env assoc_env


# ---- Redirect LOG file ----
TIME=$(date +%Y%m%d)
logfile=${LOG_DIR}/${TIME}_${PREFIX}_qc.log
# Redirect standard output and error to the log file
exec > "$logfile" 2>&1


# 1. Output a list of SNPs with call rate >0.95
plink \
    --bfile ${MT_QC_OUT_DIR}/${PREFIX}_mt_recode \
    --chr 26 \
    --geno 0.05 \
    --write-snplist \
    --out ${MT_QC_OUT_DIR}/${PREFIX}_geno05

# 2. Output a list of samples with call rate >0.95
plink2 \
    --bfile ${MT_QC_OUT_DIR}/${PREFIX}_mt_recode \
    --extract ${MT_QC_OUT_DIR}/${PREFIX}_geno05.snplist \
    --missing \
    --out ${MT_QC_OUT_DIR}/${PREFIX}_geno05_indv
awk 'NR>1{if($6<0.05){print $1,$2}}' ${MT_QC_OUT_DIR}/${PREFIX}_geno05_indv.smiss > ${MT_QC_OUT_DIR}/${PREFIX}_geno05_mind05.indvlist

# Plot the distribution of individual missing rate
run_in_env r_env Rscript 01_idv_missing.r "${MT_QC_OUT_DIR}" "${PREFIX}_geno05_indv.smiss"


# 3. Filter to samples in the previous step, output a list of SNPs with call rate >0.98
plink2 \
    --bfile ${MT_QC_OUT_DIR}/${PREFIX}_mt_recode \
    --extract ${MT_QC_OUT_DIR}/${PREFIX}_geno05.snplist \
    --keep ${MT_QC_OUT_DIR}/${PREFIX}_geno05_mind05.indvlist \
    --geno 0.02 \
    --write-snplist \
    --out ${MT_QC_OUT_DIR}/${PREFIX}_geno05_mind05_geno02


# 4. Keep samples and SNPs passing previous filters, calculate SNP-level missing rate
plink2 \
    --bfile ${MT_QC_OUT_DIR}/${PREFIX}_mt_recode \
    --extract ${MT_QC_OUT_DIR}/${PREFIX}_geno05_mind05_geno02.snplist \
    --keep ${MT_QC_OUT_DIR}/${PREFIX}_geno05_mind05.indvlist \
    --missing \
    --out ${MT_QC_OUT_DIR}/${PREFIX}_geno05_mind05_geno02_miss

plink2 \
    --bfile ${MT_QC_OUT_DIR}/${PREFIX}_mt_recode \
    --extract ${MT_QC_OUT_DIR}/${PREFIX}_geno05_mind05_geno02.snplist \
    --keep ${MT_QC_OUT_DIR}/${PREFIX}_geno05_mind05.indvlist \
    --make-bed \
    --out ${MT_QC_OUT_DIR}/${PREFIX}_geno05_mind05_geno02


