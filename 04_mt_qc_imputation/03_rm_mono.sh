#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J remove_dup
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
logfile=${LOG_DIR}/${TIME}_${PREFIX}_rm_mono.log
# Redirect standard output and error to the log file
exec > "$logfile" 2>&1


# Remove monomorhpic SNPs
plink2 \
    --bfile ${MT_QC_OUT_DIR}/${PREFIX}_geno05_mind05_geno02 \
    --not-chr 0 \
    --mac 1 \
    --make-bed \
    --out ${MT_QC_OUT_DIR}/${PREFIX}_geno05_mind05_geno02_mono

# Convert to VCF
plink2 \
    --bfile ${MT_QC_OUT_DIR}/${PREFIX}_geno05_mind05_geno02_mono \
    --ref-from-fa --fa ${MT_FASTA} \
    --recode vcf \
    --out ${MT_QC_OUT_DIR}/${PREFIX}_geno05_mind05_geno02_mono

