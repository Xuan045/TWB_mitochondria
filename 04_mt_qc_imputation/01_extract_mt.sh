#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J extract_mt
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
mkdir -p $LOG_DIR
TIME=$(date +%Y%m%d)
logfile=${LOG_DIR}/${TIME}_${PREFIX}_extractMT.log
# Redirect standard output and error to the log file
exec > "$logfile" 2>&1

mkdir -p "$MT_QC_OUT_DIR"

# 1. Extract SNVs on MT
plink2 \
    --bfile ${INPUT_PLINK} \
    --chr 26 \
    --make-bed \
    --missing \
    --snps-only \
    --out ${MT_QC_OUT_DIR}/${PREFIX}_mt

# 2. Reorder the ref, alt allele
plink2 \
    --bfile ${MT_QC_OUT_DIR}/${PREFIX}_mt \
    --normalize \
    --ref-from-fa -fa ${MT_FASTA} \
    --make-bed \
    --out ${MT_QC_OUT_DIR}/${PREFIX}_mt_recode

# 2-2. Calculate the AF
plink2 \
    --bfile ${MT_QC_OUT_DIR}/${PREFIX}_mt_recode \
    --chr 26 \
    --geno-counts \
    --freq \
    --out ${MT_QC_OUT_DIR}/${PREFIX}_mt_recode

# 3. Export VCF
plink2 \
    --bfile ${MT_QC_OUT_DIR}/${PREFIX}_mt_recode \
    --export vcf \
    --out ${MT_QC_OUT_DIR}/${PREFIX}_mt_recode
bcftools +fill-tags ${MT_QC_OUT_DIR}/${PREFIX}_mt_recode.vcf -- -t AN,AC,AF > ${MT_QC_OUT_DIR}/${PREFIX}_mt_recode_tmp.vcf
mv ${MT_QC_OUT_DIR}/${PREFIX}_mt_recode_tmp.vcf ${MT_QC_OUT_DIR}/${PREFIX}_mt_recode.vcf

# QC summary for samples and SNPs
nindv=$(cat ${INPUT_PLINK}.fam | wc -l)
nsnps=$(cat ${INPUT_PLINK}.bim | wc -l)
echo "Nindv "$nindv > ${PROJECT_DIR}/qc_${PREFIX}_qc_summary.txt
echo "Nsnps "$nsnps >> ${PROJECT_DIR}/qc_${PREFIX}_qc_summary.txt
echo -en "\n" >> ${PROJECT_DIR}/qc_${PREFIX}_qc_summary.txt

echo "Nsnps on chrM "$(cat ${MT_QC_OUT_DIR}/${PREFIX}_mt_recode.bim | wc -l) >> ${MT_QC_OUT_DIR}/${PREFIX}_qc_summary.txt
echo -en "\n" >> ${MT_QC_OUT_DIR}/${PREFIX}_qc_summary.txt
echo "snp missing rate:" >> ${MT_QC_OUT_DIR}/${PREFIX}_qc_summary.txt
awk 'NR > 1 { sum += $5; if (NR == 2 || $5 < min) min = $5; if (NR == 2 || $5 > max) max = $5; } END { print "Min:", min; print "Max:", max; print "Average:", sum / (NR - 1); }' ${MT_QC_OUT_DIR}/${PREFIX}_mt.vmiss >> ${MT_QC_OUT_DIR}/${PREFIX}_qc_summary.txt

echo -en "\nindiviual missing rate:\n" >> ${MT_QC_OUT_DIR}/${PREFIX}_qc_summary.txt
awk 'NR > 1 { sum += $6; if (NR == 2 || $6 < min) min = $6; if (NR == 2 || $6 > max) max = $6; } END { print "Min:", min; print "Max:", max; print "Average:", sum / (NR - 1); }' ${MT_QC_OUT_DIR}/${PREFIX}_mt.smiss >> ${MT_QC_OUT_DIR}/${PREFIX}_qc_summary.txt

# Plot the distribution of individual missing rate
run_in_env r_env Rscript 01_idv_missing.r "${MT_QC_OUT_DIR}" "${PREFIX}_mt.smiss"


