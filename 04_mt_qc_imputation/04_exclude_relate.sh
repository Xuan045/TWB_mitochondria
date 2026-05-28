#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J exclude_relate
#SBATCH -p ngs186G
#SBATCH -c 28
#SBATCH --mem=175G
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
logfile=${LOG_DIR}/${TIME}_${PREFIX}_king.log
# Redirect standard output and error to the log file
exec > "$logfile" 2>&1

# merge two lists
cat "$RELATED_INDV_LIST" "$NON_EAS_INDV_LIST" | sort | uniq > ${MT_QC_OUT_DIR}/${PREFIX}_remove.indv

PARA="mind05_geno02_mono_remove"
plink2 \
    --bfile ${MT_QC_OUT_DIR}/${PREFIX}_geno05_mind05_geno02_mono \
    --remove ${MT_QC_OUT_DIR}/${PREFIX}_remove.indv \
    --make-bed \
    --out ${MT_QC_OUT_DIR}/${PREFIX}_${PARA}

plink2 \
    --bfile ${MT_QC_OUT_DIR}/${PREFIX}_${PARA} \
    --recode vcf id-paste=iid \
    --out ${MT_QC_OUT_DIR}/${PREFIX}_${PARA}

# PCA on study sample
plink2 \
    --bfile ${MT_QC_OUT_DIR}/${PREFIX}_${PARA} \
	--allow-extra-chr \
	--make-bed \
	--output-chr 26 \
	--chr 26 \
	--out ${MT_QC_OUT_DIR}/${PREFIX}_${PARA}_pca

plink2 \
    --bfile ${MT_QC_OUT_DIR}/${PREFIX}_${PARA}_pca \
    --chr-set 30 \
    --pca 20 \
    --out ${MT_QC_OUT_DIR}/${PREFIX}_${PARA}

