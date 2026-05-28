#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J imputation
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

WKDIR="${PROJECT_DIR}/imputation_wgs1465.shapeit2.impute2"
mkdir -p $WKDIR
cd $WKDIR

INPUT="${MT_QC_OUT_DIR}/${PREFIX}_mind05_geno02_mono_remove"
PARA="twb2_study"
ref_para="wgs1465"

# ---- Redirect LOG file ----
TIME=$(date +%Y%m%d)
logfile=${LOG_DIR}/${TIME}_twb_phasing_imputation.log
exec > "$logfile" 2>&1

# 1. Phasing of study genotypes with SHAPEIT2 
# Filter out variants with AF deviated from WGS (AF_diff > 0.1)
# First get gcount file using PLINK
plink2 \
    --bfile ${INPUT} \
    --geno-counts \
    --make-bed \
    --out ${WKDIR}/${PARA}

# Filter out variants with AF difference with WGS greater than 0.1 using python script
run_in_env data_env python3 \
    "$REMOVE_DEVIATE_PY" "${REF_PANEL_VCF}" "${WKDIR}/${PARA}" > ${WKDIR}/non_deviate.snplist

plink2 \
    --bfile ${INPUT} \
    --output-chr chrM \
    --extract ${WKDIR}/non_deviate.snplist \
    --make-bed \
    --out ${WKDIR}/${PARA}


$SHAPEIT2 \
    --input-bed ${WKDIR}/${PARA}.bed ${WKDIR}/${PARA}.bim ${WKDIR}/${PARA}.fam \
    --rho 4.0E-12 \
    --thread 14 \
    --force \
    -O ${WKDIR}/twb2.shapeit2

# 2. Transform reference panel to OXFORD format
plink2 \
    --vcf ${REF_PANEL_VCF} \
    --recode oxford \
    --out ${WKDIR}/${ref_para}.ref.oxf
sed -i '2d' ${WKDIR}/${ref_para}.ref.oxf.sample
cut -d' ' -f2,3,4,5 < ${WKDIR}/${ref_para}.ref.oxf.gen > ${WKDIR}/${ref_para}.ref.legend
echo "rsID position a0 a1" > ${WKDIR}/header1.txt
cat ${WKDIR}/header1.txt ${WKDIR}/${ref_para}.ref.legend > ${WKDIR}/fin.${ref_para}.ref.legend

# 3. Imputation
${IMPUTE2} \
    -merge_ref_panels \
    -m $MTMAP \
    -h ${WKDIR}/${ref_para}.ref.oxf.gen \
    -l ${WKDIR}/fin.${ref_para}.ref.legend \
    -known_haps_g $WKDIR/twb2.shapeit2.haps \
    -int 1 16579 \
    -Ne 20000 \
    -o ${WKDIR}/twb2.EAS_kingMT.shapeit2.impute2.diploid

