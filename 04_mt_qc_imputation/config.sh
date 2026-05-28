#!/bin/bash

# ==============================================================================
# Mitochondrial Quality Control (mt_qc) Pipeline Configuration
# This file centralizes paths, settings, and tools for the mt_qc workflow.
# ==============================================================================

# -----------------------------------------------------------------------------
# Environment Management Setup
# -----------------------------------------------------------------------------
# Supported options: "micromamba", "conda"
export ENV_MANAGER="micromamba" # *MODIFY*

# Initialize base environment based on the selected manager
if [ "$ENV_MANAGER" == "micromamba" ]; then
    export MAMBA_EXE="/work/opt/ohpc/pkg/rcec/pkg/spack/opt/spack/linux-centos7-cascadelake/gcc-11.3.0/micromamba-1.4.2-hxbsylg3ivhyv6puis4s3l5o2b4tpwtn/bin/micromamba"
    export MAMBA_ROOT_PREFIX="/home/u4432941/micromamba"
    eval "$("$MAMBA_EXE" shell hook --shell bash --prefix "$MAMBA_ROOT_PREFIX")"
elif [ "$ENV_MANAGER" == "conda" ]; then
    ml biology
    ml Anaconda/Anaconda3
fi

# Function to activate an environment
activate_env() {
    local env_name=$1
    if [ "$ENV_MANAGER" == "micromamba" ]; then
        micromamba activate "$env_name"
    elif [ "$ENV_MANAGER" == "conda" ]; then
        conda activate "$env_name" 
    fi
}

# Function to run a command inside an environment
run_in_env() {
    local env_name=$1
    shift # Remove env_name from arguments
    if [ "$ENV_MANAGER" == "micromamba" ]; then
        micromamba run -n "$env_name" "$@"
    elif [ "$ENV_MANAGER" == "conda" ]; then
        conda run -n "$env_name" "$@"
    fi
}

# ------------------------------------------------------------------------------
# Project Global Paths
# ------------------------------------------------------------------------------
# PROJECT_DIR is where the logs and outputs are stored
export PROJECT_DIR="/staging/reserve/jacobhsu/TWB/TWBR11106-05/Phenotypes_TX/TWB1492_mt/microarray/mt_qc/test" # *MODIFY*
export MT_QC_OUT_DIR="${PROJECT_DIR}/qc_${PREFIX}"
export LOG_DIR="${PROJECT_DIR}/logs"

# ------------------------------------------------------------------------------
# Batch Selection & Genotype Data
# ------------------------------------------------------------------------------
# twb1
# export INPUT_PLINK="/staging/reserve/jacobhsu/TWB/TWB/microarray/TWB1/27719/Genotyped/TWB1.hg19"
# export PREFIX="twb1"

# twb2
export INPUT_PLINK="/staging/reserve/jacobhsu/TWB/TWB/microarray/TWB2_genotyped_120163/TWB2.hg38.v4" # *MODIFY*
export PREFIX="twb2"

# ------------------------------------------------------------------------------
# Reference Files & Genotype Maps
# ------------------------------------------------------------------------------
export MT_FASTA="/staging/reserve/jacobhsu/TWB/TWBR11106-05/Phenotypes_TX/TWB1492_mt/chrM_ref/hg38_v0_chrM_Homo_sapiens_assembly38.chrM.fasta"
export MISC_DIR="/staging/reserve/jacobhsu/TWB/TWBR11106-05/Phenotypes_TX/TWB1492_mt/microarray/misc"
export GMAP="${MISC_DIR}/chrMT.NC012920.gmap"
export MTMAP="${MISC_DIR}/mtgeneticmap.map"

# ------------------------------------------------------------------------------
# Imputation & Validation Reference Panels
# ------------------------------------------------------------------------------
export REF_DIR="/staging/reserve/jacobhsu/TWB/TWBR11106-05/Phenotypes_TX/TWB1492_mt/microarray/mt_qc/validate_imputation/wgs_split"
export REF_PANEL_VALIDATE="${REF_DIR}/wgs1000.ref.vcf.gz"
export REF_PANEL_VCF="${REF_DIR}/recode_vcf/TWB1465_wgs_recode.vcf.gz"

# ------------------------------------------------------------------------------
# Autosomal QC Exclusion Lists (used in step 04)
# ------------------------------------------------------------------------------
export AUTOSOMAL_QC_DIR="/staging/reserve/jacobhsu/TWB/TWBR11106-05/Phenotypes_TX/TWB-SV/outputs" # *MODIFY*
export NON_EAS_INDV_LIST="${AUTOSOMAL_QC_DIR}/check_ibd/twb2_ref_nonpredpop.indvlist"
export RELATED_INDV_LIST="${AUTOSOMAL_QC_DIR}/check_ibd/twb2.ldpr.2degree.king.cutoff.out.id"

# ------------------------------------------------------------------------------
# Executable Tool Paths
# ------------------------------------------------------------------------------
export PLINK="/opt/ohpc/Taiwania3/pkg/biology/PLINK/PLINK_v1.90/plink"
export PLINK2="/opt/ohpc/Taiwania3/pkg/biology/PLINK2/PLINK_v2.00a2.3_AVX2/plink2"
export IMPUTE2="/work/opt/ohpc/Taiwania3/pkg/biology/IMPUTE2/impute_v2.3.2/impute2"
export SHAPEIT2="/work/opt/ohpc/Taiwania3/pkg/biology/SHAPEIT/shapeit_v2.904.glibcv2.17/bin/shapeit"

# ------------------------------------------------------------------------------
# Helper Script Paths
# ------------------------------------------------------------------------------
export REMOVE_DEVIATE_PY="${MISC_DIR}/remove_deviate.py"
