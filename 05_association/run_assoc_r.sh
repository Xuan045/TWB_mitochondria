#!/usr/bin/sh

para="twb2_imputed_maf.01_info0.7_prob0_EAS_king_covBatch_AgeSex_effect_2"
linear_para="${para}.rmOut.rint"
logistic_para="${para}.case1000"
outdir="/staging/biology/u4432941/TWB1492_mt/microarray_association/${para}/"
linear_script="/staging/biology/u4432941/TWB1492_mt/microarray_association/linear_glm.R"
logistic_script="/staging/biology/u4432941/TWB1492_mt/microarray_association/logistic_glm.R"
geno_file="/staging/biology/u4432941/TWB1492_mt/microarray_association/geno_data/impute_for_assoc_info0.7_prob0_af0.01_EAS_kingMT_effect_2.txt"
logistic_pheno_file="/staging/biology/u4432941/TWB1492_mt/microarray_association/pheno_data/twb2_survey.final.txt"
linear_pheno_file="/staging/biology/u4432941/TWB1492_mt/microarray_association/pheno_data/lab_test_rmOut_rint.txt"
pca_file="/staging/biology/u4432941/TWB1492_mt/microarray/pca_twb2/twb2_pca.eigenvec"

# set -euo pipefail
mkdir -p $outdir
module load pkg/Anaconda3
module load compiler/gcc/9.4.0

# Log file settings
TIME=`date +%Y%m%d%H%M`
logfile=${outdir}/${TIME}_run_assoc.log

# Redirect standard output and error to the log file
exec > "$logfile" 2>&1

# Echo used files into the log file
echo "Output Directory: $outdir"
echo "Genotype File: $geno_file"
echo "Logistic Phenotype File: $logistic_pheno_file"
echo "Linear Phenotype File: $linear_pheno_file"
echo "PCA File: $pca_file"
echo "-------"


# python $py_file $input_file $output_file $outdir/exclude_var.txt
/work/opt/ohpc/Taiwania3/pkg/biology/R/R_v4.1.0/bin/Rscript $linear_script $geno_file $linear_pheno_file $pca_file $outdir ${linear_para}
/work/opt/ohpc/Taiwania3/pkg/biology/R/R_v4.1.0/bin/Rscript $logistic_script $geno_file $logistic_pheno_file $pca_file $outdir $logistic_para 1000
