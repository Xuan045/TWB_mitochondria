#!/usr/bin/bash

wkdir="/staging/biology/u4432941/TWB1492_mt/test"
#Mitochondrial reference path
refdir="/staging/biology/u4432941/TWB1492_mt/chrM_ref"
para="NGS20140611E"
input="/staging/biology/u4432941/TWB1492_mt/TWB_WGS_MAP/${para}/${para}*.cram"
input_index="${input}.crai"
HAPLOCKECK="/staging/biology/u4432941/apps/haplocheckCLI"

REF="/staging/reserve/paylong_ntu/AI_SHARE/reference/GATK_bundle/2.8/hg38/Homo_sapiens_assembly38.fasta"
CHRM_REF="${refdir}/hg38_v0_chrM_Homo_sapiens_assembly38.chrM.fasta"
CHRM_SHIFTED_REF="${refdir}/hg38_v0_chrM_Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta"
CHRM_SHIFTEDBACK_CHAIN="${refdir}/hg38_v0_chrM_ShiftBack.chain"
CHRM_INTERVAL_LIST="${refdir}/hg38_v0_chrM_non_control_region.chrM.interval_list"
CHRM_SHIFTED_INTERVAL_LIST="${refdir}/hg38_v0_chrM_control_region_shifted.chrM.interval_list"
BLACK_LIST="${refdir}/hg38_v0_chrM_blacklist_sites.hg38.chrM.bed"

PICARD="/opt/ohpc/Taiwania3/pkg/biology/Picard/picard_v2.26.0/picard.jar"
GATK4="/opt/ohpc/Taiwania3/pkg/biology/GATK/gatk_v4.2.3.0"
BWA="/opt/ohpc/Taiwania3/pkg/biology/BWA/BWA_v0.7.17"
NUM_THREAD=14

# Other settings
mkdir -p ${wkdir}
mkdir -p ${wkdir}/${para}
cd ${wkdir}/${para}
temp="${wkdir}/${para}/temp"
mkdir -p ${temp}
haplocheck="${wkdir}/${para}/check_contamination"
mkdir -p ${haplocheck}

TIME=`date +%Y%m%d%H%M`
logfile=./${TIME}_${para}_runMT.log
exec 3<&1 4<&2
exec >$logfile 2>&1
set -e

# *******************************************************
# 1. Print Reads from ChrM   
# *******************************************************
${GATK4}/gatk PrintReads \
	-R ${REF} \
	-L chrM \
	--read-filter MateOnSameContigOrNoMappedMateReadFilter \
	--read-filter MateUnmappedAndUnmappedReadFilter \
	-I ${input} \
	--read-index ${input_index} \
	-O ${para}.subset.bam

# *******************************************************
# 2. Revert Sam To unmapped BAM (uBAM)  
# *******************************************************
java -Xmx48g -jar ${PICARD} RevertSam \
  	INPUT=${para}.subset.bam \
	OUTPUT_BY_READGROUP=false \
	OUTPUT=${temp}/${para}.unmapped.bam \
	VALIDATION_STRINGENCY=LENIENT \
	ATTRIBUTE_TO_CLEAR=FT \
	ATTRIBUTE_TO_CLEAR=CO \
	SORT_ORDER=queryname \
	RESTORE_ORIGINAL_QUALITIES=false

# *******************************************************************
# 3. Align control and non-control region separately using BWA-MEM   
# *******************************************************************
### Get BWA version
bwa_version=`${BWA}/bwa 2>&1 | grep -e '^Version' | sed 's/Version: //'`
bwa_commandline="bwa mem -K 100000000 -p -v 3 -t ${NUM_THREAD} -Y ${CHRM_REF}"
output_basename=${para}.noncontrol

### Non-control region ***********************************************************************
## Convert uBAM to Fastq, then perform alignment and combine with uBAM
java -Xms48g -jar ${PICARD} SamToFastq \
  		INPUT=${temp}/${para}.unmapped.bam \
  		FASTQ=${temp}/${output_basename}.fastq \
  		INTERLEAVE=true \
  		NON_PF=true
${BWA}/bwa mem -K 100000000 -p -v 3 -t ${NUM_THREAD} -Y ${CHRM_REF} ${temp}/${output_basename}.fastq > ${temp}/${output_basename}.bwa.aligned.bam
java -Xms48g -jar ${PICARD}	MergeBamAlignment \
		VALIDATION_STRINGENCY=SILENT \
		EXPECTED_ORIENTATIONS=FR \
		ATTRIBUTES_TO_RETAIN=X0 \
		ATTRIBUTES_TO_REMOVE=NM \
		ATTRIBUTES_TO_REMOVE=MD \
		ALIGNED_BAM=${temp}/${output_basename}.bwa.aligned.bam \
		UNMAPPED_BAM=${temp}/${para}.unmapped.bam \
		OUTPUT=${temp}/${output_basename}.mba.bam \
		REFERENCE_SEQUENCE=${CHRM_REF} \
		PAIRED_RUN=true \
		SORT_ORDER="unsorted" \
		IS_BISULFITE_SEQUENCE=false \
		ALIGNED_READS_ONLY=false \
		CLIP_ADAPTERS=false \
		MAX_RECORDS_IN_RAM=2000000 \
		ADD_MATE_CIGAR=true \
		MAX_INSERTIONS_OR_DELETIONS=-1 \
		PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
		PROGRAM_RECORD_ID="bwamem" \
		PROGRAM_GROUP_VERSION=${bwa_version} \
		PROGRAM_GROUP_COMMAND_LINE="${bwa_commandline}" \
		PROGRAM_GROUP_NAME="bwamem" \
		UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
		ALIGNER_PROPER_PAIR_FLAGS=true \
		UNMAP_CONTAMINANT_READS=true \
		ADD_PG_TAG_TO_READS=false

## Mark duplicates
java -Xms48g -jar ${PICARD} MarkDuplicates \
		INPUT=${temp}/${output_basename}.mba.bam \
		OUTPUT=${temp}/${output_basename}.md.bam \
		METRICS_FILE=${temp}/${output_basename}.metrics \
		VALIDATION_STRINGENCY=SILENT \
		OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
		ASSUME_SORT_ORDER="queryname" \
		CLEAR_DT="false" \
		ADD_PG_TAG_TO_READS=false

## Sorting
java -Xms48g -jar ${PICARD} SortSam \
		INPUT=${temp}/${output_basename}.md.bam \
		OUTPUT=${temp}/${output_basename}.realigned.sorted.bam \
		SORT_ORDER="coordinate" \
		CREATE_INDEX=true \
		MAX_RECORDS_IN_RAM=300000

### Control region ***********************************************************************
shifted_output_basename=${para}.control
shifted_bwa_commandline="bwa mem -K 100000000 -p -v 3 -t ${NUM_THREAD} -Y ${CHRM_SHIFTED_REF}"

## Perform alignment using ChrM reference shifted by 8000 bp, and combine with uBAM
${BWA}/bwa mem -K 100000000 -p -v 3 -t ${NUM_THREAD} -Y ${CHRM_SHIFTED_REF} ${temp}/${output_basename}.fastq > ${temp}/${shifted_output_basename}.bwa.aligned.bam
java -Xms48g -jar ${PICARD}	MergeBamAlignment \
		VALIDATION_STRINGENCY=SILENT \
		EXPECTED_ORIENTATIONS=FR \
		ATTRIBUTES_TO_RETAIN=X0 \
		ATTRIBUTES_TO_REMOVE=NM \
		ATTRIBUTES_TO_REMOVE=MD \
		ALIGNED_BAM=${temp}/${shifted_output_basename}.bwa.aligned.bam \
		UNMAPPED_BAM=${temp}/${para}.unmapped.bam \
		OUTPUT=${temp}/${shifted_output_basename}.mba.bam \
		REFERENCE_SEQUENCE=${CHRM_SHIFTED_REF} \
		PAIRED_RUN=true \
		SORT_ORDER="unsorted" \
		IS_BISULFITE_SEQUENCE=false \
		ALIGNED_READS_ONLY=false \
		CLIP_ADAPTERS=false \
		MAX_RECORDS_IN_RAM=2000000 \
		ADD_MATE_CIGAR=true \
		MAX_INSERTIONS_OR_DELETIONS=-1 \
		PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
		PROGRAM_RECORD_ID="bwamem" \
		PROGRAM_GROUP_VERSION=${bwa_version} \
		PROGRAM_GROUP_COMMAND_LINE="${shifted_bwa_commandline}" \
		PROGRAM_GROUP_NAME="bwamem" \
		UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
		ALIGNER_PROPER_PAIR_FLAGS=true \
		UNMAP_CONTAMINANT_READS=true \
		ADD_PG_TAG_TO_READS=false

## Mark duplicates
java -Xms48g -jar ${PICARD} MarkDuplicates \
		INPUT=${temp}/${shifted_output_basename}.mba.bam \
		OUTPUT=${temp}/${shifted_output_basename}.md.bam \
		METRICS_FILE=${temp}/${shifted_output_basename}.metrics \
		VALIDATION_STRINGENCY=SILENT \
		OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
		ASSUME_SORT_ORDER="queryname" \
		CLEAR_DT="false" \
		ADD_PG_TAG_TO_READS=false

## Sorting
java -Xms48g -jar ${PICARD} SortSam \
		INPUT=${temp}/${shifted_output_basename}.md.bam \
		OUTPUT=${temp}/${shifted_output_basename}.realigned.sorted.bam \
		SORT_ORDER="coordinate" \
		CREATE_INDEX=true \
		MAX_RECORDS_IN_RAM=300000

# *********************************************************
# 4. Call variants using GATK Mutect2 (mitochondria mode)   
# *********************************************************
### Noncontrol region
${GATK4}/gatk Mutect2 \
        -R ${CHRM_REF} \
        -I ${temp}/${output_basename}.realigned.sorted.bam \
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --read-filter MateUnmappedAndUnmappedReadFilter \
        -O ${temp}/${output_basename}.raw.vcf \
        -L chrM:576-16024  \
        --annotation StrandBiasBySample \
        --mitochondria-mode \
        --max-reads-per-alignment-start 75 \
        --max-mnp-distance 0

### Control region
${GATK4}/gatk Mutect2 \
        -R ${CHRM_SHIFTED_REF} \
        -I ${temp}/${shifted_output_basename}.realigned.sorted.bam \
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --read-filter MateUnmappedAndUnmappedReadFilter \
        -O ${temp}/${shifted_output_basename}.raw.vcf \
        -L chrM:8025-9144 \
        --annotation StrandBiasBySample \
        --mitochondria-mode \
        --max-reads-per-alignment-start 75 \
        --max-mnp-distance 0

# **********************************
# 5. Lift over and combine raw vcfs   
# **********************************
# Lifts over shifted vcf of control region and combines it with the rest of the chrM calls.
### Lift over shifted vcf
java -jar ${PICARD} LiftoverVcf \
		I=${temp}/${shifted_output_basename}.raw.vcf \
		O=${temp}/${para}.raw.shifted_back.vcf \
		R=${CHRM_REF} \
		CHAIN=${CHRM_SHIFTEDBACK_CHAIN} \
		REJECT=${para}.raw.rejected.vcf

### Combine vcfs
java -jar ${PICARD} MergeVcfs \
		I=${temp}/${para}.raw.shifted_back.vcf \
		I=${temp}/${output_basename}.raw.vcf \
		O=${para}.raw.merged.vcf

# **************************************
# 6. Merge stats files Mutect2 produced   
# **************************************
${GATK4}/gatk MergeMutectStats \
	--stats ${temp}/${shifted_output_basename}.raw.vcf.stats \
	--stats ${temp}/${output_basename}.raw.vcf.stats \
	-O ${para}.raw.combined.stats

# ***********************
# 7. Initial Filter   
# ***********************
${GATK4}/gatk FilterMutectCalls \
	-V ${para}.raw.merged.vcf \
  	-R ${CHRM_REF} \
  	-O ${temp}/${para}.initial_filtered.vcf \
  	--stats ${para}.raw.combined.stats \
  	--max-alt-allele-count 4 \
  	--mitochondria-mode \
  	--min-allele-fraction 0.0 \
  	--contamination-estimate 0.0

${GATK4}/gatk VariantFiltration \
	-V ${temp}/${para}.initial_filtered.vcf \
  	-O ${temp}/${para}.initial_filtered.sorted.vcf \
  	--apply-allele-specific-filters \
  	--mask ${BLACK_LIST} \
  	--mask-name "blacklisted_site"

# Split multi-allelics and remove non-pass sites   
${GATK4}/gatk LeftAlignAndTrimVariants \
  	-R ${CHRM_REF} \
  	-V ${temp}/${para}.initial_filtered.sorted.vcf \
	-O ${temp}/${para}.initial_filtered.sorted.split.vcf \
	--split-multi-allelics \
	--dont-trim-alleles \
	--keep-original-ac

# splitAndPassOnly.vcf file for Haplocheck
${GATK4}/gatk SelectVariants \
    -V ${temp}/${para}.initial_filtered.sorted.split.vcf \
    -O ${haplocheck}/${para}.initial_filtered.splitAndPassOnly.vcf \
    --exclude-filtered

# ***************************
# 8. Get contamination   
# ***************************
# Use Haplocheck to estimate levels of contamination in mitochondria
java -jar ${HAPLOCKECK}/haplocheckCLI.jar ${haplocheck}

sed 's/\"//g' output > ${haplocheck}/${para}.output_noquotes

grep "SampleID" ${haplocheck}/${para}.output_noquotes > headers
FORMAT_ERROR="Bad contamination file format"
if [ `awk '{print $2}' headers` != "Contamination" ]; then
  echo ${FORMAT_ERROR}; exit 1
fi
if [ `awk '{print $6}' headers` != "HgMajor" ]; then
  echo ${FORMAT_ERROR}; exit 1
fi
if [ `awk '{print $8}' headers` != "HgMinor" ]; then
  echo ${FORMAT_ERROR}; exit 1
fi
if [ `awk '{print $14}' headers` != "MeanHetLevelMajor" ]; then
  echo ${FORMAT_ERROR}; exit 1
fi
if [ `awk '{print $15}' headers` != "MeanHetLevelMinor" ]; then
  echo ${FORMAT_ERROR}; exit 1
fi

grep -v "SampleID" ${haplocheck}/${para}.output_noquotes > ${haplocheck}/output-data
awk '{print $2}' ${haplocheck}/output-data > ${haplocheck}/contamination.txt
awk '{print $6}' ${haplocheck}/output-data > ${haplocheck}/major_haplogroup.txt
awk '{print $8}' ${haplocheck}/output-data > ${haplocheck}/minor_haplogroup.txt
awk '{print $14}' ${haplocheck}/output-data > ${haplocheck}/mean_het_major.txt
awk '{print $15}' ${haplocheck}/output-data > ${haplocheck}/mean_het_minor.txt

rm output* headers

# ***************************
# 9. Filter contamination   
# ***************************
# Check if there's contamintion
if [ `awk '{print $2}' ${haplocheck}/output-data` == "YES" ]; then
	contamination=`awk '{print $14}' ${haplocheck}/output-data`; else
	contamination=0.0
fi

${GATK4}/gatk FilterMutectCalls \
    -V ${temp}/${para}.initial_filtered.sorted.vcf \
    -R ${CHRM_REF} \
    -O ${temp}/${para}.contamination_filtered.vcf \
    --stats ${para}.raw.combined.stats \
    --max-alt-allele-count 4 \
    --mitochondria-mode \
    --contamination-estimate ${contamination}

${GATK4}/gatk VariantFiltration \
    -V ${temp}/${para}.contamination_filtered.vcf \
    -O ${temp}/${para}.contamination_filtered.sorted.vcf \
    --apply-allele-specific-filters \
    --mask ${BLACK_LIST} \
    --mask-name "blacklisted_site"

# **********************************
# 10. Coverage At EveryBase   
# **********************************
# *****Remove this when there's a GVCF solution*****
## Non-control region
java -jar ${PICARD} CollectHsMetrics \
    I=${temp}/${output_basename}.realigned.sorted.bam \
    R=${CHRM_REF} \
    PER_BASE_COVERAGE=non_control_region.tsv \
    O=non_control_region.metrics \
    TI=${CHRM_INTERVAL_LIST} \
    BI=${CHRM_INTERVAL_LIST} \
    COVMAX=20000 \
    SAMPLE_SIZE=1

## Control region
java -jar ${PICARD} CollectHsMetrics \
    I=${temp}/${shifted_output_basename}.realigned.sorted.bam \
    R=${CHRM_SHIFTED_REF} \
    PER_BASE_COVERAGE=control_region_shifted.tsv \
    O=control_region_shifted.metrics \
    TI=${CHRM_SHIFTED_INTERVAL_LIST} \
    BI=${CHRM_SHIFTED_INTERVAL_LIST} \
    COVMAX=20000 \
    SAMPLE_SIZE=1

R --vanilla <<CODE
    shift_back = function(x) {
    if (x < 8570) {
        return(x + 8000)
    } else {
        return (x - 8569)
    }
    }

    control_region_shifted = read.table("control_region_shifted.tsv", header=T)
    shifted_back = sapply(control_region_shifted[,"pos"], shift_back)
    control_region_shifted[,"pos"] = shifted_back

    beginning = subset(control_region_shifted, control_region_shifted[,'pos']<8000)
    end = subset(control_region_shifted, control_region_shifted[,'pos']>8000)

    non_control_region = read.table("non_control_region.tsv", header=T)
    combined_table = rbind(beginning, non_control_region, end)
    write.table(combined_table, "per_base_coverage.tsv", row.names=F, col.names=T, quote=F, sep="\t")

CODE

# Remove temporary files
rm non_control_region* control_region_shifted*

# **********************************
# 11. Split Multi-allelic Sites   
# **********************************
${GATK4}/gatk LeftAlignAndTrimVariants \
    -R ${CHRM_REF} \
    -V ${temp}/${para}.contamination_filtered.sorted.vcf \
    -O ${para}.sorted.final.split.vcf \
    --split-multi-allelics \
    --dont-trim-alleles \
    --keep-original-ac

# ***********************
# 12. Remove temporary files   
# ***********************
rm -r ${temp}
