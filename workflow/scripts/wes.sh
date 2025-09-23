#!/bin/bash
# ================================================
# Whole Exome Sequencing (WES) Mutation Calling Pipeline
# Author: Sanjay SG Krishna
# Li Lab / Purcell Lab, Monash University, Australia
# ================================================

#SBATCH -p comp
#SBATCH --job-name=WES_Pipeline_Normal_Tumor
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=19G   # 16 × 19G = 304G total
#SBATCH --account=youraccount  #youraccount [Note: replace with your actual account]
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yourmail
#SBATCH --time=2-04:00:00

#fail at any error
set -euo pipefail
# --------------------------------------------
# Load required modules
# --------------------------------------------
module load samtools
module load gatk
module load fastqc
# module load bwa
module load bwa-mem2
# --------------------------------------------
# CONFIGURATION
# --------------------------------------------
SAMPLE_ID="ORG125"
BASE_DIR="/home/$USER/xy86_scratch2/SANJAY/Organoid"
DATA_DIR="${BASE_DIR}/WES/data"
OUTPUT_DIR="${BASE_DIR}/WES/results/${SAMPLE_ID}"

# Input FASTQ files
#Tumor
TUMOR_R1="${DATA_DIR}/Tumor/${SAMPLE_ID}T_Data_L1_1.fq-002.gz"
TUMOR_R2="${DATA_DIR}/Tumor/${SAMPLE_ID}T_Data_L1_2.fq-001.gz"
#Normal
NORMAL_R1="${DATA_DIR}/Normal/${SAMPLE_ID}N_Data_L1_1.fq-002.gz"
NORMAL_R2="${DATA_DIR}/Normal/${SAMPLE_ID}N_Data_L1_2.fq-001.gz"

# --------------------------------------------
# Download reference files if not already present
# (Uncomment the following lines if you need to download them)
#wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
#wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi
# wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
# wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi
# --------------------------------------------

## Reference files
REF_GENOME="${BASE_DIR}/references/GRCh38.primary_assembly.genome.fa"
KNOWN_SITES_DBSNP="${BASE_DIR}/variant_calling_resources/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
KNOWN_SITES_MILLS="${BASE_DIR}/variant_calling_resources/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
KNOWN_SITES_indel="${BASE_DIR}/variant_calling_resources/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz"
KNOWN_SITES_small_exon="${BASE_DIR}/variant_calling_resources/small_exac_common_3.hg38.vcf.gz"
KNOWN_SITES_PON="${BASE_DIR}/variant_calling_resources/1000g_pon.hg38.vcf.gz"
KNOWN_SITES_Germline="${BASE_DIR}/variant_calling_resources/af-only-gnomad.hg38.vcf.gz"

THREADS=14
TMEMORY="19g"
MEMORY="300g"

# --------------------------------------------
# SETUP
# --------------------------------------------
mkdir -p "${OUTPUT_DIR}"/{logs,temp,qc,variants}
LOG_FILE="${OUTPUT_DIR}/logs/pipeline.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=================================================="
echo "WES Mutation Detection Pipeline"
echo "Author: Sanjay SG Krishna"
echo "Co-Author:Chen Li"
echo "Lab: Li Lab / Purcell Lab, Monash University, Australia"
echo "WES Pipeline Started: $(date)"
echo "Sample: $SAMPLE_ID"
echo "=================================================="

# --------------------------------------------
# Helper function: check file exists and >0 size
# --------------------------------------------
check_file() {
    local file=$1
    local stage=$2
    if [[ ! -f "$file" ]]; then
        echo "ERROR ($stage): File not found → $file"
        exit 1
    elif [[ ! -s "$file" ]]; then
        echo "ERROR ($stage): File is empty → $file"
        exit 1
    else
        echo "✓ [$stage] $(basename "$file") exists and size = $(du -h "$file" | cut -f1)"
    fi
}

# --------------------------------------------
# Validate input FASTQ files
# --------------------------------------------
echo "Validating input FASTQ files..."
for file in "$TUMOR_R1" "$TUMOR_R2" "$NORMAL_R1" "$NORMAL_R2"; do
    check_file "$file" "Input FASTQ"
done

# --------------------------------------------
# Validate reference genome
# --------------------------------------------
check_file "$REF_GENOME" "Reference Genome"

# Index reference genome if needed
if [[ ! -f "${REF_GENOME}.fai" ]]; then
    echo "Indexing reference genome with samtools..."
    samtools faidx "$REF_GENOME"
    check_file "${REF_GENOME}.fai" "Reference Index"
fi

if [[ ! -f "${REF_GENOME%.fa}.dict" ]]; then
    echo "Creating sequence dictionary with GATK..."
    gatk CreateSequenceDictionary -R "$REF_GENOME"
    check_file "${REF_GENOME%.fa}.dict" "Reference Dict"
fi

if [[ ! -f "${REF_GENOME}.bwt.2bit.64" ]]; then
    echo "Creating BWA index..."
    # bwa index "$REF_GENOME"
    bwa-mem2 index "$REF_GENOME"
    check_file "${REF_GENOME}.bwt.2bit.64" "bwa-mem2 Index"
fi
echo "✓ Reference genome ready"
# --------------------------------------------
# Step1-10: Function to run alignment + BQSR pipeline
# --------------------------------------------
run_sample_pipeline() {
    local sample_type=$1
    local fastq_r1=$2
    local fastq_r2=$3

    local sample_name="${SAMPLE_ID}_${sample_type}"
    local aligned_bam="${OUTPUT_DIR}/temp/${sample_name}_aligned.bam"
    local sorted_bam="${OUTPUT_DIR}/temp/${sample_name}_sorted.bam"
    local dedup_bam="${OUTPUT_DIR}/temp/${sample_name}_dedup.bam"
    local recal_table="${OUTPUT_DIR}/temp/${sample_name}_recal_data.table"
    local post_recal_table="${OUTPUT_DIR}/temp/${sample_name}_post_recal_data.table"
    local bqsr_plots="${OUTPUT_DIR}/qc/${sample_name}_BQSR_plots.pdf"

    echo "[$(date)] Aligning ${sample_type} reads..."
    bwa-mem2 mem -t $THREADS -M \
        -R "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA\tLB:${sample_name}_lib1\tPU:${sample_name}_lane1" \
        "$REF_GENOME" "$fastq_r1" "$fastq_r2" | \
        samtools view -@ $THREADS -bS -o "$aligned_bam"

    check_file "$aligned_bam" "Alignment Output"

    echo "[$(date)] Sorting BAM..."
    gatk --java-options "-Xmx${MEMORY}" SortSam \
        -I "$aligned_bam" \
        -O "$sorted_bam" \
        -SO coordinate \
        --CREATE_INDEX true
    check_file "$sorted_bam" "Sorted BAM"
    # check_file "${sorted_bam}" "Sorted BAM File"

    echo "[$(date)] Marking duplicates..."
    gatk --java-options "-Xmx${MEMORY}" MarkDuplicates \
        -I "$sorted_bam" \
        -O "$dedup_bam" \
        -M "${OUTPUT_DIR}/qc/${sample_name}_dedup_metrics.txt" \
        --CREATE_INDEX true
    check_file "$dedup_bam" "Dedup BAM"
    # check_file "${dedup_bam}" "Dedup BAM File"

    echo "[$(date)] Collecting alignment metrics..."
    gatk --java-options "-Xmx${MEMORY}" CollectAlignmentSummaryMetrics \
        -R "$REF_GENOME" \
        -I "$dedup_bam" \
        -O "${OUTPUT_DIR}/qc/${sample_name}_alignment_metrics.txt"
    check_file "${OUTPUT_DIR}/qc/${sample_name}_alignment_metrics.txt" "Alignment Metrics"

    echo "[$(date)] Base Quality Score Recalibration (BQSR)..."
    gatk --java-options "-Xmx${MEMORY}" BaseRecalibrator \
        -R "$REF_GENOME" \
        -I "$dedup_bam" \
        --known-sites "$KNOWN_SITES_DBSNP" \
        --known-sites "$KNOWN_SITES_MILLS" \
        --known-sites "$KNOWN_SITES_indel" \
        -O "$recal_table"
    # --known-sites "$KNOWN_SITES_small_exon" \

    check_file "$recal_table" "BQSR Table"

    gatk --java-options "-Xmx${MEMORY}" ApplyBQSR \
        -R "$REF_GENOME" \
        -I "$dedup_bam" \
        --bqsr-recal-file "$recal_table" \
        -O "${OUTPUT_DIR}/temp/${sample_name}_score_recalibrated.bam"

    check_file "${OUTPUT_DIR}/temp/${sample_name}_score_recalibrated.bam" "Recalibrated BAM"

    gatk --java-options "-Xmx${MEMORY}" BaseRecalibrator \
        -R "$REF_GENOME" \
        -I "${OUTPUT_DIR}/temp/${sample_name}_score_recalibrated.bam" \
        --known-sites "$KNOWN_SITES_DBSNP" \
        --known-sites "$KNOWN_SITES_MILLS" \
        --known-sites "$KNOWN_SITES_indel" \
        -O "$post_recal_table"
        # --known-sites "$KNOWN_SITES_small_exon" \

    check_file "$post_recal_table" "Post-BQSR Table"

    gatk --java-options "-Xmx${MEMORY}" AnalyzeCovariates \
        -before "$recal_table" \
        -after "$post_recal_table" \
        -plots "$bqsr_plots"
    check_file "$bqsr_plots" "BQSR Plots"

    echo "[$(date)] ${sample_type} sample pipeline completed."
}
# Tumor sample
# run_sample_pipeline "Tumor" "$TUMOR_R1" "$TUMOR_R2"

# Normal sample
# run_sample_pipeline "Normal" "$NORMAL_R1" "$NORMAL_R2"


echo "=================================================="
echo "WES Pipeline Completed Successfully: $(date)"
echo "Next Steps: BQSR, Mutect2 Variant Calling"
echo "=================================================="


# --------------------------------------------
# Function: Mutect2 variant calling with paired Tumor-Normal
# --------------------------------------------
run_mutect2_pipeline() {
    local sample_name=$1
    local tumor_bam=$2
    local normal_bam=$3

    local variant_dir="${OUTPUT_DIR}/variants"
    mkdir -p "$variant_dir"

    # local ref_genome="$Cryptic_dir/GRCh38/GRCh38.primary_assembly.genome.fa"
    # local small_exac_vcf="$Cryptic_dir/variant_calling_resources/small_exac_common_3.hg38.vcf.gz"
    #find read groups
    # extract sample name(s) (SM tag) from tumor BAM read groups

    local BASE_TUMOR=$(samtools view -H "$tumor_bam" | awk '$1~/^@RG/ {for (i=1;i<=NF;i++) { if ($i~/^SM:/) { split($i,aa,":"); print aa[2] } }}' | sort -u)
    local BASE_NORMAL=$(samtools view -H "$normal_bam" | awk '$1~/^@RG/ {for (i=1;i<=NF;i++) { if ($i~/^SM:/) { split($i,aa,":"); print aa[2] } }}' | sort -u)
    echo "Detected tumor sample(s) in RG: ${BASE_TUMOR}"
    echo "Detected normal sample(s) in RG: ${BASE_NORMAL}"
    

    # Step 17: Mutect2
    echo "[$(date)] Running Mutect2 for $sample_name..."
    gatk Mutect2 \
        -R "$REF_GENOME" \
        -I "$tumor_bam" \
        -I "$normal_bam" \
        -tumor "$BASE_TUMOR" \
        -normal "$BASE_NORMAL" \
        -O "${variant_dir}/${sample_name}_mutect.vcf.gz" \
        --f1r2-tar-gz "${variant_dir}/${sample_name}_f1r2.tar.gz" \
        &> "${variant_dir}/${sample_name}_Mutect2.log"
    
    check_file "${variant_dir}/${sample_name}_mutect.vcf.gz" "Mutect2 VCF"


    # echo "[$(date)] Running Mutect2 with panel-of-normals PON and germline-resource for $sample_name..."
    # gatk Mutect2 \
    #     -R "$REF_GENOME" \
    #     -I "$tumor_bam" \
    #     -I "$normal_bam" \
    #     -tumor "$BASE_TUMOR" \
    #     -normal "$BASE_NORMAL" \
    #     -O "${variant_dir}/${sample_name}_mutect_pon_germ.vcf.gz" \
    #     --panel-of-normals "$KNOWN_SITES_PON" \
    #     --germline-resource "$KNOWN_SITES_Germline" \
    #     --f1r2-tar-gz "${variant_dir}/${sample_name}_pon_germ_f1r2.tar.gz" \
    #     &> "${variant_dir}/${sample_name}_Mutect2.log"
    
    # check_file "${variant_dir}/${sample_name}_mutect.vcf.gz" "Mutect2 VCF"

    Step 18: Learn Read Orientation Model
    echo "[$(date)] Learning read orientation model..."
    gatk LearnReadOrientationModel \
        -I "${variant_dir}/${sample_name}_f1r2.tar.gz" \
        -O "${variant_dir}/${sample_name}_read-orientation-model.tar.gz"

    Step 19: Get Pileup Summaries
    echo "[$(date)] Getting pileup summaries..."
    gatk GetPileupSummaries \
        -I "$tumor_bam" \
        -V "$KNOWN_SITES_small_exon" \
        -L "$KNOWN_SITES_small_exon" \
        -O "${variant_dir}/${sample_name}_tumor_getpileupsummaries.table" \
        &> "${variant_dir}/${sample_name}_GPS.log"

    gatk GetPileupSummaries \
        -I "$normal_bam" \
        -V "$KNOWN_SITES_small_exon" \
        -L "$KNOWN_SITES_small_exon" \
        -O "${variant_dir}/${sample_name}_normal_getpileupsummaries.table" \
        &> "${variant_dir}/${sample_name}_normal_GPS.log"

    Step 20: Calculate Contamination
    echo "[$(date)] Calculating contamination..."
    gatk CalculateContamination \
        -I "${variant_dir}/${sample_name}_tumor_getpileupsummaries.table" \
        -matched "${variant_dir}/${sample_name}_normal_getpileupsummaries.table" \
        --tumor-segmentation "${variant_dir}/${sample_name}_segments.table" \
        -O "${variant_dir}/${sample_name}_calculatecontamination.table" \
        &> "${variant_dir}/${sample_name}_CC.log"

    # Step 21: Filter Mutect Calls
    echo "[$(date)] Filtering Mutect2 calls..."
    gatk FilterMutectCalls \
        -V "${variant_dir}/${sample_name}_mutect.vcf.gz" \
        --tumor-segmentation "${variant_dir}/${sample_name}_segments.table" \
        --contamination-table "${variant_dir}/${sample_name}_calculatecontamination.table" \
        --stats "${variant_dir}/${sample_name}_mutect.vcf.gz.stats" \
        -ob-priors "${variant_dir}/${sample_name}_read-orientation-model.tar.gz" \
        -O "${variant_dir}/${sample_name}_mutect_COV_filtered.vcf" \
        -R "$REF_GENOME" \
        &> "${variant_dir}/${sample_name}_FMC.log"

    # Step 22: Select Variants
    echo "[$(date)] Selecting high-confidence variants..."
    gatk SelectVariants \
        -V "${variant_dir}/${sample_name}_mutect_COV_filtered.vcf" \
        --exclude-filtered true \
        --sites-only-vcf-output true \
        -O "${variant_dir}/${sample_name}_mutect_COV_pass.vcf" \
        &> "${variant_dir}/${sample_name}_SV.log"

    echo "[$(date)] Mutect2 pipeline completed for $sample_name."
}
#Tumor and normal BAMs have been BQSR processed
TUMOR_BAM="${OUTPUT_DIR}/temp/${SAMPLE_ID}_Tumor_score_recalibrated.bam"
NORMAL_BAM="${OUTPUT_DIR}/temp/${SAMPLE_ID}_Normal_score_recalibrated.bam"

# run_mutect2_pipeline "$SAMPLE_ID" "$TUMOR_BAM" "$NORMAL_BAM"



#Making Database for ANNOVAR
# Step 23: Annotate Variants with ANNOVAR
# echo "[$(date)] Annotating variants with ANNOVAR..."
variant_dir="${OUTPUT_DIR}/variants"
sample_name="${SAMPLE_ID}_mutect"

## Index gatk IndexFeatureFile

# gatk IndexFeatureFile \
#     -I "${variant_dir}/${sample_name}_COV_pass.vcf"


#referance make 

# gatk FastaAlternateReferenceMaker \
#     -R "$REF_GENOME" \
#     -V "${variant_dir}/${sample_name}_COV_pass.vcf" \
#     -O "${variant_dir}/${sample_name}_COV_pass.fa"
# https://console.cloud.google.com/storage/browser/broad-public-datasets/funcotator/funcotator_dataSources.v1.8.hg38.20230908s?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))
# broad-public-datasets/funcotator/funcotator_dataSources.v1.8.hg38.20230908s


# gatk --java-options "-Xms380G -Xmx380G -XX:ParallelGCThreads=18" Funcotator \
#     -R  "$REF_GENOME" \
#     -V "${variant_dir}/${sample_name}_COV_pass.vcf" \
#     -O "${variant_dir}/${sample_name}_COV_pass_Funcotator.maf" \
#     --output-file-format MAF \
#     --ref-version hg38 \
#     --data-sources-path 


echo "Running ANNOVAR_DB"

# # # Input VCF
INPUT_VCF="${variant_dir}/${sample_name}_COV_pass.vcf"

# # ANNOVAR database path (adjust if different on your system)
ANNOVAR_DB="/home/sson0030/xy86_scratch2/SANJAY/Organoid/tool/annovar/humandb"

# Run ANNOVAR
perl "$BASE_DIR/tool/annovar/table_annovar.pl" \
    "$INPUT_VCF" \
    "$ANNOVAR_DB" \
    -buildver hg38 \
    -out "${variant_dir}/${sample_name}_annovar" \
    -remove \
    -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a \
    -operation g,r,f,f,f \
    -nastring . \
    --csvout \
    --polish
    &> "${variant_dir}/${sample_name}_ANNOVAR.log"
# ,cosmic70,clinvar_20210501
echo "[$(date)] ANNOVAR annotation completed for $sample_name."
## Annoavr 
echo "ANNOVAR VCFs with coding changes:"
# perl "$BASE_DIR/tool/annovar/coding_change.pl"\


# #download 
# https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_merged_vep_115_GRCh38.tar.gz
# ## Run ensembl VEP annotation
echo "[$(date)] Annotating variants with Ensembl VEP..."
# #activate conda environment
# source /home/$USER/xy86_scratch2/SANJAY/miniconda3/bin/activate vep
module load vep
#check if vep is installed
if ! command -v vep &> /dev/null
then
    echo "ERROR: VEP could not be found. Please ensure VEP is installed and in your PATH."
    exit 1
fi

# Set variables
# VEP="/path/to/vep"                            # Path to VEP executable
# INPUT_VCF="${variant_dir}/${sample_name}_mutect_COV_pass.vcf"
# OUTPUT_VCF="${variant_dir}/${sample_name}_mutec_VEP.vcf"                   # Desired output file
# CACHE_DIR="/home/sson0030/xy86_scratch2/SANJAY/Organoid/tool/ensembl-vep/.vep"                # Directory where you unpacked the cache
# ASSEMBLY="GRCh38"                              # Genome build
# SPECIES="homo_sapiens"                         # Species
# CACHE_TYPE="merged"                             # Use merged cache (Ensembl + RefSeq)

# vep \
#     --input_file $INPUT_VCF \
#     --output_file $OUTPUT_VCF \
#     --everything \
#     -cache \
#     --dir_cache /home/sson0030/xy86_scratch2/SANJAY/Organoid/tool/ensembl-vep/.vep/ \
#     --cache_version 90 \
#     --offline \
#     --force_overwrite \
#     --fork 8 \
#     --vcf



# # # Run VEP
# vep \
#   --input_file $INPUT_VCF \
#   --output_file $OUTPUT_VCF \
#   --vcf \
#   --assembly $ASSEMBLY \
#   --dir_cache $CACHE_DIR \
#   --merged \
#   --species $SPECIES \
#   --cache \
#   --offline \
#   --force_overwrite \
#   --fork 8 \
#   --everything


# vep --cache --offline \
#     --input_file ${INPUT_VCF} \
#     --fasta ${REFERENCE_FA} \
#     --dir_cache ${VEP_CACHE_DIR} \
#     --symbol --protein --hgvs \
#     --vcf --force_overwrite \
#     --output_file VEP_annotated.vcf
