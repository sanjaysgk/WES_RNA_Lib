#!/bin/bash
# ================================================
# RNA-seq Tumor-Normal Mutation Calling Pipeline (Adapted)
# Author: Sanjay SG Krishna
# Li Lab / Purcell Lab, Monash University, Australia
# Adapted from Cryptic pipeline for tumor-normal analysis
# ================================================

#SBATCH -p comp
#SBATCH --job-name=RNA_TumorNormal_Pipeline
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G
#SBATCH --account= #youraccount [Note: replace with your actual account]
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yourmain@monash 
#SBATCH --time=1-12:00:00

set -euo pipefail

# --------------------------------------------
# CONFIGURATION
# --------------------------------------------
SAMPLE_ID="ORG125"
BASE_DIR="/home/$USER/xy86_scratch2/SANJAY/Organoid"
DATA_DIR="${BASE_DIR}/RNA"
OUTPUT_DIR="${BASE_DIR}/RNA/results/${SAMPLE_ID}"
TEMP_DIR="${OUTPUT_DIR}/temp"

# Create output directories
mkdir -p "$OUTPUT_DIR" "$TEMP_DIR" "$OUTPUT_DIR/Tumor" "$OUTPUT_DIR/Normal" "$OUTPUT_DIR/Variant_Calling"

# Tumor/Normal FASTQ files
TUMOR_R1="${DATA_DIR}/Tumor/${SAMPLE_ID}T-LCL6191_L4_1.fq.gz"
TUMOR_R2="${DATA_DIR}/Tumor/${SAMPLE_ID}T-LCL6191_L4_2.fq.gz"
NORMAL_R1="${DATA_DIR}/Normal/${SAMPLE_ID}N-LCL6190_L4_1.fq.gz"
NORMAL_R2="${DATA_DIR}/Normal/${SAMPLE_ID}N-LCL6190_L4_2.fq.gz"

# References
REF_GENOME="${BASE_DIR}/references/GRCh38.primary_assembly.genome.fa"
REF_GTF="${BASE_DIR}/references/gencode.v44.primary_assembly.annotation.gtf"
STAR_INDEX="${BASE_DIR}/references/STAR_index"

# Variant calling resources
KNOWN_SITES_DBSNP="${BASE_DIR}/variant_calling_resources/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
KNOWN_SITES_MILLS="${BASE_DIR}/variant_calling_resources/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
KNOWN_SITES_INDEL="${BASE_DIR}/variant_calling_resources/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz"

#Variant anotations
VEP_CACHE_DIR="${BASE_DIR}/tool/ensembl-vep/.vep"
VEP_VERSION=90
PROT_SEQ_FLANK_L=15

# Parameters
THREADS=14
LOG_FILE="$OUTPUT_DIR/pipeline_status.log"

# --------------------------------------------
# Helper Functions
# --------------------------------------------
check_file() {
    local file=$1 stage=$2
    if [[ ! -f "$file" || ! -s "$file" ]]; then
        echo "ERROR ($stage): File missing or empty → $file" | tee -a "$LOG_FILE"
        exit 1
    else
        echo "✓ [$stage] $(basename "$file") exists and size = $(du -h "$file" | cut -f1)" | tee -a "$LOG_FILE"
    fi
}

run_step() {
    local step_number=$1
    local step_description=$2
    local command=$3
    
    echo "[$(date)] Running Step $step_number: $step_description..." | tee -a "$LOG_FILE"
    eval "$command"
    
    if [ $? -ne 0 ]; then
        echo "Step $step_number failed: $step_description. Exiting." | tee -a "$LOG_FILE"
        exit 1
    else
        echo "Step $step_number completed successfully." | tee -a "$LOG_FILE"
    fi
}

# Initialize log file
touch "$LOG_FILE"
echo "Pipeline started at $(date)" | tee -a "$LOG_FILE"
echo "Sample ID: $SAMPLE_ID" | tee -a "$LOG_FILE"

# Validate input files
echo "Validating input files..." | tee -a "$LOG_FILE"
for file in "$TUMOR_R1" "$TUMOR_R2" "$NORMAL_R1" "$NORMAL_R2" "$REF_GENOME" "$REF_GTF"; do
    check_file "$file" "Input"
done

# --------------------------------------------
# TUMOR SAMPLE PROCESSING
# --------------------------------------------
echo "================================================" | tee -a "$LOG_FILE"
echo "PROCESSING TUMOR SAMPLE" | tee -a "$LOG_FILE"
echo "================================================" | tee -a "$LOG_FILE"

# Step 1: STAR Alignment for Tumor
run_step "1T" "STAR Alignment - Tumor" "
module purge
module load gcc/6.1.0
module load star/2.7.9a
cd '$OUTPUT_DIR/Tumor'
STAR --genomeDir '$STAR_INDEX' \
     --readFilesIn '$TUMOR_R1' '$TUMOR_R2' \
     --runThreadN $THREADS \
     --twopassMode Basic \
     --readFilesCommand zcat \
     --outFileNamePrefix '${SAMPLE_ID}_Tumor_' \
     --outSAMtype BAM SortedByCoordinate \
     --outReadsUnmapped Fastx \
     &> '${OUTPUT_DIR}/Tumor/star_align.log'
"

# Step 2: Infer Experiment - Tumor
run_step "2T" "Infer Experiment - Tumor" "
module purge
module load rseqc
cd '$OUTPUT_DIR/Tumor'
infer_experiment.py --version
infer_experiment.py \
    -i '$OUTPUT_DIR/Tumor/${SAMPLE_ID}_Tumor_Aligned.sortedByCoord.out.bam' \
    -r '$BASE_DIR/variant_calling_resources/gencode_assembly.bed' \
    &> '$OUTPUT_DIR/Tumor/infer_tumor.txt'
"

# Step 3: StringTie - Tumor
run_step "3T" "StringTie Assembly - Tumor" "
module purge
module load gcc/5.4.0
module load stringtie/2.1.5
cd '$OUTPUT_DIR/Tumor'
stringtie -p $THREADS --rf \
    -G '$REF_GTF' \
    -l 'COV_Tumor' \
    -o '${SAMPLE_ID}_Tumor_Transcripts.gtf' \
    '${SAMPLE_ID}_Tumor_Aligned.sortedByCoord.out.bam' \
    &> 'stringtie_tumor.log'
"

# Step 4: GFFCompare - Tumor
run_step "4T" "GFFCompare - Tumor" "
module purge
module load gffcompare
cd '$OUTPUT_DIR/Tumor'
gffcompare -r '$REF_GTF' \
    -o 'tumor_prefix' \
    -R -V -C \
    '${SAMPLE_ID}_Tumor_Transcripts.gtf' \
    &> 'gffcompare_tumor.log'
"

# Step 5: Convert BAM to SAM for tumor processing
run_step "5T" "Convert BAM to SAM - Tumor" "
module purge
module load samtools
cd '$OUTPUT_DIR/Tumor'
samtools view -h -o '${SAMPLE_ID}_Tumor_Aligned.out.sam' \
    '${SAMPLE_ID}_Tumor_Aligned.sortedByCoord.out.bam'
"

# Step 6: Sort SAM by Query Name - Tumor
run_step "6T" "Sort SAM by Query Name - Tumor" "
module purge
module load gatk/4.2.5.0
cd '$OUTPUT_DIR/Tumor'
gatk SortSam \
    -I '${SAMPLE_ID}_Tumor_Aligned.out.sam' \
    -O '${SAMPLE_ID}_Tumor_Aligned_sortedqn.bam' \
    -SO queryname \
    --TMP_DIR '$TEMP_DIR' \
    &> 'sortqn_tumor.log'
"

# Step 7: FastqToSam - Tumor
run_step "7T" "FastqToSam - Tumor" "
module purge
module load gatk/4.2.5.0
cd '$OUTPUT_DIR/Tumor'
gatk FastqToSam \
    -F1 '$TUMOR_R1' \
    -F2 '$TUMOR_R2' \
    -O '${SAMPLE_ID}_Tumor_unmapped.bam' \
    -SM '${SAMPLE_ID}_Tumor' \
    -RG '${SAMPLE_ID}_Tumor_RG1' \
    -PL ILLUMINA \
    &> 'fastqtosam_tumor.log'
"

# Step 8: Merge BAM Alignments - Tumor
run_step "8T" "Merge BAM Alignments - Tumor" "
cd '$OUTPUT_DIR/Tumor'
gatk MergeBamAlignment \
    -ALIGNED '${SAMPLE_ID}_Tumor_Aligned_sortedqn.bam' \
    -UNMAPPED '${SAMPLE_ID}_Tumor_unmapped.bam' \
    -R '$REF_GENOME' \
    -O '${SAMPLE_ID}_Tumor_merged_alignments.bam' \
    &> 'mergebam_tumor.log'
"

# Step 9: Mark Duplicates - Tumor
run_step "9T" "Mark Duplicates - Tumor" "
cd '$OUTPUT_DIR/Tumor'
gatk MarkDuplicates \
    -I '${SAMPLE_ID}_Tumor_merged_alignments.bam' \
    -O '${SAMPLE_ID}_Tumor_marked_duplicates.bam' \
    -M '${SAMPLE_ID}_Tumor_marked_dup_metrics.txt' \
    -R '$REF_GENOME' \
    --READ_NAME_REGEX null \
    &> 'markdup_tumor.log'
"

# Step 10: Split NCigar Reads - Tumor
run_step "10T" "Split NCigar Reads - Tumor" "
module purge
module load gatk/4.2.5.0
cd '$OUTPUT_DIR/Tumor'
gatk SplitNCigarReads \
    -I '${SAMPLE_ID}_Tumor_marked_duplicates.bam' \
    -O '${SAMPLE_ID}_Tumor_splitNcigar.bam' \
    -R '$REF_GENOME' \
    &> 'splitncigar_tumor.log'
"

# --------------------------------------------
# NORMAL SAMPLE PROCESSING
# --------------------------------------------
echo "================================================" | tee -a "$LOG_FILE"
echo "PROCESSING NORMAL SAMPLE" | tee -a "$LOG_FILE"
echo "================================================" | tee -a "$LOG_FILE"

# Repeat steps 1-10 for Normal sample
# Step 1N: STAR Alignment for Normal
run_step "1N" "STAR Alignment - Normal" "
module purge
module load gcc/6.1.0
module load star/2.7.9a
cd '$OUTPUT_DIR/Normal'
STAR --genomeDir '$STAR_INDEX' \
     --readFilesIn '$NORMAL_R1' '$NORMAL_R2' \
     --runThreadN $THREADS \
     --twopassMode Basic \
     --readFilesCommand zcat \
     --outFileNamePrefix '${SAMPLE_ID}_Normal_' \
     --outSAMtype BAM SortedByCoordinate \
     --outReadsUnmapped Fastx \
     &> '${OUTPUT_DIR}/Normal/star_align.log'
"

# Step 2N: Infer Experiment - Normal
run_step "2N" "Infer Experiment - Normal" "
module purge
module load rseqc
cd '$OUTPUT_DIR/Normal'
infer_experiment.py --version
infer_experiment.py \
    -i '$OUTPUT_DIR/Normal/${SAMPLE_ID}_Normal_Aligned.sortedByCoord.out.bam' \
    -r '$BASE_DIR/variant_calling_resources/gencode_assembly.bed' \
    &> '$OUTPUT_DIR/Normal/infer_normal.txt'
"

# Continue with remaining normal processing steps...
# Step 3N through 10N (similar to tumor but for normal)

# Steps 3N-10N for Normal (abbreviated for space)
run_step "3N-10N" "Complete Normal Processing Pipeline" "
# StringTie
module purge
module load gcc/5.4.0
module load stringtie/2.1.5
cd '$OUTPUT_DIR/Normal'
stringtie -p $THREADS --rf -G '$REF_GTF' -l 'COV_Normal' \
    -o '${SAMPLE_ID}_Normal_Transcripts.gtf' \
    '${SAMPLE_ID}_Normal_Aligned.sortedByCoord.out.bam' &> 'stringtie_normal.log'

# GFFCompare
module purge
module load gffcompare
gffcompare -r '$REF_GTF' -o 'normal_prefix' -R -V -C \
    '${SAMPLE_ID}_Normal_Transcripts.gtf' &> 'gffcompare_normal.log'

# Convert BAM to SAM
module purge
module load samtools
samtools view -h -o '${SAMPLE_ID}_Normal_Aligned.out.sam' \
    '${SAMPLE_ID}_Normal_Aligned.sortedByCoord.out.bam'

# Sort by query name
module purge
module load gatk/4.2.5.0
gatk SortSam -I '${SAMPLE_ID}_Normal_Aligned.out.sam' \
    -O '${SAMPLE_ID}_Normal_Aligned_sortedqn.bam' -SO queryname \
    --TMP_DIR '$TEMP_DIR' &> 'sortqn_normal.log'

# FastqToSam
gatk FastqToSam -F1 '$NORMAL_R1' -F2 '$NORMAL_R2' \
    -O '${SAMPLE_ID}_Normal_unmapped.bam' -SM '${SAMPLE_ID}_Normal' \
    -RG '${SAMPLE_ID}_Normal_RG1' -PL ILLUMINA &> 'fastqtosam_normal.log'

# Merge alignments
gatk MergeBamAlignment -ALIGNED '${SAMPLE_ID}_Normal_Aligned_sortedqn.bam' \
    -UNMAPPED '${SAMPLE_ID}_Normal_unmapped.bam' -R '$REF_GENOME' \
    -O '${SAMPLE_ID}_Normal_merged_alignments.bam' &> 'mergebam_normal.log'

# Mark duplicates
gatk MarkDuplicates -I '${SAMPLE_ID}_Normal_merged_alignments.bam' \
    -O '${SAMPLE_ID}_Normal_marked_duplicates.bam' \
    -M '${SAMPLE_ID}_Normal_marked_dup_metrics.txt' -R '$REF_GENOME' \
    --READ_NAME_REGEX null &> 'markdup_normal.log'

# Split NCigar reads
gatk SplitNCigarReads -I '${SAMPLE_ID}_Normal_marked_duplicates.bam' \
    -O '${SAMPLE_ID}_Normal_splitNcigar.bam' -R '$REF_GENOME' \
    &> 'splitncigar_normal.log'
"

# --------------------------------------------
# BASE QUALITY SCORE RECALIBRATION (BQSR)
# --------------------------------------------
echo "================================================" | tee -a "$LOG_FILE"
echo "BASE QUALITY SCORE RECALIBRATION" | tee -a "$LOG_FILE"
echo "================================================" | tee -a "$LOG_FILE"

# Step 11: Base Recalibrator - Tumor
run_step "11T" "Base Recalibrator - Tumor" "
module purge
module load gatk/4.2.5.0
cd '$OUTPUT_DIR/Tumor'
gatk BaseRecalibrator \
    -I '${SAMPLE_ID}_Tumor_splitNcigar.bam' \
    -O '${SAMPLE_ID}_Tumor_recal_data.table' \
    -R '$REF_GENOME' \
    --known-sites '$KNOWN_SITES_INDEL' \
    --known-sites '$KNOWN_SITES_DBSNP' \
    --known-sites '$KNOWN_SITES_MILLS' \
    &> 'baserecal_tumor.log'
"

# Step 12: Apply BQSR - Tumor
run_step "12T" "Apply BQSR - Tumor" "
gatk ApplyBQSR \
    -R '$REF_GENOME' \
    -I '${SAMPLE_ID}_Tumor_splitNcigar.bam' \
    --bqsr-recal-file '${SAMPLE_ID}_Tumor_recal_data.table' \
    -O '${SAMPLE_ID}_Tumor_score_corrected.bam' \
    &> 'applybqsr_tumor.log'
"

# Step 13: Base Recalibrator - Normal  
run_step "13N" "Base Recalibrator - Normal" "
cd '$OUTPUT_DIR/Normal'
gatk BaseRecalibrator \
    -I '${SAMPLE_ID}_Normal_splitNcigar.bam' \
    -O '${SAMPLE_ID}_Normal_recal_data.table' \
    -R '$REF_GENOME' \
    --known-sites '$KNOWN_SITES_INDEL' \
    --known-sites '$KNOWN_SITES_DBSNP' \
    --known-sites '$KNOWN_SITES_MILLS' \
    &> 'baserecal_normal.log'
"

# Step 14: Apply BQSR - Normal
run_step "14N" "Apply BQSR - Normal" "
gatk ApplyBQSR \
    -R '$REF_GENOME' \
    -I '${SAMPLE_ID}_Normal_splitNcigar.bam' \
    --bqsr-recal-file '${SAMPLE_ID}_Normal_recal_data.table' \
    -O '${SAMPLE_ID}_Normal_score_corrected.bam' \
    &> 'applybqsr_normal.log'
"
# --------------------------------------------
# Optional: Second round of BQSR and plotting
# --------------------------------------------
echo "================================================" | tee -a "$LOG_FILE"
echo "OPTIONAL: SECOND ROUND OF BQSR AND PLOTTING" | tee -a "$LOG_FILE"
echo "================================================" | tee -a "$LOG_FILE"
#### Optinal: Second round of BQSR and plotting for both Tumor and Normal
Step 13.A: Apply BQSR - Tumor
run_step "13aT" "Second Base Recalibration - Tumor" "
module purge
module load gatk/4.2.5.0
cd '$OUTPUT_DIR/Tumor'
gatk BaseRecalibrator \
    -I '${SAMPLE_ID}_Tumor_score_corrected.bam' \
    -O '${SAMPLE_ID}_Tumor_post_recal_data.table' \
    -R '$REF_GENOME' \
    --known-sites '$KNOWN_SITES_INDEL' \
    --known-sites '$KNOWN_SITES_DBSNP' \
    --known-sites '$KNOWN_SITES_MILLS' \
    &> 'post_bqsr_tumor.log'
"
Step 13.B: Analyse Covariates - Tumor
run_step "13bT" "Analyse Covariates - Tumor" "
module purge
module load R
module load gatk/4.2.5.0
cd '$OUTPUT_DIR/Tumor'
gatk AnalyzeCovariates \
    -before '${SAMPLE_ID}_Tumor_recal_data.table' \
    -after '${SAMPLE_ID}_Tumor_post_recal_data.table' \
    -plots '${SAMPLE_ID}_Tumor_bqsr_plots.pdf' \
    &> 'analyzecovariates_tumor.log'
"
####### ##############################################

##### Optional: Second round of BQSR and plotting for Normal

#Step 14.A: Apply Second BQSR - Normal
run_step "14aN" "Second Base Recalibration - Normal" "
module purge
module load gatk/4.2.5.0
cd '$OUTPUT_DIR/Normal'
gatk BaseRecalibrator \
    -I '${SAMPLE_ID}_Normal_score_corrected.bam' \
    -O '${SAMPLE_ID}_Normal_post_recal_data.table' \
    -R '$REF_GENOME' \
    --known-sites '$KNOWN_SITES_INDEL' \
    --known-sites '$KNOWN_SITES_DBSNP' \
    --known-sites '$KNOWN_SITES_MILLS' \
    &> 'post_bqsr_normal.log'
"
# # Step 14.B: Analyse Covariates - Normal
run_step "14bN" "Analyse Covariates - Normal" "
module purge
module load R
module load gatk/4.2.5.0
cd '$OUTPUT_DIR/Normal'
gatk AnalyzeCovariates \
    -before '${SAMPLE_ID}_Normal_recal_data.table' \
    -after '${SAMPLE_ID}_Normal_post_recal_data.table' \
    -plots '${SAMPLE_ID}_Normal_bqsr_plots.pdf' \
    &> 'analyzecovariates_normal.log'
"


# --------------------------------------------
# VARIANT CALLING WITH MUTECT2
# --------------------------------------------
echo "================================================" | tee -a "$LOG_FILE"
echo "VARIANT CALLING WITH MUTECT2" | tee -a "$LOG_FILE"
echo "================================================" | tee -a "$LOG_FILE"

# Step 15: Mutect2 Tumor-Normal Variant Calling
run_step "15" "Mutect2 Variant Calling" "
module purge
module load gatk/4.2.5.0
module load samtools
cd '$OUTPUT_DIR/Variant_Calling'

# Extract sample names from BAM headers
TUMOR_SM=\$(samtools view -H '$OUTPUT_DIR/Tumor/${SAMPLE_ID}_Tumor_score_corrected.bam' | awk '\$1~/^@RG/ {for(i=1;i<=NF;i++){if(\$i~/^SM:/){split(\$i,a,\":\");print a[2]}}}' | sort -u)
NORMAL_SM=\$(samtools view -H '$OUTPUT_DIR/Normal/${SAMPLE_ID}_Normal_score_corrected.bam' | awk '\$1~/^@RG/ {for(i=1;i<=NF;i++){if(\$i~/^SM:/){split(\$i,a,\":\");print a[2]}}}' | sort -u)

echo \"Detected Tumor sample: \$TUMOR_SM\" | tee -a '$LOG_FILE'
echo \"Detected Normal sample: \$NORMAL_SM\" | tee -a '$LOG_FILE'

gatk Mutect2 \
    -R '$REF_GENOME' \
    -I '$OUTPUT_DIR/Tumor/${SAMPLE_ID}_Tumor_score_corrected.bam' \
    -I '$OUTPUT_DIR/Normal/${SAMPLE_ID}_Normal_score_corrected.bam' \
    -tumor \"\$TUMOR_SM\" \
    -normal \"\$NORMAL_SM\" \
    -O '${SAMPLE_ID}_mutect2.vcf.gz' \
    --f1r2-tar-gz '${SAMPLE_ID}_f1r2.tar.gz' \
    &> '${SAMPLE_ID}_mutect2.log'
"

echo "================================================" | tee -a "$LOG_FILE"
echo "RNA-seq Tumor-Normal Pipeline Completed Successfully!" | tee -a "$LOG_FILE"
echo "Completed at: $(date)" | tee -a "$LOG_FILE"
echo "Output directory: $OUTPUT_DIR" | tee -a "$LOG_FILE"
echo "Variants output: $OUTPUT_DIR/Variant_Calling/${SAMPLE_ID}_mutect2.vcf.gz" | tee -a "$LOG_FILE"
echo "================================================" | tee -a "$LOG_FILE"


--------------------------------------------
Step 16: Learn Model

run_step "16" "Learn Read Orientation Model" "
cd '$OUTPUT_DIR/Variant_Calling'
module purge
module load gatk/4.2.5.0
gatk LearnReadOrientationModel \
    -I '${SAMPLE_ID}_f1r2.tar.gz' \
    -O '${SAMPLE_ID}_read-orientation-model.tar.gz' \
    &> 'learn_read_orientation.log'
"
Step 17: Get Pileup summarys
run_step "17" "Get Pileup Summaries - Tumor" "
cd '$OUTPUT_DIR/Variant_Calling'
module purge
module load gatk/4.2.5.0
gatk GetPileupSummaries \
    -I '$OUTPUT_DIR/Tumor/${SAMPLE_ID}_Tumor_score_corrected.bam' \
    -V '$BASE_DIR/variant_calling_resources/small_exac_common_3.hg38.vcf.gz' \
    -L '$BASE_DIR/variant_calling_resources/small_exac_common_3.hg38.vcf.gz' \
    -O '${SAMPLE_ID}_tumor_pileups.table' \
    &> 'get_pileup_tumor.log'
"
Step 18: Get Pileup summarys - Normal
run_step "18" "Get Pileup Summaries - Normal" "
module purge
module load gatk/4.2.5.0
cd '$OUTPUT_DIR/Variant_Calling'
gatk GetPileupSummaries \
    -I '$OUTPUT_DIR/Normal/${SAMPLE_ID}_Normal_score_corrected.bam' \
    -V '$BASE_DIR/variant_calling_resources/small_exac_common_3.hg38.vcf.gz' \
    -L '$BASE_DIR/variant_calling_resources/small_exac_common_3.hg38.vcf.gz' \
    -O '${SAMPLE_ID}_normal_pileups.table' \
    &> 'get_pileup_normal.log'
"
Step 19: Calculate Contamination
run_step "19" "Calculate Contamination" "
module purge
module load gatk/4.2.5.0
cd '$OUTPUT_DIR/Variant_Calling'
gatk CalculateContamination \
    -I '${SAMPLE_ID}_tumor_pileups.table' \
    -matched '${SAMPLE_ID}_normal_pileups.table' \
    -O '${SAMPLE_ID}_contamination.table' \
    --tumor-segmentation '${SAMPLE_ID}_segments.table' \
    &> 'calculate_contamination.log'
"
Step 20: Filter Mutect Calls
run_step "20" "Filter Mutect Calls" "
module purge
module load gatk/4.2.5.0
cd '$OUTPUT_DIR/Variant_Calling'
gatk FilterMutectCalls \
    -V '${SAMPLE_ID}_mutect2.vcf.gz' \
    -R '$REF_GENOME' \
    --contamination-table '${SAMPLE_ID}_contamination.table' \
    --tumor-segmentation '${SAMPLE_ID}_segments.table' \
    --ob-priors '${SAMPLE_ID}_read-orientation-model.tar.gz' \
    -O '${SAMPLE_ID}_mutect2_filtered.vcf.gz' \
    &> 'filter_mutect_calls.log'
"
# Step 21: Select Variants
run_step "21" "Select Variants" "
module purge
module load gatk/4.2.5.0
cd '$OUTPUT_DIR/Variant_Calling'
gatk SelectVariants \
    -V '${SAMPLE_ID}_mutect2_filtered.vcf.gz' \
    --exclude-filtered true \
    -O '${SAMPLE_ID}_mutect2_filtered_pass.vcf.gz' \
    &> 'select_variants.log'
"

#    --sites-only-vcf-output true \

--------------------------------------------
Functions for key steps in the pipeline
--------------------------------------------

Step 22: VEP Annotation

run_step "22" "VEP Annotation" "
module purge
echo \"Starting VEP version ${VEP_VERSION} annotation at \$(date)\" | tee -a \"$LOG_FILE\"

module load vep/${VEP_VERSION}
cd \"${OUTPUT_DIR}/Variant_Calling\"
export VEP_PLUGINS=\"${BASE_DIR}/tool/VEP_plugins\"
vep --cache --offline \
    --input_file \"${SAMPLE_ID}_mutect2_filtered_pass.vcf.gz\" \
    --fasta \"${REF_GENOME}\" \
    --dir_cache \"${VEP_CACHE_DIR}\" \
    --symbol --protein --hgvs --tsl \
    --vcf --force_overwrite \
    --output_file \"${SAMPLE_ID}_VEP_v${VEP_VERSION}_annotated.vcf\" \
    --plugin Frameshift \
    --plugin Wildtype \
    &> \"vep_annotation.log\"
"


echo "VEP annotation completed!"
echo "Output saved to: ${OUTPUT_VCF}"


## Step 23: Make FASTA Database from VCF
run_step "23" "Make FASTA Database from VCF" "
module purge
echo \"Starting FASTA generation from VEP-annotated VCF at \$(date)\" | tee -a \"$LOG_FILE\"

# Activate conda environment
source /home/$USER/xy86_scratch2/SANJAY/miniconda3/etc/profile.d/conda.sh
conda activate immunopep-tools

# Create output directory if it doesn't exist
mkdir -p \"${OUTPUT_DIR}/Variant_Calling/FASTA\"
cd \"${OUTPUT_DIR}/Variant_Calling/FASTA\"

# Run fasta_generator script
python \"${BASE_DIR}/tool/fasta_generator/generate_protein_fasta.py\" \
    \"${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_VEP_v${VEP_VERSION}_annotated.vcf\" \
    ${PROT_SEQ_FLANK_L} \
    \"${SAMPLE_ID}_Tumor_mutant_only_full_flank_${PROT_SEQ_FLANK_L}.fasta\" \
    --mutant-only \
    --pass-only \
    -s \"${SAMPLE_ID}_Tumor\" \
    -d full \
    &> \"${OUTPUT_DIR}/Variant_Calling/FASTA/${SAMPLE_ID}_fasta_generation.log\"

echo \"FASTA generation completed!\"
echo \"Output saved to: ${OUTPUT_DIR}/Variant_Calling/FASTA/${SAMPLE_ID}_Tumor_mutant_only_full_flank_${PROT_SEQ_FLANK_L}.fasta\"
"

