#!/bin/bash
# ================================================
# RNA-seq Tumor-Normal Mutation Calling Pipeline (Adapted)
# Author: Sanjay SG Krishna
# Li Lab / Purcell Lab, Monash University, Australia
# Adapted from Cryptic pipeline for tumor-normal analysis
# ================================================

#SBATCH -p comp
#SBATCH --job-name=RNA_Tumor_Mutect_Pipeline
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --account=xy86
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sanjay.sondekoppagopalakrishna@monash.edu
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

#for IPG
IPG_programs="/home/$USER/xy86_scratch/sanjay/ATLANTIS/RNAseq/Analysis/Cryptic/IPG_slurm"


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

####Mutation Calling with Mutect2
echo "================================================" | tee -a "$LOG_FILE"
echo "MUTATION CALLING WITH MUTECT2" | tee -a "$LOG_FILE"
echo "================================================" | tee -a "$LOG_FILE"
# Step 15: Mutect2 Variant Calling - Tumor only
run_step "15" "Mutect2 Variant Calling - Tumor" "
module purge
module load gatk/4.2.5.0
cd '$OUTPUT_DIR/Tumor'
gatk Mutect2 \
    -R '$REF_GENOME' \
    -I '${SAMPLE_ID}_Tumor_score_corrected.bam' \
    -O '$OUTPUT_DIR/Variant_Calling/${SAMPLE_ID}_Tumor_mutect.vcf.gz' \
    --f1r2-tar-gz '$OUTPUT_DIR/Variant_Calling/${SAMPLE_ID}_Tumor_f1r2.tar.gz' \
    &> 'mutect2_tumor.log'
"
# Step 16: learn Model - Tumor
run_step "16" "Learn Model - Tumor" "
module purge
module load gatk/4.2.5.0
cd '$OUTPUT_DIR/Tumor'
gatk LearnReadOrientationModel \
    -I '$OUTPUT_DIR/Variant_Calling/${SAMPLE_ID}_Tumor_f1r2.tar.gz' \
    -O '$OUTPUT_DIR/Variant_Calling/${SAMPLE_ID}_Tumor_read-orientation-model.tar.gz' \
    &> 'learnmodel_tumor.log'
"   

Step 17: Get the pileup summaries for tumor
run_step "17" "Get Pileup Summaries - Tumor" "
module purge
module load gatk/4.2.5.0
cd '$OUTPUT_DIR/Tumor'
gatk GetPileupSummaries \
    -I '${SAMPLE_ID}_Tumor_score_corrected.bam' \
    -V '${BASE_DIR}/variant_calling_resources/small_exac_common_3.hg38.vcf.gz' \
    -L '${BASE_DIR}/variant_calling_resources/small_exac_common_3.hg38.vcf.gz' \
    -O '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_pileups.table' \
    &> 'getpileup_tumor.log'
"
Step 18: Calculate contamination - Tumor
run_step "18" "Calculate Contamination - Tumor" "
module purge
module load gatk/4.2.5.0
cd '${OUTPUT_DIR}/Tumor'
gatk CalculateContamination \
    -I '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_pileups.table' \
    -O '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_contamination.table' \
    --tumor-segmentation '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_segments.table' \
    &> 'calculatecontamination_tumor.log'
"
Step 19: Filter Mutect2 Calls - Tumor
run_step "19" "Filter Mutect2 Calls - Tumor" "
module purge
module load gatk/4.2.5.0
cd '${OUTPUT_DIR}/Tumor'
gatk FilterMutectCalls \
    -V '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_mutect.vcf.gz' \
    -R '$REF_GENOME' \
    --contamination-table '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_contamination.table' \
    --ob-priors '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_read-orientation-model.tar.gz' \
    -O '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_mutect_filtered.vcf.gz' \
    &> 'filtermutect_tumor.log'
"
Step 20: Select Variants - Tumor
run_step "20" "Select Variants - Tumor" "
module purge
module load gatk/4.2.5.0
cd '${OUTPUT_DIR}/Tumor'
gatk SelectVariants \
    -V '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_mutect_filtered.vcf.gz' \
    --exclude-filtered true \
    --sites-only-vcf-output true \
    -O '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_mutect_filtered_pass.vcf.gz' \
    &> 'selectvariants_tumor.log'
"


Step 21: Curate Variants with IPG
gzip -d '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_mutect_filtered_pass.vcf.gz'

run_step "21" "Curate Variants with IPG" "
module purge
cd '${OUTPUT_DIR}/Tumor'
\"$IPG_programs/curate_vcf.o\" \
    '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_mutect_filtered_pass.vcf' \
    -d \
    &> 'curate_variants_tumor.log'
"


Step 22: Make Index for Curated VCF
run_step "22" "Make Index for Curated VCF" "
module purge
module load gatk/4.2.5.0
cd '${OUTPUT_DIR}/Tumor'

# Index unmasked VCF
gatk IndexFeatureFile -I "${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_mutect_filtered_pass_unmasked.vcf" &> index_curated_vcf_tumor_unmasked.log

# Index indel VCF
gatk IndexFeatureFile -I "${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_mutect_filtered_pass_indel.vcf" &> index_curated_vcf_tumor_indel.log

"

Step 23: Make alternative Genome with Curated Variants
run_step "23" "Make alternative Genome with Curated Variants" "
module purge
module load gatk/4.2.5.0
cd '${OUTPUT_DIR}/Tumor'
gatk FastaAlternateReferenceMaker \
    -R '$REF_GENOME' \
    -V '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_mutect_filtered_pass_unmasked.vcf' \
    -O '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_unmasked_alternate_genome.fa' \
    &> 'fasta_unmasked_alternate_genome_tumor.log'

gatk FastaAlternateReferenceMaker \
    -R '$REF_GENOME' \
    -V '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_mutect_filtered_pass_indel.vcf' \
    -O '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_indel_alternate_genome.fa' \
    &> 'fasta_indel_alternate_genome_tumor.log'    
"

#Step 24: Revert Header of Curated VCF
run_step "24" "Revert Header of Curated VCF" "
module purge
cd '${OUTPUT_DIR}/Tumor'
\"$IPG_programs/revert_headers.o\" \
    '$REF_GENOME' \
    '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_unmasked_alternate_genome.fa' \
    &> 'revert_header_unmasked_tumor.log'

mv tmpc.fasta tmpc_unmasked.fasta #renaming for clarity 

\"$IPG_programs/revert_headers.o\" \
    '$REF_GENOME' \
    '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_indel_alternate_genome.fa' \
    &> 'revert_header_indel_tumor.log'
mv tmpc.fasta tmpc_indel.fasta #renaming for clarity
"

Step 25: Sort Assembly gtf
run_step "25" "Sort Assembly gtf" "
module purge
cd '${OUTPUT_DIR}/Tumor'
"$IPG_programs/gff3sort-master/gff3sort.pl" \
    --chr_order original \
    'tumor_prefix.combined.gtf' > '${SAMPLE_ID}_Tumor_Transcripts_sorted.gtf'
"

Step 26: LiftOver
run_step "26" "LiftOver" "
module purge
cd '${OUTPUT_DIR}/Tumor'

"$IPG_programs/alt_liftover.o" \
-v "$OUTPUT_DIR/Variant_Calling/${SAMPLE_ID}_Tumor_mutect_filtered_pass_unmasked.vcf" \
-g '${SAMPLE_ID}_Tumor_Transcripts_sorted.gtf' \
-s _alt_liftover_unmasked \
 &> 'liftover_unmasked_tumor.log'

"$IPG_programs/alt_liftover.o" \
-v "$OUTPUT_DIR/Variant_Calling/${SAMPLE_ID}_Tumor_mutect_filtered_pass_indel.vcf" \
-g '${SAMPLE_ID}_Tumor_Transcripts_sorted.gtf' \
-s _alt_liftover_indel \
 &> 'liftover_indel_tumor.log'
"

Step 27: gffread for lifted over VCF
run_step "27" "gffread for lifted over VCF" "
cd \"${OUTPUT_DIR}/Tumor\"
module purge

\"${IPG_programs}/gffread\" \
  -F -W \
  -w \"${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_unmasked_alternate_genome.fa\" \
  -g \"${OUTPUT_DIR}/Tumor/tmpc_unmasked.fasta\" \
  \"${SAMPLE_ID}_Tumor_Transcripts_sorted_alt_liftover_unmasked.gtf\" \
  &> gffread_unmasked_tumor.log

\"${IPG_programs}/gffread\" \
  -F -W \
  -w \"${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_indel_alternate_genome.fa\" \
  -g \"${OUTPUT_DIR}/Tumor/tmpc_indel.fasta\" \
  \"${SAMPLE_ID}_Tumor_Transcripts_sorted_alt_liftover_indel.gtf\" \
  &> gffread_indel_tumor.log

\"${IPG_programs}/gffread\" \
  -F -W \
  -w \"${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_ref_transcriptome.fa\" \
  -g '$REF_GENOME' \
  '${SAMPLE_ID}_Tumor_Transcripts_sorted.gtf' \
  &> gffread_reference_tumor.log
"

Step 28: Triple Translate
run_step "28" "Triple Translate" "
module purge
cd \"${OUTPUT_DIR}/Tumor\"

\"${IPG_programs}/triple_translate.o\" \
  -c tumor_prefix.tracking \
  \"${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_unmasked_alternate_genome.fa\" \
  &> triple_translate_unmasked_tumor.log

\"${IPG_programs}/triple_translate.o\" \
  -c tumor_prefix.tracking \
  \"${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_indel_alternate_genome.fa\" \
  &> triple_translate_indel_tumor.log

\"${IPG_programs}/triple_translate.o\" \
  -c tumor_prefix.tracking \
  \"${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_ref_transcriptome.fa\" \
  &> triple_translate_reference_tumor.log
"

Step 29: Squish
run_step "29" "Squish" "
module purge
cd \"${OUTPUT_DIR}/Tumor\"

\"${IPG_programs}/squish.o\" \
  -d \"${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_ref_transcriptome_3translate.fasta\" \
  -d \"${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_unmasked_alternate_genome_3translate.fasta\" \
  -d \"${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_indel_alternate_genome_3translate.fasta\" \
  -o \"${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Tumor_Cryptic_DB.fasta\" \
  &> squish_tumor.log
"



###### Normal processing with Mutect2 can be added similarly if needed
# --------------------------------------------

# Step 15: Mutect2 Variant Calling - Normal only
run_step "15N" "Mutect2 Variant Calling - Normal" "
module purge
module load gatk/4.2.5.0
cd '${OUTPUT_DIR}/Normal'
gatk Mutect2 \
    -R '$REF_GENOME' \
    -I '${SAMPLE_ID}_Normal_score_corrected.bam' \
    -O '$OUTPUT_DIR/Variant_Calling/${SAMPLE_ID}_Normal_mutect.vcf.gz' \
    --f1r2-tar-gz '$OUTPUT_DIR/Variant_Calling/${SAMPLE_ID}_Normal_f1r2.tar.gz' \
    &> 'mutect2_normal.log'
"

# Step 16: Learn Model - Normal
run_step "16N" "Learn Model - Normal" "
module purge
module load gatk/4.2.5.0
cd '$OUTPUT_DIR/Normal'
gatk LearnReadOrientationModel \
    -I '$OUTPUT_DIR/Variant_Calling/${SAMPLE_ID}_Normal_f1r2.tar.gz' \
    -O '$OUTPUT_DIR/Variant_Calling/${SAMPLE_ID}_Normal_read-orientation-model.tar.gz' \
    &> 'learnmodel_normal.log'
"

# Step 17: Get the pileup summaries for Normal
run_step "17N" "Get Pileup Summaries - Normal" "
module purge
module load gatk/4.2.5.0
cd '$OUTPUT_DIR/Normal'
gatk GetPileupSummaries \
    -I '${SAMPLE_ID}_Normal_score_corrected.bam' \
    -V '${BASE_DIR}/variant_calling_resources/small_exac_common_3.hg38.vcf.gz' \
    -L '${BASE_DIR}/variant_calling_resources/small_exac_common_3.hg38.vcf.gz' \
    -O '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_pileups.table' \
    &> 'getpileup_normal.log'
"

# Step 18: Calculate contamination - Normal
run_step "18N" "Calculate Contamination - Normal" "
module purge
module load gatk/4.2.5.0
cd '${OUTPUT_DIR}/Normal'
gatk CalculateContamination \
    -I '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_pileups.table' \
    -O '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_contamination.table' \
    --tumor-segmentation '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_segments.table' \
    &> 'calculatecontamination_normal.log'
"

# Step 19: Filter Mutect2 Calls - Normal
run_step "19N" "Filter Mutect2 Calls - Normal" "
module purge
module load gatk/4.2.5.0
cd '${OUTPUT_DIR}/Normal'
gatk FilterMutectCalls \
    -V '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_mutect.vcf.gz' \
    -R '$REF_GENOME' \
    --contamination-table '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_contamination.table' \
    --ob-priors '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_read-orientation-model.tar.gz' \
    -O '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_mutect_filtered.vcf.gz' \
    &> 'filtermutect_normal.log'
"

# Step 20: Select Variants - Normal
run_step "20N" "Select Variants - Normal" "
module purge
module load gatk/4.2.5.0
cd '${OUTPUT_DIR}/Normal'
gatk SelectVariants \
    -V '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_mutect_filtered.vcf.gz' \
    --exclude-filtered true \
    --sites-only-vcf-output true \
    -O '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_mutect_filtered_pass.vcf.gz' \
    &> 'selectvariants_normal.log'
"

# Step 21: Curate Variants with IPG - Normal

gzip -d "${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_mutect_filtered_pass.vcf.gz"

run_step "21N" "Curate Variants with IPG - Normal" "
module purge
cd '${OUTPUT_DIR}/Normal'
\"$IPG_programs/curate_vcf.o\" \
    '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_mutect_filtered_pass.vcf' \
    -d \
    &> 'curate_variants_normal_unmasked.log'

cd '${OUTPUT_DIR}/Normal'
\"$IPG_programs/curate_vcf.o\" \
    '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_mutect_filtered_pass.vcf' \
    &> 'curate_variants_normal_indel.log'
"

# Step 22: Make Index for Curated VCF - Normal
run_step "22N" "Make Index for Curated VCF - Normal" "
module purge
module load gatk/4.2.5.0
cd '${OUTPUT_DIR}/Normal'

gatk IndexFeatureFile -I '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_mutect_filtered_pass_unmasked.vcf' &> index_curated_vcf_normal_unmasked.log

gatk IndexFeatureFile -I '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_mutect_filtered_pass_indel.vcf' &> index_curated_vcf_normal_indel.log
"

# Step 23: Make alternative Genome with Curated Variants - Normal
run_step "23N" "Make alternative Genome with Curated Variants - Normal" "
module purge
module load gatk/4.2.5.0
cd '${OUTPUT_DIR}/Normal'
gatk FastaAlternateReferenceMaker \
    -R '$REF_GENOME' \
    -V '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_mutect_filtered_pass_unmasked.vcf' \
    -O '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_unmasked_alternate_genome.fa' \
    &> 'fasta_unmasked_alternate_genome_normal.log'

gatk FastaAlternateReferenceMaker \
    -R '$REF_GENOME' \
    -V '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_mutect_filtered_pass_indel.vcf' \
    -O '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_indel_alternate_genome.fa' \
    &> 'fasta_indel_alternate_genome_normal.log'
"

Step 24: Revert Header of Curated VCF - Normal
run_step "24N" "Revert Header of Curated VCF - Normal" "
module purge
cd '${OUTPUT_DIR}/Normal'
\"$IPG_programs/revert_headers.o\" \
    '$REF_GENOME' \
    '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_unmasked_alternate_genome.fa' \
    &> 'revert_header_unmasked_normal.log'
mv tmpc.fasta tmpc_unmasked_normal.fasta

\"$IPG_programs/revert_headers.o\" \
    '$REF_GENOME' \
    '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_indel_alternate_genome.fa' \
    &> 'revert_header_indel_normal.log'
mv tmpc.fasta tmpc_indel_normal.fasta
"

# Step 25: Sort Assembly gtf - Normal
run_step "25N" "Sort Assembly gtf - Normal" "
module purge
cd '${OUTPUT_DIR}/Normal'
\"$IPG_programs/gff3sort-master/gff3sort.pl\" \
    --chr_order original \
    'normal_prefix.combined.gtf' > '${SAMPLE_ID}_Normal_Transcripts_sorted.gtf'
"

# Step 26: LiftOver - Normal
run_step "26N" "LiftOver - Normal" "
module purge
cd '${OUTPUT_DIR}/Normal'
\"$IPG_programs/alt_liftover.o\" \
  -v '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_mutect_filtered_pass_unmasked.vcf' \
  -g '${SAMPLE_ID}_Normal_Transcripts_sorted.gtf' \
  -s _alt_liftover_unmasked &> 'liftover_unmasked_normal.log'

\"$IPG_programs/alt_liftover.o\" \
  -v '${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_mutect_filtered_pass_indel.vcf' \
  -g '${SAMPLE_ID}_Normal_Transcripts_sorted.gtf' \
  -s _alt_liftover_indel &> 'liftover_indel_normal.log'
"

Step 27: gffread for lifted over VCF - Normal
run_step "27N" "gffread for lifted over VCF - Normal" "
cd '${OUTPUT_DIR}/Normal'
module purge

\"${IPG_programs}/gffread\" -F -W \
  -w \"${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_unmasked_alternate_genome.fa\" \
  -g \"${OUTPUT_DIR}/Normal/tmpc_unmasked_normal.fasta\" \
  \"${SAMPLE_ID}_Normal_Transcripts_sorted_alt_liftover_unmasked.gtf\" \
  &> gffread_unmasked_normal.log

\"${IPG_programs}/gffread\" -F -W \
  -w \"${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_indel_alternate_genome.fa\" \
  -g \"${OUTPUT_DIR}/Normal/tmpc_indel_normal.fasta\" \
  \"${SAMPLE_ID}_Normal_Transcripts_sorted_alt_liftover_indel.gtf\" \
  &> gffread_indel_normal.log

\"${IPG_programs}/gffread\" -F -W \
  -w \"${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_ref_transcriptome.fa\" \
  -g '$REF_GENOME' \
  '${SAMPLE_ID}_Normal_Transcripts_sorted.gtf' \
  &> gffread_reference_normal.log
"

Step 28: Triple Translate - Normal
run_step "28N" "Triple Translate - Normal" "
module purge
cd \"${OUTPUT_DIR}/Normal\"

\"${IPG_programs}/triple_translate.o\" \
  -c normal_prefix.tracking \
  \"${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_unmasked_alternate_genome.fa\" \
  &> triple_translate_unmasked_normal.log

\"${IPG_programs}/triple_translate.o\" \
  -c normal_prefix.tracking \
  \"${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_indel_alternate_genome.fa\" \
  &> triple_translate_indel_normal.log

\"${IPG_programs}/triple_translate.o\" \
  -c normal_prefix.tracking \
  \"${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_ref_transcriptome.fa\" \
  &> triple_translate_reference_normal.log
"

# Step 29: Squish - Normal
run_step "29N" "Squish - Normal" "
module purge
cd \"${OUTPUT_DIR}/Normal\"

\"${IPG_programs}/squish.o\" \
  -d \"${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_ref_transcriptome_3translate.fasta\" \
  -d \"${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_unmasked_alternate_genome_3translate.fasta\" \
  -d \"${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_indel_alternate_genome_3translate.fasta\" \
  -o \"${OUTPUT_DIR}/Variant_Calling/${SAMPLE_ID}_Normal_Cryptic_DB.fasta\" \
  &> squish_normal.log
"
