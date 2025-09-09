## Workflow overview

This workflow is a best-practice workflow for `<detailed description>`.
The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

1. Download genome reference from NCBI
2. Validate downloaded genome (`python` script)
3. Simulate short read sequencing data on the fly (`dwgsim`)
4. Check quality of input read data (`FastQC`)
5. Collect statistics from tool output (`MultiQC`)

## Running the workflow

### Input data

This template workflow creates artificial sequencing data in `*.fastq.gz` format.
It does not contain actual input data.
The simulated input files are nevertheless created based on a mandatory table linked in the `config.yaml` file (default: `.test/samples.tsv`).
The sample sheet has the following layout:

| sample  | condition | replicate | read1                      | read2                      |
| ------- | --------- | --------- | -------------------------- | -------------------------- |
| sample1 | wild_type | 1         | sample1.bwa.read1.fastq.gz | sample1.bwa.read2.fastq.gz |
| sample2 | wild_type | 2         | sample2.bwa.read1.fastq.gz | sample2.bwa.read2.fastq.gz |

### Parameters

This table lists all parameters that can be used to run the workflow.

| parameter          | type | details                               | default                        |
| ------------------ | ---- | ------------------------------------- | ------------------------------ |
| **samplesheet**    |      |                                       |                                |
| path               | str  | path to samplesheet, mandatory        | "config/samples.tsv"           |
| **get_genome**     |      |                                       |                                |
| ncbi_ftp           | str  | link to a genome on NCBI's FTP server | link to _S. cerevisiae_ genome |
| **simulate_reads** |      |                                       |                                |
| read_length        | num  | length of target reads in bp          | 100                            |
| read_number        | num  | number of total reads to be simulated | 10000                          |
Got it 👍 — we’ll adapt this README/workflow overview to **your WES Mutation Calling Pipeline** instead of the toy demo (genome download, dwgsim, etc.).

Here’s a rewritten version that reflects your pipeline steps (Tumor/Normal WES, alignment, QC, variant calling, etc.):

---

# WES Mutation Calling Pipeline

This workflow is a **whole exome sequencing (WES) best-practice pipeline** for detecting somatic variants from paired tumor/normal samples.
The workflow is built using [Snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

1. **Prepare references**

   * Validate and index the reference genome (FASTA, dict, fai, BWA index)
   * Index known variant databases (dbSNP, etc.)

2. **Quality control of raw data**

   * Run FastQC on raw FASTQ reads
   * Summarize results with MultiQC

3. **Alignment and BAM processing**

   * Align tumor and normal samples with BWA-MEM2
   * Sort, index, and mark duplicates with GATK
   * Collect alignment QC metrics

4. **Base quality score recalibration (BQSR)**

   * Recalibrate BAMs using known sites
   * Apply recalibration

5. **Somatic variant calling**

   * Run GATK Mutect2 on paired tumor/normal
   * Filter calls and create final VCF

6. **Variant QC & reporting**

   * Generate summary statistics
   * Produce MultiQC report of all steps

---

## Running the workflow

### Input data

The workflow expects **paired-end FASTQ files** for both **tumor** and **normal** samples.
The input samples are described in a mandatory table defined in `config/samples.tsv`.

Example layout:

| sample\_id | type   | read1                      | read2                      |
| ---------- | ------ | -------------------------- | -------------------------- |
| ORG125     | tumor  | ORG125T\_Data\_L1\_1.fq.gz | ORG125T\_Data\_L1\_2.fq.gz |
| ORG125     | normal | ORG125N\_Data\_L1\_1.fq.gz | ORG125N\_Data\_L1\_2.fq.gz |

### Parameters

This table lists all configurable parameters.

| parameter        | type | details                         | default                        |
| ---------------- | ---- | ------------------------------- | ------------------------------ |
| **samplesheet**  | str  | path to sample sheet, mandatory | `config/samples.tsv`           |
| **ref\_genome**  | str  | path to reference FASTA         | `references/GRCh38.fa`         |
| **known\_sites** | str  | dbSNP / known indels VCF        | `references/dbsnp.hg38.vcf.gz` |
| **threads**      | int  | number of CPU threads per job   | 16                             |
| **memory**       | str  | memory per process              | 19G                            |

---

## Example run

Dry run (see planned steps):

```bash
snakemake -n --use-conda
```

Run on 16 cores with conda:

```bash
snakemake --cores 16 --use-conda
```

Run on SLURM cluster:

```bash
snakemake --profile slurm
```

---

🔥 This replaces the toy workflow with a **real WES mutation calling pipeline** while keeping the same structure and readability.

Do you want me to also **rewrite the `Snakefile` rules** (alignment, BQSR, Mutect2, etc.) to match your bash pipeline? That way you’d have a **fully automated Snakemake version** of your script.
