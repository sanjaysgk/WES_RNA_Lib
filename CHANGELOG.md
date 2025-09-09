# Changelog

## [1.0.0] (2025-09-09)

### Features

* Initial commit: WES & RNA-seq pipeline for tumor/normal samples.
* Integrated BWA-mem2 alignment with read groups.
* Added GATK SortSam, MarkDuplicates, and CollectAlignmentSummaryMetrics steps.
* Automatic reference genome indexing if not present.
* Configurable sample IDs, input fastq files, reference genome, and known sites.
* Automated directory creation for logs, temp files, QC, and variants.
* Setup for SLURM cluster execution with customizable resources (threads, memory, time, partition).
* Included modular Snakemake workflow structure (`rules`, `scripts`, `config`, `workflow/Snakefile`).
* Added environment setup script (`install_conda.sh`) for automatic Miniconda installation and environment creation.
* Added README.md with usage instructions, deployment options, authorship, and references.
* Added GitHub Actions workflows for linting, tests, and CI integration.

### Bug Fixes / Improvements

* Fixed alignment and sorting commands to ensure proper BAM generation.
* Removed temporary Miniconda installer after environment setup.
* Corrected EOF and quoting issues in `install_conda.sh`.
* Verified SLURM headers and resource allocation for alignment and variant calling steps.
* Added detailed logging for input validation, reference genome preparation, and pipeline steps.
* Updated README badges and references to point to your repository (`https://github.com/sanjaysgk/WES_RNA_Pipeline.git`).

