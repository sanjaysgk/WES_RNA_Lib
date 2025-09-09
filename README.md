# Snakemake workflow: `WES_RNA_Pipeline`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/sanjaysgk/WES_RNA_Pipeline/workflows/Tests/badge.svg?branch=main)](https://github.com/sanjaysgk/WES_RNA_Pipeline/actions?query=branch%3Amain+workflow%3ATests)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/sanjaysgk/WES_RNA_Pipeline)

A Snakemake workflow for **Whole Exome Sequencing (WES) and RNA-seq processing**, including alignment, BAM processing, quality control, and somatic mutation calling with MuTect2.

- [Snakemake workflow: `WES_RNA_Pipeline`](#snakemake-workflow-wes_rna_pipeline)
  - [Usage](#usage)
  - [Deployment options](#deployment-options)
  - [Authors](#authors)
  - [References](#references)
  - [TODO](#todo)

## Usage

The workflow is designed to process WES and RNA-seq data starting from raw FASTQ files and produce processed BAM files, quality control reports, and variant calls.

Detailed information about input data and workflow configuration can be found in the [`config/README.md`](config/README.md).

If you use this workflow in a publication, please give credit by citing this repository URL: `https://github.com/sanjaysgk/WES_RNA_Pipeline`

## Deployment options

To run the workflow from the command line:

```bash
cd path/to/WES_RNA_Pipeline
```