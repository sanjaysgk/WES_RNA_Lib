# Reference genome and variant sites
REF_GENOME = config["references"]["genome_fa"]
KNOWN_SITES_DBSNP = config["references"]["dbsnp"]
KNOWN_SITES_MILLS = config["references"]["mills"]
KNOWN_SITES_INDEL = config["references"]["indels"]
KNOWN_SITES_SMALL_EXON = config["references"]["small_exon"]
KNOWN_SITES_PON = config["references"]["pon"]
KNOWN_SITES_GERMLINE = config["references"]["germline"]

# Sample and directories
SAMPLE_ID = config["sample_id"]
BASE_DIR = config["base_dir"]
DATA_DIR = config["data"]["wes_dir"]
OUTPUT_DIR = config["data"]["output_dir"]

# Input FASTQs
TUMOR_R1 = config["fastq"]["tumor_r1"]
TUMOR_R2 = config["fastq"]["tumor_r2"]
NORMAL_R1 = config["fastq"]["normal_r1"]
NORMAL_R2 = config["fastq"]["normal_r2"]

THREADS = config["threads"]

rule index_reference:
    input:
        fa=config["references"]["genome_fa"]
    output:
        fai=config["references"]["genome_fa"] + ".fai",
        dict=config["references"]["genome_fa"].replace(".fa", ".dict")
    threads: config["threads"]
    shell:
        """
        samtools faidx {input.fa}
        gatk CreateSequenceDictionary -R {input.fa} -O {output.dict}
        """

#https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling/blob/main/workflow/rules/ref.smk
#snakemake wrapper for common tasks https://snakemake-wrappers.readthedocs.io/en/stable/
#
#Writen by: johanneskoester [https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling/commits?author=johanneskoester]

# Rule: prepare reference genome
rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "logs/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    wrapper:
        "0.74.0/bio/reference/ensembl-sequence"

checkpoint genome_faidx:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    wrapper:
        "0.74.0/bio/samtools/faidx"


rule genome_dict:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.dict",
    log:
        "logs/samtools/create_dict.log",
    conda:
        "../envs/samtools.yaml"
    cache: True
    shell:
        "samtools dict {input} > {output} 2> {log} "


rule bwa_index:
    input:
        "resources/genome.fasta",
    output:
        multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "0.74.0/bio/bwa/index"


rule get_vep_cache:
    output:
        directory("resources/vep/cache"),
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    log:
        "logs/vep/cache.log",
    wrapper:
        "0.74.0/bio/vep/cache"


rule get_vep_plugins:
    output:
        directory("resources/vep/plugins"),
    log:
        "logs/vep/plugins.log",
    params:
        release=config["ref"]["release"],
    wrapper:
        "0.74.0/bio/vep/plugins"
