rule bwa_mem2:
    input:
        r1=config["data"]["tumor_r1"],
        r2=config["data"]["tumor_r2"],
        fa=config["references"]["genome_fa"]
    output:
        bam="results/bam/{sample}.sorted.bam"
    threads: config["threads"]
    shell:
        """
        bwa-mem2 mem -t {threads} {input.fa} {input.r1} {input.r2} | \
        samtools sort -o {output.bam}
        samtools index {output.bam}
        """
