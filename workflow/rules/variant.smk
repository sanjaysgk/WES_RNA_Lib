rule mutect2:
    input:
        tumor="results/bam/{sample}.sorted.bam",
        normal="results/bam/{sample}_normal.sorted.bam",
        fa=config["references"]["genome_fa"]
    output:
        vcf="results/variants/{sample}.final.vcf.gz"
    threads: config["threads"]
    shell:
        """
        gatk Mutect2 \
            -R {input.fa} \
            -I {input.tumor} \
            -I {input.normal} \
            --germline-resource {config[references][af_gnomad]} \
            --panel-of-normals {config[references][mills]} \
            -O {output.vcf}
        """
