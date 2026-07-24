"""
bowtie2 比对到去冗余基因集 + samtools 排序
"""


rule quant_gene_cat_bowtie2_build:
    input:
        rules.cd_hit.output,
    output:
        idx1="gene_prediction/gene_catalogue.1.bt2",
        idx2="gene_prediction/gene_catalogue.2.bt2",
        idx3="gene_prediction/gene_catalogue.3.bt2",
        idx4="gene_prediction/gene_catalogue.4.bt2",
        rev1="gene_prediction/gene_catalogue.rev.1.bt2",
        rev2="gene_prediction/gene_catalogue.rev.2.bt2",
    benchmark:
        ".log/gene_quantification/bowtie2_build.bm"
    log:
        ".log/gene_quantification/bowtie2_build.log",
    conda:
        config["conda"]["bowtie2"]
    threads: config["threads"]["high"]
    shell:
        "bowtie2-build --threads {threads} {input} gene_prediction/gene_catalogue 2> {log}"


rule quant_bowtie2_map:
    input:
        r1=rules.bowtie2_host_removal.output.r1,
        r2=rules.bowtie2_host_removal.output.r2,
        idx=rules.quant_gene_cat_bowtie2_build.output,
    output:
        temp("gene_quantification/{sample}.gene_catalogue.sam"),
    benchmark:
        ".log/gene_quantification/bowtie2_map/{sample}.bowtie2_map.bm"
    log:
        ".log/gene_quantification/bowtie2_map/{sample}.bowtie2_map.log",
    conda:
        config["conda"]["bowtie2"]
    threads: config["threads"]["high"]
    params:
        extra=config["params"]["bowtie2"]["mapping"],
    shell:
        "bowtie2 -x gene_prediction/gene_catalogue {params.extra} --threads {threads} -1 {input.r1} -2 {input.r2} -S {output} 2> {log}"


rule samtools_sort_index:
    input:
        rules.quant_bowtie2_map.output,
    output:
        bam="{sample}.gene_catalogue.sorted.bam",
        bai="{sample}.gene_catalogue.sorted.bam.bai",
    benchmark:
        ".log/quantification/samtools/{sample}.samtools_sort.bm"
    log:
        ".log/quantification/samtools/{sample}.samtools_sort.log",
    conda:
        config["conda"]["samtools"]
    threads: config["threads"]["medium"]
    shell:
        """
        samtools sort -o {output.bam} {input} 2> {log}
        samtools index {output.bam}
        """
