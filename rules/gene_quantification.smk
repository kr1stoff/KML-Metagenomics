"""
bowtie2 比对到去冗余基因集 + samtools 排序
"""


rule bowtie2_build:
    input:
        "{sample}.cd-hit.fna",
    output:
        idx1="{sample}.cd-hit.1.bt2",
        idx2="{sample}.cd-hit.2.bt2",
        idx3="{sample}.cd-hit.3.bt2",
        idx4="{sample}.cd-hit.4.bt2",
        rev1="{sample}.cd-hit.rev.1.bt2",
        rev2="{sample}.cd-hit.rev.2.bt2",
    benchmark:
        ".log/quantification/bowtie2_build/{sample}.bowtie2_build.bm"
    log:
        ".log/quantification/bowtie2_build/{sample}.bowtie2_build.log",
    conda:
        config["conda"]["basic2"]
    threads: config["threads"]["high"]
    shell:
        "bowtie2-build --threads {threads} "
        "{input} {wildcards.sample}.cd-hit "
        "2> {log}"


rule bowtie2_map:
    input:
        r1="{sample}_host_removed.1.fastq.gz",
        r2="{sample}_host_removed.2.fastq.gz",
        idx=rules.bowtie2_build.output,
    output:
        "{sample}.cd-hit.bowtie2.sam",
    benchmark:
        ".log/quantification/bowtie2_map/{sample}.bowtie2_map.bm"
    log:
        ".log/quantification/bowtie2_map/{sample}.bowtie2_map.log",
    conda:
        config["conda"]["basic2"]
    threads: config["threads"]["high"]
    params:
        extra=config["bowtie2"]["mapping_extra"],
    shell:
        "bowtie2 -x {wildcards.sample}.cd-hit "
        "{params.extra} "
        "--threads {threads} "
        "-1 {input.r1} -2 {input.r2} "
        "-S {output} "
        "2> {log}"


rule samtools_sort:
    input:
        "{sample}.cd-hit.bowtie2.sam",
    output:
        "{sample}.cd-hit.bowtie2.sorted.bam",
    benchmark:
        ".log/quantification/samtools/{sample}.samtools_sort.bm"
    log:
        ".log/quantification/samtools/{sample}.samtools_sort.log",
    conda:
        config["conda"]["basic2"]
    threads: config["threads"]["medium"]
    shell:
        "samtools view -bS {input} | "
        "samtools sort -o {output} - "
        "2> {log}"
