"""
bowtie2 去宿主 + 解压 (供 megahit 使用)
"""


rule bowtie2_host_removal:
    input:
        r1="{sample}.cleaned.1.fastq.gz",
        r2="{sample}.cleaned.2.fastq.gz",
    output:
        r1="{sample}_host_removed.1.fastq.gz",
        r2="{sample}_host_removed.2.fastq.gz",
        sam=temp("{sample}.bowtie2.sam"),
    benchmark:
        ".log/host_removal/bowtie2/{sample}.bowtie2_host_removal.bm"
    log:
        ".log/host_removal/bowtie2/{sample}.bowtie2_host_removal.log",
    conda:
        config["conda"]["basic2"]
    threads: config["threads"]["high"]
    params:
        extra=config["bowtie2"]["host_removal_extra"],
        host_ref=config["database"]["host_reference"],
    shell:
        "bowtie2 -x {params.host_ref} "
        "-1 {input.r1} -2 {input.r2} "
        "{params.extra} "
        "--threads {threads} "
        "--un-conc-gz {wildcards.sample}_host_removed.%.fastq.gz "
        "-S {output.sam} "
        "2> {log}"


rule gunzip_host_removed:
    """解压去宿主后的 fastq，供 megahit 组装使用"""
    input:
        r1="{sample}_host_removed.1.fastq.gz",
        r2="{sample}_host_removed.2.fastq.gz",
    output:
        r1=temp("{sample}_host_removed.1.fastq"),
        r2=temp("{sample}_host_removed.2.fastq"),
    benchmark:
        ".log/host_removal/gunzip/{sample}.gunzip.bm"
    log:
        ".log/host_removal/gunzip/{sample}.gunzip.log",
    shell:
        "gunzip -c {input.r1} > {output.r1} 2>> {log}; "
        "gunzip -c {input.r2} > {output.r2} 2>> {log}"
