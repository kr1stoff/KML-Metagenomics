"""
megahit 宏基因组组装 + seqtk 过滤 (>500bp)
"""


rule megahit:
    input:
        r1="{sample}_host_removed.1.fastq",
        r2="{sample}_host_removed.2.fastq",
    output:
        directory("{sample}_megahit"),
    benchmark:
        ".log/assembly/megahit/{sample}.megahit.bm"
    log:
        ".log/assembly/megahit/{sample}.megahit.log",
    conda:
        config["conda"]["meta"]
    threads: config["threads"]["high"]
    params:
        extra=config["megahit"]["extra"],
    shell:
        "megahit {params.extra} -t {threads} "
        "-1 {input.r1} -2 {input.r2} "
        "-o {output} "
        "2> {log}"


rule seqtk_filter_contigs:
    """保留长度 > 500bp 的 contig"""
    input:
        "{sample}_megahit/final.contigs.fa",
    output:
        "{sample}.contigs.gt500.fa",
    benchmark:
        ".log/assembly/seqtk/{sample}.seqtk_gt500.bm"
    log:
        ".log/assembly/seqtk/{sample}.seqtk_gt500.log",
    conda:
        config["conda"]["qc"]
    shell:
        "seqtk seq -L 500 {input} > {output} 2> {log}"
