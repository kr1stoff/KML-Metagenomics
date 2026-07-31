"""
megahit 宏基因组组装 + seqtk 过滤 (>500bp)
"""


rule megahit:
    input:
        r1=rules.gunzip_host_removed.output.r1,
        r2=rules.gunzip_host_removed.output.r2,
    output:
        dir=directory("assembly/megahit/{sample}"),
        fa="assembly/megahit/{sample}/final.contigs.fa",
        # todo 删掉中间文件
        inter=directory("assembly/megahit/{sample}/intermediate_contigs"),
    benchmark:
        ".log/assembly/{sample}.megahit.bm"
    log:
        ".log/assembly/{sample}.megahit.log",
    conda:
        config["conda"]["megahit"]
    threads: config["threads"]["high"]
    params:
        extra=config["params"]["megahit"],
    shell:
        "megahit {params.extra} -t {threads} -1 {input.r1} -2 {input.r2} -o {output.dir} --force 2> {log}"


# 保留长度 > 500bp 的 contig
rule seqtk_filter_contigs:
    input:
        rules.megahit.output.fa,
    output:
        "assembly/contigs_gt500/{sample}.fa",
    benchmark:
        ".log/assembly/seqtk/{sample}.seqtk_gt500.bm"
    log:
        ".log/assembly/seqtk/{sample}.seqtk_gt500.log",
    conda:
        config["conda"]["seqtk"]
    params:
        "-L 500",
    shell:
        "seqtk seq {params} {input} > {output} 2> {log}"


rule quast:
    input:
        rules.seqtk_filter_contigs.output,
    output:
        "assembly/quast/{sample}/report.tsv",
    benchmark:
        ".log/assembly/{sample}.quast.bm"
    log:
        ".log/assembly/{sample}.quast.log",
    conda:
        config["conda"]["quast"]
    shell:
        # 不需要多线程
        "quast {input} -o assembly/quast/{wildcards.sample} 2> {log}"


rule assembly_stats:
    input:
        expand("assembly/quast/{sample}/report.tsv", sample=samples),
    output:
        "upload/assembly_stats.xlsx",
    benchmark:
        ".log/assembly/stats.bm"
    log:
        ".log/assembly/stats.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/assembly_stats.py"
