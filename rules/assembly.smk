"""
megahit 宏基因组组装 + seqtk 过滤 (>500bp)
"""


rule megahit:
    input:
        r1=rules.gunzip_host_removed.output.r1,
        r2=rules.gunzip_host_removed.output.r2,
    output:
        dir=directory("megahit/{sample}"),
        fa="megahit/{sample}/final.contigs.fa",
        # 删掉中间文件
        inter=temp(directory("megahit/{sample}/intermediate_contigs")),
    benchmark:
        ".log/assembly/megahit/{sample}.megahit.bm"
    log:
        ".log/assembly/megahit/{sample}.megahit.log",
    conda:
        config["conda"]["megahit"]
    threads: config["threads"]["high"]
    params:
        extra=config["params"]["megahit"],
    shell:
        "megahit {params.extra} -t {threads} -1 {input.r1} -2 {input.r2} -o {output.dir} 2> {log}"


# 保留长度 > 500bp 的 contig
rule seqtk_filter_contigs:
    input:
        rules.megahit.output.fa,
    output:
        "megahit/contigs_gt500/{sample}.fa",
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
