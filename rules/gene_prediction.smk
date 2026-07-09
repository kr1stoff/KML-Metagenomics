"""
MetaGeneMark 基因预测 + seqtk 过滤 (>100bp)
"""


rule metagenemark:
    input:
        "{sample}.contigs.gt500.fa",
    output:
        gff="{sample}.gm.gff",
        faa="{sample}.gm.faa",
        fna="{sample}.gm.fna",
    benchmark:
        ".log/gene_prediction/metagenemark/{sample}.metagenemark.bm"
    log:
        ".log/gene_prediction/metagenemark/{sample}.metagenemark.log",
    params:
        bin=config["software"]["metagenemark_bin"],
        model=config["software"]["metagenemark_model"],
    threads: 1
    shell:
        "{params.bin} -m {params.model} "
        "-a -d -f G -p {threads} "
        "{input} "
        "-o {output.gff} "
        "-A {output.faa} "
        "-D {output.fna} "
        "-L {wildcards.sample}.gm.log "
        "2> {log}"


rule seqtk_filter_genes:
    """保留长度 > 100bp 的基因序列"""
    input:
        "{sample}.gm.fna",
    output:
        "{sample}.gm.gt100bp.fna",
    benchmark:
        ".log/gene_prediction/seqtk/{sample}.seqtk_gt100bp.bm"
    log:
        ".log/gene_prediction/seqtk/{sample}.seqtk_gt100bp.log",
    conda:
        config["conda"]["qc"]
    shell:
        "seqtk seq -L 100 {input} > {output} 2> {log}"
