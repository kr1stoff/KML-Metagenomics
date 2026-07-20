"""
MetaGeneMark 基因预测 + seqtk 过滤 (>100bp)
"""


rule metagenemark:
    input:
        rules.seqtk_filter_contigs.output,
    output:
        gff="gene_prediction/{sample}.gm.gff",
        faa="gene_prediction/{sample}.gm.faa",
        fna="gene_prediction/{sample}.gm.fna",
    benchmark:
        ".log/gene_prediction/{sample}.metagenemark.bm"
    log:
        ".log/gene_prediction/{sample}.metagenemark.log",
    params:
        gmhmmp=config["software"]["gmhmmp"],
        model=config["database"]["metagenemark_model"],
        extra=config["params"]["gmhmmp"],
    shell:
        "{params.gmhmmp} -m {params.model} {params.extra} {input} -o {output.gff} -A {output.faa} -D {output.fna} -L {log}"


# 保留长度 > 100bp 的基因序列
rule seqtk_filter_genes:
    input:
        rules.metagenemark.output.fna,
    output:
        "gene_prediction/{sample}.gm.gt100bp.fna",
    benchmark:
        ".log/gene_prediction/{sample}.seqtk_gt100bp.bm"
    log:
        ".log/gene_prediction/{sample}.seqtk_gt100bp.log",
    conda:
        config["conda"]["seqtk"]
    shell:
        "seqtk seq -L 100 {input} > {output} 2> {log}"
