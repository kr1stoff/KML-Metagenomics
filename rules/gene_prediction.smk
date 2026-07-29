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


rule format_fasta_head:
    input:
        rules.metagenemark.output.fna,
    output:
        # todo 删除中间文件
        "gene_prediction/{sample}.format_head.fna",
    benchmark:
        ".log/gene_prediction/{sample}.format_fasta_head.bm"
    log:
        ".log/gene_prediction/{sample}.format_fasta_head.log",
    shell:
        "sed 's/|.*/_'{wildcards.sample}'/g' {input} > {output} 2> {log}"


# 保留长度 > 100bp 的基因序列
rule seqtk_filter_genes:
    input:
        rules.format_fasta_head.output,
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


# 合并所有样本的基因序列
rule merge_all_genes:
    input:
        expand("gene_prediction/{sample}.gm.gt100bp.fna", sample=samples),
    output:
        "gene_prediction/all_samples.orf.fna",
    benchmark:
        ".log/gene_prediction/all_samples.orf.bm"
    log:
        ".log/gene_prediction/all_samples.orf.log",
    shell:
        "cat {input} > {output} 2> {log}"


rule cd_hit:
    input:
        rules.merge_all_genes.output,
    output:
        "gene_prediction/gene_catalogue.raw.fna",
    benchmark:
        ".log/gene_prediction/all_samples.cd-hit.bm"
    log:
        ".log/gene_prediction/all_samples.cd-hit.log",
    threads: config["threads"]["high"]
    conda:
        config["conda"]["cd_hit"]
    params:
        extra=config["params"]["cd_hit"],
    shell:
        # protein 用 cd-hit, nucleotide 用 cd-hit-est
        "cd-hit-est -T {threads} {params.extra} -i {input} -o {output} 2> {log}"
