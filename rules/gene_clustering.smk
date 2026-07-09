"""
cd-hit 基因去冗余 (95% 相似度)
"""


rule cd_hit:
    input:
        "{sample}.gm.gt100bp.fna",
    output:
        fna="{sample}.cd-hit.fna",
        clstr="{sample}.cd-hit.fna.clstr",
    benchmark:
        ".log/gene_clustering/cd-hit/{sample}.cd_hit.bm"
    log:
        ".log/gene_clustering/cd-hit/{sample}.cd_hit.log",
    conda:
        config["conda"]["meta"]
    threads: config["threads"]["high"]
    params:
        extra=config["cd_hit"]["extra"],
    shell:
        "cd-hit {params.extra} "
        "-T {threads} "
        "-i {input} "
        "-o {output.fna} "
        "2> {log}"
