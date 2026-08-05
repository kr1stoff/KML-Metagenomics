# Gene catalogue 基本信息统计表 + 长度分布统计图
rule gene_stats:
    input:
        rules.gene_catalogue_unigene.output.unigene,
    output:
        xlsx="upload/gene_catalogue_stats.xlsx",
        png="upload/gene_catalogue_length_distribution.png",
    benchmark:
        ".log/gene_statistics/gene_stats.bm"
    log:
        ".log/gene_statistics/gene_stats.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/gene_stats.py"


# Core-pan 分析
