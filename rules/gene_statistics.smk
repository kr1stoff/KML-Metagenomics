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


rule core_pan_stats:
    input:
        rules.gene_abundance_table.output,
    output:
        "upload/core_pan_stats.png",
    benchmark:
        ".log/gene_statistics/core_pan_stats.bm"
    log:
        ".log/gene_statistics/core_pan_stats.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/core_pan_stats.py"


rule gene_abundance_heatmap:
    input:
        rules.gene_abundance_table.output[0],
    output:
        png="upload/gene_abundance_heatmap.png",
        csv="upload/gene_abundance_heatmap.csv",
    benchmark:
        ".log/gene_statistics/gene_abundance_heatmap.bm"
    log:
        ".log/gene_statistics/gene_abundance_heatmap.log",
    conda:
        config["conda"]["python"]
    params:
        # 支持 "spearman", "pearson" 两种算法
        method="spearman"
    script:
        "../scripts/gene_abundance_heatmap.py"


rule gene_count_boxplot:
    input:
        abund=rules.gene_abundance_table.output[0],
        meta=config["metadata"]
    output:
        png="upload/groups_gene_count_boxplot.png",
        csv="upload/groups_gene_count_boxplot.csv",
    benchmark:
        ".log/gene_statistics/gene_count_boxplot.bm"
    log:
        ".log/gene_statistics/gene_count_boxplot.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/gene_count_boxplot.py"
