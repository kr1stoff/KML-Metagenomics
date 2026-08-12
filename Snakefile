"""
KML-Metagenomics 宏基因组分析流程
步骤: fastp → fastqc → multiqc → 去宿主 → 组装 → 基因预测 → 基因聚类 → 丰度定量
"""

import pandas as pd


configfile: workflow.source_path("config.yaml")


config_schema: workflow.source_path("config.schema.yaml")

# 读取样本表
samples_df = pd.read_csv(config["samples_tsv"], sep="\t", header=None, dtype=str)
samples_df.columns = ["sample", "fq1", "fq2"]
samples = samples_df["sample"].tolist()


rule all:
    input:
        # 基因丰度
        "upload/gene_catalogue_stats.xlsx",
        "gene_quantification/gene_reads_table.tsv",
        "upload/core_pan_stats.png",
        "upload/gene_abundance_heatmap.png",
        # 组装
        "upload/assembly_stats.xlsx",
        # 质控
        "qc/multiqc/multiqc_report.html",
        "upload/data_process_stats.xlsx",


include: "rules/rawdata.smk"
include: "rules/fastqc.smk"
include: "rules/host_removal.smk"
include: "rules/assembly.smk"
include: "rules/gene_prediction.smk"
include: "rules/gene_quantification.smk"
include: "rules/gene_statistics.smk"
