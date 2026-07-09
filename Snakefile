"""
KML-Metagenomics 宏基因组分析流程
步骤: fastp → fastqc → multiqc → 去宿主 → 组装 → 基因预测 → 基因聚类 → 丰度定量
"""

import pandas as pd


configfile: workflow.source_path("config.yaml")
config_schema: workflow.source_path("config.schema.yaml")

# 读取样本表
samples = pd.read_csv(config["samples_tsv"], sep="\t", header=None, dtype=str)
samples.columns = ["sample_id", "fq1", "fq2"]


rule all:
    input:
        expand("{sample}.cd-hit.bowtie2.sorted.bam",
               sample=samples["sample_id"].tolist()),
        "qc/multiqc_data/multiqc_report.html",


include: "rules/qc.smk"
include: "rules/host_removal.smk"
include: "rules/assembly.smk"
include: "rules/gene_prediction.smk"
include: "rules/gene_clustering.smk"
include: "rules/quantification.smk"
