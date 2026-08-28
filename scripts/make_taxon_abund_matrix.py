#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""按 snakemake 规则生成物种丰度矩阵 / 物种基因数矩阵。

规则: rules/taxon_classification.smk 中的 make_taxon_abund_matrix
输入: 所有样本的 taxon_classification/{sample}.linage_abund.txt 明细文件
输出: taxon_classification/matrix/{rank}.{kind}.txt (行=物种, 列=样本, 缺失补 0)
通配符:
    rank: K/P/C/O/F/G/S (分类级别)
    kind: abund (丰度矩阵, 取 abund_cum 列) / gene (基因数矩阵, 取 gene_cum 列)
"""

import os
import sys
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

# 分类级别前缀
LEVEL_PREFIX = {"K": "k__", "P": "p__", "C": "c__", "O": "o__", "F": "f__", "G": "g__", "S": "s__"}
# 矩阵类型 -> 明细文件中的数值列
KIND_COLUMN = {"abund": "abund_cum", "gene": "gene_cum"}


def format_taxon(lineage: str, level: str) -> str:
    """把 lineage 字符串截取到指定级别, 格式化为 k__xxx 形式"""
    # 明细中该级别的 lineage 形如 "Bacteria; Firmicutes", 末位即该级别的物种名
    parts = lineage.split("; ")
    return LEVEL_PREFIX[level] + parts[-1]


def main() -> None:
    """读取各样本明细, 按通配符 rank/kind 生成对应矩阵"""
    rank = snakemake.wildcards.rank
    kind = snakemake.wildcards.kind
    out_path = str(snakemake.output)
    value_col = KIND_COLUMN[kind]

    # 读取所有样本明细, 只保留目标级别, 按格式化物种名聚合为 序列
    series = {}
    for fp in snakemake.input:
        sample = os.path.basename(fp).split(".")[0]   # 样本名 (文件名前缀)
        df = pd.read_csv(fp, sep="\t")
        sub = df[df["level"] == rank]                 # 只取该级别的行
        taxa = sub["lineage"].map(lambda s: format_taxon(s, rank))
        series[sample] = pd.Series(sub[value_col].values, index=taxa).groupby(level=0).sum()

    # 合并为 行=物种, 列=样本 的矩阵, 缺失补 0, 物种名排序
    matrix = pd.DataFrame(series).fillna(0).sort_index()
    matrix.to_csv(out_path, sep="\t", index_label="taxon")

    sys.stderr.write(f"级别 {rank} / 类型 {kind}: {len(matrix)} 个物种, {len(series)} 个样本\n")
    sys.stderr.write(f"矩阵已写出: {out_path}\n")


if __name__ == "__main__":
    main()
