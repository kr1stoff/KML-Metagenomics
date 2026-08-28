#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""通用热图: 对当前目录下所有 X.abund_cum.tsv / X.gene_cum.tsv (X=K/P/C/O/F/G/S) 分别绘图.
取全局总和前30的物种, x=样本, y=物种, 颜色=log10(值+1); 不限制分类级别与丰度/基因计数"""

import glob
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

sys.stderr = open(snakemake.log[0], "w")

# 各级别首字母对应的名称 (用于图标题)
LEVEL_NAME = {"K": "Kingdom", "P": "Phylum", "C": "Class", "O": "Order", "F": "Family", "G": "Genus", "S": "Species"}

TOP_N = int(snakemake.params.get("top_n", 30))          # 展示的物种数上限
CMAP = "YlOrRd"                                         # 热图色阶


def main() -> None:
    input_file = snakemake.input[0]
    out_png = snakemake.output[0]
    # 读取矩阵: 行=物种, 列=样本
    df = pd.read_csv(input_file, sep="\t", index_col="taxon")
    # 按全部样本总和降序取前 TOP_N 个物种
    top = df.sum(axis=1).nlargest(TOP_N).index
    sub = df.loc[top]
    # 去掉物种名前缀, 如 g__Actinomyces -> Actinomyces
    sub.index = sub.index.map(lambda s: s.split("__", 1)[-1])
    # log10(x+1) 变换, 缓解数值跨度大导致的色阶失真
    log_df = np.log10(sub + 1)
    # 画热图
    plt.figure(figsize=(8, 10))
    sns.heatmap(log_df, cmap=CMAP, cbar_kws={"label": "log10(value + 1)"})
    plt.xlabel("Sample")                         # x 轴标签
    plt.ylabel("Taxon")                          # y 轴标签
    level = os.path.basename(input_file)[0]              # 级别首字母
    kind = "gene count" if "gene" in input_file else "abundance"  # 丰度或基因计数
    plt.title(f"Top {TOP_N} taxa ({LEVEL_NAME[level]}, {kind})")  # 图标题
    plt.tight_layout()                           # 自动调整布局
    plt.savefig(out_png, dpi=300)                # 保存高清图
    plt.close()                                  # 关闭画布
    print(out_png)


if __name__ == "__main__":
    main()
