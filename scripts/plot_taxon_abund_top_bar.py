#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""通用堆叠柱形图: 对当前目录下所有 X.abund_cum.tsv (X=K/P/C/O/F/G/S) 分别绘图.
x=样本, y=物种百分比; 物种数>10 时仅展示前10, 其余(含 *__unknown)归入 Others"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import sys

sys.stderr = open(snakemake.log[0], "w")

# 各级别首字母对应的名称 (用于图标题)
LEVEL_NAME = {"K": "Kingdom", "P": "Phylum", "C": "Class", "O": "Order", "F": "Family", "G": "Genus", "S": "Species"}

TOP_N = int(snakemake.params.get("top_n", 10))          # 展示的物种数上限
OTHERS_COLOR = "#999999"                                # Others 使用的灰色


def make_plot_df(df: pd.DataFrame) -> pd.DataFrame:
    """输入丰度矩阵(行=物种, 列=样本), 返回绘图矩阵: 前10物种+Others, 或全部物种(物种数<=10时)"""
    df = df.copy()
    # 将 *__unknown 改名为 Others, 使其必然归入 Others 分组
    df.index = df.index.map(lambda s: "Others" if s.endswith("__unknown") else s)
    total = df.sum(axis=1)                       # 每个物种在所有样本中的总丰度
    if len(total) > TOP_N:                       # 物种数超过上限时截断
        top = total.drop("Others", errors="ignore").nlargest(TOP_N).index  # 全局前 TOP_N 的物种
        others = total.index.difference(top)     # 剩余物种名 (含 Others)
        plot_df = df.loc[top].copy()             # 前 TOP_N 物种的丰度矩阵
        plot_df.loc["Others"] = df.loc[others].sum(axis=0)  # 其余物种按样本加和成一行
    else:                                        # 物种数少时全部保留
        plot_df = df
    # 每列(样本)转换为百分比, 除零补 0
    return plot_df.div(plot_df.sum(axis=0), axis=1).fillna(0) * 100


def plot_one(input_tsv: str, out_png: str) -> None:
    """对单个丰度表绘图并保存 png, 输出文件名为 <前缀>.abund_cum.top<TOP_N>_bar.png"""
    level = os.path.basename(input_tsv)[0]       # 文件前缀首字符即级别 (K/P/C/O/F/G/S)
    df = pd.read_csv(input_tsv, sep="\t", index_col="taxon")  # 读取丰度矩阵
    plot_df = make_plot_df(df)                   # 组装绘图矩阵 (前10+Others 或全部)
    n = len(plot_df.index)                       # 绘图物种数 (含 Others)
    # 颜色: 各物种用 tab20 循环取色, 最后一行 Others 用灰色
    colors = [plt.cm.tab20(i % 20) for i in range(n - 1)] + [OTHERS_COLOR]
    ax = plot_df.T.plot(kind="bar", stacked=True, figsize=(10, 7), color=colors, width=0.7)
    ax.set_xlabel("Sample")                      # x 轴标签
    ax.set_ylabel("Relative abundance (%)")      # y 轴标签
    ax.set_title(f"{LEVEL_NAME[level]} relative abundance (top {TOP_N} + Others)")  # 图标题
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", title="Taxon")  # 图例放右侧
    plt.tight_layout()                           # 自动调整布局
    plt.savefig(out_png, dpi=300)                # 保存高清图
    plt.close()                                  # 关闭画布
    print(out_png)


def main() -> None:
    input_tsv = snakemake.input[0]
    output_png = snakemake.output[0]
    plot_one(input_tsv, output_png)


if __name__ == "__main__":
    main()
