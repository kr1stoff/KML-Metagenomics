# -*- coding: utf-8 -*-
"""
Core-Pan 基因累积曲线分析
基于基因存在/缺失矩阵，通过随机抽样估计不同样本数量下的核心基因数和泛基因数，
绘制核心基因和泛基因的累积曲线。

用法（Snakemake rule）：
    rule core_pan_gene_analysis:
        input:  "gene_presence_absence.tsv"
        output: "core_pan_curve.pdf"
        log:    "core_pan_gene_analysis.log"
        script: "scripts/core_pan_gene_analysis.py"

输入格式：TSV 文件，行为基因，列为样本，值为 0/1（存在/缺失）。
"""

import random
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---- Snakemake 日志重定向 ----
sys.stderr = open(snakemake.log[0], "w")

# ---- 输入输出 ----
input_file = snakemake.input[0]
output_file = snakemake.output[0]

# ---- 读取基因存在/缺失矩阵 ----
# 每一行是一个基因，每一列是一个样本，值为 True（存在）或 False（缺失）
presence = pd.read_csv(input_file, sep="\t", index_col=0).astype(bool)

samples = presence.columns
n_samples = len(samples)

# ---- 累积曲线参数 ----
repeat = 100  # 每个样本数下的随机抽样次数，用于稳定估计

sample_sizes = list(range(1, n_samples + 1))  # 样本数量序列 [1, 2, ..., N]
core_curve = []  # 每个样本数下的平均核心基因数
pan_curve = []   # 每个样本数下的平均泛基因数

# ---- 随机抽样估计累积曲线 ----
for n in sample_sizes:
    core_list = []  # 当前样本数下 repeat 次抽样的核心基因数
    pan_list = []   # 当前样本数下 repeat 次抽样的泛基因数

    for _ in range(repeat):
        # 从所有样本中随机抽取 n 个
        selected = random.sample(list(samples), n)
        subset = presence[selected]

        # 核心基因：在被选样本中都存在的基因
        core = subset.all(axis=1).sum()
        # 泛基因：在被选样本中至少存在一次的基因
        pan = subset.any(axis=1).sum()

        core_list.append(core)
        pan_list.append(pan)

    # 取 repeat 次随机抽样的均值，作为该样本数下的估计
    core_curve.append(np.mean(core_list))
    pan_curve.append(np.mean(pan_list))

# ---- 绘图 ----
fig, axes = plt.subplots(
    nrows=1,
    ncols=2,
    figsize=(10, 4),
    dpi=300,
)

# 左图：核心基因累积曲线
axes[0].plot(sample_sizes, core_curve, marker="o")
axes[0].set_xlabel("Number of samples")
axes[0].set_ylabel("Number of core genes")
axes[0].set_title("Core gene accumulation curve")

# 右图：泛基因累积曲线
axes[1].plot(sample_sizes, pan_curve, marker="o")
axes[1].set_xlabel("Number of samples")
axes[1].set_ylabel("Number of pan genes")
axes[1].set_title("Pan gene accumulation curve")

plt.tight_layout()
plt.savefig(output_file, bbox_inches="tight")
