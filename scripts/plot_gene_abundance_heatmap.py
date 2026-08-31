#!/home/mengxf/miniforge3/envs/python3.12/bin/python
"""基于基因丰度表绘制样本间相关热图。

输入: gene_abundance_table.tsv (行=基因, 列=样本, tab 分隔)
输出: 样本 x 样本 相关热图 CSV 和 PNG (横纵坐标均为样本)
"""

from pathlib import Path
import pandas as pd
import seaborn as sns
import scipy.spatial.distance as ssd
from scipy.cluster.hierarchy import linkage as scipy_linkage
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys

sys.stderr = open(snakemake.log[0], "w")

# IO
in_path = Path(snakemake.input[0])
out_png = snakemake.output['png']
out_csv = snakemake.output['csv']
# 支持 "spearman", "pearson" 两种算法
method = snakemake.params['method']

# MAIN
# 1. 读取基因丰度表: 行=基因, 列=样本
df = pd.read_csv(in_path, sep="\t", index_col=0)
df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
sys.stderr.write(f"基因数: {df.shape[0]}, 样本数: {df.shape[1]}")
sys.stderr.write(f"样本: {list(df.columns)}")

# 2. 计算样本间相关性矩阵 (横纵坐标均为样本)
#    关键: 不转置! df 行=基因 列=样本, df.corr() 按列(样本)计算相关,
#    若用 df.T.corr() 会退化成 基因x基因 级别的计算, 慢几个数量级.
corr = df.corr(method=method)
sys.stderr.write(f"\n{method} 相关性矩阵 ({corr.shape[0]} x {corr.shape[1]}):\n")
sys.stderr.write(corr.round(4).to_string())

# 3. 用 scipy 预计算聚类: pdist 得到 condensed distance, linkage 聚类,
#    行/列复用同一份 linkage (相关矩阵对称), 避免 seaborn 重复计算, 显著加速
dist = ssd.pdist(corr.values, metric="euclidean")
link = scipy_linkage(dist, method="average", optimal_ordering=False)

n = corr.shape[0]
show_annot = n <= 30  # 样本过多时省略格内标注, 避免拥挤与绘图变慢
g = sns.clustermap(
    corr,
    row_linkage=link, col_linkage=link,
    cmap="RdBu_r", center=0,
    vmin=-1, vmax=1,
    annot=show_annot, fmt=".2f",
    annot_kws={"size": 11},
    linewidths=0.8, linecolor="white",
    figsize=(8, 8),
    xticklabels=True, yticklabels=True,
)
g.ax_heatmap.set_xlabel("Samples")
g.ax_heatmap.set_ylabel("Samples")
g.fig.suptitle(f"Sample-to-Sample {method.title()} Correlation of Gene Abundance",
                y=1.02, fontsize=14)
g.savefig(out_png, dpi=300, bbox_inches="tight")
sys.stderr.write(f"\n热图已保存至: {out_png}")

# 4. 同时输出相关性矩阵 CSV
csv_path = Path(out_csv)
corr.to_csv(csv_path, sep=",")
sys.stderr.write(f"相关性矩阵已保存至: {csv_path}")
