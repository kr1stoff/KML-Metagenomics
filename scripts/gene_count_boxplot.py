0#!/home/mengxf/miniforge3/envs/python3.12/bin/python
"""统计各样本非零基因数量, 并按分组绘制箱线图。

输入:
    gene_abundance_table.tsv (行=基因, 列=样本, tab 分隔, 0 表示缺失)
    metadata.tsv             (两列: 样本ID<TAB>分组)
输出:
    gene_count_boxplot.png (分组箱线图 + 散点)
    gene_count.csv         (每样本非零基因数量明细)
"""

import pandas as pd
import seaborn as sns
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys

sys.stderr = open(snakemake.log[0], "w")


# IO
in_file = snakemake.input["abund"]
metadata_file = snakemake.input["meta"]
out_path = snakemake.output["png"]
csv_path = snakemake.output["csv"]

# MAIN
# 1. 读取基因丰度表, 统计每个样本非零基因数量
df = pd.read_csv(in_file, sep="\t", index_col=0)
df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
sys.stderr.write(f"基因数: {df.shape[0]}, 样本数: {df.shape[1]}")

# 每样本非零(>0)基因数量
nonzero = (df > 0).sum(axis=0)
counts = nonzero.rename("gene_count").to_frame()
counts["sample"] = counts.index

# 2. 读取 metadata 合并分组
meta = pd.read_csv(metadata_file, sep="\t")
meta = meta[meta["sample"].notna()]  # 去掉末尾空行
data = counts.merge(meta, on="sample", how="left")
missing = data["group"].isna()
if missing.any():
    sys.stderr.write(f"警告: 以下样本在 metadata 中无分组: {list(data.loc[missing, 'sample'])}")
data = data.dropna(subset=["group"]).sort_values("group")
sys.stderr.write("\n各分组样本数:")
sys.stderr.write(data.groupby("group")["sample"].count().to_string())

# 输出明细 CSV
data[["sample", "group", "gene_count"]].to_csv(csv_path, index=False, sep="\t")
sys.stderr.write(f"\n明细已保存至: {csv_path}")

# 3. 绘制分组箱线图
groups = data["group"].unique()
n_groups = len(groups)
sns.set_theme(style="whitegrid", font_scale=1.1)
plt.figure(figsize=(max(6, 2.4 * n_groups), 5.5))

bp = sns.boxplot(data=data, x="group", y="gene_count",
                    hue="group", palette="Set2", width=0.5,
                    legend=False, fliersize=0)
sns.stripplot(data=data, x="group", y="gene_count",
                color="0.25", size=6, jitter=0.15, alpha=0.85)

bp.set_xlabel("Group")
bp.set_ylabel("Number of nonzero genes")
bp.set_title("Nonzero Gene Count per Sample by Group", fontsize=13, y=1.20)

# 在箱体上方标注各分组均值 ± 标准差
y_span = data["gene_count"].max() - data["gene_count"].min()
# 设置标注的 y 偏移量, 以避免与箱体重叠
y_offset = max(y_span * 0.20, 30)
for xpos, grp in enumerate(groups):
    vals = data.loc[data["group"] == grp, "gene_count"]
    mean, std = vals.mean(), vals.std()
    ymax = vals.max()
    bp.text(xpos, ymax + y_offset,
            f"n={len(vals)}\n{mean:.0f}±{std:.0f}",
            ha="center", va="bottom", fontsize=9, color="dimgray")

plt.tight_layout()
plt.savefig(out_path, dpi=300, bbox_inches="tight")
sys.stderr.write(f"\n箱线图已保存至: {out_path}")

# 打印每组统计摘要
sys.stderr.write("\n各组非零基因数量统计:")
sys.stderr.write(data.groupby("group")["gene_count"].describe().round(1).to_string())
