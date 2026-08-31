#!/home/mengxf/miniforge3/envs/python3.12/bin/python
"""统计各分组的非零基因集合, 绘制分组间基因重叠韦恩图。

输入:
    gene_abundance_table.tsv (行=基因, 列=样本, tab 分隔, 0 表示缺失)
    metadata.tsv             (两列: 样本ID<TAB>分组)
输出:
    gene_venn.png (分组基因集合韦恩图, 支持 2 或 3 个分组)
    gene_venn.csv (各组基因数及各交集区域计数明细)
依赖:
    matplotlib-venn

出错时不会中断 snakemake: 图片与 CSV 仍会生成, 其中记录报错信息。
"""

import sys
import traceback

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn2, venn3

sys.stderr = open(snakemake.log[0], "w")


# IO
in_file = snakemake.input["abund"]
metadata_file = snakemake.input["meta"]
out_path = snakemake.output["png"]
csv_path = snakemake.output["csv"]


def fail_gracefully(msg):
    """出错时也生成图片与 CSV, 以文字记录报错信息, 保证 snakemake 不中断。"""
    sys.stderr.write("\n" + msg)
    plt.close("all")
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.axis("off")
    ax.set_title("Gene Venn - Error", fontsize=14, y=1.02)
    ax.text(0.02, 0.98, msg, transform=ax.transAxes,
            fontsize=9, va="top", ha="left", family="monospace")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    # 同样生成 CSV, 确保 output 文件全部存在
    pd.DataFrame({"region": ["ERROR"], "gene_count": [msg]}).to_csv(
        csv_path, index=False, sep="\t")


try:
    # MAIN
    # 1. 读取基因丰度表, 每个样本取其非零(>0)基因作为基因集合
    df = pd.read_csv(in_file, sep="\t", index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    sys.stderr.write(f"基因数: {df.shape[0]}, 样本数: {df.shape[1]}")

    sample_genes = {col: set(df.index[df[col] > 0]) for col in df.columns}

    # 2. 读取 metadata, 按分组聚合组内样本基因并集
    meta = pd.read_csv(metadata_file, sep="\t")
    meta = meta[meta["sample"].notna()]  # 去掉末尾空行
    missing = meta[~meta["sample"].isin(sample_genes)]
    if not missing.empty:
        sys.stderr.write(f"警告: metadata 中以下样本在丰度表中不存在: {list(missing['sample'])}")
    valid = meta[meta["sample"].isin(sample_genes)]
    sys.stderr.write("\n各分组样本数:")
    sys.stderr.write(valid.groupby("group")["sample"].count().to_string())

    group_genes = {}
    for grp, sub in valid.groupby("group"):
        genes = set()
        for s in sub["sample"]:
            # genes = genes | sample_genes[s], 新的基因更新到集合, 保持非重复
            genes |= sample_genes[s]
        group_genes[grp] = genes

    # 3. 参与绘制的分组: metadata 中的全部分组 (最多 3 组)
    sel = list(group_genes.keys())
    if len(sel) > 3:
        sys.stderr.write(f"\n警告: 分组数 {len(sel)} > 3, 韦恩图仅绘制前 3 组: {sel[:3]}")
        sel = sel[:3]
    if len(sel) < 2:
        fail_gracefully(f"错误: 韦恩图至少需要 2 个分组, 当前仅 {len(sel)} 组!")
        sys.exit(0)  # 不中断 snakemake

    sys.stderr.write("\n参与绘制的分组及其非零基因数:")
    sys.stderr.write(pd.Series({g: len(group_genes[g]) for g in sel}).to_string())

    # 4. 绘制韦恩图
    sets = [group_genes[g] for g in sel]
    fig, ax = plt.subplots(figsize=(8, 6))
    if len(sel) == 2:
        venn2(subsets=sets, set_labels=list(sel), ax=ax)
    else:
        venn3(subsets=sets, set_labels=list(sel), ax=ax)
    ax.set_title("Gene Set Overlap between Groups", fontsize=13, y=1.05)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    sys.stderr.write(f"\n韦恩图已保存至: {out_path}")

    # 5. 输出各区域交集计数 CSV
    a, b = sets
    if len(sel) == 2:
        regions = {
            f"{sel[0]} only": len(a - b),
            f"{sel[1]} only": len(b - a),
            f"{sel[0]}∩{sel[1]}": len(a & b),
        }
    else:
        c = sets[2]
        regions = {
            f"{sel[0]} only": len(a - b - c),
            f"{sel[1]} only": len(b - a - c),
            f"{sel[2]} only": len(c - a - b),
            f"{sel[0]}∩{sel[1]}": len((a & b) - c),
            f"{sel[0]}∩{sel[2]}": len((a & c) - b),
            f"{sel[1]}∩{sel[2]}": len((b & c) - a),
            f"{sel[0]}∩{sel[1]}∩{sel[2]}": len(a & b & c),
        }

    out = pd.DataFrame({
        "region": list(regions.keys()),
        "gene_count": list(regions.values()),
    })
    out.to_csv(csv_path, index=False, sep="\t")
    sys.stderr.write(f"\n区域统计已保存至: {csv_path}")
    sys.stderr.write("\n各区域基因数:")
    sys.stderr.write(out.to_string(index=False))
except Exception:
    # 任何报错: 生成记录错误信息的图片与 CSV, 不中断 snakemake
    fail_gracefully(traceback.format_exc())
    sys.exit(0)
