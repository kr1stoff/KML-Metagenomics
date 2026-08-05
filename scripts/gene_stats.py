#!/home/mengxf/miniforge3/envs/python3.12/bin/python
"""统计 FNA 基因序列的基本指标。"""

import sys
import pandas as pd
import numpy as np
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

sys.stderr = open(snakemake.log[0], "w")

# 标准细菌启动子密码子 (start codon) 和终止密码子 (stop codon)
START_CODONS = {"ATG", "GTG", "TTG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}


def parse_fna(filepath):
    """解析 FNA 文件，返回 (header, seq) 列表。"""
    seqs = []
    current_header = None
    current_seq = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    seqs.append((current_header, "".join(current_seq)))
                current_header = line[1:]  # 去掉 >
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_header is not None:
            seqs.append((current_header, "".join(current_seq)))
    return seqs


def has_start(seq):
    """序列前3位是否为起始密码子"""
    return seq[:3] in START_CODONS


def has_stop(seq):
    """序列末3位是否为终止密码子"""
    return seq[-3:] in STOP_CODONS


def main():
    fna_path = snakemake.input[0]
    seqs = parse_fna(fna_path)

    total_len = 0
    gc_count = 0
    count_start_only = 0
    count_stop_only = 0
    count_none = 0
    count_all = 0
    lengths = []

    for _, seq in seqs:
        total_len += len(seq)
        gc_count += seq.count("G") + seq.count("C")
        lengths.append(len(seq))

        s = has_start(seq)
        e = has_stop(seq)

        if s and e:
            count_all += 1
        elif s and not e:
            count_start_only += 1
        elif not s and e:
            count_stop_only += 1
        else:
            count_none += 1

    n = len(seqs)
    avg_len = total_len / n if n > 0 else 0
    gc_percent = (gc_count / total_len * 100) if total_len > 0 else 0

    pct_start = count_start_only / n * 100 if n > 0 else 0
    pct_end   = count_stop_only / n * 100 if n > 0 else 0
    pct_none  = count_none / n * 100 if n > 0 else 0
    pct_all   = count_all / n * 100 if n > 0 else 0

    # 输出 Excel
    sample_name = Path(fna_path).stem
    data = {
        "Sample": [sample_name],
        "ORFs_NO.": [n],
        "Integrity_start": [f"{count_start_only} ({pct_start:.2f}%)"],
        "Integrity_end": [f"{count_stop_only} ({pct_end:.2f}%)"],
        "Integrity_none": [f"{count_none} ({pct_none:.2f}%)"],
        "Integrity_all": [f"{count_all} ({pct_all:.2f}%)"],
        "Total_Len.(Mbp)": [round(total_len / 1e6, 4)],
        "Average_Len.(bp)": [round(avg_len, 2)],
        "GC_Percent": [round(gc_percent, 2)],
    }
    df_out = pd.DataFrame(data)
    xlsx_path = snakemake.output.xlsx
    df_out.to_excel(xlsx_path, index=False)
    sys.stderr.write(f"\nExcel 已保存至: {xlsx_path}")

    sys.stderr.write(f"{'指标':<30} {'数值':>25}")
    sys.stderr.write("-" * 57)
    sys.stderr.write(f"{'ORFs_NO.':<30} {n:>25,}")
    sys.stderr.write(f"{'Integrity_start':<30} {count_start_only:>15,} ({pct_start:.2f}%)")
    sys.stderr.write(f"{'Integrity_end':<30}   {count_stop_only:>15,} ({pct_end:.2f}%)")
    sys.stderr.write(f"{'Integrity_none':<30}  {count_none:>15,} ({pct_none:.2f}%)")
    sys.stderr.write(f"{'Integrity_all':<30}   {count_all:>15,} ({pct_all:.2f}%)")
    sys.stderr.write(f"{'Total_Len.(Mbp)':<30} {total_len / 1e6:>15.4f}")
    sys.stderr.write(f"{'Average_Len.(bp)':<30} {avg_len:>15.2f}")
    sys.stderr.write(f"{'GC_Percent':<30}       {gc_percent:>15.2f}%")

    # ============================================================
    # 基因长度分布图
    # ============================================================
    # x 轴上限定为 75 分位长度，超出部分归入最后一个 bin
    # cap = np.percentile(lengths, 75)
    # 直接设置 3000 吧
    cap = 3000
    clipped = np.clip(lengths, None, cap)

    n_bins = 50
    counts, bin_edges = np.histogram(clipped, bins=n_bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_width = bin_edges[1] - bin_edges[0]
    pct = counts / counts.sum() * 100

    fig, ax1 = plt.subplots(figsize=(12, 6))

    # 主 y 轴: 条形图 (序列数目)
    ax1.bar(bin_centers, counts, width=bin_width * 0.9,
            color="steelblue", alpha=0.75, edgecolor="white")
    ax1.set_xlabel("Gene Length (bp)", fontsize=12)
    ax1.set_ylabel("Sequence Count", fontsize=12, color="steelblue")
    ax1.tick_params(axis="y", labelcolor="steelblue")

    # x 轴刻度: 取 11 个均匀刻度，最后一个标签标注为 "≥{cap}"
    x_ticks = np.linspace(0, cap, 11)
    x_labels = [f"{int(v):,}" for v in x_ticks]
    x_labels[-1] = f"≥ {int(cap):,}"
    ax1.set_xticks(x_ticks)
    ax1.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=9)
    ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

    # 副 y 轴: 折线图 (百分比)
    ax2 = ax1.twinx()
    ax2.plot(bin_centers, pct, color="red", marker="o", linewidth=1.5,
             markersize=3)
    ax2.set_ylabel("Percentage (%)", fontsize=12, color="red")
    ax2.tick_params(axis="y", labelcolor="red")

    # 主标题
    ax1.set_title(f"Gene Length Distribution — {sample_name}", fontsize=14,
                  fontweight="bold")

    fig.tight_layout()
    chart_path = snakemake.output.png
    fig.savefig(chart_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    sys.stderr.write(f"长度分布图已保存至: {chart_path}")


if __name__ == "__main__":
    main()
