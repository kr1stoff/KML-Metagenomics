#!/home/mengxf/miniforge3/envs/python3.12/bin/python
"""统计 FNA 基因序列的基本指标。"""

import sys
import pandas as pd
from pathlib import Path

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
    fna_path = sys.argv[1] if len(sys.argv) > 1 else "gene_catalogue.unigene.fna.demo"
    seqs = parse_fna(fna_path)

    total_len = 0
    gc_count = 0
    count_start_only = 0
    count_stop_only = 0
    count_none = 0
    count_all = 0

    for _, seq in seqs:
        total_len += len(seq)
        gc_count += seq.count("G") + seq.count("C")

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
    xlsx_path = Path(fna_path).with_suffix(".stats.xlsx")
    df_out.to_excel(xlsx_path, index=False)
    print(f"\nExcel 已保存至: {xlsx_path}")

    print(f"{'指标':<30} {'数值':>25}")
    print("-" * 57)
    print(f"{'ORFs_NO.':<30} {n:>25,}")
    print(f"{'Integrity_start':<30} {count_start_only:>15,} ({pct_start:.2f}%)")
    print(f"{'Integrity_end':<30}   {count_stop_only:>15,} ({pct_end:.2f}%)")
    print(f"{'Integrity_none':<30}  {count_none:>15,} ({pct_none:.2f}%)")
    print(f"{'Integrity_all':<30}   {count_all:>15,} ({pct_all:.2f}%)")
    print(f"{'Total_Len.(Mbp)':<30} {total_len / 1e6:>15.4f}")
    print(f"{'Average_Len.(bp)':<30} {avg_len:>15.2f}")
    print(f"{'GC_Percent':<30}       {gc_percent:>15.2f}%")


if __name__ == "__main__":
    main()
