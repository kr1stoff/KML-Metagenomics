#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""根据 megan6 分类结果和样本基因丰度表生成 krona 输入及基因分类明细文件。

输入:
    1. megan6 分类结果 (test_out): 跳过 # 注释行,
       第1列基因名, 第2列注释级别, 第3列物种 lineage (去掉前面 NCBI; 等前缀)
    2. 样本基因丰度表: 第1列 gene_id, 第2列 gene_len, 第3列 reads, 第4列 norm_abund
输出 (均写入当前工作目录, 不写入输入文件所在目录):
    1. krona 输入文件 (ktImportText 格式):
       每行 数值(TAB)层级1(TAB)层级2..., 相同 lineage 的 norm_abund 加和合并
    2. 基因分类明细: lineage(去前缀)  level  gene_len  reads  norm_abund
"""

import pandas as pd
import math
import sys

sys.stderr = open(snakemake.log[0], "w")


# MEGAN_FILE = "/data/mengxf/Develop/KML260617-MetaGenomics/results/260722/test_out"
# ABUND_FILE = "/data/mengxf/Develop/KML260617-MetaGenomics/results/260722/gene_quantification/SRR23604277.sample_gene_abundance.tsv"
# KRONA_OUT = "/data/mengxf/Develop/KML260617-MetaGenomics/work/260819-classification/SRR23604277.krona_input.txt"
# DETAIL_OUT = "/data/mengxf/Develop/KML260617-MetaGenomics/work/260819-classification/SRR23604277.gene_taxonomy.tsv"
MEGAN_FILE = snakemake.input['megan']
ABUND_FILE = snakemake.input['abund']
KRONA_OUT = snakemake.output['krona']
DETAIL_OUT = snakemake.output['detail']


def clean_lineage(lineage: str) -> list:
    """去掉 NCBI; cellular organisms; 等固定前缀, 返回层级列表"""
    parts = [p.strip() for p in str(lineage).split(';') if p.strip()]
    # 过滤 NCBI 分类体系的固定根节点
    parts = [p for p in parts if p not in ('NCBI', 'cellular organisms')]
    # 无有效分类时统一标记为 Unclassified
    return parts if parts else ['Unclassified']


def load_megan(path: str) -> pd.DataFrame:
    """读取 megan6 分类结果, 跳过 # 注释行, 生成 gene_id/level/lineage_clean"""
    megan = pd.read_csv(path, sep='\t', comment='#', header=None,
                        names=['gene_id', 'level', 'lineage_raw'], dtype=str)
    # 对 lineage 做前缀清理, 得到层级列表
    megan['lineage_clean'] = megan['lineage_raw'].apply(clean_lineage)
    return megan


def load_abundance(path: str) -> pd.DataFrame:
    """读取样本基因丰度表, 返回 gene_id/gene_len/reads/norm_abund"""
    abund = pd.read_csv(path, sep='\t', dtype={'gene_id': str})
    return abund


def main() -> None:
    """主流程: 合并分类与丰度, 聚合后写出 krona 输入和分类明细文件"""
    # 读取两个输入表
    megan = load_megan(MEGAN_FILE)
    abund = load_abundance(ABUND_FILE)
    # 以丰度表为基准关联 megan6 分类信息
    merged = abund.merge(megan, on='gene_id', how='left')
    # 未关联到分类的基因补默认值 (NaN 统一替换为 ['Unclassified'])
    merged['lineage_clean'] = merged['lineage_clean'].apply(
        lambda x: x if isinstance(x, list) else ['Unclassified'])

    # ---- krona 输入: 相同 lineage 的 norm_abund 加和合并 ----
    # 用 lineage 层级元组作为分组键, 聚合 norm_abund 求和
    grouped = merged.groupby(merged['lineage_clean'].apply(tuple))['norm_abund'].sum()
    krona_lines = ['\t'.join([f'{math.ceil(val)}', *lineage])
                   for lineage, val in grouped.items()]
    with open(KRONA_OUT, 'w') as f:
        f.write('\n'.join(krona_lines) + '\n')

    # ---- 分类明细: lineage(去前缀)  level  gene_len  reads  norm_abund ----
    merged['lineage_text'] = merged['lineage_clean'].apply(lambda x: '; '.join(x))
    detail_df: pd.DataFrame = merged[['lineage_text', 'level',
                                      'gene_len', 'reads', 'norm_abund']]
    detail_df.to_csv(DETAIL_OUT, sep='\t',
                     header=['lineage', 'level', 'gene_len', 'reads', 'norm_abund'],
                     index=False)

    # 打印统计信息
    n_classified = int((merged['lineage_clean'].apply(lambda x: x != ['Unclassified'])).sum())
    print(f"丰度表基因总数: {len(abund)}")
    print(f"关联到有效分类的基因数: {n_classified}")
    print(f"未分类基因数: {len(abund) - n_classified}")
    print(f"聚合后 krona 层级数: {len(grouped)}")
    print(f"krona 输入已写出: {KRONA_OUT}")
    print(f"分类明细已写出: {DETAIL_OUT}")


if __name__ == '__main__':
    main()
