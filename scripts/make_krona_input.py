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
    2. 基因分类明细: lineage(去前缀)  level  gene_count  norm_abund
       按 lineage 合并, norm_abund 求和, gene_count 为该 lineage 包含的预测基因数
"""

import pandas as pd
import math
import re
import sys

sys.stderr = open(snakemake.log[0], "w")


# MEGAN_FILE = "/data/mengxf/Develop/KML260617-MetaGenomics/results/260722/taxon_classification/diamond_nr_miro.daa.meganized.info"
# ABUND_FILE = "/data/mengxf/Develop/KML260617-MetaGenomics/results/260722/gene_quantification/SRR23604277.sample_gene_abundance.tsv"
# KRONA_OUT = "/data/mengxf/Develop/KML260617-MetaGenomics/work/260819-classification/SRR23604277.krona_input.txt"
# DETAIL_OUT = "/data/mengxf/Develop/KML260617-MetaGenomics/work/260819-classification/SRR23604277.gene_taxonomy.tsv"
MEGAN_FILE = snakemake.input['megan']
ABUND_FILE = snakemake.input['abund']
KRONA_OUT = snakemake.output['krona']
DETAIL_OUT = snakemake.output['detail']


# 病毒界 (kingdom) 名称, 统一替换为 Viruses
VIRUS_KINGDOMS = {
    'Abadenavirae', 'Bamfordvirae', 'Helvetiavirae', 'Heunggongvirae',
    'Loebvirae', 'Orthornavirae', 'Pararnavirae', 'Sangervirae',
    'Shotokuvirae', 'Trapavirae', 'Zilligvirae',
}
# 域 [D] 归为界 [K] 的生物类群 (细菌/古菌)
DOMAIN_TO_KINGDOM = {'Bacteria', 'Archaea'}
# 七级分类级别顺序
LEVEL_ORDER = ['K', 'P', 'C', 'O', 'F', 'G', 'S']


def parse_lineage(raw: str) -> list:
    """解析 '[D] Bacteria; [P] Firmicutes; ...' 格式, 返回 [(级别, 名称), ...]"""
    items = []
    for part in str(raw).split(';'):
        # 去掉首尾空白, 跳过空段 (行尾分号产生的空段)
        part = part.strip()
        if not part:
            continue
        # 匹配 '[X] name' 形式, 级别取 D/K/P/C/O/F/G/S
        m = re.match(r'^\[([DKPCOFGS])\]\s*(.+)$', part)
        if m:
            items.append((m.group(1), m.group(2)))
    return items


def format_lineage(parsed: list) -> list:
    """统一格式化 lineage: 域[D]归界[K], 病毒界替换为 Viruses, 返回 [(级别, 名称), ...]"""
    # 无有效分类 (未注释或无法识别的 Others) 统一标记为 Others
    if not parsed:
        return [('K', 'Others')]
    levels = list(parsed)
    top_lv, top_nm = levels[0]
    if top_lv == 'D' and top_nm in DOMAIN_TO_KINGDOM:
        # 细菌/古菌: 域 [D] 归为界 [K]
        levels[0] = ('K', top_nm)
    elif len(levels) > 1 and levels[1][0] == 'K' and levels[1][1] in VIRUS_KINGDOMS:
        # 病毒: '[D] Others; [K] <virae>; [P] ...' 第二层界替换为 Viruses
        levels = [('K', 'Viruses')] + levels[2:]
    elif top_lv == 'K' and top_nm in VIRUS_KINGDOMS:
        # 病毒直接以 [K] <virae> 开头时同样替换为 Viruses
        levels = [('K', 'Viruses')] + levels[1:]
    elif top_lv == 'D' and top_nm == 'Others':
        # 其他以 Others 为顶层的行归为未分类
        return [('K', 'Others')]
    return levels


def load_megan(path: str) -> pd.DataFrame:
    """读取 megan6 分类结果, 生成 gene_id / lineage_clean(统一格式化后的列表)"""
    megan = pd.read_csv(path, sep='\t', comment='#', header=None,
                        names=['gene_id', 'level_raw', 'lineage_raw'], dtype=str)
    # 解析并统一格式化 lineage, 得到 [(级别, 名称), ...]
    megan['lineage_clean'] = megan['lineage_raw'].apply(
        lambda x: format_lineage(parse_lineage(x)) if pd.notna(x) else [('K', 'Others')])
    return megan


def load_abundance(path: str) -> pd.DataFrame:
    """读取样本基因丰度表, 返回 gene_id/gene_len/reads/norm_abund"""
    abund = pd.read_csv(path, sep='\t', dtype={'gene_id': str})
    return abund


def main() -> None:
    """主流程: 合并分类与丰度, 聚合后写出 krona 输入和分类明细文件"""
    # 读取两个输入表
    # megan6 分类信息是所有样本
    megan = load_megan(MEGAN_FILE)
    # 是当前样本经过过滤后的, 基因数目差很多
    abund = load_abundance(ABUND_FILE)
    # 以丰度表为基准关联 megan6 分类信息
    merged = abund.merge(megan, on='gene_id', how='left')
    # 未关联到分类的基因补默认值 (NaN 统一替换为 Others)
    merged['lineage_clean'] = merged['lineage_clean'].apply(
        lambda x: x if isinstance(x, list) else [('K', 'Others')])

    # ---- krona 输入: 相同 lineage 的 norm_abund 加和合并 ----
    # 取统一格式化后的名称序列作为 krona 层级
    merged['lineage_names'] = merged['lineage_clean'].apply(
        lambda x: tuple(nm for _, nm in x))
    # 按层级元组分组, 聚合 norm_abund 求和
    grouped = merged.groupby(merged['lineage_names'])['norm_abund'].sum()
    krona_lines = ['\t'.join([f'{math.ceil(val)}', *names])
                    for names, val in grouped.items()]
    with open(KRONA_OUT, 'w') as f:
        f.write('\n'.join(krona_lines) + '\n')

    # ---- 分类明细: 每个层级节点一行, 含直接/累加基因数与丰度 ----
    merged['detail_path'] = merged['lineage_clean'].apply(
        lambda x: x if x == [('K', 'Others')] else x)
    merged['detail_text'] = merged['detail_path'].apply(
        lambda x: '; '.join(nm for _, nm in x))
    # 直接注释: 每个基因只在最深注释节点计一次 (按补全后的完整路径分组)
    direct = merged.groupby('detail_text', as_index=False).agg(
        gene_direct=('gene_id', 'count'),
        abund_direct=('norm_abund', 'sum'),
    )
    # 累加: 每个基因展开到其全部祖先节点, 节点累加值 = 该节点及其所有后代直接值之和
    node_rows = []  # 每元素: (节点完整路径, 节点级别, 该基因的 norm_abund)
    for _, row in merged.iterrows():
        path = []  # 累积前缀名称
        for lv, nm in row['detail_path']:
            path.append(nm)
            node_rows.append(('; '.join(path), lv, row['norm_abund']))
    node_df = pd.DataFrame(node_rows, columns=['lineage', 'level', 'norm_abund'])
    # 按节点聚合: 基因数用 count, 丰度用 sum, level 取首个 (同一节点级别一致)
    cum = node_df.groupby('lineage', as_index=False).agg(
        level=('level', 'first'),
        gene_cum=('norm_abund', 'count'),
        abund_cum=('norm_abund', 'sum'),
    )
    # 合并直接值与累加值, 补全缺失的直接值 (纯祖先节点无直接注释) 为 0
    detail = cum.merge(direct.rename(columns={'detail_text': 'lineage'}),
                        on='lineage', how='left')
    detail['gene_direct'] = detail['gene_direct'].fillna(0).astype(int)
    detail['abund_direct'] = detail['abund_direct'].fillna(0.0)
    detail.to_csv(DETAIL_OUT, sep='\t', index=False)

    # 打印统计信息
    n_classified = int(merged['lineage_names'].apply(lambda x: x != ('Others',)).sum())
    sys.stderr.write(f"丰度表基因总数: {len(abund)}\n")
    sys.stderr.write(f"关联到有效分类的基因数: {n_classified}\n")
    sys.stderr.write(f"未分类基因数: {len(abund) - n_classified}\n")
    sys.stderr.write(f"聚合后 krona 层级数: {len(grouped)}\n")
    sys.stderr.write(f"明细层级节点数: {len(detail)}\n")
    sys.stderr.write(f"krona 输入已写出: {KRONA_OUT}\n")
    sys.stderr.write(f"分类明细已写出: {DETAIL_OUT}\n")


if __name__ == '__main__':
    main()
