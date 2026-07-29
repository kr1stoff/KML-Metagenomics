from pathlib import Path
import pandas as pd
import json
import sys

# 设置日志输出文件
sys.stderr = open(snakemake.log[0], "w")

# 读取输入文件
input_files = snakemake.input
output_reads_table = snakemake.output[0]

# 分析
sample_gene_reads = {}
gene_len_dict = {}

for idx_file in input_files:
    with open(idx_file, 'r') as f:
        for line in f:
            sample = Path(idx_file).stem.split('.')[0]
            gene_id, gene_len, reads = line.strip().split('\t')[:3]
            sample_gene_reads.setdefault(sample, {})[gene_id] = int(reads)
            gene_len_dict[gene_id] = int(gene_len)

df = pd.DataFrame.from_dict(sample_gene_reads, orient='index').fillna(0).astype(int).T
df.insert(0, "gene_length", df.index.map(gene_len_dict))
df.reset_index(inplace=True, names=['gene_id'])
df.to_csv(output_reads_table, sep='\t', index=False)
