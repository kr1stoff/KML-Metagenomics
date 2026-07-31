from pathlib import Path
import pandas as pd
import json
import sys

# 设置日志输出文件
sys.stderr = open(snakemake.log[0], "w")

# 读取输入文件
input_files = snakemake.input
output_reads_table = snakemake.output["reads"]
output_gene_length = snakemake.output["glen"]

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

# 输出样本-基因reads计数表
df = pd.DataFrame.from_dict(sample_gene_reads, orient='index').fillna(0).astype(int).T
df.to_csv(output_reads_table, sep='\t')

# 输出基因长度统计
with open(output_gene_length, 'w') as f:
    json.dump(gene_len_dict, f, indent=4)
