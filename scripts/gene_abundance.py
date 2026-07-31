import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

# 输入输出文件
input_file = snakemake.input[0]
output_file = snakemake.output[0]

# 分析
df = pd.read_csv(input_file, sep="\t", usecols=range(3))
df.columns = ['gene_id', 'gene_len', 'reads']
# 计算归一化丰度. 基因reads数/基因长度/总reads数
df['norm_abund'] = df['reads'] / df['gene_len'] / df['reads'].sum() * 10**9
df.to_csv(output_file, sep="\t", index=False)
