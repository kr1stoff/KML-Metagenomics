import pandas as pd
from functools import reduce
from pathlib import Path
import sys

sys.stderr = open(snakemake.log[0], "w")

# 输入输出文件
input_files = snakemake.input
output_file = snakemake.output[0]

# 分析
dfs = []
for input_file in input_files:
    sample = Path(input_file).stem.split('.')[0]
    df = pd.read_csv(input_file, sep="\t", usecols=['gene_id', 'norm_abund'])
    df.rename(columns={'norm_abund': sample}, inplace=True)
    dfs.append(df)
df_all = reduce(lambda left, right: pd.merge(left, right, on='gene_id', how='outer'), dfs).fillna(0)
df_all.to_csv(output_file, sep="\t", index=False)
