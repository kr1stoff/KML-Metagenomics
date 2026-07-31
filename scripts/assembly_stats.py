import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

# 输入输出
# reports = glob('/data/mengxf/Develop/KML260617-MetaGenomics/results/260722/assembly/quast/*/report.tsv')
reports = snakemake.input
output = snakemake.output[0]

# 分析
dfs = []

for report in reports:
    # Total length, # contigs (>= 0 bp), N50, N90, Largest contig
    df = pd.read_csv(report, sep='\t', index_col=0)
    indices = ['Total length', '# contigs (>= 0 bp)', 'N50', 'N90', 'Largest contig']
    df = df.loc[indices]
    # quast 不计算平均长度
    df.loc['Average len'] = df.loc['Total length'] / df.loc['# contigs (>= 0 bp)']
    df = df.astype(int).T
    df.rename(columns={'# contigs (>= 0 bp)': 'Contigs num'}, inplace=True)
    dfs.append(df)

df_concat = pd.concat(dfs, axis=0)
df_concat.to_excel(output, index=True)
