import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

df_fastp = pd.read_table(snakemake.input.fastp, sep="\t")
df_nonhost = pd.read_table(snakemake.input.nonhost, sep="\t")
df_merged = pd.merge(df_fastp, df_nonhost, on="Sample", how="left")
df_merged.to_excel(snakemake.output[0], index=False)
