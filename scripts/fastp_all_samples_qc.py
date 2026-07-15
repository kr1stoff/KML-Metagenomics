# fastp 数据转表格

import numpy as np
import pandas as pd
import json
from pathlib import Path
import sys

# logging
sys.stderr = open(snakemake.log[0], "w")


def fastp_all_samples_qc(files_fastp_json) -> pd.DataFrame:
    title = ["Sample", "Clean_Reads", "Total_Base", "Q20",
             "Q30", "Q20_Rate", "Q30_Rate", "Average_Length", "GC"]
    df = pd.DataFrame(columns=title)
    for js_path in files_fastp_json:
        js_data = json.loads(open(js_path, "r").read())
        sample = Path(js_path).stem
        mean_lengths = np.array(
            [v for k, v in js_data["summary"]["after_filtering"].items() if k.endswith("mean_length")])
        out = [
            sample,
            js_data["summary"]["after_filtering"]["total_reads"],
            js_data["summary"]["after_filtering"]["total_bases"],
            js_data["summary"]["after_filtering"]["q20_bases"],
            js_data["summary"]["after_filtering"]["q30_bases"],
            js_data["summary"]["after_filtering"]["q20_rate"],
            js_data["summary"]["after_filtering"]["q30_rate"],
            mean_lengths.mean(),
            js_data["summary"]["after_filtering"]["gc_content"],
        ]
        df.loc[len(df)] = out
    return df


# main
df = fastp_all_samples_qc(snakemake.input)
df.to_csv(snakemake.output['tsv'], index=False, sep="\t")
df.to_excel(snakemake.output['excel'], index=False)
