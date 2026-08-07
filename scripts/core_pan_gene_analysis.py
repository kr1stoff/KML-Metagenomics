# ChatGPT 20260805 Core-Pan 基因分析

import random
import numpy as np
import pandas as pd

presence = pd.read_csv(
    "gene_presence.tsv",
    sep="\t",
    index_col=0
).astype(bool)

samples = presence.columns

repeat = 100

core_curve = []
pan_curve = []

for n in range(1, len(samples)+1):

    core_list = []
    pan_list = []

    for _ in range(repeat):

        selected = random.sample(list(samples), n)

        subset = presence[selected]

        core = subset.all(axis=1).sum()

        pan = subset.any(axis=1).sum()

        core_list.append(core)
        pan_list.append(pan)

    core_curve.append(np.mean(core_list))
    pan_curve.append(np.mean(pan_list))


import matplotlib.pyplot as plt


fig, axes = plt.subplots(
    nrows=1,
    ncols=2,
    figsize=(10,4),
    dpi=300
)


# Core gene
axes[0].plot(
    x,
    core_curve,
    marker="o"
)

axes[0].set_xlabel(
    "Number of samples"
)

axes[0].set_ylabel(
    "Number of core genes"
)

axes[0].set_title(
    "Core gene accumulation curve"
)


# Pan gene
axes[1].plot(
    x,
    pan_curve,
    marker="o"
)

axes[1].set_xlabel(
    "Number of samples"
)

axes[1].set_ylabel(
    "Number of pan genes"
)

axes[1].set_title(
    "Pan gene accumulation curve"
)


plt.tight_layout()

plt.savefig(
    "core_pan_curve.pdf",
    bbox_inches="tight"
)

plt.show()
