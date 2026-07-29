# 三篇文献基因丰度计算公式汇总

## 1. Buchfink et al. 2014 — DIAMOND 论文

**文献信息**: *Fast and sensitive protein alignment using DIAMOND*, Nature Methods, 2015

**与丰度计算的关系**: 本文主要介绍 DIAMOND 序列比对算法（替代 BLASTX），**本身不直接定义基因丰度计算公式**。它的核心贡献是提供上游比对工具，将测序 reads 比对到蛋白参考数据库（如 NCBI-nr、KEGG），为后续丰度定量提供比对结果。该工具比对速度是 BLASTX 的 ~20,000 倍，是宏基因组丰度计算的**上游比对步骤**的基础工具。

---

## 2. Cotillard et al. 2013 — 肠道微生物基因丰度

**文献信息**: *Dietary intervention impact on gut microbial gene richness*, Nature, 2013

**使用的工具**: METEOR 软件（定量宏基因组分析平台）

**基因丰度计算分两步**：

**第1步 — 基因长度归一化**：
```
Normalized_Abundance_i = Uniquely mapped reads_i / Gene length_i
```

**第2步 — 转化为频率**（分母为总 reads 数，非标准化丰度之和）：
```
Gene_Frequency_{i,s} = Normalized_Abundance_{i,s} / Total_uniquely_mapped_reads_s
                      = (reads_{i,s} / length_i) / Total_reads_s
                      = reads_{i,s} / (length_i * Total_reads_s)
```
> 其中 s 表示样本，Total_reads_s = Σ_j reads_{j,s}（该样本内所有基因的 reads 之和）。
> **分母是样本维度的总 reads 数**（"for a given sample"），而非跨样本求和。每个样本独立归一化，消除样本间测序深度差异。
> 因此这些"频率"严格意义上**不保证总和为 1**，与 RPKM/10⁹ 等价。

**物种/基因簇丰度**（后聚合）：
```
Cluster_Abundance_cluster = sum_{i in cluster} Gene_Frequency_i
```

原文依据（Methods 部分）：
> "Abundance of each gene in an individual was normalized with METEOR by dividing the number of reads that uniquely mapped to a gene by its nucleotide length. After that, normalized gene abundances were transformed in frequencies by dividing them with the total number of uniquely mapped reads for a given sample."

> "The group abundance of each cluster was computed as the sum of the frequencies of its genes."

**过滤标准**: 仅保留 unique mapping reads，最多允许 3 个 mismatch。

---

## 3. Villar et al. 2015 — 海洋浮游宏基因组

**文献信息**: *Environmental characteristics of Agulhas rings affect interocean plankton transport*, Science, 2015

**基因丰度计算公式**: 标准 **RPKM**（Reads Per Kilobase per Million mapped reads）

```
RPKM_i = (Reads mapped to gene_i * 10^9) / (Gene_length_i * Total mapped reads)
```

**KO 功能组丰度**（后聚合）：
```
KO_Abundance_KO = sum_{i in KO group} RPKM_i
```

原文依据（Materials and methods）：
> "Gene abundances were computed for the set of genes... by counting the number of reads from each sample that mapped to each KO-associated gene. Abundances were normalized as reads per kilobase per million mapped reads (RPKM). Gene abundances were then aggregated (summed) for each KO group."

**方法细节**: 基于 OM-RGC（Ocean Microbial Reference Gene Catalog）参考基因集，使用 BLAST 进行比对鉴定。

---

## 三文公式对比总结

| 特征 | Cotillard 2013 | Villar 2015 |
|------|----------------|-------------|
| **方法名称** | METEOR（两步归一化） | RPKM |
| **核心公式** | `reads / (length * total_reads)` | `reads * 10^9 / (length * total_reads)` |
| **最终量纲** | 与 RPKM/10⁹ 等价 | RPKM 值 |
| **基因簇聚合** | 基因频率求和 | RPKM 求和 |
| **比对工具** | SOLiD / corona_lite | BLAST |
| **参考基因集** | IGC 基因目录 (3.3M genes) | OM-RGC 基因目录 |
| **过滤策略** | 仅 unique mapped (<=3 mismatches) | - |
| **Buchfink 2014** 角色 | —（DIAMOND 为上游比对工具，可替代 BLAST/BLASTX） | — |

两种方法**数学上仅相差常数因子 10⁹**：Cotillard 的 Gene Frequency = RPKM / 10⁹。均为 reads_i/(length_i × total_reads) 的线性变换，仅系数不同。
