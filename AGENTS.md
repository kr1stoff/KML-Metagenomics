# KML-Metagenomics

宏基因组分析流程：双端测序数据的质控过滤、去宿主、组装、基因预测、基因去冗余与丰度定量。

## 分析流程

```
fastp → fastqc ─┐
                 ├→ multiqc
                 │
fastp ───────────┘
  ↓
bowtie2 (去宿主) → gunzip → megahit → seqtk (>500bp)
                                          ↓
                  MetaGeneMark → seqtk (>100bp) → cd-hit
                                                      ↓
                                    bowtie2-build → bowtie2 (比对) → samtools sort
```

| 步骤 | 工具 | 说明 |
|------|------|------|
| 质控过滤 | fastp | 接头/质量过滤，输出清洗后 FASTQ 及 HTML/JSON 报告 |
| 质控报告 | fastqc | 对清洗后 FASTQ 生成质量报告 |
| 汇总报告 | multiqc | 汇总 fastp + fastqc 报告 |
| 去宿主 | bowtie2 | 比对到宿主参考 (hg19)，提取未比对 reads |
| 组装 | megahit | 宏基因组组装 (`meta-large` 预设) |
| 过滤 | seqtk | 保留 > 500bp 的 contig |
| 基因预测 | MetaGeneMark | 预测基因，输出 GFF / 蛋白 / 核酸序列 |
| 过滤 | seqtk | 保留 > 100bp 的基因序列 |
| 去冗余 | cd-hit | 95% 相似度聚簇去冗余 |
| 定量 | bowtie2 + samtools | 比对 reads 到基因集，输出排序 BAM |

## 目录结构

```
KML-Metagenomics/
├── Snakefile              # 主规则文件
├── config.yaml            # 运行配置
├── config.schema.yaml     # 配置校验 schema
├── rules/
│   ├── qc.smk             # fastp / fastqc / multiqc
│   ├── host_removal.smk   # bowtie2 去宿主 + gunzip
│   ├── assembly.smk       # megahit + seqtk
│   ├── gene_prediction.smk  # MetaGeneMark + seqtk
│   ├── gene_clustering.smk   # cd-hit
│   └── quantification.smk    # bowtie2 比对 + samtools
└── tests/
    ├── input.tsv          # 样本输入表 (样本ID, R1, R2)
    └── pipe.sh            # 原始分析脚本
```

## 环境要求

需要预装 [Snakemake](https://snakemake.readthedocs.io/) (≥ 7.0) 及三个 conda 环境：

| conda 环境 | 包含工具 |
|------------|----------|
| `qc` | fastp, seqtk, fastqc, multiqc |
| `basic2` | bowtie2, samtools |
| `meta` | megahit, cd-hit |

环境名可在 `config.yaml` 的 `conda` 分组下修改。

其他外部依赖（需自行配置路径）：
- MetaGeneMark：`software.metagenemark_bin` / `software.metagenemark_model`
- 宿主参考基因组：`database.host_reference`（bowtie2 index 前缀）

## 使用方式

1. 编辑 `input.tsv`，填写样本信息（三列，无表头，制表符分隔）：

   ```
   <样本ID>	<R1.fastq.gz路径>	<R2.fastq.gz路径>
   ```

2. 按需修改 `config.yaml` 中的参考路径、软件路径、线程数、参数等。

3. 试运行（检查 DAG 与命令，不实际执行）：

   ```bash
   snakemake --cores 32 --use-conda --dry-run
   ```

4. 正式运行：

   ```bash
   snakemake --cores 32 --use-conda --rerun-incomplete --scheduler greedy --config samples=$PWD/tests/input.tsv --directory /data/mengxf/Develop/KML260617-MetaGenomics/results/260714
   ```

   集群提交可追加 `--cluster` / `--profile` 等参数。

## 配置说明 (`config.yaml`)

| 分组 | 字段 | 说明 |
|------|------|------|
| samples | — | 样本输入表路径 |
| workflow | steps | 启用步骤列表 |
| conda | qc / basic2 / meta | 各工具 conda 环境名 |
| threads | low / medium / high | 线程分级 |
| database | host_reference | 宿主参考基因组索引前缀 |
| software | metagenemark_bin / metagenemark_model | MetaGeneMark 可执行文件及模型路径 |
| fastp | extra | fastp 额外参数 |
| bowtie2 | host_removal_extra / mapping_extra | bowtie2 比对参数 |
| megahit | extra | megahit 组装参数 |
| cd_hit | extra | cd-hit 聚类参数 |

`config.schema.yaml` 在运行时对 `config.yaml` 做结构与类型校验，配置项缺失或类型错误会直接报错。

## 输出文件

每个样本（以 `<sample>` 为前缀）生成：

- `<sample>.cleaned.{1,2}.fastq.gz` —— 清洗后 FASTQ
- `<sample>.fastp.html` / `.fastp.json` —— fastp 报告
- `<sample>.cleaned.{1,2}_fastqc.html` / `.zip` —— fastqc 报告
- `<sample>_host_removed.{1,2}.fastq.gz` —— 去宿主后 FASTQ
- `<sample>_megahit/final.contigs.fa` —— 组装 contig
- `<sample>.contigs.gt500.fa` —— > 500bp contig
- `<sample>.gm.{gff,faa,fna}` —— 基因预测结果
- `<sample>.gm.gt100bp.fna` —— > 100bp 基因序列
- `<sample>.cd-hit.fna` / `.clstr` —— 去冗余基因集
- `<sample>.cd-hit.bowtie2.sorted.bam` —— 基因定量 BAM
- `qc/multiqc_data/multiqc_report.html` —— 汇总报告

日志与耗时统计位于 `.log/` 目录。

## 注意事项

- `gunzip` 为原始样本数据存在压缩问题时设置的解压步骤（已设为 `temp()` 自动清理）。若 megahit 支持直接读取 gz 输入，可跳过此步。
- 原始脚本中 megahit 输出目录为 `<sample>_megahit`，但下游 seqtk 读取路径为 `megahit/`，本项目已统一修正为 `<sample>_megahit`。
- 去宿主的 SAM 文件仅作中间产物，已设为 `temp()`，不影响后续分析。

## AI 协作准则

遵循 `E:\xiangfu.meng\AGENTS\AI_RULES` 中的个人准则，要点如下。

### 编程准则

1. `Python` 代码遵循 `PEP-8` 规范，不要写抽象代码。
2. 不要设计抽象复杂的代码单元，代码简单直接。
3. 不要造无用的轮子，效率优先。
4. 除非明确要求，否则不要写 `argparse` 或 `click` 之类的命令行参数。
5. 测试运行的输入不要写到输入文件的目录，所有的增改删都在当前工作目录执行。
6. 任何删除、修改现有文件的操作都要先经过我同意。

### 代码风格

1. 除基础元素赋值或简单逻辑外，逐行写注释。
2. 类和函数输入输出参数要标明类型。

### 分析习惯

1. 每次执行完分析都总结后追加更新到工作目录下的 `AGENTS.md` 文件，如果不存在则新建。
2. `Python` 表格处理优先用 `pandas`，作图用 `seaborn` 和 `matplotlib`，无法实现或提示词有明确指定库的情况除外。
3. 输入数据文件过大时，不要都读进缓存，读部分看大概数据结构。

### 解释器地址

- Windows 个人办公电脑：`C:\Users\xiangfu.meng\AppData\Local\miniconda3\envs\python3.12\python.exe`
- Ubuntu 服务器：`/home/mengxf/miniforge3/envs/python3.12/bin/python`

## 更新记录

### 2026-08-27 分类明细 abund_cum 保留小数

- 问题：`scripts/make_krona_input.py` 输出的基因分类明细表中出现 `gene_cum > abund_cum`，看似不合理。
- 原因：`gene_cum` 是基因个数（count），`abund_cum` 是 `norm_abund` 之和（sum），量纲不同，不可直接比较。`norm_abund = reads / gene_len / total_reads * 10^9`，单个基因的值可为 <1 的小数，低丰度基因聚集的分类节点天然满足 `gene_cum > abund_cum`，属正常现象。
- 修改：`abund_cum` 不再 `astype(int)` 向下取整，保留小数，避免丰度 <1 的节点被截断为 0 造成观感异常；`gene_cum` 仍为整数。
