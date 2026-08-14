# KML-Metagenomics

宏基因组分析流程：双端测序数据的质控过滤、去宿主、组装、基因预测、基因去冗余与丰度定量。

## 流程介绍
看 `AGENTS.md` 文件

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
   snakemake --cores 32 --use-conda --rerun-incomplete --scheduler greedy \
      --config samples=$PWD/tests/input.tsv \
      metadata=$PWD/tests/metadata.tsv \
      --directory /data/mengxf/Develop/KML260617-MetaGenomics/results/260714
   ```

## 注意事项

- `gunzip` 为原始样本数据存在压缩问题时设置的解压步骤（已设为 `temp()` 自动清理）。若 megahit 支持直接读取 gz 输入，可跳过此步。
- 原始脚本中 megahit 输出目录为 `<sample>_megahit`，但下游 seqtk 读取路径为 `megahit/`，本项目已统一修正为 `<sample>_megahit`。
- 去宿主的 SAM 文件仅作中间产物，已设为 `temp()`，不影响后续分析。

## 开发
- 20260806 diamond nr_lite 构建
   - 在 `/data/mengxf/Database/NCBI/blast/db/nr` 目录， 等待测试完成后，按照细菌、古菌、真菌、病毒进行拆分.
   - `docs\make_diamond_nr_database.md` 记录了构建diamond小nr数据库的思路和方法.
   - `diamond` 在 `rnaseq` conda 环境，taxonkit 在 `basic`
    ```bash
   /home/mengxf/miniforge3/envs/basic/bin/blastx \
      -query test.fna \
      -db /data/mengxf/Database/NCBI/blast/db/nr/nr \
      -out result_blastn.tsv \
      -outfmt "6 qseqid sseqid ssciname staxid pident qcovs length" \
      -evalue 1e-5 \
      -num_threads 32
   ```
- 20260806 core-pan 分析在 `scripts\core_pan_gene_analysis.py` 待整合
