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
- ⌛️20260824 - eggNOG  
   在 `meta` 环境

- ⌛️20260824 - krona
   - 准备 krona 输入, 顺便输出样本注释后物种, 层级, 丰富度信息表. `/data/mengxf/Develop/KML260617-MetaGenomics/work/260819-classification/make_krona_input.py`
   - 运行多样本 krona
      ```bash
      mamba -n meta run ktImportText SRR23604277.krona_input.txt SRR23604277.krona_input.txt2
      ```
- 🔄20260825 - 统一格式化 megan6 物种输出, 门[P]是可以对齐的, 域[D]/界[K] 都归为界[K]
   - 细菌(Bacteria),古菌(Archaea): 细菌古菌是域. "[D] Bacteria; [P] Firmicutes; [C] Clostridia; [O] Eubacteriales;" 门[P]上面是域[D]
   - 真菌(Fungi): 真菌是界, 门上面是界[K]
   - 病毒(Viruses)所有界
      ```text
      Abadenavirae
      Bamfordvirae
      Helvetiavirae
      Heunggongvirae
      Loebvirae
      Orthornavirae
      Pararnavirae
      Sangervirae
      Shotokuvirae
      Trapavirae
      Zilligvirae
      ```

- 20260821 - diamond(daa) + megan6
   ```bash
   # daa-meganizer 就地修改, 直接修改输入文件
   mamba -n meta run daa-meganizer -i test.daa -mdb /data/mengxf/Database/MEGAN6/megan-map-Feb2022.db --longReads --threads 32
   # 会在输入的原.daa文件上修改
   # -i 是diamond比对后输出的结果文件
   # -mdb 是适应megan的映射文件
   # --longReads对较长的contigs/gene开启该模式

   # 提取NCBI分类信息
   mamba -n meta run daa2info -i test.daa -o test_out -l -m -r2c Taxonomy -p true -r true
   # -o输出的文件中即每条序列的物种注释信息

   # 若要提取物种注释和功能注释信息
   # mamba -n meta run daa2info -i taxon_classification/diamond_nr_miro.daa -o test_out -l -m -r2c Taxonomy GTDB KEGG EC EGGNOG INTERPRO2GO SEED
   mamba -n meta run daa2info \
      -i taxon_classification/diamond_nr_miro.daa \
      -o test_out \
      --list true\
      --listMore true \
      --names true \
      --paths true \
      --prefixRank true \
      --majorRanksOnly true \
      --read2class Taxonomy
   ```

   加运行完成标记, 就不用复制了. 如果不行可以用复制的方式
   ```snakemake
   rule meganize_daa:
      input:
         daa = "results/{sample}.daa",
         a2t = "db/prot_acc2tax.bin"
      output:
         flag = "results/{sample}.daa.meganized"  # 标记文件
      shell:
         """
         daa-meganizer -i {input.daa} -a2t {input.a2t}
         touch {output.flag}
         """
   ```

- ✅20260806 diamond nr_lite 构建
   - 在 `/data/mengxf/Database/NCBI/blast/db/nr` 目录， 等待测试完成后，按照细菌、古菌、真菌、病毒进行拆分.
   - `docs\make_diamond_nr_database.md` 记录了构建diamond小nr数据库的思路和方法.
   - `diamond` 在 `rnaseq` conda 环境，taxonkit 在 `basic`

- ✅20260806 core-pan 分析在 `scripts\core_pan_gene_analysis.py` 待整合
