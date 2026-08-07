# ✴️Claude 202607313 - 建库方法
提取思路取决于你的 NR 数据库是 **BLAST格式数据库**还是**纯fasta文件**,两种方式如下:

## 方案一:NR是BLAST数据库格式(推荐,效率最高)

如果你的NR是通过 `makeblastdb` 建好的库(有.pin/.phr/.psq等文件),可以直接用 `blastdbcmd` 按taxid批量抽取,不需要自己解析fasta头。

**第一步:获取taxonomy分类信息,拿到四大类群的完整taxid列表**

用 `taxonkit`(需要先下载NCBI taxdump):

```bash
# 下载并解压taxdump
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
mkdir -p ~/.taxonkit && tar -zxvf taxdump.tar.gz -C ~/.taxonkit

# 获取四个类群下所有子taxid(包括各级子分类)
# Bacteria=2, Archaea=2157, Fungi=4751, Viruses=10239
taxonkit list --ids 2 > bacteria_taxids.txt
taxonkit list --ids 2157 > archaea_taxids.txt
taxonkit list --ids 4751 > fungi_taxids.txt
taxonkit list --ids 10239 > viruses_taxids.txt
```

**第二步:用blastdbcmd按taxid列表抽取序列**

```bash
blastdbcmd -db nr -taxidlist bacteria_taxids.txt -out bacteria.fasta
blastdbcmd -db nr -taxidlist archaea_taxids.txt -out archaea.fasta
blastdbcmd -db nr -taxidlist fungi_taxids.txt -out fungi.fasta
blastdbcmd -db nr -taxidlist viruses_taxids.txt -out viruses.fasta

# 合并成一个新库
cat bacteria.fasta archaea.fasta fungi.fasta viruses.fasta > microbe_nr.fasta

# 重新建库(可选,如果后续要用blast比对)
makeblastdb -in microbe_nr.fasta -dbtype prot -out microbe_nr
```

## 方案二:NR是纯fasta文件

这种情况需要先建立 accession → taxid 的映射,再筛选。

```bash
# 下载accession2taxid映射(比较大,几个GB)
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz

# 用taxonkit获取四大类群的taxid列表(同上)

# 从accession2taxid中筛出属于这四类的accession
zcat prot.accession2taxid.gz | awk 'NR==FNR{a[$1];next} $3 in a' bacteria_taxids.txt - > bacteria_acc.txt
# (archaea/fungi/viruses同理)

# 用seqkit根据accession列表从fasta中抽取序列
seqkit grep -f bacteria_acc.txt nr.fasta -o bacteria.fasta
# (其余三类同理,最后cat合并)
```

**两种方案的选择建议:**
- 如果你本地有 `nr` 的BLAST库文件(常见于用blast+比对宏基因组的流程),**方案一**又快又省事,强烈推荐
- 如果只有一个巨大的 `nr.fasta` 平文本文件,只能走**方案二**,但 accession2taxid 文件很大(未压缩十几GB),对磁盘和内存有一定要求

你现在NR是哪种形式(BLAST库还是fasta文件)?我可以帮你把对应的命令细化成可直接跑的脚本。

---

# ✴️Claude 20260806 - 合并库还是分开库
这个问题的答案跟"建库"那次不同——**跑样本比对时,合并库通常更快**,原因和上次相反侧重点不一样。

## 为什么合并库更快

跑BLAST(或DIAMOND)比对时,决定总耗时的主要是:
- **查询序列(query)要扫描比对的次数**
- **数据库索引加载开销**
- **多线程能否充分利用**

如果分成4个子库分别跑:
- 你的样本query文件要被**读取加载4次**(每次跑一个子库都要重新load query)
- 数据库句柄、索引也要**初始化4次**
- 4次运行的进程调度、I/O开销都是独立的,叠加起来比合并跑一次更费时
- 而且4个子库总的序列比对计算量,跟合并库的计算量本质上是一样的(该比的序列都要比,没有减少)

所以:**分开跑 = 合并跑的计算量 + 3次额外的query加载与索引初始化开销**,肯定更慢。

## 合并库不会丢失分类群信息——只要建库时保留taxid映射

很多人分开建库,就是怕合并后不知道每条比对结果属于哪个类群。其实不用担心,只要建库时加上taxid映射:

```bash
makeblastdb -in microbe_nr.fasta -dbtype prot \
  -parse_seqids \
  -taxid_map accession_taxid.txt \
  -out microbe_nr
```

跑比对时指定输出格式带上 `staxids`:

```bash
blastx -query sample.fasta -db microbe_nr \
  -outfmt "6 qseqid sseqid pident length evalue staxids" \
  -num_threads 16 \
  -out sample_vs_microbe.tsv
```

结果里每一行都带 `staxid`,后续用taxonkit把taxid映射回"细菌/真菌/古菌/病毒"四大类,直接按需拆分统计即可,不影响你后面分类群分析。

## 一句话总结

| | 建库阶段 | 比对样本阶段 |
|---|---|---|
| **推荐** | 分开建4个子库,最后cat合并成总库 | 直接用合并好的总库跑一次比对 |
| 原因 | 需要保留各类群fasta文件时才分开跑,cat合并本身不耗时 | 分开跑=多付3次query加载+索引初始化的开销,总比对计算量不会减少 |

所以最佳流程是:**建库时按类群分开抽取(方便管理/统计各类群序列数)→ 建库时合并 → 跑比对时用合并后的总库一次性跑,靠staxids字段区分类群**。