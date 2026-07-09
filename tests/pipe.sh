# 原始数据
fastp -g -q 5 -u 50 -n 15 -l 150 --overlap_diff_limit 1 --overlap_diff_percent_limit 10 -w 16 -i SRR23604269.raw.1.fastq.gz -I SRR23604269.raw.2.fastq.gz -o SRR23604269.cleaned.1.fastq.gz -O SRR23604269.cleaned.2.fastq.gz -h SRR23604269.fastp.html -j SRR23604269.fastp.json
# 去宿主, --un-conc-gz 用 % 匹配 1/2
mamba -n basic2 run bowtie2 -x /data/mengxf/Database/reference/hg19/hg19.fa -1 SRR23604269.cleaned.1.fastq.gz -2 SRR23604269.cleaned.2.fastq.gz --end-to-end --sensitive --no-hd --no-sq -I 200 -X 400 --threads 16 --un-conc-gz SRR23604269_host_removed.%.fastq.gz -S SRR23604269.bowtie2.sam
# todo fastp 加上求宿主后的数据量

# 组装
# ! SRR 样本压缩有问题，解压后可以用，真实样本不需要这一步
gunzip -k SRR23604269_host_removed.1.fastq.gz &
gunzip -k SRR23604269_host_removed.2.fastq.gz
# megahit -m memory参数宏基因组推荐 64G 以上, 按照当前服务器的线程数, 内存数进行调整. 
# 当前服务器为 32 线程, 256G 内存. -m 64000000000 -t 16 (64G内存 + 16线程)
mamba -n meta run megahit --presets meta-large -m 64000000000 -t 16 -1 SRR23604269_host_removed.1.fastq -2 SRR23604269_host_removed.2.fastq -o SRR23604269_megahit
# 保留大于 500bp 的序列
seqtk seq -L 500 megahit/final.contigs.fa > SRR23604269.contigs.gt500.fa

# 基因预测
/data/mengxf/Software/MetaGeneMark_linux_64/mgm/gmhmmp -m /data/mengxf/Software/MetaGeneMark_linux_64/mgm/MetaGeneMark_v1.mod -a -d -f G -p 1 SRR23604269.contigs.gt500.fa -o SRR23604269.gm.gff -A SRR23604269.gm.faa -D SRR23604269.gm.fna -L SRR23604269.gm.log
# 保留大于 100bp 的序列
seqtk seq -L 100 SRR23604269.gm.fna > SRR23604269.gm.gt100bp.fna
# 基因去冗余. -T 线程, -M 内存
mamba -n meta run cd-hit -T 16 -G 0 -aS 0.9 -g 1 -d 0 -c 0.95 -n 5 -M 8000 -i SRR23604269.gm.gt100bp.fna -o SRR23604269.cd-hit.fna
# bowtie2 比对 fastq 到去冗余后的基因
mamba -n basic2 run bowtie2-build --threads 16 SRR23604269.cd-hit.fna SRR23604269.cd-hit
mamba -n basic2 run bowtie2 -x SRR23604269.cd-hit --end-to-end --sensitive --no-hd --no-sq -I 200 -X 400 --threads 16 -1 SRR23604269_host_removed.1.fastq.gz -2 SRR23604269_host_removed.2.fastq.gz -S SRR23604269.cd-hit.bowtie2.sam
samtools view -bS SRR23604269.cd-hit.bowtie2.sam | samtools sort -o SRR23604269.cd-hit.bowtie2.sorted.bam
# 过滤掉reads数小于2的基因

# 统计起始/终止密码子数量等信息
