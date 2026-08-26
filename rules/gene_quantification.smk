"""
bowtie2 比对到去冗余基因集 + samtools 排序
"""


rule quant_gene_cat_bowtie2_build:
    input:
        rules.cd_hit.output,
    output:
        idx1="gene_prediction/gene_catalogue.raw.1.bt2",
        idx2="gene_prediction/gene_catalogue.raw.2.bt2",
        idx3="gene_prediction/gene_catalogue.raw.3.bt2",
        idx4="gene_prediction/gene_catalogue.raw.4.bt2",
        rev1="gene_prediction/gene_catalogue.raw.rev.1.bt2",
        rev2="gene_prediction/gene_catalogue.raw.rev.2.bt2",
    benchmark:
        ".log/gene_quantification/bowtie2_build.bm"
    log:
        ".log/gene_quantification/bowtie2_build.log",
    conda:
        config["conda"]["bowtie2"]
    threads: config["threads"]["max"]
    shell:
        "bowtie2-build --threads {threads} {input} gene_prediction/gene_catalogue.raw 2> {log}"


rule quant_bowtie2_map:
    input:
        r1=rules.bowtie2_host_removal.output.r1,
        r2=rules.bowtie2_host_removal.output.r2,
        idx=rules.quant_gene_cat_bowtie2_build.output,
    output:
        # todo 删除中间文件
        "gene_quantification/{sample}.gene_catalogue.sam",
    benchmark:
        ".log/gene_quantification/bowtie2_map/{sample}.bowtie2_map.bm"
    log:
        ".log/gene_quantification/bowtie2_map/{sample}.bowtie2_map.log",
    conda:
        config["conda"]["bowtie2"]
    threads: config["threads"]["medium"]
    params:
        extra=config["params"]["bowtie2"]["mapping"],
        rgid="{sample}",
        rgsm="{sample}",
        rgspl="Illumina",
    shell:
        "bowtie2 "
        "{params.extra} "
        "--rg-id {params.rgid} --rg SM:{params.rgsm} --rg PL:{params.rgspl} "
        "--threads {threads} "
        "-x gene_prediction/gene_catalogue.raw "
        "-1 {input.r1} "
        "-2 {input.r2} "
        "-S {output} "
        "2> {log}"


rule samtools_sort_index:
    input:
        rules.quant_bowtie2_map.output,
    output:
        bam="gene_quantification/{sample}.gene_catalogue.sorted.bam",
        bai="gene_quantification/{sample}.gene_catalogue.sorted.bam.bai",
    benchmark:
        ".log/quantification/{sample}.samtools_sort.bm"
    log:
        ".log/quantification/{sample}.samtools_sort.log",
    conda:
        config["conda"]["samtools"]
    shell:
        """
        samtools sort -o {output.bam} {input} 2> {log}
        samtools index {output.bam} 2> {log}
        """


rule samtools_idxstats:
    message:
        "samtools idxstats 统计各样本中各基因的reads数量, 过滤掉样本中 reads≤2 的基因"
    input:
        bam=rules.samtools_sort_index.output.bam,
        bai=rules.samtools_sort_index.output.bai,
    output:
        "gene_quantification/{sample}.reads_gt2.idxstats",
    benchmark:
        ".log/gene_quantification/{sample}.samtools_idxstats.bm"
    log:
        ".log/gene_quantification/{sample}.samtools_idxstats.log",
    conda:
        config["conda"]["samtools"]
    shell:
        "samtools idxstats {input.bam} | awk '$3>2' > {output} 2> {log}"


rule gene_reads_table:
    message:
        "样本-基因reads计数表, 第一列基因长度"
    input:
        expand("gene_quantification/{sample}.reads_gt2.idxstats", sample=samples),
    output:
        "gene_quantification/gene_reads_table.tsv",
    benchmark:
        ".log/gene_quantification/gene_reads_table.bm"
    log:
        ".log/gene_quantification/gene_reads_table.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/gene_reads_table.py"


rule sample_gene_abundance:
    message:
        "基因丰度, 就是 RPKM. 算法参考 docs/gene_abundance.md."
        "样本内, 基因丰度 = (基因reads/基因长度) / sum(所有基因reads数量) * 10**9"
    input:
        rules.samtools_idxstats.output,
    output:
        "gene_quantification/{sample}.sample_gene_abundance.tsv",
    benchmark:
        ".log/gene_quantification/{sample}.sample_gene_abundance.bm"
    log:
        ".log/gene_quantification/{sample}.sample_gene_abundance.log",
    conda:
        config["conda"]["python"]
    # "run:" 不能进入环境, 所以要用回 "script:"
    script:
        "../scripts/sample_gene_abundance.py"


rule gene_abundance_table:
    message:
        "样本-基因丰度表"
    input:
        expand("gene_quantification/{sample}.sample_gene_abundance.tsv", sample=samples),
    output:
        "gene_quantification/gene_abundance_table.tsv",
    benchmark:
        ".log/gene_quantification/gene_abundance_table.bm"
    log:
        ".log/gene_quantification/gene_abundance_table.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/gene_abundance_table.py"


rule gene_catalogue_unigene:
    input:
        abund=rules.gene_abundance_table.output,
        fa=rules.cd_hit.output,
    output:
        glist="gene_quantification/gene_catalogue.unigene.list",
        unigene="gene_quantification/gene_catalogue.unigene.fna",
    benchmark:
        ".log/gene_quantification/gene_catalogue_unigene.bm"
    log:
        ".log/gene_quantification/gene_catalogue_unigene.log",
    conda:
        config["conda"]["seqtk"]
    shell:
        """
        sed '1d' {input.abund} | cut -f1 > {output.glist} 2> {log}
        seqtk subseq {input.fa} {output.glist} > {output.unigene} 2>> {log}
        """
