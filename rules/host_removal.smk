"""
bowtie2 去宿主 + 解压 (供 megahit 使用)
"""


rule bowtie2_host_removal:
    input:
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2,
    output:
        # todo 删除中间文件
        r1="host_removal/{sample}_host_removed.1.fastq.gz",
        r2="host_removal/{sample}_host_removed.2.fastq.gz",
        sam="host_removal/{sample}.bowtie2.sam",
    benchmark:
        ".log/host_removal/bowtie2/{sample}.bowtie2_host_removal.bm"
    log:
        ".log/host_removal/bowtie2/{sample}.bowtie2_host_removal.log",
    conda:
        config["conda"]["bowtie2"]
    threads: config["threads"]["medium"]
    params:
        extra=config["params"]["bowtie2"]["host_removal"],
        host_ref=config["database"]["host_reference"],
    shell:
        "bowtie2 -x {params.host_ref} "
        "-1 {input.r1} -2 {input.r2} "
        "{params.extra} "
        "--threads {threads} "
        "--un-conc-gz host_removal/{wildcards.sample}_host_removed.%.fastq.gz "
        "-S {output.sam} "
        "2> {log}"


rule gunzip_host_removed:
    message: "解压去宿主后的 fastq, 供 megahit 组装使用, 真实样本不需要这一步",
    input:
        r1=rules.bowtie2_host_removal.output.r1,
        r2=rules.bowtie2_host_removal.output.r2,
    output:
        # todo 删除中间文件
        r1="host_removal/{sample}_host_removed.1.fastq",
        r2="host_removal/{sample}_host_removed.2.fastq",
    benchmark:
        ".log/host_removal/gunzip/{sample}.gunzip.bm"
    log:
        ".log/host_removal/gunzip/{sample}.gunzip.log",
    shell:
        "gunzip -c {input.r1} > {output.r1} 2>> {log}; "
        "gunzip -c {input.r2} > {output.r2} 2>> {log}"


rule nonhost_seqtk_size:
    input:
        rules.bowtie2_host_removal.output.r1,
    output:
        "host_removal/{sample}.host_removal.size.txt",
    benchmark:
        ".log/host_removal/size/{sample}.size.bm"
    log:
        ".log/host_removal/size/{sample}.size.log",
    conda:
        config["conda"]["seqtk"]
    shell:
        "seqtk size {input} > {output} 2> {log}"


rule merge_nonhost_file:
    input:
        expand("host_removal/{sample}.host_removal.size.txt", sample=samples),
    output:
        "host_removal/merge_nonhost_size.txt",
    benchmark:
        ".log/host_removal/merge_nonhost_file.bm"
    run:
        f = open(output[0], "w")
        f.write("Sample\tNonHostBases\n")
        for size_file in input:
            se_nonhost_bases = open(size_file, "r").read().strip().split("\t")[1]
            # 双端
            nonhost_bases = str(int(se_nonhost_bases) * 2)
            sample = size_file.split("/")[-1].split(".")[0]
            print(f"{sample}\t{nonhost_bases}\n", file=f)
        f.close()


rule data_process_stats:
    input:
        fastp=rules.fastp_stats_summary.output[0],
        nonhost=rules.merge_nonhost_file.output[0],
    output:
        "upload/data_process_stats.xlsx",
    log:
        ".log/host_removal/data_process_stats.log",
    benchmark:
        ".log/host_removal/data_process_stats.bm"
    conda:
        config["conda"]["python"]
    # 使用 run 没法使用 conda 环境和 log, 用 script 可以
    script:
        "../scripts/data_process_stats.py"
