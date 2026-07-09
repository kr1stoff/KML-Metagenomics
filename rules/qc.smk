"""
fastp 质控过滤 + fastqc 质控报告 + multiqc 汇总
"""


rule fastp:
    input:
        r1=lambda wc: samples.set_index("sample_id").loc[wc.sample, "fq1"],
        r2=lambda wc: samples.set_index("sample_id").loc[wc.sample, "fq2"],
    output:
        r1="{sample}.cleaned.1.fastq.gz",
        r2="{sample}.cleaned.2.fastq.gz",
        html="{sample}.fastp.html",
        json="{sample}.fastp.json",
    benchmark:
        ".log/qc/fastp/{sample}.fastp.bm"
    log:
        ".log/qc/fastp/{sample}.fastp.log",
    conda:
        config["conda"]["qc"]
    threads: config["threads"]["high"]
    shell:
        "fastp {config[fastp][extra]} -w {threads} "
        "-i {input.r1} -I {input.r2} "
        "-o {output.r1} -O {output.r2} "
        "-h {output.html} -j {output.json} "
        "2> {log}"


rule fastqc:
    input:
        r1="{sample}.cleaned.1.fastq.gz",
        r2="{sample}.cleaned.2.fastq.gz",
    output:
        html_r1="{sample}.cleaned.1_fastqc.html",
        html_r2="{sample}.cleaned.2_fastqc.html",
        zip_r1="{sample}.cleaned.1_fastqc.zip",
        zip_r2="{sample}.cleaned.2_fastqc.zip",
    benchmark:
        ".log/qc/fastqc/{sample}.fastqc.bm"
    log:
        ".log/qc/fastqc/{sample}.fastqc.log",
    conda:
        config["conda"]["qc"]
    threads: config["threads"]["medium"]
    shell:
        "fastqc -t {threads} "
        "{input.r1} {input.r2} "
        "2> {log}"


rule multiqc:
    input:
        expand("{sample}.fastp.json", sample=samples["sample_id"].tolist()),
        expand("{sample}.cleaned.1_fastqc.zip", sample=samples["sample_id"].tolist()),
        expand("{sample}.cleaned.2_fastqc.zip", sample=samples["sample_id"].tolist()),
    output:
        "qc/multiqc_data/multiqc_report.html",
    benchmark:
        ".log/qc/multiqc/multiqc.bm"
    log:
        ".log/qc/multiqc/multiqc.log",
    conda:
        config["conda"]["qc"]
    threads: config["threads"]["low"]
    shell:
        "multiqc -f -o qc/multiqc_data . "
        "2> {log}"
