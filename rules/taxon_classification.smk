rule diamond_nr_miro:
    input:
        rules.gene_catalogue_unigene.output.unigene,
    output:
        # todo 测试后 temp() 中间文件
        "taxon_classification/diamond_nr_miro.daa",
    benchmark:
        ".log/taxon_classification/diamond_nr_miro.bm"
    log:
        ".log/taxon_classification/diamond_nr_miro.log",
    params:
        db=config["database"]["nr_micro"],
        extra=config["params"]["diamond"],
    conda:
        config["conda"]["diamond"]
    threads: config["threads"]["max"]
    shell:
        # --outfmt 100 输出 DIAMOND alignment archive (DAA) 格式, 用于后续 MEGAN6 做 LCA 和功能分析
        """
        diamond blastx \
            {params.extra} \
            --query {input} \
            --db {params.db} \
            --outfmt 100 \
            --out {output} \
            --threads {threads} \
            --log \
            2> {log}
        """


rule daa_meganizer:
    message:
        "MEGAN6 LCA 物种分类"
    input:
        rules.diamond_nr_miro.output,
    output:
        flag="taxon_classification/diamond_nr_miro.daa.meganized",
    benchmark:
        ".log/taxon_classification/daa_meganizer.bm"
    log:
        ".log/taxon_classification/daa_meganizer.log",
    params:
        mdb=config["database"]["megan_map"],
    conda:
        config["conda"]["megan6"]
    threads: config["threads"]["max"]
    shell:
        # daa-meganizer 就地修改, 直接修改输入文件
        # -i 是 diamond 比对后输出的结果文件
        # -mdb 是适应 megan 的映射文件
        # --longReads对较长的 contigs/gene 开启该模式
        """
        daa-meganizer -i {input} -mdb {params.mdb} --longReads --threads {threads} 2> {log}
        touch {output.flag}
        """


rule daa2info:
    message:
        "提取 NCBI 分类信息"
    input:
        flag=rules.daa_meganizer.output.flag,
        daa="taxon_classification/diamond_nr_miro.daa",
    output:
        "taxon_classification/diamond_nr_miro.daa.meganized.info",
    benchmark:
        ".log/taxon_classification/daa2info.bm"
    log:
        ".log/taxon_classification/daa2info.log",
    conda:
        config["conda"]["megan6"]
    params:
        extra=config["params"]["daa2info"],
    threads: config["threads"]["max"]
    shell:
        "daa2info -i {input.daa} -o {output} {params.extra} 2> {log}"


rule make_krona_input:
    message:
        "生成 krona 输入文件和 lineage 丰度明细文件"
    input:
        megan=rules.daa2info.output[0],
        abund=rules.sample_gene_abundance.output[0],
    output:
        krona="taxon_classification/krona_input/{sample}.txt",
        detail="taxon_classification/linage_abund/{sample}.txt",
    benchmark:
        ".log/taxon_classification/{sample}.make_krona_input.bm"
    log:
        ".log/taxon_classification/{sample}.make_krona_input.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/make_krona_input.py"


rule krona_ktImportText:
    message:
        "运行 krona 所有样本合并输入"
    input:
        expand("taxon_classification/krona_input/{sample}.txt", sample=samples),
    output:
        "upload/all.krona.html",
    benchmark:
        ".log/taxon_classification/krona_ktImportText.bm"
    log:
        ".log/taxon_classification/krona_ktImportText.log",
    conda:
        config["conda"]["krona"]
    shell:
        "ktImportText {input} -o {output} 2> {log}"


rule make_taxon_abund_matrix:
    message:
        "生成物种丰度矩阵和物种基因数矩阵"
    input:
        expand("taxon_classification/linage_abund/{sample}.txt", sample=samples),
    output:
        # common.smk 中的 RANKS,KINDS 用于下游流程输入使用
        "taxon_classification/taxon_abund_matrix/{rank}.{kind}.txt",
    wildcard_constraints:
        rank="K|P|C|O|F|G|S",
        kind="gene|abund",
    log:
        ".log/taxon_classification/make_taxon_abund_matrix-{rank}-{kind}.log",
    benchmark:
        ".log/taxon_classification/make_taxon_abund_matrix-{rank}-{kind}.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/make_taxon_abund_matrix.py"


rule plot_taxon_abund_top_bar:
    message:
        "绘制物种丰度 Top10 柱状图"
    input:
        "taxon_classification/taxon_abund_matrix/{rank}.abund.txt",
    output:
        "upload/taxon_plot/{rank}.abund.top_bar.png",
    log:
        ".log/taxon_classification/taxon_plot/plot_taxon_abund_top_bar-{rank}.log",
    benchmark:
        ".log/taxon_classification/taxon_plot/plot_taxon_abund_top_bar-{rank}.bm"
    params:
        top=10,
    conda:
        config["conda"]["python"]
    script:
        "../scripts/plot_taxon_abund_top_bar.py"


rule plot_taxon_top_heatmap:
    message:
        "绘制物种丰度 Top30 热图"
    input:
        rules.make_taxon_abund_matrix.output,
    output:
        "upload/taxon_plot/{rank}.{kind}.top_heatmap.png",
    log:
        ".log/taxon_classification/taxon_plot/plot_taxon_top_heatmap-{rank}-{kind}.log"
    benchmark:
        ".log/taxon_classification/taxon_plot/plot_taxon_top_heatmap-{rank}-{kind}.bm"
    params:
        top=30,
    conda:
        config["conda"]["python"]
    script:
        "../scripts/plot_taxon_top_heatmap.py"
