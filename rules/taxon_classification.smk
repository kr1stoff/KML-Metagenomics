rule diamond_nr_miro:
    input:
        rules.gene_catalogue_unigene.output.unigene,
    output:
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
