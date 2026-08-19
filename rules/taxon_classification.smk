rule diamond_nr_miro:
    input:
        rules.gene_catalogue_unigene.output.unigene,
    output:
        "taxon_classification/diamond_nr_miro.tsv"
    benchmark:
        ".log/taxon_classification/diamond_nr_miro.bm"
    log:
        ".log/taxon_classification/diamond_nr_miro.log"
    params:
        db=config["database"]["nr_micro"],
        evalue=config["params"]["diamond"],
    conda:
        config["conda"]["diamond"]
    threads:
        config["threads"]["max"]
    shell:
        """
        diamond blastx --query {input} \
            --db {params.db} \
            --outfmt 6 qseqid sseqid sscinames staxids pident qcovhsp length \
            --out {output} \
            --threads {threads} \
            --log \
            2> {log}
        """
