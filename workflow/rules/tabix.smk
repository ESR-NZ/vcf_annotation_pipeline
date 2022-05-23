rule tabix:
    input:
        "../results/annotated/{sample}_filtered_dbnsfp.vcf.gz"
    output:
        temp("../results/annotated/{sample}_filtered_dbnsfp.vcf.gz.tbi")
    singularity:
        "docker://staphb/htslib:1.15"
    message:
        "Tabix indexing {input}"
    shell:
        "tabix {input}"