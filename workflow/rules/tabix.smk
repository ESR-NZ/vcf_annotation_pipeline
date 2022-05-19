rule tabix:
    input:
        "../results/annotated/{sample}_filtered_dbnsfp.vcf.gz"
    output:
        temp("../results/annotated/{sample}_filtered_dbnsfp.vcf.gz.tbi")
    singularity:
        "docker://staphb/samtools:1.15"
    message:
        "Tabix indexing {input}"
    shell:
        "tabix {input}"