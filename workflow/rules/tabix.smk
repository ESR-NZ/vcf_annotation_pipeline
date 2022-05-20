rule tabix:
    input:
        "../results/annotated/{sample}_filtered_dbnsfp.vcf.gz"
    output:
        temp("../results/annotated/{sample}_filtered_dbnsfp.vcf.gz.tbi")
    singularity:
        "docker://biocontainers/tabix:v1.9-11-deb_cv1"
    message:
        "Tabix indexing {input}"
    shell:
        "tabix {input}"