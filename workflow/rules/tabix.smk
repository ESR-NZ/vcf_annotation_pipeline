rule tabix:
    input:
        "../results/annotated/{sample}_filtered_scoutfiltered_dbnsfp.vcf.gz"
    output:
        temp("../results/annotated/{sample}_filtered_scoutfiltered_dbnsfp.vcf.gz.tbi")
    conda:
        "../envs/bgzip.yaml"
    message:
        "Tabix indexing {input}"
    shell:
        "tabix {input}"