rule bgzip:
    input:
        "../results/annotated/{sample}_filtered_scoutfiltered_dbnsfp.vcf"
    output:
        temp("../results/annotated/{sample}_filtered_scoutfiltered_dbnsfp.vcf.gz")
    conda:
        "../envs/bgzip.yaml"
    message:
        "Bgzipping {input}"
    shell:
        "bgzip < {input} > {output}"