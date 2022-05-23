rule bgzip:
    input:
        "../results/annotated/{sample}_filtered_dbnsfp.vcf"
    output:
        temp("../results/annotated/{sample}_filtered_dbnsfp.vcf.gz")
    singularity:
        "docker://staphb/htslib:1.15"
    message:
        "Bgzipping {input}"
    shell:
        "bgzip < {input} > {output}"