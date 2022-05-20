rule bgzip:
    input:
        "../results/annotated/{sample}_filtered_dbnsfp.vcf"
    output:
        temp("../results/annotated/{sample}_filtered_dbnsfp.vcf.gz")
    singularity:
        "docker://biocontainers/samtools:v1.9-4-deb_cv1"
    message:
        "Bgzipping {input}"
    shell:
        "bgzip < {input} > {output}"