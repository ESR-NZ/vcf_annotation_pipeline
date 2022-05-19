rule bcftools_view_multiallelicsites:
    input:
        "../results/annotated/{sample}_filtered_annotated.vcf"
    output:
        temp("../results/readyforscout/{sample}_filtered_annotated_multiallelicsites.vcf.gz")
    params:
        "-O z --max-alleles 2 --exclude-types indels"
    log:
        "logs/bcftools_view_multiallelicsites/{sample}.log"
    benchmark:
        "benchmarks/bcftools_view_multiallelicsites/{sample}.tsv"
    singularity:
        "docker://biocontainers/bcftools:v1.9-1-deb_cv1"
    message:
        "Filtering out multiallelic sites in {input}"
    shell:
        "bcftools view {params} {input} > {output} 2> {log}"