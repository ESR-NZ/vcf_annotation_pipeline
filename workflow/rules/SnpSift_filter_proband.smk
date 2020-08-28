rule SnpSift_filter_proband:
    input:
        "../results/filtered/{sample}_filtered_multiallelicsites.vcf"
    output:
        temp("../results/filtered/{sample}_filtered_scoutfiltered.vcf")
    params:
        " ' ( isVariant ( GEN[0] ) ) ' " # Proband must be in the first column of the vcf
    log: 
        "logs/SnpSift_filter_proband/{sample}.log"
    benchmark:
        "benchmarks/SnpSift_filter_proband/{sample}.tsv"
    conda:
        "../envs/SnpSift.yaml"
    message:
        "Filtering for variants only found in the proband"
    shell:
        "cat {input} | SnpSift filter {params} > {output}"