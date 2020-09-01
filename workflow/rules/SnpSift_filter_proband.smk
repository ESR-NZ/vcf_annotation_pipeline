rule SnpSift_filter_proband:
    input:
        "../results/readyforscout/{sample}_filtered_annotated_multiallelicsites.vcf"
    output:
        temp("../results/readyforscout/{sample}_filtered_annotated_multiallelicsites_probandonly.vcf")
    params:
        " '( isVariant ( GEN[{sample}] ) ) ' "
    log: 
        "logs/SnpSift_filter_proband/{sample}.log"
    benchmark:
        "benchmarks/SnpSift_filter_proband/{sample}.tsv"
    conda:
        "../envs/SnpSift.yaml"
    message:
        "Filtering for variants only found in the proband in {input}"
    shell:
        "cat {input} | SnpSift filter {params} > {output}"