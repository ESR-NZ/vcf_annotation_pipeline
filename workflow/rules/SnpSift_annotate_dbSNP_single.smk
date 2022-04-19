rule SnpSift_annotate_dbSNP_single:
    input:
        vcf = "../results/annotated/{sample}_filtered_dbnsfp_vep_cadd.vcf",
        dbsnp = expand("{dbsnp}", dbsnp = config['dbSNP'])
    output:
        "../results/annotated/{sample}_filtered_annotated.vcf"
    params:
        maxmemory = expand('"-Xmx{maxmemory}"', maxmemory = config['MAXMEMORY'])
    log: 
        "logs/SnpSift_annotate_dbSNP/{sample}.log"
    benchmark:
        "benchmarks/SnpSift_annotate_dbSNP/{sample}.tsv"
    conda:
        "../envs/SnpSift.yaml"
    message:
        "Using SnpSift to annotate {input.vcf} with dbSNP"
    shell:
        "SnpSift {params.maxmemory} annotate {input.dbsnp} {input.vcf} > {output} 2> {log}"