rule SnpSift_annotate_dbSNP:
    input:
        vcf = "../results/annotated/{sample}_filtered_scoutfiltered_dbnsfp_vep_cadd.vcf",
        dbsnp = expand("{dbsnp}", dbsnp = config['dbSNP'])
    output:
        protected("../results/annotated/{sample}_filtered_scoutfiltered_annotated.vcf")
    params:
        "--regions"
    log: 
        "logs/SnpSift_annotate_dbSNP/{sample}.log"
    benchmark:
        "benchmarks/SnpSift_annotate_dbSNP/{sample}.tsv"
    conda:
        "../envs/SnpSift.yaml"
    message:
        "Using SnpSift to annotate {input.vcf} with dbSNP"
    shell:
        "SnpSift -Xmx16g annotate {input.dbsnp} {input.vcf} > {output}"