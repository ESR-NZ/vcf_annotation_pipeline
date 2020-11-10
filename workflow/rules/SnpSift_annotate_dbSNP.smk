if config['DATA'] == "Single" or config['DATA'] == 'single':
    outfile = "../results/annotated/{sample}_filtered_annotated.vcf"

elif config['DATA'] == "Cohort" or config['DATA'] == 'cohort':
    outfile = "../results/annotated/{sample}_filtered_dbnsfp_vep_cadd_dbsnp.vcf"

rule SnpSift_annotate_dbSNP:
    input:
        vcf = "../results/annotated/{sample}_filtered_dbnsfp_vep_cadd.vcf",
        dbsnp = expand("{dbsnp}", dbsnp = config['dbSNP'])
    output:
        protected(outfile)
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