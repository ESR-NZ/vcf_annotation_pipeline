rule SnpSift_dbnsfp:
    input:
        vcf = "../results/filtered/{sample}_filtered.vcf",
        dbnsfp = expand("{dbnsfp}", dbnsfp = config['dbNSFP'])
    output:
        temp("../results/annotated/{sample}_filtered_dbnsfp.vcf")
    params:
        maxmemory = expand('"-Xmx{maxmemory}"', maxmemory = config['MAXMEMORY']),
        other = "-v"
    log: 
        "logs/SnpSift_dbnsfp/{sample}.log"
    benchmark:
        "benchmarks/SnpSift_dbnsfp/{sample}.tsv"
    conda:
        "../envs/SnpSift.yaml"
    message:
        "Using the dbNSFP database to annotate {input.vcf} with functional predictions from multiple algorithms (SIFT, Polyphen2, LRT and MutationTaster, PhyloP and GERP++, etc.)"
    shell:
        "SnpSift {params.maxmemory} dbnsfp {input.vcf} > {output} -db {input.dbnsfp} {params.other} 2> {log}"