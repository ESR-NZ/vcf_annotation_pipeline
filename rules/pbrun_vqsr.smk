rule pbrun_vqsr:
    input:
        vcf = "../vcf/{sample}_raw_snps_indels_AS_g.vcf"
    output:
        vcf = "filtered/{sample}_filtered.vcf",
        recal = "recalibrated/{sample}.recal",
        tranches = "recalibrated/{sample}.tranches"
    params:
        resources_indels = expand("{resources}", resources = config['FILTERING']['COHORT']['INDELS']),
        resources_snps = expand("{resources}", resources = config['FILTERING']['COHORT']['SNPS'])
    log: 
        "logs/pbrun_vqsr/{sample}.log"
    benchmark:
        "benchmarks/pbrun_vqsr/{sample}.pbrun.vqsr"
    message:
        "Building a recalibration model to score variant quality and apply a score cutoff to filter variants"
    shell:
        "pbrun vqsr --in-vcf {input.vcf} --out-vcf {output.vcf} --out-recal {output.recal} --out-tranches {output.tranches} {params.resources_indels} {params.resources_snps} --annotation QD --annotation MQ --annotation MQRankSum -annotation ReadPosRankSum"