rule pbrun_vqsr_snp:
    input:
        "../results/filtered/{sample}_tmp_vqsr_recal_indels.vcf"
    output:
        vcf = protected("../results/filtered/{sample}_filtered.vcf"),
        recal = temp("../results/filtered/{sample}_recal_snp"),
        tranches = temp("../results/filtered/{sample}_tranches_snp")
    resources:
        gpu = 1
    params:
        snptranche = expand("--truth-sensitivity-level {snptranche}", snptranche = config['FILTERING']['TRANCHE']['SNPS']),
        resources = expand("{resources}", resources = config['FILTERING']['COHORT']['SNPS']),
        other = "--mode SNP --annotation QD --annotation MQ --annotation MQRankSum --annotation ReadPosRankSum --annotation FS --annotation SOR --annotation DP"
    log: 
        "logs/pbrun_vqsr_snp/{sample}.log"
    benchmark:
        "benchmarks/pbrun_vqsr_snp/{sample}.tsv"
    message:
        "Building a recalibration model to score variant quality and using machine learning to filter out probable artifacts from the variant callset (snps) for {input}"
    shell:
        "pbrun vqsr --in-vcf {input} --out-vcf {output.vcf} --out-recal {output.recal} --out-tranches {output.tranches} {params.snptranche} {params.resources} {params.other} &> {log}"