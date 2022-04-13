rule pbrun_vqsr_snp:
    input:
        vcf = "../results/filtered/{sample}_tmp_vqsr_recal_indels.vcf"
    output:
        vcf = protected("../results/filtered/{sample}_filtered.vcf"),
        recal = temp("../results/filtered/{sample}_recal_snps"),
        tranches = temp("../results/filtered/{sample}_tranches_snps")
    params:
        snptranche = expand("--truth-sensitivity-level {snptranche}", snptranche = config['TRANCHE']['SNPS']),
        resources = get_cohort_snp_filtering_command,
        other = "--mode SNP --annotation QD --annotation MQ --annotation MQRankSum --annotation ReadPosRankSum"
    log:
        "logs/pbrun_vqsr_snp/{sample}.log"
    benchmark:
        "benchmarks/pbrun_vqsr_snp/{sample}.tsv"
    message:
        "Building a recalibration model to score variant quality in {input.vcf} and apply a score cutoff to filter variants (snp's)."
    shell:
        "/opt/parabricks/3.6.1/parabricks/pbrun vqsr --in-vcf {input.vcf} --out-vcf {output.vcf} --out-recal {output.recal} --out-tranches {output.tranches} {params.snptranche} {params.resources} {params.other} &> {log}"