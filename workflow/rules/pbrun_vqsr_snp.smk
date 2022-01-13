rule pbrun_vqsr_snp:
    input:
        vcf = "../../human_genomics_pipeline/results/called/{sample}_raw_snps_indels.g.vcf"    
    output:
        vcf = protected("../results/filtered/{sample}_filtered.vcf"),
        recal = temp("../results/filtered/{sample}_recal_snps"),
        index = temp("../results/filtered/{sample}_recal_snps.idx"),
        tranches = temp("../results/filtered/{sample}_tranches_snps")
    params:
        snptranche = expand("--truth-sensitivity-level {snptranche}", snptranche = config['TRANCHE']['SNPS']),
        resources = get_cohort_snp_filtering_command,
        other = "--mode SNP --annotation QD --annotation MQ --annotation MQRankSum --annotation ReadPosRankSum"
    resources:
        gpu = config['GPU']
    log:
        "logs/pbrun_vqsr_snp/{sample}.log"
    benchmark:
        "benchmarks/pbrun_vqsr_snp/{sample}.tsv"
    message:
        "Building a recalibration model to score variant quality in {input.vcf} and apply a score cutoff to filter variants (snp's)."
    shell:
        "pbrun vqsr --in-vcf {input.vcf} --out-vcf {output.vcf} --out-recal {output.recal} --out-tranches {output.recal} {params.snptranche} {params.resources} {params.other} &> {log}"