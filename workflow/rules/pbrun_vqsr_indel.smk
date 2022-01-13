rule pbrun_vqsr_indel:
    input:
        vcf = "../results/filtered/{sample}_tmp_vqsr_recal_indels.vcf"
    output:
        vcf = temp("../results/filtered/{sample}_tmp_vqsr_recal_indels.vcf"),
        vcfindex = temp("../results/filtered/{sample}_tmp_vqsr_recal_indels.vcf.idx"),
        recal = temp("../results/filtered/{sample}_recal_indels"),
        index = temp("../results/filtered/{sample}_recal_indels.idx"),
        tranches = temp("../results/filtered/{sample}_tranches_indels")
    params:
        indeltranche = expand("--truth-sensitivity-level {indeltranche}", indeltranche = config['TRANCHE']['INDELS']),
        resources = get_cohort_indel_filtering_command,
        other = "--mode INDEL -an QD -an MQ -an MQRankSum -an ReadPosRankSum"
    resources:
        gpu = config['GPU']
    log:
        "logs/pbrun_vqsr_indel/{sample}.log"
    benchmark:
        "benchmarks/pbrun_vqsr_indel/{sample}.tsv"
    message:
        "Building a recalibration model to score variant quality in {input.vcf} and apply a score cutoff to filter variants (indels)."
    shell:
        "pbrun vqsr --in-vcf {input.vcf} --out-vcf {output.vcf} --out-recal {output.recal} --out-tranches {output.recal} {params.indeltranche} {params.resources} {params.other} &> {log}"