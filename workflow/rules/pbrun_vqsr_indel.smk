rule pbrun_vqsr_indel:
    input:
        "../../vcf/{sample}_raw_snps_indels_g.vcf"
    output:
        vcf = temp("../results/filtered/{sample}_tmp_vqsr_recal_indels.vcf"),
        recal = temp("../results/filtered/{sample}_recal_indels"),
        tranches = temp("../results/filtered/{sample}_tranches_indels")
    resources:
        gpu = 1
    params:
        indeltranche = expand("--truth-sensitivity-level {indeltranche}", indeltranche = config['FILTERING']['TRANCHE']['INDELS']),
        resources = expand("{resources}", resources = config['FILTERING']['COHORT']['INDELS']),
        other = "--mode INDEL --annotation QD --annotation MQ --annotation MQRankSum --annotation ReadPosRankSum --annotation FS --annotation SOR --annotation DP"
    log: 
        "logs/pbrun_vqsr_indel/{sample}.log"
    benchmark:
        "benchmarks/pbrun_vqsr_indel/{sample}.tsv"
    message:
        "Building a recalibration model to score variant quality and using machine learning to filter out probable artifacts from the variant callset (indels) for {input}"
    shell:
        "pbrun vqsr --in-vcf {input} --out-vcf {output.vcf} --out-recal {output.recal} --out-tranches {output.tranches} {params.indeltranche} {params.resources} {params.other} &> {log}"