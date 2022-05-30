rule gatk_ApplyVQSR_indel:
    input:
        vcf = "../../human_genomics_pipeline/results/called/{sample}_raw_snps_indels.vcf",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME']),
        recal = "../results/filtered/{sample}_recal_indels",
        recalindex = "../results/filtered/{sample}_recal_indels.idx",
        tranches = "../results/filtered/{sample}_tranches_indels"
    output:
        vcf = temp("../results/filtered/{sample}_tmp_vqsr_recal_indels.vcf"),
        index = temp("../results/filtered/{sample}_tmp_vqsr_recal_indels.vcf.idx")
    params:
        padding = config['WES']['PADDING'],
        intervals = config['WES']['INTERVALS'],
        other = "-mode INDEL"
    log:
        "logs/gatk_VQSR_indel/{sample}.log"
    benchmark:
        "benchmarks/gatk_VQSR_indel/{sample}.tsv"
    singularity:
        "docker://broadinstitute/gatk:4.2.6.1"
    message:
        "Using machine learning to filter out probable artifacts from the variant callset (indels)"
    shell: 
        "gatk ApplyVQSR -V {input.vcf} -R {input.refgenome} --recal-file {input.recal} --tranches-file {input.tranches} -O {output.vcf} {params.padding} {params.intervals} {params.other} &> {log}"