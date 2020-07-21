rule gatk_VQSR_indel:
    input:
        vcf = "../vcf/{sample}_raw_snps_indels_AS_g.vcf",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME']),
        recal = "recalibrated/{sample}.recal.indels",
        recalindex = "recalibrated/{sample}.recal.indels.idx",
        tranches = "recalibrated/{sample}.tranches.indels"
    output:
        vcf = temp("filtered/{sample}_tmp_vqsr_recal_indels.vcf"),
        index = temp("filtered/{sample}_tmp_vqsr_recal_indels.vcf.idx")
    params:
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS']),
        other = "-mode INDEL -ts-filter-level 99.0 -OVI true"
    log:
        "logs/gatk_vqsr_indels/{sample}.log"
    benchmark:
        "benchmarks/gatk_vqsr_indels/{sample}.recal.indels"
    conda:
        "../envs/gatk4.yaml"
    message:
        "Using machine learning to filter out probable artifacts from the variant callset (indels)"
    threads: 4
    shell: 
        "gatk ApplyVQSR -V {input.vcf} -R {input.refgenome} --recal-file {input.recal} --tranches-file {input.tranches} -O {output.vcf} {params.padding} {params.intervals} {params.other}"