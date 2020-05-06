rule gatk4_VQSR_SNP:
    input:
        vcf = "filtered/{sample}_tmp_vqsr_recal_indels.vcf",
        genome = expand("{genome}", genome = config['FILEDIR']['GENOME']),
        recal = "filtered/{sample}_recal_snps",
        recalindex = "filtered/{sample}_recal_snps.idx",
        tranches = "filtered/{sample}_tranches_snps"
    output:
        "filtered/{sample}_filtered.vcf"
    params:
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS']),
        other = "-mode SNP -ts-filter-level 99.0 -OVI true"
    log:
        "logs/gatk_vqsr_snps/{sample}.log"
    benchmark:
        "benchmarks/gatk_vqsr_snps/{sample}.recal.snps"
    conda:
        "../envs/gatk4.yaml"
    message:
        "Using machine learning to filter out probable artifacts from the variant callset (snps)"
    threads: 4
    shell:
        "gatk ApplyVQSR -V {input.vcf} -R {input.genome} --recal-file {input.recal} --tranches-file {input.tranches} -O {output} {params.padding} {params.intervals} {params.other}"