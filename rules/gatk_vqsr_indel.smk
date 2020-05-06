rule gatk4_VQSR_indel:
    input:
        vcf = "../vcf/{sample}_raw_snps_indels_AS_g.vcf",
        genome = expand("{genome}", genome = config['FILEDIR']['GENOME']),
        recal = "filtered/{sample}_recal_indels",
        recalindex = "filtered/{sample}_recal_indels.idx",
        tranches = "filtered/{sample}_tranches_indels"
    output:
        vcf = temp("filtered/{sample}_tmp_vqsr_recal_indels.vcf"),
        index = temp("filtered/{sample}_tmp_vqsr_recal_indels.vcf.idx")
    params:
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS']),
        other = "-mode SNP -ts-filter-level 99.0 -OVI true"
    log:
        "logs/gatk_vqsr_indels/{sample}.log"
    benchmark:
        report("benchmarks/gatk_vqsr_indels/{sample}.recal.indels", caption = benchmarking.rst, category = "Benchmarking")
    conda:
        "../envs/gatk4.yaml"
    message:
        "Using machine learning to filter out probable artifacts from the variant callset (indels)"
    threads: 4
    shell: 
        "gatk ApplyVQSR -V {input.vcf} -R {input.genome} --recal-file {input.recal} --tranches-file {input.tranches} -O {output.vcf} {params.padding} {params.intervals} {params.other}"