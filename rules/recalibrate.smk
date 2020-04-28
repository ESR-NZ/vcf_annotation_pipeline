rule gatk4_VQSR_indel:
    input:
        vcf = "genotyped/{sample}.genotype.vcf",
        recal = "recalibrated/{sample}.recal.indels",
        tranches = "recalibrated/{sample}.tranches.indels"
    output:
        vcf = temp("recalibrated/{sample}.tmp.vqsr.recal.indels.vcf"),
        index = temp("recalibrated/{sample}.tmp.vqsr.recal.indels.vcf.idx")
    params:
        genome = expand("{genome}", genome=config["GENOME"]),
        mode = "INDEL",
        sensitivity = expand("{sensitivity}", sensitivity = config["SENSITIVITY"]),
        varindex = "true"
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
        "gatk ApplyVQSR -V {input.vcf} --recal-file {input.recal} -O {output.vcf} --tranches-file {input.tranches} -R {params.genome} -mode {params.mode} -ts-filter-level {params.sensitivity} -OVI {params.varindex}"

rule gatk4_VQSR_SNP:
    input:
        vcf = "recalibrated/{sample}.tmp.vqsr.recal.indels.vcf",
        recal = "recalibrated/{sample}.recal.snps",
        tranches = "recalibrated/{sample}.tranches.snps"
    output:
        vcf = "recalibrated/{sample}.vqsr.recal.vcf"
    params:
        genome = expand("{genome}", genome = config["GENOME"]),
        mode = "SNP",
        sensitivity = expand("{sensitivity}", sensitivity = config["SENSITIVITY"]),
        varindex = "true"
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
        "gatk ApplyVQSR -V {input.vcf} --recal-file {input.recal} --tranches-file {input.tranches} -O {output.vcf} -R {params.genome} -mode {params.mode} -ts-filter-level {params.sensitivity} -OVI {params.varindex}"