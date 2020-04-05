rule gatk4_VariantRecalibrator_indel:
    input:
        vcf="genotyped/{sample}.genotype.vcf"
    output:
        report="recalibrated/{sample}.recal.indels",
        tranches="recalibrated/{sample}.tranches.indels",
        rscript="recalibrated/{sample}.plots.indels.R"
    log: 
        "logs/gatk_recal_indels/{sample}.log"
    benchmark:
        "benchmarks/gatk_recal_indels/{sample}.recal.indels"
    conda:
        "../envs/gatk4.yaml"
    shell:
        """
        gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
            -R {GENOME} \
            -V {input.vcf} \
            -O {output.report} \
            --tranches-file {output.tranches} \
            --rscript-file {output.rscript} \
            --trust-all-polymorphic \
            --max-gaussians 4 \
            -mode INDEL \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
            -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
            -resource:mills,known=false,training=true,truth=true,prior=12.0 {MILLS} \
            -resource:1000G,known=false,training=true,truth=false,prior=10.0 {INDEL1000G} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {DBSNP}
        """

rule gatk4_VariantRecalibrator_SNP:
    input:
        vcf="genotyped/{sample}.genotype.vcf"
    output:
        report="recalibrated/{sample}.recal.snps",
        tranches="recalibrated/{sample}.tranches.snps",
        rscript="recalibrated/{sample}.plots.snps.R"
    log: 
        "logs/gatk_recal_snps/{sample}.log"
    benchmark:
        "benchmarks/gatk_recal_snps/{sample}.snp.recal"
    conda:
        "../envs/gatk4.yaml"
    shell:
        """
        gatk --java-options "-Xmx3g -Xms3g" VariantRecalibrator \
            -R {GENOME} \
            -V {input.vcf} \
            -O {output.report} \
            --tranches-file {output.tranches} \
            --rscript-file {output.rscript} \
            --trust-all-polymorphic \
            --max-gaussians 6 \
            -mode SNP \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
            -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {HAPMAP} \
            -resource:omni,known=false,training=true,truth=false,prior=12.0 {OMNI} \
            -resource:1000G,known=false,training=true,truth=false,prior=10.0 {SNP1000G} \
            -resource:resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {DBSNP}
        """

rule gatk4_VQSR_indel:
    input:
        vcf="genotyped/{sample}.genotype.vcf",
        recal="recalibrated/{sample}.recal.indels",
        tranches="recalibrated/{sample}.tranches.indels"
    output:
        vcf=temp("recalibrated/{sample}.tmp.vqsr.recal.indels.vcf")
    log:
        "logs/gatk_vqsr_indels/{sample}.log"
    benchmark:
        "benchmarks/gatk_vqsr_indels/{sample}.recal.indels"
    conda:
        "../envs/gatk4.yaml"
    threads: 4
    shell:
        "gatk ApplyVQSR -R {GENOME} -V {input.vcf} -O {output.vcf} --recal-file {input.recal} --tranches-file {input.tranches} --truth-sensitivity-filter-level 99.0 --create-output-variant-index true -mode INDEL"

rule gatk4_VQSR_SNP:
    input:
        vcf="recalibrated/{sample}.tmp.vqsr.recal.indels.vcf",
        recal="recalibrated/{sample}.recal.snps",
        tranches="recalibrated/{sample}.tranches.snps"
    output:
        vcf="recalibrated/{sample}.vqsr.recal.vcf"
    log:
        "logs/gatk_vqsr_snps/{sample}.log"
    benchmark:
        "benchmarks/gatk_vqsr_snps/{sample}.recal.snps"
    conda:
        "../envs/gatk4.yaml"
    threads: 4
    shell:
        "gatk ApplyVQSR -R {GENOME} -V {input.vcf} -O {output.vcf} --recal-file {input.recal} --tranches-file {input.tranches} --truth-sensitivity-filter-level 99.0 --create-output-variant-index true -mode SNP"