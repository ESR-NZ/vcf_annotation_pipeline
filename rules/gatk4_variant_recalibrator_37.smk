rule gatk4_VariantRecalibrator_indel:
    input:
        vcf = "genotyped/{sample}.genotype.vcf"
    output:
        report("recalibrated/{sample}.plots.indels.R.pdf", caption = "../report/recalibration.rst", category = "Recalibration - Indels"),
        recal = temp("recalibrated/{sample}.recal.indels"),
        index = temp("recalibrated/{sample}.recal.indels.idx"),
        tranches = "recalibrated/{sample}.tranches.indels",
        rscript = "recalibrated/{sample}.plots.indels.R"

    params:
        genome = expand("{genome}", genome = config["GENOME"]),
        mills = expand("{mills}", mills = config["MILLS"]),
        indel1000g = expand("{indel1000g}", indel1000g = config["INDEL1000G"]),
        dbsnp = expand("{dbsnp}", dbsnp = config["dbSNP"]),
        mode = "INDEL",
        gaussians = "4"
    log: 
        "logs/gatk_recal_indels/{sample}.log"
    benchmark:
        "benchmarks/gatk_recal_indels/{sample}.recal.indels"
    conda:
        "../envs/gatk4.yaml"
    message:
        "Building a recalibration model to score variant quality for indels"
    shell:
        """
        gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
            -V {input.vcf} \
            -O {output.recal} \
            --tranches-file {output.tranches} \
            --rscript-file {output.rscript} \
            -R {params.genome} \
            -resource:mills,known=false,training=true,truth=true,prior=12.0 {params.mills} \
            -resource:1000G,known=false,training=true,truth=false,prior=10.0 {params.indel1000g} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} \
            -mode {params.mode} \
            --max-gaussians {params.gaussians} \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
            -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP
        """

rule gatk4_VariantRecalibrator_SNP:
    input:
        vcf = "genotyped/{sample}.genotype.vcf"
    output:
        report("recalibrated/{sample}.plots.snps.R.pdf", caption = "../report/recalibration.rst", category = "Recalibration - SNP's"),
        report("recalibrated/{sample}.tranches.snps.pdf", caption = "../report/recalibration.rst", category = "Recalibration - SNP's"),
        recal = temp("recalibrated/{sample}.recal.snps"),
        index = temp("recalibrated/{sample}.recal.snps.idx"),
        tranches = "recalibrated/{sample}.tranches.snps",
        rscript = "recalibrated/{sample}.plots.snps.R"

    params:
        genome = expand("{genome}", genome = config["GENOME"]),
        hapmap = expand("{hapmap}", hapmap = config["HAPMAP"]),
        omni = expand("{omni}", omni = config["OMNI"]),
        snp1000g = expand("{snp1000g}", snp1000g = config["SNP1000G"]),
        dbsnp = expand("{dbsnp}", dbsnp = config["dbSNP"]),
        mode = "SNP",
        gaussians = "4"
    log: 
        "logs/gatk_recal_snps/{sample}.log"
    benchmark:
        "benchmarks/gatk_recal_snps/{sample}.snp.recal"
    conda:
        "../envs/gatk4.yaml"
    message:
        "Building a recalibration model to score variant quality for snps"
    shell:
        """
        gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
            -V {input.vcf} \
            -O {output.recal} \
            --tranches-file {output.tranches} \
            --rscript-file {output.rscript} \
            -R {params.genome} \
            -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.hapmap} \
            -resource:omni,known=false,training=true,truth=false,prior=12.0 {params.omni} \
            -resource:1000G,known=false,training=true,truth=false,prior=10.0 {params.snp1000g} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} \
            -mode {params.mode} \
            --max-gaussians {params.gaussians} \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP
        """