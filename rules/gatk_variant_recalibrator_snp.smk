rule gatk4_VariantRecalibrator_SNP:
    input:
        "../vcf/{sample}_raw_snps_indels_AS_g.vcf"
    output:
        report("recalibrated/{sample}.plots.snps.R.pdf", caption = "../report/recalibration.rst", category = "Recalibration - SNP's"),
        report("recalibrated/{sample}.tranches.snps.pdf", caption = "../report/recalibration.rst", category = "Recalibration - SNP's"),
        recal = temp("recalibrated/{sample}.recal.snps"),
        index = temp("recalibrated/{sample}.recal.snps.idx"),
        tranches = "recalibrated/{sample}.tranches.snps",
        rscript = "recalibrated/{sample}.plots.snps.R"
    params:
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME']),
        resources = expand("{resources}", resources = config['FILTERING']['COHORT']['SNPS']),
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS']),
        other = "-mode SNP --max-gaussians 4 --trust-all-polymorphic"
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
            -V {input} \
            -O {output.recal} \
            --tranches-file {output.tranches} \
            --rscript-file {output.rscript} \
            -R {params.refgenome} \
            {params.resources} \
            {params.padding} \
            {params.intervals} \
            {params.other} \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP
        """