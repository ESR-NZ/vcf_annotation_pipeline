rule gatk4_VariantRecalibrator_indel:
    input:
        "../vcf/{sample}_raw_snps_indels_AS_g.vcf",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME']),
        resources = expand("{resources}", resources = config['FILTERING']['COHORT']['INDELS'])
    output:
        report("recalibrated/{sample}.plots.indels.R.pdf", caption = "../report/recalibration.rst", category = "Recalibration - Indels"),
        recal = temp("recalibrated/{sample}.recal.indels"),
        index = temp("recalibrated/{sample}.recal.indels.idx"),
        tranches = "recalibrated/{sample}.tranches.indels",
        rscript = "recalibrated/{sample}.plots.indels.R"
    params:
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS']),
        other = "-mode INDEL --max-gaussians 4 --trust-all-polymorphic"
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
            -V {input} \
            -O {output.recal} \
            --tranches-file {output.tranches} \
            --rscript-file {output.rscript} \
            -R {params.refgenome} \
            {params.resources} \
            {params.padding} \
            {params.intervals} \
            {params.other} \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
            -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP
        """