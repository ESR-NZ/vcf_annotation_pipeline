rule gatk_VariantRecalibrator_snp:
    input:
        vcf = "../../human_genomics_pipeline/results/called/{sample}_raw_snps_indels.vcf",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        report("../results/filtered/{sample}_plots_snps.R.pdf", caption = "../report/recalibration.rst", category = "Recalibration - SNP's"),
        report("../results/filtered/{sample}_tranches_snps.pdf", caption = "../report/recalibration.rst", category = "Recalibration - SNP's"),
        recal = temp("../results/filtered/{sample}_recal_snps"),
        index = temp("../results/filtered/{sample}_recal_snps.idx"),
        tranches = temp("../results/filtered/{sample}_tranches_snps"),
        rscript = "../results/filtered/{sample}_plots_snps.R"
    params:
        maxmemory = expand('"-Xmx{maxmemory}"', maxmemory = config['MAXMEMORY']),
        snptranche = expand("-tranche {snptranche}", snptranche = config['TRANCHE']['SNPS']),
        padding = config['WES']['PADDING'],
        intervals = config['WES']['INTERVALS'],
        resources = get_cohort_snp_filtering_command,
        other = "-mode SNP -an QD -an MQ -an MQRankSum -an ReadPosRankSum"
    log: 
        "logs/gatk_VariantRecalibrator_snp/{sample}.log"
    benchmark:
        "benchmarks/gatk_VariantRecalibrator_snp/{sample}.tsv"
    conda:
        "../envs/gatk4.yaml"
    message:
        "Building a recalibration model to score variant quality (snps)"
    shell:
        """
        gatk --java-options {params.maxmemory} VariantRecalibrator \
            -V {input.vcf} \
            -O {output.recal} \
            --tranches-file {output.tranches} \
            --rscript-file {output.rscript} \
            -R {input.refgenome} \
            {params.snptranche} \
            {params.resources} \
            {params.padding} \
            {params.intervals} \
            {params.other} \
            &> {log}
        """
