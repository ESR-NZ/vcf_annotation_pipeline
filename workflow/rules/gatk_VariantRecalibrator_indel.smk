rule gatk_VariantRecalibrator_indel:
    input:
        vcf = "../../human_genomics_pipeline/results/called/{sample}_raw_snps_indels.g.vcf",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        report("../results/filtered/{sample}_plots_indels.R.pdf", caption = "../report/recalibration.rst", category = "Recalibration - Indels"),
        recal = temp("../results/filtered/{sample}_recal_indels"),
        index = temp("../results/filtered/{sample}_recal_indels.idx"),
        tranches = temp("../results/filtered/{sample}_tranches_indels"),
        rscript = "../results/filtered/{sample}_plots_indels.R"
    params:
        maxmemory = expand('"-Xmx{maxmemory}"', maxmemory = config['MAXMEMORY']),
        indeltranche = expand("-tranche {indeltranche}", indeltranche = config['TRANCHE']['INDELS']),
        padding = config['WES']['PADDING'],
        intervals = config['WES']['INTERVALS'],
        resources = get_cohort_indel_filtering_command,
        other = "-mode INDEL -an QD -an MQ -an MQRankSum -an ReadPosRankSum"
    log: 
        "logs/gatk_VariantRecalibrator_indel/{sample}.log"
    benchmark:
        "benchmarks/gatk_VariantRecalibrator_indel/{sample}.tsv"
    conda:
        "../envs/gatk4.yaml"
    message:
        "Building a recalibration model to score variant quality (indels)"
    shell:
        """
        gatk --java-options {params.maxmemory} VariantRecalibrator \
            -V {input.vcf} \
            -O {output.recal} \
            --tranches-file {output.tranches} \
            --rscript-file {output.rscript} \
            -R {input.refgenome} \
            {params.indeltranche} \
            {params.resources} \
            {params.padding} \
            {params.intervals} \
            {params.other} \
            &> {log}
        """