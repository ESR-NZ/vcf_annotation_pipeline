rule gatk_VQSR_snp:
    input:
        vcf = "../results/filtered/{sample}_tmp_vqsr_recal_indels.vcf",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME']),
        recal = "../results/filtered/{sample}_recal_snps",
        recalindex = "../results/filtered/{sample}_recal_snps.idx",
        tranches = "../results/filtered/{sample}_tranches_snps"
    output:
        protected("../results/filtered/{sample}_filtered.vcf")
    params:
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS']),
        other = "-mode SNP -ts-filter-level 99.0 -OVI true"
    log:
        "logs/gatk_VQSR_snp/{sample}.log"
    benchmark:
        "benchmarks/gatk_VQSR_snp/{sample}.tsv"
    conda:
        "../envs/gatk4.yaml"
    message:
        "Using machine learning to filter out probable artifacts from the variant callset (snps)"
    threads: 4
    shell:
        "gatk ApplyVQSR -V {input.vcf} -R {input.refgenome} --recal-file {input.recal} --tranches-file {input.tranches} -O {output} {params.padding} {params.intervals} {params.other} &> {log}"