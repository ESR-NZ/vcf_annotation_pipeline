rule gatk_CNNScoreVariants:
    input:
        vcf = "../../human_genomics_pipeline/results/called/{sample}_raw_snps_indels.vcf",
        bams = "../../human_genomics_pipeline/results/mapped/{sample}_recalibrated.bam",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        vcf = temp("../results/filtered/{sample}_scored.vcf"),
        index = temp("../results/filtered/{sample}_scored.vcf.idx")
    params:
        maxmemory = expand('"-Xmx{maxmemory}"', maxmemory = config['MAXMEMORY']),
        padding = config['WES']['PADDING'],
        intervals = config['WES']['INTERVALS'])
        other = "-tensor-type read_tensor"
    log:
        "logs/gatk_CNNScoreVariants/{sample}.log"
    benchmark:
        "benchmarks/gatk_CNNScoreVariants/{sample}.tsv"
    singularity:
        "docker://broadinstitute/gatk:4.1.7.0"
    threads: config['THREADS']
    message:
        "Annotating {input.vcf} with scores from a Convolutional Neural Network (CNN) (2D model with pre-trained architecture)"
    shell:
        "gatk --java-options {params.maxmemory} CNNScoreVariants -V {input.vcf} -I {input.bams} -R {input.refgenome} -O {output.vcf} --intra-op-threads {threads} {params.padding} {params.intervals} {params.other} &> {log}"