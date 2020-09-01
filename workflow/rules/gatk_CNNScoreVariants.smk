rule gatk_CNNScoreVariants:
    input:
        vcf = "../../vcf/{sample}_raw_snps_indels.vcf",
        bams = "../../bams/{sample}_recalibrated.bam",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        vcf = temp("../results/filtered/{sample}_scored.vcf"),
        index = temp("../results/filtered/{sample}_scored.vcf.idx")
    params:
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS']),
        other = "-tensor-type read_tensor"
    log:
        "logs/gatk_CNNScoreVariants/{sample}.log"
    benchmark:
        "benchmarks/gatk_CNNScoreVariants/{sample}.tsv"
    singularity:
        "docker://broadinstitute/gatk:4.1.7.0"
    threads: 12
    message:
        "Annotating {input.vcf} with scores from a Convolutional Neural Network (CNN) (2D model with pre-trained architecture)"
    shell:
        """
        gatk --java-options "-Xmx64g -Xms64g" CNNScoreVariants \
            -V {input.vcf} -I {input.bams} -R {input.refgenome} -O {output.vcf} --intra-op-threads {threads} {params.padding} {params.intervals} {params.other} &> {log}
        """