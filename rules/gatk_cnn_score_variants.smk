rule gatk4_CNNScoreVariants:
    input:
        vcf = "../vcf/{sample}_raw_snps_indels_AS_g.vcf",
        bams = "../bams/{sample}_recal.bam",
        genome = expand("{genome}", genome = config['FILEDIR']['GENOME'])
    output:
        "filtered/{sample}_scored.vcf"
    params:
        tdir = expand("{tdir}", tdir = config['TEMPDIR']),
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS']),
        other = "-tensor-type read_tensor"
    log:
        "logs/gatk_score_variants/{sample}.log"
    benchmark:
        "benchmarks/gatk_score_variants/{sample}.gatkscorevariants"
    singularity:
        "docker://broadinstitute/gatk:4.1.7.0"
    threads: 12
    message:
        "Annotating vcf with scores from a Convolutional Neural Network (CNN) (2D model with pre-trained architecture)"
    shell:
        """
        gatk --java-options "-Xmx64g -Xms64g" CNNScoreVariants \
            -V {input.vcf} -I {input.bams} -R {input.genome} -O {output} --inter-op-threads {threads} --intra-op-threads {threads} --tmp-dir {params.tdir} {params.padding} {params.intervals} {params.other}
        """