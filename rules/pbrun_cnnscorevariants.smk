rule pbrun_cnnscorevariants:
    input:
        vcf = "../vcf/{sample}_raw_snps_indels_AS_g.vcf",
        bams = "../bams/{sample}_bwa_recal.bam",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        "filtered/{sample}_scored.vcf"
    log:
        "logs/pbrun_score_variants/{sample}.log"
    benchmark:
        "benchmarks/pbrun_score_variants/{sample}.pbrunscorevariants"
    message:
        "Annotating vcf with scores from a Convolutional Neural Network (CNN) (2D model with pre-trained architecture)"
    shell:
        "pbrun cnnscorevariants --ref {input.refgenome} --in-bam {input.bams} --in-vcf {input.vcf} --out-vcf {output}"