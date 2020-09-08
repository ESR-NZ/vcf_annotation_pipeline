rule pbrun_cnnscorevariants:
    input:
        vcf = "../../vcf/{sample}_raw_snps_indels.vcf",
        bams = "../../bams/{sample}_recalibrated.bam",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        temp("../results/filtered/{sample}_scored.vcf")
    resources:
        gpu = 1
    log:
        "logs/pbrun_cnnscorevariants/{sample}.log"
    benchmark:
        "benchmarks/pbrun_cnnscorevariants/{sample}.tsv"
    message:
        "Annotating {input.vcf} with scores from a Convolutional Neural Network (CNN) (2D model with pre-trained architecture)"
    shell:
        "pbrun cnnscorevariants --in-vcf {input.vcf} --in-bam {input.bams} --ref {input.refgenome} --out-vcf {output} &> {log}"