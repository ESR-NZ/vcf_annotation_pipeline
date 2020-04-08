rule gatk4_GenotypeGVCFs:
    input:
        vcf = "../vcf/{sample}.raw.snps.indels.AS.g.vcf"
    output:
        vcf = "genotyped/{sample}.genotype.vcf"
    params:
        genome = expand("{genome}", genome = config["GENOME"]),
        dbsnp = expand("{dbsnp}", dbsnp = config["dbSNP"])
    log:
        "logs/gatk_genotype/{sample}.log"
    benchmark:
        "benchmarks/gatk_genotype/{sample}.genotype"
    conda:
        "../envs/gatk4.yaml"
    message:
        "Running initial genotyping on vcf data from this directory: {input}"
    shell:
        """
        gatk --java-options "-Xmx64g -Xms64g" GenotypeGVCFs \
            -V {input.vcf} \
            -O {output.vcf} \
            -R {params.genome} \
            -D {params.dbsnp} \
            -G StandardAnnotation \
            -G AS_StandardAnnotation
        """