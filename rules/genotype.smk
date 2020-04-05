rule gatk4_GenotypeGVCFs:
    input:
        vcf=expand("{sampledir}{sample}.raw.snps.indels.AS.g.vcf", sampledir=SAMPLEDIR, sample=SAMPLES)
    output:
        vcf="genotyped/{sample}.genotype.vcf"
    log:
        "logs/gatk_genotype/{sample}.log"
    benchmark:
        "benchmarks/gatk_genotype/{sample}.genotype"
    conda:
        "../envs/gatk4.yaml"
    shell:
        """
        gatk --java-options "-Xmx64g -Xms64g" GenotypeGVCFs \
            -R {GENOME} \
            -V {input.vcf} \
            -O {output.vcf} \
            -D {DBSNP} \
            -G StandardAnnotation -G AS_StandardAnnotation
        """