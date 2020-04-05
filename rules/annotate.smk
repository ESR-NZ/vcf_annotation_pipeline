rule SnpSift:
    input:
        vcf="recalibrated/{sample}.vqsr.recal.vcf"
    output:
        vcf="annotated/{sample}.vqsr.recal.dbnsfp.vcf"
    log: 
        "logs/snpsift/{sample}.log"
    benchmark:
        "benchmarks/snipsift/{sample}.dbnsfp"
    conda:
        "../envs/dbnsfp.yaml"
    shell:
        "SnpSift -Xmx16g dbnsfp -v -db {DBNSFP} {input.vcf} > {output.vcf}"

rule VEP:
    input:
        vcf="annotated/{sample}.vqsr.recal.dbnsfp.vcf"
    output:
        vcf="annotated/{sample}.vqsr.recal.dbnsfp.vep.vcf"
    log: 
        "logs/vep/{sample}.log"
    benchmark:
        "benchmarks/vep/{sample}.vep"
    conda:
        "../envs/vep.yaml"
    threads: 4
    shell:
        "vep -v --assembly GRCh37 --cache --dir {VEP} --fasta {GENOME} -i {input.vcf} -o {output.vcf} --stats_text --everything --vcf --force_overwrite --offline"

rule GENMOD:
    input:
        vcf="annotated/{sample}.vqsr.recal.dbnsfp.vep.vcf"
    output:
        vcf="annotated/{sample}.vqsr.recal.dbnsfp.vep.genmod.vcf"
    log: 
        "logs/genmod/{sample}.log"
    benchmark:
        "benchmarks/genmod/{sample}.genmod"
    singularity:
        "shub://sirselim/singularity-genmod:latest"
    threads: 4
    shell:
        "genmod annotate {input.vcf} --regions -c {CADD} -o {output.vcf}"