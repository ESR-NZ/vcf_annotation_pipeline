rule SnpSift:
    input:
        vcf="recalibrated/{sample}.vqsr.recal.vcf"
    output:
        vcf="annotated/{sample}.vqsr.recal.dbnsfp.vcf"
    params:
        dbnsfp=expand("{dbnsfp}", dbnsfp=config["dbNSFP"])
    log: 
        "logs/snpsift/{sample}.log"
    benchmark:
        "benchmarks/snipsift/{sample}.dbnsfp"
    conda:
        "../envs/dbnsfp.yaml"
    message:
        "Using the dbNSFP database to annotate variants with functional predictions from multiple algorithms (SIFT, Polyphen2, LRT and MutationTaster, PhyloP and GERP++, etc.)"
    shell:
        "SnpSift -Xmx16g dbnsfp -v -db {params.dbnsfp} {input.vcf} > {output.vcf}"

rule VEP:
    input:
        vcf="annotated/{sample}.vqsr.recal.dbnsfp.vcf"
    output:
        vcf="annotated/{sample}.vqsr.recal.dbnsfp.vep.vcf"
    params:
        genome=expand("{genome}", genome=config["GENOME"]),
        vep=expand("{vep}", vep=config["VEP"])
    log: 
        "logs/vep/{sample}.log"
    benchmark:
        "benchmarks/vep/{sample}.vep"
    conda:
        "../envs/vep.yaml"
    threads: 4
    message:
        "Using the VEP database to determine the effect of the variants"
    shell:
        "vep -v --assembly GRCh37 --cache --dir {params.vep} --fasta {params.genome} -i {input.vcf} -o {output.vcf} --stats_text --everything --vcf --force_overwrite --offline"

rule GENMOD:
    input:
        vcf="annotated/{sample}.vqsr.recal.dbnsfp.vep.vcf"
    output:
        vcf="annotated/{sample}.vqsr.recal.dbnsfp.vep.genmod.vcf"
    params:
        cadd=expand("{cadd}", cadd=config["CADD"])
    log: 
        "logs/genmod/{sample}.log"
    benchmark:
        "benchmarks/genmod/{sample}.genmod"
    singularity:
        "shub://sirselim/singularity-genmod:latest"
    threads: 4
    message:
        "Using the CADD database to annotate the variants with deleteriousness scores"
    shell:
        "genmod annotate {input.vcf} --regions -c {params.cadd} -o {output.vcf}"
