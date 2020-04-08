rule SnpSift:
    input:
        vcf = "recalibrated/{sample}.vqsr.recal.vcf"
    output:
        vcf = "annotated/{sample}.vqsr.recal.dbnsfp.vcf"
    params:
        dbnsfp = expand("{dbnsfp}", dbnsfp = config["dbNSFP"])
    log: 
        "logs/snpsift/{sample}.log"
    benchmark:
        "benchmarks/snipsift/{sample}.dbnsfp"
    conda:
        "../envs/dbnsfp.yaml"
    message:
        "Using the dbNSFP database to annotate variants with functional predictions from multiple algorithms (SIFT, Polyphen2, LRT and MutationTaster, PhyloP and GERP++, etc.)"
    shell:
        "SnpSift -Xmx16g dbnsfp {input.vcf} > {output.vcf} -db {params.dbnsfp} -v"

rule VEP:
    input:
        vcf = "annotated/{sample}.vqsr.recal.dbnsfp.vcf"
    output:
        vcf = "annotated/{sample}.vqsr.recal.dbnsfp.vep.vcf"
    params:
        genome = expand("{genome}", genome = config["GENOME"]),
        vep = expand("{vep}", vep = config["VEP"]),
        assembly = "GRCh37",
        extra = " --format=vcf --cache --stats_text --everything --vcf --force_overwrite --offline"
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
        "vep -i {input.vcf} -o {output.vcf} --fasta {params.genome} -v --dir {params.vep} --assembly {params.assembly} {params.extra}"

rule GENMOD:
    input:
        vcf = "annotated/{sample}.vqsr.recal.dbnsfp.vep.vcf"
    output:
        vcf = "annotated/{sample}.vqsr.recal.dbnsfp.vep.genmod.vcf"
    params:
        cadd = expand("{cadd}", cadd = config["CADD"])
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
        "genmod annotate {input.vcf} -o {output.vcf} -c {params.cadd} --regions"
