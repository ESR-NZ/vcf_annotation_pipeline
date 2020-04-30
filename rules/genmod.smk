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