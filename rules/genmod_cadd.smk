rule genmod_CADD:
    input:
        vcf = "annotated/{sample}_filtered_dbnsfp_vep.vcf",
        cadd = expand("{cadd}", cadd = config['FILEDIR']['CADD'])
    output:
        vcf = "annotated/{sample}_filtered_annotated.vcf"
    params:
        "--regions"
    log: 
        "logs/genmod_cadd/{sample}.log"
    benchmark:
        report("benchmarks/genmod_cadd/{sample}.genmodcadd", caption = benchmarking.rst, category = "Benchmarking")
    singularity:
        "shub://sirselim/singularity-genmod:latest"
    threads: 4
    message:
        "Using the CADD database to annotate the variants with deleteriousness scores"
    shell:
        "genmod annotate {input.vcf} -c {input.cadd} -o {output.vcf} {params}"