rule genmod_annotate_CADD:
    input:
        vcf = "../results/annotated/{sample}_filtered_dbnsfp_vep.vcf.gz",
        cadd = expand("{cadd}", cadd = config['CADD'])
    output:
        temp("../results/annotated/{sample}_filtered_dbnsfp_vep_cadd.vcf")
    log: 
        "logs/genmod_annotate_CADD/{sample}.log"
    benchmark:
        "benchmarks/genmod_annotate_CADD/{sample}.tsv"
    singularity:
        "docker://clinicalgenomics/genmod:3.7.4"
    message:
        "Using genmod to annotate {input.vcf} with CADD"
    shell:
        "genmod annotate {input.vcf} -c {input.cadd} -o {output} &> {log}"
