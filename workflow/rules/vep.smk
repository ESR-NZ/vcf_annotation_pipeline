rule vep:
    input:
        vcf = "../results/annotated/{sample}_filtered_dbnsfp.vcf.gz",
        index = "../results/annotated/{sample}_filtered_dbnsfp.vcf.gz.tbi",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME']),
        vep = expand("{vep}", vep = config['VEP'])
    output:
        report("../results/annotated/{sample}_filtered_dbnsfp_vep.vcf.gz_summary.txt", caption = "../report/vep.rst", category = "Variant effect predictor"),
        vcf = temp("../results/annotated/{sample}_filtered_dbnsfp_vep.vcf.gz")
    params:
        build = config['BUILD'],
        other = "--compress_output bgzip --offline --stats_text --everything --vcf --force_overwrite"
    log: 
        "logs/vep/{sample}.log"
    benchmark:
        "benchmarks/vep/{sample}.tsv"
    singularity:
        "docker://ensemblorg/ensembl-vep:release_106.1"
    threads: config['THREADS']
    message:
        "Using the VEP database to determine the effect of the variants in {input.vcf}"
    shell:
        "vep -i {input.vcf} --fasta {input.refgenome} --dir {input.vep} -o {output.vcf} --assembly {params.build} {params.other} --fork {threads} &> {log}"