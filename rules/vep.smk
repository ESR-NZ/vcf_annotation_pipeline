rule vep:
    input:
        vcf = "annotated/{sample}_filtered_dbnsfp.vcf",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME']),
        vep = expand("{vep}", vep = config['VEP'])
    output:
        report("annotated/{sample}_filtered_dbnsfp_vep.vcf_summary.txt", caption = "../report/vep.rst", category = "Variant effect predictor"),
        vcf = temp("annotated/{sample}_filtered_dbnsfp_vep.vcf")
    params:
        build = expand("{build}", build = config['BUILD']),
        other = "--format=vcf --cache --stats_text --everything --vcf --force_overwrite --offline -v"
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
        "vep -i {input.vcf} --fasta {input.refgenome} --dir {input.vep} -o {output.vcf} --assembly {params.build} {params.other} --fork {threads}"