rule vep:
    input:
        vcf = "annotated/{sample}_filtered_dbnsfp.vcf",
        genome = expand("{genome}", genome = config['FILEDIR']['GENOME']),
        vep = expand("{vep}", vep = config['FILEDIR']['VEP']),
    output:
        report("annotated/{sample}_filtered_dbnsfp_vep.vcf_summary.txt", caption = "../report/vep.rst", category = "Variant effect predictor"),
        report("annotated/{sample}_filtered_dbnsfp_vep.vcf_warnings.txt", caption = "../report/vep.rst", category = "Variant effect predictor"),
        vcf = temp("annotated/{sample}_filtered_dbnsfp_vep.vcf")
    params:
        assembly = expand("{build}", build = config['BUILD']),
        extra = " --format=vcf --cache --stats_text --everything --vcf --force_overwrite --offline -v"
    log: 
        "logs/vep/{sample}.log"
    benchmark:
        report("benchmarks/vep/{sample}.vep", caption = "../report/benchmarking.rst", category = "Benchmarking")
    conda:
        "../envs/vep.yaml"
    threads: 4
    message:
        "Using the VEP database to determine the effect of the variants"
    shell:
        "vep -i {input.vcf} --fasta {input.genome} --dir {input.vep} -o {output.vcf} --assembly {params.assembly} {params.extra}"