rule VEP:
    input:
        vcf = "annotated/{sample}.vqsr.recal.dbnsfp.vcf"
    output:
        report("annotated/{sample}.vqsr.recal.dbnsfp.vep.vcf_summary.txt", caption = "../report/vep.rst", category = "Variant effect predictor"),
        report("annotated/{sample}.vqsr.recal.dbnsfp.vep.vcf_warnings.txt", caption = "../report/vep.rst", category = "Variant effect predictor"),
        vcf = temp("annotated/{sample}.vqsr.recal.dbnsfp.vep.vcf")
    params:
        genome = expand("{genome}", genome = config["GENOME"]),
        vep = expand("{vep}", vep = config["VEP"]),
        assembly = expand("{build}", build = config["BUILD"]),
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