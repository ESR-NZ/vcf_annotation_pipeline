rule SnpSift:
    input:
        vcf = "recalibrated/{sample}.vqsr.recal.vcf"
    output:
        vcf = temp("annotated/{sample}.vqsr.recal.dbnsfp.vcf")
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