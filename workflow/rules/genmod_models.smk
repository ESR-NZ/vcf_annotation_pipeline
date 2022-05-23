rule genmod_models:
    input:
        vcf = "../results/annotated/{sample}_filtered_dbnsfp_vep_cadd_dbsnp_posteriors_denovo.vcf",
        pedigree = "../../pedigrees/{sample}_pedigree.ped",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        protected("../results/annotated/{sample}_filtered_annotated.vcf")
    params:
        "--vep"
    log: 
        "logs/genmod_models/{sample}.log"
    benchmark:
        "benchmarks/genmod_models/{sample}.tsv"
    conda:
        "../envs/genmod.yaml"
    threads: config['THREADS']
    message:
        "Annotating {input.vcf} with patterns of inheritance"
    shell:
        "genmod models {input.vcf} -f {input.pedigree} -o {output} {params} -p {threads} &> {log}"
