rule genmod_score_cohort:
    input:
        vcf = get_genmod_score_input,
        refgenome = config['REFGENOME']
    params:
        get_genmod_score_params
    output:
        "../results/readyforscout/{sample}_filtered_annotated_readyforscout.vcf"
    log: 
        "logs/genmod_score/{sample}.log"
    benchmark:
        "benchmarks/genmod_score/{sample}.tsv"
    singularity:
        "docker://clinicalgenomics/genmod:3.7.4"
    message:
        "Scoring the variants in {input.vcf} based on several annotations"
    shell:
        "genmod score {input.vcf} {params} -o {output} &> {log}"