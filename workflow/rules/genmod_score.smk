rule genmod_score:
    input:
        vcf = get_genmod_score_input,
        refgenome = config['REFGENOME']
    params:
        get_genmod_score_params
    output:
        protected(get_genmod_score_output)
    log: 
        "logs/genmod_score/{sample}.log"
    benchmark:
        "benchmarks/genmod_score/{sample}.tsv"
    conda:
        "../envs/genmod.yaml"
    message:
        "Scoring the variants in {input.vcf} based on several annotations"
    shell:
        "genmod score {input.vcf} {params} -o {output} &> {log}"