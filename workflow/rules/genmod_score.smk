if config['DATA'] == "Single" or config['DATA'] == 'single':
    infile = "../results/readyforscout/{sample}_filtered_annotated_multiallelicsites.vcf.gz"
    outfile = "../results/readyforscout/{sample}_filtered_annotated_readyforscout.vcf.gz"
    params = "--score_config scripts/score_single.ini"
elif config['DATA'] == "Cohort" or config['DATA'] == 'cohort':
    infile = "../results/readyforscout/{sample}_filtered_annotated_multiallelicsites_probandonly.vcf"
    outfile = "../results/readyforscout/{sample}_filtered_annotated_readyforscout.vcf"
    params = "--score_config scripts/score_cohort.ini --family_file ../../pedigrees/{sample}_pedigree.ped"

rule genmod_score:
    input:
        vcf = infile,
        refgenome = config['REFGENOME']
    params:
        params
    output:
        protected(outfile)
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