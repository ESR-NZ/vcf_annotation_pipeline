rule gatk_CalculateGenotypePosteriors:
    input:
        vcf = "../results/annotated/{sample}_filtered_dbnsfp_vep_cadd_dbsnp.vcf",
        pedigree = "../../pedigrees/{sample}_pedigree.ped",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        vcf = temp("../results/annotated/{sample}_filtered_dbnsfp_vep_cadd_dbsnp_posteriors.vcf"),
        index = temp("../results/annotated/{sample}_filtered_dbnsfp_vep_cadd_dbsnp_posteriors.vcf.idx")
    params:
        maxmemory = expand('"-Xmx{maxmemory}"', maxmemory = config['MAXMEMORY']),
        other = "--skip-population-priors"
    log: 
        "logs/gatk_CalculateGenotypePosteriors/{sample}.log"
    benchmark:
        "benchmarks/gatk_CalculateGenotypePosteriors/{sample}.tsv"
    singularity:
        "docker://broadinstitute/gatk:4.2.6.1"
    message:
        "Calculating genotype posterior probabilities given family and/or known population genotypes for {input.vcf}"
    shell:
        "gatk --java-options {params.maxmemory} CalculateGenotypePosteriors -R {input.refgenome} -V {input.vcf} -O {output.vcf} --pedigree {input.pedigree} {params.other} &> {log}"