rule gatk_VariantAnnotator_PossibleDeNovo:
    input:
        vcf = "../results/annotated/{sample}_filtered_dbnsfp_vep_cadd_dbsnp_posteriors.vcf",
        pedigree = "../../pedigrees/{sample}_pedigree.ped",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        vcf = temp("../results/annotated/{sample}_filtered_dbnsfp_vep_cadd_dbsnp_posteriors_denovo.vcf"),
        index = temp("../results/annotated/{sample}_filtered_dbnsfp_vep_cadd_dbsnp_posteriors_denovo.vcf.idx")
    params:
        maxmemory = expand('"-Xmx{maxmemory}"', maxmemory = config['MAXMEMORY']),
        other = "-A PossibleDeNovo"
    log: 
        "logs/gatk_VariantAnnotator_PossibleDeNovo/{sample}.log"
    benchmark:
        "benchmarks/gatk_VariantAnnotator_PossibleDeNovo/{sample}.tsv"
    singularity:
        "docker://broadinstitute/gatk:4.2.6.1"
    message:
        "Marking variants in {input.vcf} that are possible denovo mutations"
    shell:
        "gatk --java-options {params.maxmemory} VariantAnnotator -R {input.refgenome} -V {input.vcf} -O {output.vcf} --pedigree {input.pedigree} {params.other} &> {log}"