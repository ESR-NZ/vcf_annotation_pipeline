rule gatk_FilterVariantTranches:
    input:
        "../results/filtered/{sample}_scored.vcf" 
    output:
        protected("../results/filtered/{sample}_filtered.vcf")
    params:
        snptranche = expand("--snp-tranche {snptranche}", snptranche = config['TRANCHE']['SNPS']),
        indeltranche = expand("--indel-tranche {indeltranche}", indeltranche = config['TRANCHE']['INDELS']),
        tdir = config['TEMPDIR'],
        padding = config['WES']['PADDING'],
        intervals = config['WES']['INTERVALS'],
        resources = get_single_filtering_command,
        other = "--info-key CNN_2D"
    log:
        "logs/gatk_FilterVariantTranches/{sample}.log"
    benchmark:
        "benchmarks/gatk_FilterVariantTranches/{sample}.tsv"
    singularity:
        "docker://broadinstitute/gatk:4.2.6.1"
    message:
        "Applying tranche filtering to variant calls in {input}"
    shell:
        "gatk FilterVariantTranches -V {input} -O {output} {params.snptranche} {params.indeltranche} --tmp-dir {params.tdir} {params.padding} {params.intervals} {params.resources} {params.other} &> {log}"