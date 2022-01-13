rule gatk_FilterVariantTranches:
    input:
        "../results/filtered/{sample}_scored.vcf" 
    output:
        protected("../results/filtered/{sample}_filtered.vcf")
    params:
        snptranche = expand("--snp-tranche {snptranche}", snptranche = config['FILTERING']['TRANCHE']['SNPS']),
        indeltranche = expand("--indel-tranche {indeltranche}", indeltranche = config['FILTERING']['TRANCHE']['INDELS']),
        tdir = config['TEMPDIR'],
        padding = config['WES']['PADDING'],
        intervals = config['WES']['INTERVALS'],
        resources = config['FILTERING']['SINGLE'],
        other = "--info-key CNN_2D"
    log:
        "logs/gatk_FilterVariantTranches/{sample}.log"
    benchmark:
        "benchmarks/gatk_FilterVariantTranches/{sample}.tsv"
    singularity:
        "docker://broadinstitute/gatk:4.1.7.0"
    message:
        "Applying tranche filtering to variant calls in {input}"
    shell:
        "gatk FilterVariantTranches -V {input} -O {output} {params.snptranche} {params.indeltranche} --tmp-dir {params.tdir} {params.padding} {params.intervals} {params.resources} {params.other} &> {log}"