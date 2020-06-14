rule gatk4_FilterVariantTranches:
    input:
        vcf = "filtered/{sample}_scored.vcf" 
    output:
        "filtered/{sample}_filtered.vcf"
    params:
        tdir = expand("{tdir}", tdir = config['TEMPDIR']),
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS']),
        resources = expand ("{resources}", resources = config['FILTERING']['SINGLE']),
        other = "--info-key CNN_2D --snp-tranche 99.95 --indel-tranche 99.4"
    log:
        "logs/gatk_filter_variant_tranches/{sample}.log"
    benchmark:
        "benchmarks/gatk_filter_variant_tranches/{sample}.gatkfiltervartranch"
    singularity:
        "docker://broadinstitute/gatk:4.1.7.0"
    message:
        "Applying tranche filtering to variant calls"
    shell:
        "gatk FilterVariantTranches -V {input.vcf} -O {output} --tmp-dir {params.tdir} {params.padding} {params.intervals} {params.resources} {params.other}"