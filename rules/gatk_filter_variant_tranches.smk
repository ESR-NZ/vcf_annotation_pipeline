# Add more resources
rule gatk4_FilterVariantTranches:
    input:
        vcf = "filtered/{sample}_scored.vcf",
        hapmap = expand("{hapmap}", hapmap = config['FILEDIR']['HAPMAP']),
        mills = expand("{mills}", mills = config['FILEDIR']['MILLS'])
    output:
        "filtered/{sample}_filtered.vcf"
    params:
        tdir = expand("{tdir}", tdir = config['TEMPDIR']),
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS']),
        other = "--info-key CNN_2D --snp-tranche 99.95 --indel-tranche 99.4"
    log:
        "logs/gatk_filter_variant_tranches/{sample}.log"
    benchmark:
        report("benchmarks/gatk_filter_variant_tranches/{sample}.gatkfiltervartranch", caption = "../report/benchmarking.rst", category = "Benchmarking")
    singularity:
        "docker://broadinstitute/gatk:4.1.7.0"
    message:
        "Applying tranche filtering to variant calls"
    shell:
        "gatk FilterVariantTranches -V {input.vcf} --resource {input.hapmap} --resource {input.mills} -O {output} --tmp-dir {params.tdir} {params.padding} {params.intervals} {params.other}"