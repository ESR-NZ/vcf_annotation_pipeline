rule setup:
    input:
        "../../human_genomics_pipeline/results/called/{sample}_raw_snps_indels.g.vcf"
    output:
        "../../human_genomics_pipeline/results/called/{sample}_raw_snps_indels.vcf"
    message:
        "Renaming {input} file extension for compatability of vcf_annotation_pipeline with human_genomics_pipeline v2.0.0 and earlier"
    shell:
        "mv -f {input} {output}"