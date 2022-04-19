rule setup:
    input:
        vcf = "../../human_genomics_pipeline/results/called/{sample}_raw_snps_indels.g.vcf",
        index = "../../human_genomics_pipeline/results/called/{sample}_raw_snps_indels.g.vcf.idx"
    output:
        vcf = "../../human_genomics_pipeline/results/called/{sample}_raw_snps_indels.vcf",
        index = "../../human_genomics_pipeline/results/called/{sample}_raw_snps_indels.vcf.idx"
    message:
        "Renaming {input.vcf} and {input.index} file extensions for compatability with vcf_annotation_pipeline version 2.0.0 onwards"
    shell:
        """
        mv -f {input.vcf} "`echo {input.vcf} | sed 's/\.g//'`"
        mv -f {input.index} "`echo {input.index} | sed 's/\.g//'`"
        """